"""The batched-extinction kernel, and the reduction that feeds off it.

solve_extinction_batch replaces a loop over solve_grid, so the thing worth
pinning down is that it *is* the same calculation: the same six contractions,
reassociated so the A_v-independent halves are computed once and the dependent
ones become single large products. These tests hold it against solve_grid at
every extinction, and hold the block reduction against a full sort.
"""

import numpy as np
import pytest

from superfit.SF_functions import (
    _top_k_per_extinction,
    assemble_results,
    build_tasks,
    merge_blocks,
    solve_extinction_batch,
    solve_grid,
)


def problem(rng, n_sn=23, n_gal=5, n_lam=140, n_av=7, ragged=True):
    """A small fit with realistic ragged coverage and a NaN or two."""

    S = np.abs(rng.normal(1.0, 0.3, (n_sn, n_lam)))
    G = np.abs(rng.normal(1.0, 0.3, (n_gal, n_lam)))

    if ragged:
        for s in range(n_sn):
            S[s, : rng.integers(0, n_lam // 4)] = np.nan
            S[s, n_lam - rng.integers(0, n_lam // 4) :] = np.nan
        for g in range(0, n_gal, 2):
            G[g, : rng.integers(0, n_lam // 8)] = np.nan

    obj = np.abs(rng.normal(1.0, 0.2, n_lam))
    obj[: n_lam // 40] = np.nan
    sig = np.full(n_lam, 0.05)

    alam = np.linspace(1.0, 3.0, n_lam)
    avs = np.linspace(-2.0, 2.0, n_av)
    red = 10.0 ** (-0.4 * avs[:, None] * alam[None, :])
    return S, G, obj, sig, red


class TestAgreesWithSolveGrid:
    @pytest.mark.parametrize("weighted", [False, True])
    def test_every_extinction_matches_one_at_a_time(self, rng, weighted):
        S, G, obj, sig, red = problem(rng)

        b, d, chi2, times = solve_extinction_batch(
            S, G, obj, sig, red, weighted=weighted
        )

        for a in range(red.shape[0]):
            b1, d1, c1, t1 = solve_grid(
                (S * red[a])[None, :, :], G[:, None, :], obj, sig, weighted=weighted
            )
            np.testing.assert_array_equal(t1, times[a])
            for new, old in ((b[a], b1), (d[a], d1), (chi2[a], c1)):
                # The same sums in a different association order.
                np.testing.assert_allclose(new, old, rtol=1e-11, atol=1e-11)
                assert (np.isfinite(new) == np.isfinite(old)).all()

    def test_the_shapes_carry_the_extinction_axis_first(self, rng):
        S, G, obj, sig, red = problem(rng, n_sn=9, n_gal=4, n_av=3)

        b, d, chi2, times = solve_extinction_batch(S, G, obj, sig, red)

        for array in (b, d, chi2, times):
            assert array.shape == (3, 4, 9)

    def test_a_single_extinction_is_not_a_special_case(self, rng):
        S, G, obj, sig, red = problem(rng, n_av=1)

        _, _, chi2, _ = solve_extinction_batch(S, G, obj, sig, red)
        _, _, c1, _ = solve_grid((S * red[0])[None, :, :], G[:, None, :], obj, sig)

        np.testing.assert_allclose(chi2[0], c1, rtol=1e-11, atol=1e-11)

    def test_no_extinction_at_all_reproduces_the_bare_templates(self, rng):
        S, G, obj, sig, _ = problem(rng)
        red = np.ones((1, S.shape[1]))

        _, _, chi2, _ = solve_extinction_batch(S, G, obj, sig, red)
        _, _, c1, _ = solve_grid(S[None, :, :], G[:, None, :], obj, sig)

        np.testing.assert_allclose(chi2[0], c1, rtol=1e-11, atol=1e-11)

    def test_a_template_with_no_overlap_is_rejected_the_same_way(self, rng):
        S, G, obj, sig, red = problem(rng)
        S[3] = np.nan

        _, _, chi2, times = solve_extinction_batch(S, G, obj, sig, red)

        assert (times[:, :, 3] == 0).all()
        assert (chi2[:, :, 3] == 0).all()


class TestTopK:
    def test_it_picks_the_same_candidates_a_full_sort_would(self, rng):
        scores = rng.random((4, 6, 30))

        index, values = _top_k_per_extinction(scores, 5, 0, 30)

        for a in range(4):
            expected = np.sort(scores[a].reshape(-1))[:5]
            np.testing.assert_allclose(values[a], expected)

    def test_indices_point_at_the_scores_they_came_from(self, rng):
        scores = rng.random((3, 4, 20))

        index, values = _top_k_per_extinction(scores, 4, 0, 20)

        for a in range(3):
            for rank in range(4):
                g, s = divmod(int(index[a, rank]), 20)
                assert scores[a, g, s] == pytest.approx(values[a, rank])

    def test_a_block_offset_maps_back_to_the_whole_bank(self, rng):
        scores = rng.random((2, 3, 10))

        index, _ = _top_k_per_extinction(scores, 3, 40, 100)

        sn = index % 100
        assert ((sn >= 40) & (sn < 50)).all()

    def test_ties_are_broken_by_position_not_by_luck(self):
        scores = np.zeros((1, 2, 3))

        index, values = _top_k_per_extinction(scores, 4, 0, 3)

        assert (values == 0).all()
        assert list(index[0]) == sorted(index[0])

    def test_asking_for_more_than_exist_is_not_an_error(self, rng):
        scores = rng.random((2, 1, 3))

        index, values = _top_k_per_extinction(scores, 99, 0, 3)

        assert values.shape == (2, 3)


class TestMergeBlocks:
    def _record(self, z, av, gal, sn, score):
        n = len(gal)
        return {
            "z": z, "extcon": av,
            "gal_index": np.asarray(gal), "sn_index": np.asarray(sn),
            "b": np.ones(n), "d": np.ones(n),
            "reduchi2": np.asarray(score, dtype=float),
            "reduchi2_once": np.asarray(score, dtype=float),
            "sn_mean": np.ones(n), "gal_mean": np.ones(n),
        }

    def test_the_best_across_blocks_wins(self):
        """Each block's top-k is only its own; the grid point's is the union."""

        records = [
            self._record(0.1, 0.0, [0, 0], [1, 2], [5.0, 9.0]),
            self._record(0.1, 0.0, [1, 1], [7, 8], [1.0, 3.0]),
        ]

        merged = merge_blocks(records, 3)

        assert len(merged) == 1
        np.testing.assert_allclose(merged[0]["reduchi2"], [1.0, 3.0, 5.0])
        assert list(merged[0]["sn_index"]) == [7, 8, 1]

    def test_each_grid_point_is_merged_on_its_own(self):
        records = [
            self._record(0.1, 0.0, [0], [1], [5.0]),
            self._record(0.1, 0.5, [0], [2], [1.0]),
            self._record(0.2, 0.0, [0], [3], [2.0]),
        ]

        merged = merge_blocks(records, 10)

        assert {(r["z"], r["extcon"]) for r in merged} == {
            (0.1, 0.0), (0.1, 0.5), (0.2, 0.0)
        }

    def test_it_keeps_no_more_than_asked_for(self):
        records = [
            self._record(0.1, 0.0, [0] * 5, list(range(5)), [5, 4, 3, 2, 1]),
            self._record(0.1, 0.0, [1] * 5, list(range(5, 10)), [9, 8, 7, 6, 0]),
        ]

        merged = merge_blocks(records, 4)

        assert len(merged[0]["reduchi2"]) == 4
        np.testing.assert_allclose(merged[0]["reduchi2"], [0, 1, 2, 3])

    def test_one_block_gets_the_same_guarantee_as_many(self):
        """The cap holds however the work was split, not just when it was."""

        records = [self._record(0.1, 0.0, [0] * 4, [0, 1, 2, 3], [4, 3, 2, 1])]

        merged = merge_blocks(records, 2)

        np.testing.assert_allclose(merged[0]["reduchi2"], [1.0, 2.0])

    def test_the_result_does_not_depend_on_how_the_bank_was_split(self):
        whole = [self._record(0.1, 0.0, [0] * 6, list(range(6)),
                              [3.0, 1.0, 5.0, 2.0, 6.0, 4.0])]
        split = [
            self._record(0.1, 0.0, [0] * 3, [0, 1, 2], [3.0, 1.0, 5.0]),
            self._record(0.1, 0.0, [0] * 3, [3, 4, 5], [2.0, 6.0, 4.0]),
        ]

        a = merge_blocks(whole, 3)[0]
        b = merge_blocks(split, 3)[0]

        np.testing.assert_allclose(a["reduchi2"], b["reduchi2"])
        assert list(a["sn_index"]) == list(b["sn_index"])


class TestAssembleResults:
    def _record(self):
        return {
            "z": 0.1, "extcon": 0.5,
            "gal_index": np.array([1, 0]), "sn_index": np.array([2, 0]),
            "b": np.array([2.0, 1.0]), "d": np.array([0.5, 0.25]),
            "reduchi2": np.array([1.0, 2.0]),
            "reduchi2_once": np.array([3.0, 4.0]),
            "sn_mean": np.array([1.0, 1.0]),
            "gal_mean": np.array([1.0, 1.0]),
        }

    def test_names_are_resolved_from_the_indices(self):
        result = assemble_results(
            [self._record()],
            ["Ia-norm/x/KAST phase-band : 1.0B", "b", "Ic/y/LRIS phase-band : -2.5V"],
            ["/bank/gal/E", "/bank/gal/Sa"],
            "spectrum.flm",
            10,
        )

        assert list(result["SN"].astype(str)) == [
            "Ic/y/LRIS phase-band : -2.5V", "Ia-norm/x/KAST phase-band : 1.0B"
        ]
        assert list(result["GALAXY"].astype(str)) == ["Sa", "E"]

    def test_phase_and_band_are_split_off_the_shorthand(self):
        result = assemble_results(
            [self._record()],
            ["a", "b", "Ic/y/LRIS phase-band : -2.5V"],
            ["/bank/gal/E", "/bank/gal/Sa"],
            "spectrum.flm",
            10,
        )

        assert result["Phase"].astype(str)[0] == " -2.5"
        assert result["Band"].astype(str)[0] == "V"

    def test_the_flux_fractions_sum_to_one(self):
        result = assemble_results(
            [self._record()], ["a", "b", "c"], ["/g/E", "/g/Sa"], "s.flm", 10
        )

        total = np.asarray(result["Frac(SN)"]) + np.asarray(result["Frac(gal)"])
        np.testing.assert_allclose(total, 1.0, rtol=1e-6)

    def test_the_columns_are_the_ones_the_csv_promises(self):
        result = assemble_results(
            [self._record()], ["a", "b", "c"], ["/g/E", "/g/Sa"], "s.flm", 10
        )

        assert result.colnames == [
            "SPECTRUM", "GALAXY", "SN", "CONST_SN", "CONST_GAL", "Z", "A_v",
            "Phase", "Band", "Frac(SN)", "Frac(gal)", "CHI2/dof", "CHI2/dof2",
        ]

    def test_no_candidates_at_all_is_an_error_with_a_reason(self):
        with pytest.raises(RuntimeError, match="no candidates"):
            assemble_results([], ["a"], ["/g/E"], "s.flm", 10)

    def test_degenerate_sentinel_scores_are_dropped(self):
        record = self._record()
        record["reduchi2"] = np.array([1.0, 1e10])

        result = assemble_results(
            [record], ["a", "b", "c"], ["/g/E", "/g/Sa"], "s.flm", 10
        )

        assert len(result) == 1
