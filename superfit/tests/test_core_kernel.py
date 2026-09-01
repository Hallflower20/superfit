"""Pin the fitting maths in ``core`` against an independent implementation.

These are the tests that make a refactor safe. ``core`` computes the whole
(galaxy x supernova) chi2 grid with broadcast nansums, which is fast but
opaque; the reference below does the same thing with plain loops so that any
change in the vectorised version -- caching, reordering, dtype changes --
shows up as a mismatch rather than as quietly different science.

The reference deliberately reproduces *current* behaviour, including the
choices that are statistically questionable (an unweighted linear solve
feeding a weighted chi2, and dividing by ``(dof - 2) ** 2``). If one of those
is changed on purpose, the reference has to be updated in the same commit --
which is exactly the signal we want.
"""

import numpy as np
import pytest

from superfit.SF_functions import Alam, GridSolver, redshifted_models, solve_grid
from superfit.loggrid import LogGrid, RedshiftableTemplates, velocity_to_dlnlam


def reference_solution(int_obj, sn, gal, sigma, minimum_overlap=0.7):
    """Loop-based restatement of the linear algebra inside ``core``.

    Parameters mirror the internals of ``core``: ``sn`` is (1, n_sn, n_lam),
    ``gal`` is (n_gal, 1, n_lam), ``int_obj`` and ``sigma`` are (n_lam,).

    Returns dict of (n_gal, n_sn) arrays.
    """

    n_gal = gal.shape[0]
    n_sn = sn.shape[1]
    n_lam = sn.shape[2]

    b = np.zeros((n_gal, n_sn))
    d = np.zeros((n_gal, n_sn))
    chi2 = np.zeros((n_gal, n_sn))
    times = np.zeros((n_gal, n_sn))

    for i in range(n_gal):
        for j in range(n_sn):
            s = sn[0, j]
            g = gal[i, 0]

            sum_ss = np.nansum(s * s)
            sum_gg = np.nansum(g * g)
            sum_gs = np.nansum(g * s)
            sum_so = np.nansum(s * int_obj)
            sum_go = np.nansum(g * int_obj)

            denom = sum_ss * sum_gg - sum_gs**2
            c = 1.0 / denom

            b_ij = c * (sum_gg * sum_so - sum_gs * sum_go)
            d_ij = c * (sum_ss * sum_go - sum_gs * sum_so)

            # Negative scalings are rejected outright, not clipped.
            if b_ij < 0:
                b_ij = np.nan
            if d_ij < 0:
                d_ij = np.nan

            resid = (int_obj - (b_ij * s + d_ij * g)) / sigma
            n_valid = n_lam - np.isnan(resid**2).sum()

            b[i, j] = b_ij
            d[i, j] = d_ij
            times[i, j] = n_valid
            chi2[i, j] = np.nansum(
                (int_obj - (b_ij * s + d_ij * g)) ** 2 / sigma**2
            )

    chi2 = np.where(times / n_lam > minimum_overlap, chi2, np.inf)

    reduchi2 = chi2 / (times - 2) ** 2
    reduchi2 = np.where(reduchi2 == 0, 1e10, reduchi2)

    reduchi2_once = chi2 / (times - 2)
    reduchi2_once = np.where(reduchi2_once == 0, 1e10, reduchi2_once)

    return {
        "b": b,
        "d": d,
        "chi2": chi2,
        "times": times,
        "reduchi2": reduchi2,
        "reduchi2_once": reduchi2_once,
    }


def vectorised_solution(int_obj, sn, gal, sigma, minimum_overlap=0.7):
    """The expressions as they appear in ``core`` today, lifted out verbatim."""

    n_lam = sn.shape[2]

    c = 1 / (
        np.nansum(sn**2, 2) * np.nansum(gal**2, 2) - np.nansum(gal * sn, 2) ** 2
    )
    b = c * (
        np.nansum(gal**2, 2) * np.nansum(sn * int_obj, 2)
        - np.nansum(gal * sn, 2) * np.nansum(gal * int_obj, 2)
    )
    d = c * (
        np.nansum(sn**2, 2) * np.nansum(gal * int_obj, 2)
        - np.nansum(gal * sn, 2) * np.nansum(sn * int_obj, 2)
    )

    b[b < 0] = np.nan
    d[d < 0] = np.nan

    sn_b = b[:, :, np.newaxis]
    gal_d = d[:, :, np.newaxis]

    a = ((int_obj - (sn_b * sn + gal_d * gal)) / sigma) ** 2
    times = n_lam - np.nansum(np.isnan(a), 2)

    chi2 = np.nansum(((int_obj - (sn_b * sn + gal_d * gal)) ** 2 / sigma**2), 2)
    chi2[~(times / n_lam > minimum_overlap)] = np.inf

    reduchi2 = np.where(chi2 / (times - 2) ** 2 == 0, 1e10, chi2 / (times - 2) ** 2)
    reduchi2_once = np.where(
        chi2 / (times - 2) == 0, 1e10, chi2 / (times - 2)
    )

    return {
        "b": b,
        "d": d,
        "chi2": chi2,
        "times": times,
        "reduchi2": reduchi2,
        "reduchi2_once": reduchi2_once,
    }


class TestLinearAlgebra:
    def test_matches_reference(self, toy_fit_problem):
        p = toy_fit_problem
        ref = reference_solution(p["int_obj"], p["sn"], p["gal"], p["sigma"])
        vec = vectorised_solution(p["int_obj"], p["sn"], p["gal"], p["sigma"])

        for key in ("b", "d", "chi2", "times", "reduchi2", "reduchi2_once"):
            np.testing.assert_allclose(
                vec[key], ref[key], rtol=1e-10, atol=0, err_msg=key
            )

    def test_recovers_the_injected_truth(self, toy_fit_problem):
        """The best-ranked cell should be the pair the observation was built from."""

        p = toy_fit_problem
        ref = reference_solution(p["int_obj"], p["sn"], p["gal"], p["sigma"])

        best = np.unravel_index(np.nanargmin(ref["reduchi2"]), ref["reduchi2"].shape)
        assert best == (p["truth_gal"], p["truth_sn"])

        np.testing.assert_allclose(ref["b"][best], p["b_true"], rtol=0.05)
        np.testing.assert_allclose(ref["d"][best], p["d_true"], rtol=0.05)

    def test_nan_padding_does_not_shift_the_answer(self, toy_fit_problem):
        """Templates that do not span the full grid are handled by nansum.

        Blanking the edges of every template must not change which cell wins,
        only how many pixels each fit is scored over.
        """

        p = toy_fit_problem
        sn = p["sn"].copy()
        gal = p["gal"].copy()
        sn[..., :5] = np.nan
        sn[..., -5:] = np.nan

        ref = reference_solution(p["int_obj"], sn, gal, p["sigma"])
        vec = vectorised_solution(p["int_obj"], sn, gal, p["sigma"])

        np.testing.assert_allclose(vec["chi2"], ref["chi2"], rtol=1e-10)
        np.testing.assert_array_equal(vec["times"], ref["times"])
        assert (ref["times"] <= p["lam"].size - 10).all()

    def test_minimum_overlap_rejects_thin_templates(self, toy_fit_problem):
        """A template covering too little of the grid must score inf, not a good chi2."""

        p = toy_fit_problem
        sn = p["sn"].copy()
        # Template 0 only covers 20% of the wavelength grid.
        sn[0, 0, 24:] = np.nan

        vec = vectorised_solution(
            p["int_obj"], sn, p["gal"], p["sigma"], minimum_overlap=0.7
        )
        assert np.isinf(vec["chi2"][:, 0]).all()
        assert np.isfinite(vec["chi2"][:, 1:]).all()


class TestSolveGrid:
    """``solve_grid`` replaced the broadcast residual cube with matrix products.

    The algebra is identical, the floating-point summation order is not, so
    these compare against the direct form rather than against stored numbers.
    Ragged template coverage and NaN gaps are the cases where the two could
    plausibly diverge, so they are the cases exercised.
    """

    @staticmethod
    def _ragged_problem(rng, n_lam=200, n_sn=12, n_gal=5, gap_obj=False, gap_sig=False):
        sn = np.abs(rng.normal(1.0, 0.4, (1, n_sn, n_lam)))
        gal = np.abs(rng.normal(1.0, 0.4, (n_gal, 1, n_lam)))
        obj = np.abs(rng.normal(1.0, 0.3, n_lam))
        sigma = np.abs(rng.normal(0.05, 0.01, n_lam)) + 1e-3

        # Real templates do not all span the observed grid.
        for s in range(n_sn):
            sn[0, s, : rng.integers(0, n_lam // 4)] = np.nan
        for g in range(n_gal):
            gal[g, 0, : rng.integers(0, n_lam // 5)] = np.nan

        if gap_obj:
            obj[rng.choice(n_lam, size=n_lam // 20, replace=False)] = np.nan
        if gap_sig:
            sigma[rng.choice(n_lam, size=5, replace=False)] = np.nan

        return sn, gal, obj, sigma

    @staticmethod
    def _direct(sn, gal, int_obj, sigma):
        """The pre-optimisation broadcast implementation."""

        c = 1 / (
            np.nansum(sn**2, 2) * np.nansum(gal**2, 2)
            - np.nansum(gal * sn, 2) ** 2
        )
        b = c * (
            np.nansum(gal**2, 2) * np.nansum(sn * int_obj, 2)
            - np.nansum(gal * sn, 2) * np.nansum(gal * int_obj, 2)
        )
        d = c * (
            np.nansum(sn**2, 2) * np.nansum(gal * int_obj, 2)
            - np.nansum(gal * sn, 2) * np.nansum(sn * int_obj, 2)
        )
        b, d = b.copy(), d.copy()
        b[b < 0] = np.nan
        d[d < 0] = np.nan

        sn_b, gal_d = b[:, :, np.newaxis], d[:, :, np.newaxis]
        resid = int_obj - (sn_b * sn + gal_d * gal)
        times = int_obj.size - np.nansum(np.isnan((resid / sigma) ** 2), 2)
        chi2 = np.nansum(resid**2 / sigma**2, 2)
        return b, d, chi2, times

    @pytest.mark.parametrize(
        "gap_obj,gap_sig",
        [(False, False), (True, False), (False, True), (True, True)],
        ids=["clean", "nan-in-object", "nan-in-sigma", "nan-in-both"],
    )
    def test_agrees_with_the_direct_form(self, rng, gap_obj, gap_sig):
        sn, gal, obj, sigma = self._ragged_problem(
            rng, gap_obj=gap_obj, gap_sig=gap_sig
        )

        b0, d0, chi0, t0 = self._direct(sn, gal, obj, sigma)
        b1, d1, chi1, t1 = solve_grid(sn, gal, obj, sigma)

        # Pixel counts are integers and must agree exactly.
        np.testing.assert_array_equal(t0, t1)

        # Rejected (negative-amplitude) cells must be rejected identically.
        np.testing.assert_array_equal(np.isnan(b0), np.isnan(b1))
        np.testing.assert_array_equal(np.isnan(d0), np.isnan(d1))

        good = np.isfinite(b0) & np.isfinite(d0)
        np.testing.assert_allclose(b1[good], b0[good], rtol=1e-10)
        np.testing.assert_allclose(d1[good], d0[good], rtol=1e-10)
        np.testing.assert_allclose(chi1[good], chi0[good], rtol=1e-9)

    def test_ranking_is_preserved(self, rng):
        """What actually matters: the same templates come out on top."""

        for _ in range(20):
            sn, gal, obj, sigma = self._ragged_problem(rng)
            _, _, chi0, _ = self._direct(sn, gal, obj, sigma)
            _, _, chi1, _ = solve_grid(sn, gal, obj, sigma)

            order0 = np.argsort(np.where(np.isfinite(chi0), chi0, np.inf), axis=None)
            order1 = np.argsort(np.where(np.isfinite(chi1), chi1, np.inf), axis=None)
            np.testing.assert_array_equal(order0[:10], order1[:10])

    def test_chi2_is_never_negative(self, rng):
        """The expanded form is a difference of large terms; guard the floor."""

        sn, gal, obj, sigma = self._ragged_problem(rng)
        _, _, chi2, _ = solve_grid(sn, gal, obj, sigma)
        assert (chi2[np.isfinite(chi2)] >= 0).all()

    def test_perfect_fit_scores_zero(self):
        """An observation built exactly from one pair must score ~0 there.

        "Zero" has to be measured against the size of the terms it cancels
        between, not against a fixed epsilon. solve_grid evaluates chi2 as a
        difference of six weighted sums, each of order sum(w * obj**2) -- about
        1.5e6 here, since sigma is 0.01 over 150 pixels -- so float64 leaves a
        residue around 1e-10 and the exact value depends on the order BLAS
        accumulated in. An absolute `< 1e-12` held on Linux and failed on
        Windows and macOS at ~1.3e-9, which is a difference between BLAS
        builds rather than a difference in the answer.
        """

        n_lam = 150
        rng = np.random.default_rng(3)
        sn = np.abs(rng.normal(1.0, 0.3, (1, 4, n_lam)))
        gal = np.abs(rng.normal(1.0, 0.3, (3, 1, n_lam)))
        obj = 0.8 * sn[0, 2] + 0.5 * gal[1, 0]
        sigma = np.full(n_lam, 0.01)

        b, d, chi2, _ = solve_grid(sn, gal, obj, sigma)

        # The scale the cancellation happens at: chi2's leading term.
        scale = float(np.sum(obj**2 / sigma**2))
        assert chi2[1, 2] < 1e-12 * scale
        np.testing.assert_allclose(b[1, 2], 0.8, rtol=1e-9)
        np.testing.assert_allclose(d[1, 2], 0.5, rtol=1e-9)

    def test_weighted_solve_is_off_by_default(self, rng):
        sn, gal, obj, sigma = self._ragged_problem(rng)

        default = solve_grid(sn, gal, obj, sigma)
        explicit = solve_grid(sn, gal, obj, sigma, weighted=False)

        for a, b in zip(default, explicit):
            np.testing.assert_array_equal(a, b)

    def test_weighted_solve_never_scores_worse(self, rng):
        """It minimises the chi2 that is reported, so it cannot lose to a
        solution that minimised something else."""

        for _ in range(10):
            sn, gal, obj, sigma = self._ragged_problem(rng)
            # A strongly varying sigma is what makes the two differ at all.
            sigma = sigma * (1 + 4 * np.abs(np.sin(np.arange(sigma.size) / 20.0)))

            _, _, chi_unweighted, _ = solve_grid(sn, gal, obj, sigma)
            _, _, chi_weighted, _ = solve_grid(sn, gal, obj, sigma, weighted=True)

            good = np.isfinite(chi_unweighted) & np.isfinite(chi_weighted)
            assert (chi_weighted[good] <= chi_unweighted[good] * (1 + 1e-9)).all()

    def test_weighted_solve_matches_a_direct_least_squares(self, rng):
        """Check one cell against an explicit weighted lstsq."""

        sn, gal, obj, sigma = self._ragged_problem(rng, n_sn=4, n_gal=3)
        b, d, _, _ = solve_grid(sn, gal, obj, sigma, weighted=True)

        for g in range(3):
            for s in range(4):
                if not (np.isfinite(b[g, s]) and np.isfinite(d[g, s])):
                    continue
                ok = (
                    np.isfinite(sn[0, s])
                    & np.isfinite(gal[g, 0])
                    & np.isfinite(obj)
                    & np.isfinite(sigma)
                )
                design = np.column_stack([sn[0, s][ok], gal[g, 0][ok]])
                w = 1.0 / sigma[ok]
                coeffs, *_ = np.linalg.lstsq(
                    design * w[:, None], obj[ok] * w, rcond=None
                )
                if (coeffs < 0).any():
                    continue  # solve_grid rejects negatives; lstsq does not
                np.testing.assert_allclose(
                    [b[g, s], d[g, s]], coeffs, rtol=1e-8
                )

    def test_does_not_modify_its_inputs(self, rng):
        sn, gal, obj, sigma = self._ragged_problem(rng)
        before = [a.copy() for a in (sn, gal, obj, sigma)]

        solve_grid(sn, gal, obj, sigma)

        for got, want in zip((sn, gal, obj, sigma), before):
            np.testing.assert_array_equal(got, want)




class TestRedshiftedModels:
    """Assembling the model grid at one (redshift, extinction) point.

    These are the properties the old linear-grid ``sn_hg_arrays`` guaranteed,
    restated against the log-grid path that replaced it.
    """

    @staticmethod
    def _banks(observed_grid, max_z=0.5):
        lam_t = np.linspace(2500.0, 12000.0, 3000)
        sn_flux = 1.0 + 0.2 * np.sin(lam_t / 300.0)
        gal_flux = 2.0 + 0.1 * np.cos(lam_t / 400.0)

        sn_bank = RedshiftableTemplates.from_templates(
            [lam_t], [sn_flux], observed_grid, max_z
        )
        gal_bank = RedshiftableTemplates.from_templates(
            [lam_t], [gal_flux], observed_grid, max_z
        )
        return sn_bank, gal_bank, lam_t, sn_flux, gal_flux

    @staticmethod
    def _grid():
        return LogGrid.spanning(3500.0, 9000.0, velocity_to_dlnlam(400.0))

    def test_shapes_broadcast_into_a_grid(self):
        grid = self._grid()
        sn_bank, gal_bank, *_ = self._banks(grid)

        sn, gal = redshifted_models(sn_bank, gal_bank, grid.wavelength, 0.0, 0.0)

        n = len(grid)
        assert sn.shape == (1, 1, n)
        assert gal.shape == (1, 1, n)
        assert (sn * gal).shape == (1, 1, n)

    def test_zero_redshift_zero_extinction_is_plain_resampling(self):
        grid = self._grid()
        sn_bank, gal_bank, lam_t, sn_flux, _ = self._banks(grid)

        sn, _ = redshifted_models(sn_bank, gal_bank, grid.wavelength, 0.0, 0.0)

        expected = np.interp(
            grid.wavelength, lam_t, sn_flux, left=np.nan, right=np.nan
        )
        np.testing.assert_allclose(sn[0, 0], expected, rtol=1e-3, atol=1e-3)

    def test_redshift_stretches_and_dims(self):
        """Flux is divided by (1 + z) and features move to longer wavelength."""

        grid = self._grid()
        sn_bank, gal_bank, lam_t, sn_flux, _ = self._banks(grid)
        z = 0.1

        sn, _ = redshifted_models(sn_bank, gal_bank, grid.wavelength, z, 0.0)

        expected = (
            np.interp(
                grid.wavelength / (1 + z), lam_t, sn_flux, left=np.nan, right=np.nan
            )
            / (1 + z)
        )
        np.testing.assert_allclose(sn[0, 0], expected, rtol=5e-3, atol=5e-3)

    def test_extinction_reddens_the_supernova(self):
        grid = self._grid()
        sn_bank, gal_bank, *_ = self._banks(grid)

        clear, _ = redshifted_models(sn_bank, gal_bank, grid.wavelength, 0.0, 0.0)
        reddened, _ = redshifted_models(sn_bank, gal_bank, grid.wavelength, 0.0, 1.0)

        ratio = reddened[0, 0] / clear[0, 0]
        expected = 10 ** (-0.4 * Alam(grid.wavelength))
        np.testing.assert_allclose(ratio, expected, rtol=1e-12)

        assert (ratio > 0).all()
        assert (ratio <= 1).all()

        # And it points the right way: positive A_v takes away more blue
        # light than red. See TestExtinctionLaw below.
        assert ratio[0] < ratio[-1]

    def test_r_v_changes_the_shape_of_the_reddening(self):
        """R_v is threaded through from the config, not pinned at 3.1."""

        grid = self._grid()
        sn_bank, gal_bank, *_ = self._banks(grid)

        default, _ = redshifted_models(sn_bank, gal_bank, grid.wavelength, 0.0, 1.0)
        steeper, _ = redshifted_models(
            sn_bank, gal_bank, grid.wavelength, 0.0, 1.0, R_v=2.1
        )

        np.testing.assert_allclose(
            steeper[0, 0] / default[0, 0],
            10 ** (-0.4 * (Alam(grid.wavelength, R_v=2.1) - Alam(grid.wavelength))),
            rtol=1e-12,
        )
        assert not np.allclose(steeper[0, 0], default[0, 0])

    def test_galaxy_templates_ignore_extinction(self):
        grid = self._grid()
        sn_bank, gal_bank, *_ = self._banks(grid)

        _, gal_a = redshifted_models(sn_bank, gal_bank, grid.wavelength, 0.05, 0.0)
        _, gal_b = redshifted_models(sn_bank, gal_bank, grid.wavelength, 0.05, 2.0)

        np.testing.assert_array_equal(gal_a, gal_b)

    def test_galaxy_is_still_dimmed_by_redshift(self):
        grid = self._grid()
        sn_bank, gal_bank, *_ = self._banks(grid)

        _, at_rest = redshifted_models(sn_bank, gal_bank, grid.wavelength, 0.0, 0.0)
        _, at_z = redshifted_models(sn_bank, gal_bank, grid.wavelength, 0.5, 0.0)

        finite = np.isfinite(at_rest) & np.isfinite(at_z)
        assert at_z[finite].max() < at_rest[finite].max()

    @pytest.mark.parametrize("z", [0.0, 0.05, 0.2])
    def test_templates_off_the_grid_become_nan_not_zero(self, z):
        """Extrapolation must be NaN so nansum drops it; 0 would score well."""

        grid = LogGrid.spanning(4000.0, 8000.0, velocity_to_dlnlam(400.0))
        lam_t = np.linspace(5000.0, 6000.0, 500)
        bank = RedshiftableTemplates.from_templates(
            [lam_t], [np.ones_like(lam_t)], grid, 0.5
        )

        sn, gal = redshifted_models(bank, bank, grid.wavelength, z, 0.0)

        assert np.isnan(sn[0, 0, 0])
        assert np.isnan(sn[0, 0, -1])
        assert np.isnan(gal[0, 0, 0])
        assert np.isfinite(sn[0, 0]).any()


class TestExtinctionLaw:
    """``Alam`` returns A_lambda in magnitudes, and every caller relies on it.

    Callers form the transmission themselves, as
    ``10 ** (-0.4 * A_v * Alam(lam))``. That is only correct if ``Alam``
    hands back magnitudes. It used to return
    ``extinction.apply(ccm89(...), ones)``, the *transmission*
    ``10 ** (-0.4 * A_lambda)``, so the callers exponentiated a second time
    and the resulting curve removed more red light than blue -- reddening ran
    backwards, and negative A_v was doing the work extinction should.

    ``test_doubling_av_squares_the_transmission`` is the sharpest test here:
    the law is linear in A_v in magnitudes, and the double exponentiation
    broke exactly that. Anything that reintroduces an extra ``10 ** ...``
    fails it.
    """

    WAVELENGTHS = np.array([3500.0, 5500.0, 9000.0])

    @pytest.mark.parametrize("A_v", [0.5, 1.0, 2.5])
    @pytest.mark.parametrize("R_v", [2.1, 3.1, 4.5])
    def test_alam_is_exactly_ccm89_in_magnitudes(self, A_v, R_v):
        import extinction

        np.testing.assert_allclose(
            Alam(self.WAVELENGTHS, A_v, R_v),
            extinction.ccm89(self.WAVELENGTHS, A_v, R_v),
            rtol=1e-12,
        )

    def test_positive_extinction_removes_more_blue_than_red(self):
        transmission = 10 ** (-0.4 * 1.0 * Alam(self.WAVELENGTHS))
        assert transmission[0] < transmission[-1]

    def test_reddening_is_monotonic_across_the_optical(self):
        lam = np.linspace(3500.0, 9000.0, 2000)
        transmission = 10 ** (-0.4 * 1.0 * Alam(lam))

        # More extinction at every step towards the blue, no local reversals.
        assert (np.diff(transmission) > 0).all()

    def test_zero_extinction_is_exactly_the_identity(self):
        transmission = 10 ** (-0.4 * 0.0 * Alam(self.WAVELENGTHS))
        np.testing.assert_array_equal(transmission, np.ones(len(self.WAVELENGTHS)))

    def test_doubling_av_squares_the_transmission(self):
        """The law is linear in A_v; double exponentiation broke this."""

        once = 10 ** (-0.4 * 1.0 * Alam(self.WAVELENGTHS))
        twice = 10 ** (-0.4 * 2.0 * Alam(self.WAVELENGTHS))
        np.testing.assert_allclose(twice, once**2, rtol=1e-12)

    def test_av_scales_out_of_the_curve(self):
        """``Alam(lam, A_v)`` is ``A_v * Alam(lam)``, which is why the
        callers may pass A_v as a multiplier instead of an argument."""

        np.testing.assert_allclose(
            Alam(self.WAVELENGTHS, 1.7), 1.7 * Alam(self.WAVELENGTHS), rtol=1e-12
        )

    def test_accepts_a_python_list(self):
        """Callers pass plain sequences; ccm89 itself demands float64."""

        np.testing.assert_allclose(
            Alam([3500.0, 9000.0]), Alam(np.array([3500.0, 9000.0])), rtol=1e-12
        )


class TestPlottingReconstructsTheFittedModel:
    """The output plots rebuild the best-fit model independently of the fit.

    ``sf_class._plot_best_fits`` and ``sf_class.any_result`` redden the
    template themselves rather than reusing what the fitter produced. If the
    two disagree, the published plot shows a different model than the one that
    was scored. This pins them together.
    """

    @staticmethod
    def _template():
        lam_t = np.linspace(2500.0, 12000.0, 3000)
        return lam_t, 1.0 + 0.2 * np.sin(lam_t / 300.0)

    @pytest.mark.parametrize("z", [0.0, 0.127, 0.3])
    @pytest.mark.parametrize("A_v", [-1.0, 0.0, 1.2])
    def test_the_two_reddenings_agree(self, z, A_v):
        from scipy import interpolate

        grid = LogGrid.spanning(3500.0, 9000.0, velocity_to_dlnlam(400.0))
        lam_t, flux_t = self._template()

        # What the fitter builds: redden at lam / (1 + z), on the log grid.
        bank = RedshiftableTemplates.from_templates([lam_t], [flux_t], grid, 0.5)
        sn, _ = redshifted_models(bank, bank, grid.wavelength, z, A_v)
        fitted = sn[0, 0]

        # What the plotting path builds: redden at the template's own
        # wavelength, then shift and interpolate onto the same grid.
        extinct = flux_t * 10 ** (-0.4 * A_v * Alam(lam_t)) / (1 + z)
        plotted = interpolate.interp1d(
            lam_t * (1 + z), extinct, bounds_error=False, fill_value=np.nan
        )(grid.wavelength)

        finite = np.isfinite(fitted) & np.isfinite(plotted)
        assert finite.sum() > 100
        np.testing.assert_allclose(
            fitted[finite], plotted[finite], rtol=5e-3, atol=5e-3
        )


class TestGridSolver:
    """GridSolver hoists the A_v-independent half of solve_grid.

    The contract is stronger than "close": the per-A_v terms are evaluated by
    the same operations in the same order, so the results must be
    bit-identical to reddening the bank first and calling solve_grid on it.
    """

    @pytest.mark.parametrize("weighted", [False, True], ids=["legacy", "weighted"])
    @pytest.mark.parametrize("extcon", [-2.0, 0.0, 0.6, 2.0])
    def test_bit_identical_to_reddening_the_bank_first(self, rng, weighted, extcon):
        sn, gal, obj, sigma = TestSolveGrid._ragged_problem(
            rng, gap_obj=True, gap_sig=True
        )
        lam = np.linspace(3500.0, 9000.0, obj.size)
        reddening = 10 ** (-0.4 * extcon * Alam(lam))

        solver = GridSolver(sn, gal, obj, sigma, weighted=weighted)
        got = solver.solve(reddening)
        want = solve_grid(sn * reddening, gal, obj, sigma, weighted=weighted)

        for g, w, name in zip(got, want, ("b", "d", "chi2", "times")):
            np.testing.assert_array_equal(g, w, err_msg=name)

    def test_solve_with_no_reddening_is_solve_grid(self, rng):
        sn, gal, obj, sigma = TestSolveGrid._ragged_problem(rng)

        got = GridSolver(sn, gal, obj, sigma).solve()
        want = solve_grid(sn, gal, obj, sigma)

        for g, w in zip(got, want):
            np.testing.assert_array_equal(g, w)

    def test_one_solver_serves_a_whole_extinction_grid(self, rng):
        """Reusing the solver across A_v values must not leak state."""

        sn, gal, obj, sigma = TestSolveGrid._ragged_problem(rng)
        lam = np.linspace(3500.0, 9000.0, obj.size)
        solver = GridSolver(sn, gal, obj, sigma)

        grid = [-1.0, 0.0, 1.0]
        first_pass = [solver.solve(10 ** (-0.4 * a * Alam(lam))) for a in grid]
        # Same values again, out of order, after the solver has been used.
        again = solver.solve(10 ** (-0.4 * grid[0] * Alam(lam)))

        for g, w in zip(again, first_pass[0]):
            np.testing.assert_array_equal(g, w)


class TestSplitPhaseBand:
    """The shorthand tail ": {phase}{band}" split into its two columns."""

    @pytest.mark.parametrize(
        "shorthand,phase,band",
        [
            ("Ic/1994I/KAST phase-band : -2.57B", "-2.57", "B"),
            ("II/1999em/X phase-band : 12.0-", "12.0", "-"),
            ("Ia/2011fe/X phase-band : uB", "u", "B"),
            # The cases the old last-character split got wrong: an object the
            # phase table does not cover has no band, and the split handed the
            # Band column the last digit of the phase (or the "u").
            ("Ia/2099xx/X phase-band : 12.5", "12.5", ""),
            ("Ia/2099xx/X phase-band : u", "u", ""),
        ],
    )
    def test_split(self, shorthand, phase, band):
        from superfit.SF_functions import split_phase_band

        assert split_phase_band(shorthand) == (phase, band)
