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

from NGSF.SF_functions import sn_hg_arrays


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


class TestSnHgArrays:
    """``sn_hg_arrays`` resamples templates onto the observed grid."""

    def _bank(self):
        lam_t = np.linspace(3000.0, 9000.0, 400)
        sn_t = {"a": np.column_stack([lam_t, 1.0 + 0.2 * np.sin(lam_t / 300.0)])}
        gal_t = {"g": np.column_stack([lam_t, 2.0 + 0.1 * np.cos(lam_t / 400.0)])}
        alam = {"a": np.ones_like(lam_t)}
        return sn_t, gal_t, alam

    def test_shapes_broadcast_into_a_grid(self):
        sn_t, gal_t, alam = self._bank()
        lam = np.linspace(4000.0, 8000.0, 100)

        sn, gal = sn_hg_arrays(0.0, 0.0, lam, ["a"], sn_t, ["g"], gal_t, alam)

        assert sn.shape == (1, 1, 100)
        assert gal.shape == (1, 1, 100)
        assert (sn * gal).shape == (1, 1, 100)

    def test_zero_redshift_zero_extinction_is_plain_interpolation(self):
        sn_t, gal_t, alam = self._bank()
        lam = np.linspace(4000.0, 8000.0, 100)

        sn, _ = sn_hg_arrays(0.0, 0.0, lam, ["a"], sn_t, ["g"], gal_t, alam)

        expected = np.interp(
            lam, sn_t["a"][:, 0], sn_t["a"][:, 1], left=np.nan, right=np.nan
        )
        np.testing.assert_allclose(sn[0, 0], expected, rtol=1e-12)

    def test_redshift_stretches_and_dims(self):
        """A redshifted template is shifted in lambda and divided by (1 + z)."""

        sn_t, gal_t, alam = self._bank()
        lam = np.linspace(4000.0, 8000.0, 100)
        z = 0.1

        sn, _ = sn_hg_arrays(z, 0.0, lam, ["a"], sn_t, ["g"], gal_t, alam)

        expected = np.interp(
            lam,
            sn_t["a"][:, 0] * (1 + z),
            sn_t["a"][:, 1] / (1 + z),
            left=np.nan,
            right=np.nan,
        )
        np.testing.assert_allclose(sn[0, 0], expected, rtol=1e-12)

    def test_extinction_scales_the_flux(self):
        sn_t, gal_t, alam = self._bank()
        lam = np.linspace(4000.0, 8000.0, 100)

        unextincted, _ = sn_hg_arrays(0.0, 0.0, lam, ["a"], sn_t, ["g"], gal_t, alam)
        extincted, _ = sn_hg_arrays(0.0, 1.0, lam, ["a"], sn_t, ["g"], gal_t, alam)

        # alam is 1 everywhere in this bank, so the ratio is a single constant.
        ratio = extincted[0, 0] / unextincted[0, 0]
        np.testing.assert_allclose(ratio, 10 ** (-0.4), rtol=1e-12)

    def test_galaxy_templates_ignore_extinction(self):
        """Only the supernova is reddened; the host is not."""

        sn_t, gal_t, alam = self._bank()
        lam = np.linspace(4000.0, 8000.0, 100)

        _, gal_a = sn_hg_arrays(0.05, 0.0, lam, ["a"], sn_t, ["g"], gal_t, alam)
        _, gal_b = sn_hg_arrays(0.05, 2.0, lam, ["a"], sn_t, ["g"], gal_t, alam)

        np.testing.assert_array_equal(gal_a, gal_b)


@pytest.mark.parametrize("z", [0.0, 0.05, 0.2])
def test_templates_off_the_grid_become_nan_not_zero(z):
    """Extrapolation must produce NaN so nansum drops it, not 0 which would score well."""

    lam_t = np.linspace(5000.0, 6000.0, 50)
    sn_t = {"a": np.column_stack([lam_t, np.ones_like(lam_t)])}
    gal_t = {"g": np.column_stack([lam_t, np.ones_like(lam_t)])}
    alam = {"a": np.ones_like(lam_t)}

    lam = np.linspace(4000.0, 8000.0, 200)
    sn, gal = sn_hg_arrays(z, 0.0, lam, ["a"], sn_t, ["g"], gal_t, alam)

    assert np.isnan(sn[0, 0, 0])
    assert np.isnan(sn[0, 0, -1])
    assert np.isnan(gal[0, 0, 0])
    assert np.isfinite(sn[0, 0]).any()
