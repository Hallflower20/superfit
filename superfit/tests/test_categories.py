"""Stars and QSOs: standalone categories, only where the bank supplies them.

The modern banks carry stellar templates (``sne/star-*`` type directories)
and QSO templates (``gal/*QSO*`` files). Each is a category of its own:

* a star or QSO is fit ALONE -- no host galaxy underneath;
* a QSO is never offered as a host, so no transient is ever fit on top of one;
* stars are foreground, so their redshift is pinned to zero;
* the legacy bank has neither, so a legacy fit is byte-for-byte unchanged
  (the golden regression in test_regression.py is the proof of that).

Everything here builds small banks in a tmp_path; none of it needs a real one.
"""

import os

import numpy as np
import pytest

from superfit.config import ConfigError, as_tristate, load_config
from superfit.SF_functions import SingleTemplateSolver
from superfit.params import Parameters, is_qso_name, is_star_type


# --- the synthetic bank ----------------------------------------------------

LAM = np.arange(3000.0, 9000.0, 5.0)

WISEREP_HEADER = (
    "Obj. ID,IAU name,Obj. Type,Redshift,JD,Obs-date,Instrument,"
    "Ascii file,WL Medium,Obj. RA,Obj. DEC\n"
)


def star_flux(lam):
    return 2.0 + np.sin(lam / 300.0)


def sn_flux(lam):
    return 2.0 + np.cos(lam / 150.0)


def gal_flux(lam):
    return 1.0 + (lam - lam[0]) / (lam[-1] - lam[0])


def qso_flux(lam):
    return 2.0 + np.exp(-0.5 * ((lam - 6000.0) / 200.0) ** 2)


def write_template(path, flux):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.column_stack([LAM, flux]))


def add_object(root, kind, sn_type, name, flux):
    """One object directory with its template and wiserep metadata."""

    filename = name + ".ascii"
    for res in ("original_resolution", "binnings/10A"):
        write_template(root / res / "sne" / sn_type / name / filename, flux)

    for res in ("original_resolution", "binnings/10A"):
        csv = root / res / "sne" / sn_type / name / "wiserep_spectra.csv"
        csv.write_text(
            WISEREP_HEADER
            + "{0},,{1},0.0,,,TEST,{2},vacuum,,\n".format(name, kind, filename)
        )


def make_category_bank(root, with_stars=True, with_qsos=True):
    """A bank with one SN type, one galaxy, and optionally stars and QSOs."""

    add_object(root, "SN Ic", "Ic", "SN2099zz", sn_flux(LAM))
    if with_stars:
        add_object(root, "STAR M", "star-M", "DESI_STAR_7", star_flux(LAM))

    for res in ("original_resolution", "binnings/10A"):
        write_template(root / res / "gal" / "E", gal_flux(LAM))
        if with_qsos:
            write_template(root / res / "gal" / "DESI_QSO_7", qso_flux(LAM))

    (root / "mjd_of_maximum_brightness.csv").write_text(
        "Name,mjd_peak,band_peak,isupperlimit\nSN2099zz,55814.0,B,0\n"
    )
    return root


@pytest.fixture
def bank(tmp_path, monkeypatch):
    """A category-carrying bank, installed as the bank a fit would find."""

    from superfit import paths

    root = make_category_bank(tmp_path / "bank")
    monkeypatch.setenv("SUPERFIT_BANK_DIR", str(root))
    monkeypatch.setattr(paths, "_cached_bank_dir", None)
    return root


@pytest.fixture
def legacy_shaped_bank(tmp_path, monkeypatch):
    """The same bank with no stars and no QSOs, like the legacy one."""

    from superfit import paths

    root = make_category_bank(tmp_path / "bank", with_stars=False, with_qsos=False)
    monkeypatch.setenv("SUPERFIT_BANK_DIR", str(root))
    monkeypatch.setattr(paths, "_cached_bank_dir", None)
    return root


def observed_spectrum(tmp_path, flux_of, scale=3.0, seed=0):
    """A spectrum file that is one of the templates, rescaled, plus noise."""

    lam = np.arange(4000.0, 8000.0, 2.0)
    rng = np.random.default_rng(seed)
    flux = scale * np.interp(lam, LAM, flux_of(LAM))
    flux = flux * (1.0 + 0.01 * rng.normal(size=lam.size))
    path = tmp_path / "observed.flm"
    np.savetxt(path, np.column_stack([lam, flux]))
    return str(path)


BASE = {
    "z": 0.0,
    "temp_sn_tr": ["Ic"],
    "temp_gal_tr": ["E"],
    "n_cores": 1,
    "n_plots": 0,
    "overwrite": True,
}


def run_fit(tmp_path, flux_of, **overrides):
    from superfit import Superfit

    settings = dict(BASE, output_dir=str(tmp_path / "out"), **overrides)
    return Superfit(observed_spectrum(tmp_path, flux_of), **settings).run().results


# --- the solver ------------------------------------------------------------

class TestSingleTemplateSolver:
    """One amplitude, same conventions as the two-component GridSolver."""

    @staticmethod
    def _problem(rng, n_t=8, n_lam=300):
        T = np.abs(rng.normal(1.0, 0.4, (n_t, n_lam)))
        obj = np.abs(rng.normal(1.0, 0.3, n_lam))
        sigma = np.abs(rng.normal(0.05, 0.01, n_lam)) + 1e-3
        for i in range(n_t):
            T[i, : rng.integers(0, n_lam // 4)] = np.nan
        obj[rng.choice(n_lam, size=n_lam // 20, replace=False)] = np.nan
        sigma[rng.choice(n_lam, size=5, replace=False)] = np.nan
        return T, obj, sigma

    @staticmethod
    def _direct(T, obj, sigma, weighted):
        """Loop-based restatement, one template at a time."""

        n_t = T.shape[0]
        b = np.full(n_t, np.nan)
        chi2 = np.zeros(n_t)
        times = np.zeros(n_t)

        for i in range(n_t):
            t = T[i]
            valid = np.isfinite(t) & np.isfinite(obj) & np.isfinite(sigma)
            w = 1.0 / sigma[valid] ** 2

            if weighted:
                amp = (w * t[valid] * obj[valid]).sum() / (w * t[valid] ** 2).sum()
            else:
                # Legacy: unweighted normals; sum(t^2) over the template's
                # own coverage, sum(t*obj) over its overlap with the object.
                own = np.isfinite(t)
                amp = (
                    np.nansum(t[own] * np.where(np.isfinite(obj[own]), obj[own], 0.0))
                    / (t[own] ** 2).sum()
                )

            if amp < 0:
                continue

            b[i] = amp
            chi2[i] = (w * (obj[valid] - amp * t[valid]) ** 2).sum()
            times[i] = valid.sum()

        return b, chi2, times

    @pytest.mark.parametrize("weighted", [False, True], ids=["legacy", "weighted"])
    def test_matches_the_direct_form(self, rng, weighted):
        T, obj, sigma = self._problem(rng)

        solver = SingleTemplateSolver(T, obj, sigma, weighted=weighted)
        b, chi2, times = solver.solve()
        b0, chi0, t0 = self._direct(T, obj, sigma, weighted)

        np.testing.assert_allclose(b, b0, rtol=1e-10, equal_nan=True)
        keep = np.isfinite(b)
        np.testing.assert_allclose(chi2[keep], chi0[keep], rtol=1e-8)
        np.testing.assert_array_equal(times, t0)

    def test_reddening_is_the_same_as_dimming_the_templates(self, rng):
        T, obj, sigma = self._problem(rng)
        reddening = np.linspace(0.5, 1.5, T.shape[1])

        got = SingleTemplateSolver(T, obj, sigma).solve(reddening)
        want = SingleTemplateSolver(T * reddening, obj, sigma).solve()

        for g, w in zip(got, want):
            np.testing.assert_allclose(g, w, rtol=1e-12, equal_nan=True)

    def test_negative_amplitude_is_rejected(self, rng):
        T = np.ones((1, 100))
        obj = -np.ones(100)
        sigma = np.full(100, 0.1)

        b, chi2, times = SingleTemplateSolver(T, obj, sigma).solve()

        assert np.isnan(b[0])
        assert chi2[0] == 0.0 and times[0] == 0.0


# --- category discovery and gating -----------------------------------------

# A wavelength grid stated outright, so Parameters can be built without an
# observation in hand.
GRID = dict(load_config(), lower_lam=3800, upper_lam=8300)


class TestDiscovery:
    def test_name_predicates(self):
        assert is_qso_name("/bank/gal/DESI_QSO_39627")
        assert not is_qso_name("/bank/gal/DESI_39627")
        assert not is_qso_name("/bank/gal/E")
        assert is_star_type("star-M") and is_star_type("star-WD")
        assert not is_star_type("Ia-norm")

    def test_bank_with_categories_resolves_auto_on(self, bank):
        p = Parameters(dict(GRID, temp_sn_tr=["Ic"], temp_gal_tr=["E"]))

        assert p.star_types == ["star-M"]
        assert [os.path.basename(x) for x in p.templates_qso] == ["DESI_QSO_7"]
        assert p.fit_stars and p.fit_qsos

    def test_legacy_shaped_bank_resolves_auto_off(self, legacy_shaped_bank):
        p = Parameters(dict(GRID, temp_sn_tr=["Ic"], temp_gal_tr=["E"]))

        assert p.star_types == [] and p.templates_qso == []
        assert not p.fit_stars and not p.fit_qsos

    def test_requiring_categories_of_a_bank_without_them_is_an_error(
        self, legacy_shaped_bank
    ):
        with pytest.raises(ConfigError, match="no star templates"):
            Parameters(
                dict(GRID, temp_sn_tr=["Ic"], temp_gal_tr=["E"],
                     fit_stars=True)
            )
        with pytest.raises(ConfigError, match="no QSO templates"):
            Parameters(
                dict(GRID, temp_sn_tr=["Ic"], temp_gal_tr=["E"],
                     fit_qsos=True)
            )

    def test_a_qso_can_never_be_selected_as_a_host(self, bank):
        """Even naming the QSO in temp_gal_tr does not make it a host."""

        with pytest.raises(ConfigError, match="never hosts"):
            Parameters(
                dict(GRID, temp_sn_tr=["Ic"],
                     temp_gal_tr=["DESI_QSO_7"])
            )

    def test_tristate_parsing(self):
        assert as_tristate("auto") == "auto"
        assert as_tristate(True) is True and as_tristate(0) is False
        assert as_tristate("yes") is True
        with pytest.raises(ConfigError):
            as_tristate("maybe", "fit_stars")


# --- the fit itself ---------------------------------------------------------

class TestCategoryFit:
    def test_a_star_wins_for_a_stellar_spectrum(self, bank, tmp_path):
        results = run_fit(tmp_path, star_flux)

        best = results.iloc[0]
        assert str(best["SN"]).startswith("star-M/DESI_STAR_7")
        assert str(best["GALAXY"]) == "none"
        assert best["Z"] == 0.0
        assert best["Frac(SN)"] == 1.0 and best["Frac(gal)"] == 0.0
        assert best["CONST_GAL"] == 0.0

    def test_a_qso_wins_for_a_qso_spectrum(self, bank, tmp_path):
        results = run_fit(tmp_path, qso_flux)

        best = results.iloc[0]
        assert str(best["SN"]) == "QSO/DESI_QSO_7"
        assert str(best["GALAXY"]) == "none"

    def test_no_transient_is_ever_fit_on_a_qso(self, bank, tmp_path):
        results = run_fit(tmp_path, qso_flux)

        hosts = results["GALAXY"].astype(str)
        assert not hosts.str.contains("QSO").any()

    def test_supernovae_still_win_for_a_supernova_spectrum(self, bank, tmp_path):
        results = run_fit(tmp_path, sn_flux)

        assert str(results.iloc[0]["SN"]).startswith("Ic/SN2099zz")
        # The categories were still tried and ranked below.
        assert (results["GALAXY"].astype(str) == "none").any()

    def test_categories_can_be_switched_off(self, bank, tmp_path):
        results = run_fit(tmp_path, star_flux, stars=False, qsos=False)

        names = results["SN"].astype(str)
        assert not names.str.startswith("star").any()
        assert not names.str.startswith("QSO/").any()

    def test_a_legacy_shaped_bank_fits_exactly_as_before(
        self, legacy_shaped_bank, tmp_path
    ):
        results = run_fit(tmp_path, sn_flux)

        assert len(results) > 0
        assert not (results["GALAXY"].astype(str) == "none").any()

    def test_star_redshift_is_pinned_to_zero(self, bank, tmp_path):
        """Even fitting at z = 0.05, the star rows report the z they were fit at."""

        results = run_fit(tmp_path, star_flux, z=0.05, mask_galaxy_lines=False)

        stars = results[results["SN"].astype(str).str.startswith("star")]
        assert len(stars) > 0
        assert (stars["Z"] == 0.0).all()

    def test_standalone_rows_can_be_plotted(self, bank, tmp_path, monkeypatch):
        import matplotlib

        matplotlib.use("Agg")
        from superfit import Superfit

        fit = Superfit(
            observed_spectrum(tmp_path, star_flux),
            output_dir=str(tmp_path / "out"),
            **{k: v for k, v in BASE.items() if k != "overwrite"}
        )
        results = fit.run().results

        for prefix in ("star-M/", "QSO/"):
            matches = results.index[
                results["SN"].astype(str).str.startswith(prefix)
            ]
            assert len(matches) > 0
            path = fit.plot_rank(int(matches[0]))
            assert os.path.getsize(path) > 0
