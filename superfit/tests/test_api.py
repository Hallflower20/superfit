"""The public API: the ways a fit can be set up, and that they agree.

The point of accepting arrays is that a spectrum already in memory does not
have to be written to a file first. That is only worth anything if it gives
the same answer as the file would have, so the central test here fits the
same data both ways and compares.
"""

import json
import os

import numpy as np
import pytest

from superfit import Spectrum, Superfit
from superfit.config import ConfigError
from superfit.tests.conftest import TEST_SPECTRUM, needs_bank

pytestmark = [needs_bank, pytest.mark.endtoend]

SCIENCE_COLUMNS = [
    "GALAXY",
    "SN",
    "Z",
    "A_v",
    "CONST_SN",
    "CONST_GAL",
    "Frac(SN)",
    "Frac(gal)",
    "CHI2/dof",
    "CHI2/dof2",
]

# Settings shared by every fit here, minus the redshift.
BASE_KWARGS = dict(resolution=10, n_plots=0, show_plot=0)

FIT_KWARGS = dict(z=0.127, **BASE_KWARGS)


def outdir(tmp_path, name):
    """An existing, empty directory to write one fit's results into."""

    path = tmp_path / name
    path.mkdir(parents=True, exist_ok=True)
    return str(path) + "/"


def assert_same_science(a, b):
    """Two result tables agree on everything except what the input was called."""

    import pandas as pd

    assert len(a) == len(b)
    for column in SCIENCE_COLUMNS:
        # pandas >= 3 gives text columns a str dtype rather than object, so
        # ask whether the column is numeric rather than whether it is not.
        if pd.api.types.is_numeric_dtype(a[column]):
            np.testing.assert_allclose(
                a[column].to_numpy(), b[column].to_numpy(), rtol=1e-10, err_msg=column
            )
        else:
            assert list(a[column]) == list(b[column]), column


class TestInputForms:
    @pytest.fixture(scope="class")
    def from_file(self, tmp_path_factory):
        out = str(tmp_path_factory.mktemp("file")) + "/"
        return Superfit(TEST_SPECTRUM, output_dir=out, **FIT_KWARGS).run()

    def test_file_path_positional(self, from_file):
        assert len(from_file) > 0
        assert from_file["CHI2/dof"].notna().all()

    def test_arrays_match_the_file(self, from_file, tmp_path):
        """The headline guarantee: in-memory input is not a different fit."""

        data = np.loadtxt(TEST_SPECTRUM)
        out = str(tmp_path) + "/"

        from_arrays = Superfit(
            wavelength=data[:, 0],
            flux=data[:, 1],
            name="in_memory",
            output_dir=out,
            **FIT_KWARGS
        ).run()

        assert_same_science(from_file, from_arrays)

    def test_spectrum_object_matches_the_file(self, from_file, tmp_path):
        spectrum = Spectrum.from_file(TEST_SPECTRUM)
        out = str(tmp_path) + "/"

        result = Superfit(spectrum, output_dir=out, **FIT_KWARGS).run()

        assert_same_science(from_file, result)

    def test_raw_array_matches_the_file(self, from_file, tmp_path):
        data = np.loadtxt(TEST_SPECTRUM)
        out = str(tmp_path) + "/"

        result = Superfit(data, name="raw", output_dir=out, **FIT_KWARGS).run()

        assert_same_science(from_file, result)

    def test_one_spectrum_reused_for_several_redshifts(self, tmp_path):
        """Spectrum is immutable enough to fit repeatedly without reloading."""

        spectrum = Spectrum.from_file(TEST_SPECTRUM)
        best = {}
        for z in (0.12, 0.127):
            out = str(tmp_path / "z{}".format(z)) + "/"
            os.makedirs(out, exist_ok=True)
            results = Superfit(
                spectrum, output_dir=out, z=z, resolution=10, n_plots=0, show_plot=0
            ).run()
            best[z] = results["CHI2/dof2"].iloc[0]

        assert set(best) == {0.12, 0.127}
        assert all(np.isfinite(v) for v in best.values())


class TestInstanceIsolation:
    """Two fits in one process must not be able to retune each other.

    Settings used to be installed into a module-level global by the
    constructor and read back by ``run()``, so the fit that ran was
    configured by whichever Superfit was built *last*. Silently fitting at
    someone else's redshift is about the worst failure mode this package
    could have, and it also ruled out any concurrency.
    """

    def test_building_a_second_fit_does_not_disturb_the_first(self, tmp_path):
        first = Superfit(
            TEST_SPECTRUM, output_dir=outdir(tmp_path, "a"), z=0.10, **BASE_KWARGS
        )
        second = Superfit(
            TEST_SPECTRUM, output_dir=outdir(tmp_path, "b"), z=0.20, **BASE_KWARGS
        )

        assert first.parameters.redshift == pytest.approx([0.10])
        assert second.parameters.redshift == pytest.approx([0.20])

        # The bug: `first.run()` after `second` was constructed fitted at 0.20.
        results = first.run()
        assert results["Z"].to_numpy() == pytest.approx(0.10)

    def test_each_fit_keeps_its_own_output_location(self, tmp_path):
        a, b = outdir(tmp_path, "a"), outdir(tmp_path, "b")
        first = Superfit(TEST_SPECTRUM, output_dir=a, z=0.1, **BASE_KWARGS)
        Superfit(TEST_SPECTRUM, output_dir=b, z=0.2, **BASE_KWARGS)

        assert first.parameters.save_results_path == a

    def test_parameters_cannot_be_mutated_after_construction(self, tmp_path):
        fit = Superfit(
            TEST_SPECTRUM, output_dir=str(tmp_path) + "/", z=0.1, **BASE_KWARGS
        )
        with pytest.raises(AttributeError, match="immutable"):
            fit.parameters.resolution = 30


class TestConfigForms:
    def test_dict_config_still_works_positionally(self, tmp_path, base_parameters):
        """Superfit(config_dict) predates the spectrum argument."""

        config = dict(base_parameters)
        config.update(
            object_to_fit=TEST_SPECTRUM,
            saving_results_path=str(tmp_path) + "/",
            how_many_plots=0,
            show_plot=0,
        )

        results = Superfit(config).run()
        assert len(results) > 0

    def test_json_file_config_still_works_positionally(self, tmp_path, base_parameters):
        config = dict(base_parameters)
        config.update(
            object_to_fit=TEST_SPECTRUM,
            saving_results_path=str(tmp_path) + "/",
            how_many_plots=0,
            show_plot=0,
        )
        path = tmp_path / "params.json"
        path.write_text(json.dumps(config))

        assert len(Superfit(str(path)).run()) > 0

    def test_config_and_keyword_overrides_combine(self, tmp_path):
        results = Superfit(
            TEST_SPECTRUM,
            config={"resolution": 10, "minimum_overlap": 0.7},
            output_dir=str(tmp_path) + "/",
            **FIT_KWARGS
        ).run()

        assert len(results) > 0

    def test_used_config_is_written_and_complete(self, tmp_path):
        out = str(tmp_path) + "/"
        Superfit(TEST_SPECTRUM, output_dir=out, **FIT_KWARGS).run()

        from superfit.config import DEFAULT_CONFIG

        written = [f for f in os.listdir(out) if f.endswith("_used.json")]
        assert written, "no *_used.json written"

        with open(os.path.join(out, written[0])) as handle:
            recorded = json.load(handle)

        # A three-keyword call must still record every effective setting.
        assert set(recorded) == set(DEFAULT_CONFIG)
        assert recorded["z_exact"] == 0.127


class TestErrors:
    def test_no_spectrum_at_all_is_a_clear_error(self):
        with pytest.raises(ConfigError, match="No spectrum to fit"):
            Superfit(z=0.1)

    def test_spectrum_and_arrays_together_is_rejected(self):
        data = np.loadtxt(TEST_SPECTRUM)
        with pytest.raises(TypeError, match="not both"):
            Superfit(TEST_SPECTRUM, wavelength=data[:, 0], flux=data[:, 1])

    def test_wavelength_without_flux_is_rejected(self):
        with pytest.raises(TypeError, match="together"):
            Superfit(wavelength=np.linspace(4000, 5000, 10))

    def test_included_errors_without_an_error_column(self, tmp_path):
        with pytest.raises(ValueError, match="needs an uncertainty column"):
            Superfit(
                TEST_SPECTRUM,
                output_dir=str(tmp_path) + "/",
                error_model="included",
                **FIT_KWARGS
            ).run()

    def test_included_errors_work_when_supplied(self, tmp_path):
        data = np.loadtxt(TEST_SPECTRUM)
        results = Superfit(
            wavelength=data[:, 0],
            flux=data[:, 1],
            error=np.abs(data[:, 1]) * 0.05 + 1e-3,
            name="with_error",
            output_dir=str(tmp_path) + "/",
            error_model="included",
            **FIT_KWARGS
        ).run()

        assert len(results) > 0
        assert np.isfinite(results["CHI2/dof"]).all()


def test_run_returns_the_results_and_sets_the_attribute(tmp_path):
    fit = Superfit(TEST_SPECTRUM, output_dir=str(tmp_path) + "/", **FIT_KWARGS)
    returned = fit.run()

    assert returned is fit.results
    assert list(returned.columns)[:3] == ["SPECTRUM", "GALAXY", "SN"]


def test_superfit_method_is_an_alias_for_run(tmp_path):
    fit = Superfit(TEST_SPECTRUM, output_dir=str(tmp_path) + "/", **FIT_KWARGS)
    assert len(fit.superfit()) > 0
