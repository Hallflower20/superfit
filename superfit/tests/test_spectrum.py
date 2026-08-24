"""The Spectrum type: construction, coercion, and preprocessing.

Spectrum exists so that a fit does not require a file on disk. These tests
cover the ways one can be built and the guarantees the fitter relies on --
sorted wavelengths, a usable error column, and binning that is bit-for-bit
what the previous file-based path produced.
"""

import numpy as np
import pytest

from superfit.Header_Binnings import bin_spectrum_bank
from superfit.spectrum import Spectrum
from superfit.tests.conftest import TEST_SPECTRUM


class TestConstruction:
    def test_from_arrays(self):
        lam = np.linspace(4000.0, 8000.0, 100)
        flux = np.ones(100)

        spectrum = Spectrum(lam, flux, name="synthetic")

        assert len(spectrum) == 100
        assert spectrum.name == "synthetic"
        assert spectrum.error is None
        np.testing.assert_array_equal(spectrum.wavelength, lam)

    def test_from_arrays_with_error(self):
        lam = np.linspace(4000.0, 8000.0, 100)
        spectrum = Spectrum(lam, np.ones(100), error=np.full(100, 0.1))

        np.testing.assert_allclose(spectrum.error, 0.1)
        assert spectrum.as_array().shape == (100, 3)

    def test_from_file(self):
        spectrum = Spectrum.from_file(TEST_SPECTRUM)

        assert len(spectrum) > 1000
        assert spectrum.name.startswith("SN2021urb")
        assert spectrum.wavelength[0] < spectrum.wavelength[-1]

    def test_accepts_python_lists(self):
        spectrum = Spectrum([4000, 4010, 4020], [1.0, 2.0, 3.0])
        assert isinstance(spectrum.wavelength, np.ndarray)
        assert spectrum.wavelength.dtype == float

    def test_unsorted_input_is_sorted_with_columns_kept_together(self):
        """Descending wavelength is a common export quirk, not an error."""

        lam = np.array([5000.0, 4000.0, 6000.0])
        flux = np.array([2.0, 1.0, 3.0])
        error = np.array([0.2, 0.1, 0.3])

        spectrum = Spectrum(lam, flux, error=error)

        np.testing.assert_array_equal(spectrum.wavelength, [4000, 5000, 6000])
        np.testing.assert_array_equal(spectrum.flux, [1.0, 2.0, 3.0])
        np.testing.assert_array_equal(spectrum.error, [0.1, 0.2, 0.3])

    @pytest.mark.parametrize(
        "kwargs,match",
        [
            (dict(wavelength=[1, 2, 3], flux=[1, 2]), "different lengths"),
            (dict(wavelength=[], flux=[]), "empty"),
            (dict(wavelength=[1, 2], flux=[1, 2], error=[1]), "error has length"),
            (dict(wavelength=[1, np.nan], flux=[1, 2]), "NaN or inf"),
            (dict(wavelength=[1, 1], flux=[1, 2]), "duplicate"),
            (dict(wavelength=[1, 2], flux=[np.nan, np.nan]), "entirely NaN"),
        ],
    )
    def test_rejects_unusable_input(self, kwargs, match):
        with pytest.raises(ValueError, match=match):
            Spectrum(**kwargs)


class TestCoerce:
    def test_passes_a_spectrum_through(self):
        spectrum = Spectrum([4000, 4010], [1.0, 2.0])
        assert Spectrum.coerce(spectrum) is spectrum

    def test_reads_a_path(self):
        assert isinstance(Spectrum.coerce(TEST_SPECTRUM), Spectrum)

    def test_wraps_a_two_column_array(self):
        array = np.column_stack([np.linspace(4000, 5000, 50), np.ones(50)])
        spectrum = Spectrum.coerce(array)
        assert len(spectrum) == 50
        assert spectrum.error is None

    def test_wraps_a_three_column_array_as_wavelength_flux_error(self):
        array = np.column_stack(
            [np.linspace(4000, 5000, 50), np.ones(50), np.full(50, 0.1)]
        )
        spectrum = Spectrum.coerce(array)
        np.testing.assert_allclose(spectrum.error, 0.1)

    def test_rejects_something_that_is_not_a_spectrum(self):
        with pytest.raises(TypeError, match="cannot interpret"):
            Spectrum.coerce(42)


class TestPreprocessing:
    def test_normalised_has_unit_median(self, noisy_spectrum):
        spectrum = Spectrum(noisy_spectrum[:, 0], noisy_spectrum[:, 1] * 1e-17)
        np.testing.assert_allclose(np.median(spectrum.normalised().flux), 1.0)

    def test_without_telluric_blanks_the_band_and_keeps_the_rest(self):
        lam = np.linspace(7000.0, 8000.0, 400)
        spectrum = Spectrum(lam, np.ones_like(lam))

        out = spectrum.without_telluric()

        in_band = (lam >= 7594) & (lam <= 7680)
        assert np.isnan(out.flux[in_band]).all()
        assert np.isfinite(out.flux[~in_band]).all()

    def test_without_host_lines_drops_pixels(self):
        lam = np.linspace(4000.0, 8000.0, 4000)
        spectrum = Spectrum(lam, np.ones_like(lam))

        out = spectrum.without_host_lines(0.0)

        assert len(out) < len(spectrum)

    def test_preprocessing_does_not_mutate_the_original(self, noisy_spectrum):
        spectrum = Spectrum(noisy_spectrum[:, 0], noisy_spectrum[:, 1])
        before = spectrum.flux.copy()

        spectrum.normalised()
        spectrum.without_telluric()
        spectrum.binned(20)

        np.testing.assert_array_equal(spectrum.flux, before)

    def test_interpolated_onto_is_nan_outside_coverage(self):
        spectrum = Spectrum(np.linspace(5000.0, 6000.0, 100), np.ones(100))

        out = spectrum.interpolated_onto(np.linspace(4000.0, 7000.0, 50))

        assert np.isnan(out[0]) and np.isnan(out[-1])
        assert np.isfinite(out).any()


class TestBinning:
    def test_flux_matches_the_original_binning_routine(self):
        """The file-based path used bin_spectrum_bank; results must not move."""

        spectrum = Spectrum.from_file(TEST_SPECTRUM)
        reference = bin_spectrum_bank(
            np.column_stack([spectrum.wavelength, spectrum.flux]), 10
        )

        binned = spectrum.binned(10)

        np.testing.assert_array_equal(binned.wavelength, reference[:, 0])
        np.testing.assert_array_equal(binned.flux, reference[:, 1])

    def test_an_error_column_does_not_perturb_the_flux(self):
        spectrum = Spectrum.from_file(TEST_SPECTRUM)
        with_error = spectrum.copy(error=np.abs(spectrum.flux) * 0.05 + 1e-3)

        np.testing.assert_array_equal(
            spectrum.binned(10).flux, with_error.binned(10).flux
        )

    def test_error_is_carried_through_binning(self):
        spectrum = Spectrum.from_file(TEST_SPECTRUM)
        spectrum = spectrum.copy(error=np.abs(spectrum.flux) * 0.05 + 1e-3)

        binned = spectrum.binned(10)

        assert binned.error is not None
        assert binned.error.shape == binned.flux.shape
        assert np.isfinite(binned.error).all()
        assert (binned.error > 0).all()

    def test_binning_reduces_the_sample_count(self):
        spectrum = Spectrum.from_file(TEST_SPECTRUM)
        assert len(spectrum.binned(10)) < len(spectrum)

    @pytest.mark.parametrize("sampling", [10.0, 12.0, 25.0])
    def test_a_spectrum_already_coarser_than_the_resolution_keeps_its_error(
        self, sampling
    ):
        """bin_spectrum_bank passes such a spectrum through unbinned. The
        error has to take the same branch, or the two come back different
        lengths and Spectrum rejects the pair -- which made a three-column
        spectrum sampled at exactly the default 10 A unfittable."""

        lam = np.arange(4000.0, 9000.0, sampling)
        flux = 1.0 + 0.1 * np.sin(lam / 200.0)
        spectrum = Spectrum(lam, flux, error=np.full_like(flux, 0.02))

        binned = spectrum.binned(10)

        assert binned.error is not None
        assert binned.error.shape == binned.flux.shape
        assert np.isfinite(binned.error).all()

    def test_the_passthrough_error_is_on_the_same_scale_as_the_flux(self):
        """The flux is divided by its median on the way through; so is the error."""

        lam = np.arange(4000.0, 9000.0, 10.0)
        flux = np.full_like(lam, 4.0)
        flux[::2] = 6.0                        # median 5, well away from zero
        spectrum = Spectrum(lam, flux, error=np.full_like(flux, 0.5))

        binned = spectrum.binned(10)

        np.testing.assert_allclose(binned.error, 0.5 / 5.0)


class TestLengthCheck:
    def test_accepts_a_normal_spectrum(self):
        Spectrum.from_file(TEST_SPECTRUM).check_long_enough(10)

    def test_rejects_a_stub(self):
        lam = np.linspace(5000.0, 5050.0, 40)
        spectrum = Spectrum(lam, np.ones_like(lam), name="stub")

        with pytest.raises(ValueError, match="At least"):
            spectrum.check_long_enough(10)

    def test_the_message_names_the_spectrum(self):
        lam = np.linspace(5000.0, 5050.0, 40)
        spectrum = Spectrum(lam, np.ones_like(lam), name="my_object")

        with pytest.raises(ValueError, match="my_object"):
            spectrum.check_long_enough(10)
