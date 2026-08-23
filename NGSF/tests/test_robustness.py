"""How the fit behaves when a spectrum contains extreme points.

Real spectra arrive with cosmic-ray hits, dead columns, negative sky
subtraction residuals and NaN gaps. Each test here drives one of those
through the code and asserts what currently happens.

Tests marked ``xfail(strict=True)`` describe behaviour that is wrong and
not yet fixed: they document the defect and will start failing loudly the
moment it is repaired, at which point the marker comes off and the
assertion becomes a guarantee.
"""

import numpy as np
import pytest

from NGSF.error_routines import linear_error, savitzky_golay
from NGSF.Header_Binnings import bin_spectrum, bin_spectrum_bank
from NGSF.SF_functions import remove_telluric


def spike(spectrum, index, height):
    """Return a copy of ``spectrum`` with one pixel multiplied up."""

    out = spectrum.copy()
    out[index, 1] = height
    return out


class TestSavitzkyGolayError:
    """The default error model (``"error_spectrum": "sg"``)."""

    def test_clean_spectrum_gives_positive_finite_errors(self, noisy_spectrum):
        err = savitzky_golay(noisy_spectrum)[:, 1]
        assert np.isfinite(err).all()
        assert (err > 0).all()

    def test_error_tracks_the_injected_noise_level(self, noisy_spectrum):
        """A 0.02 white-noise spectrum should get an error estimate near 0.02."""

        flux = noisy_spectrum[:, 1]
        err = savitzky_golay(noisy_spectrum)[:, 1]
        # savitzky_golay normalises by the mean before measuring the scatter.
        recovered = np.median(err) * np.mean(flux)
        assert 0.01 < recovered < 0.04

    @pytest.mark.xfail(
        strict=True,
        reason="savitzky_golay normalises with .mean() and savgol_filter has no "
        "NaN handling, so a single NaN pixel makes every error NaN",
    )
    def test_a_single_nan_does_not_poison_the_whole_error_spectrum(
        self, noisy_spectrum
    ):
        contaminated = noisy_spectrum.copy()
        contaminated[300, 1] = np.nan

        err = savitzky_golay(contaminated)[:, 1]

        # At most the neighbourhood of the bad pixel should be unusable.
        assert np.isfinite(err).sum() > 0.9 * err.size

    @pytest.mark.xfail(
        strict=True,
        reason="the .mean() normalisation is not outlier resistant, so one "
        "cosmic ray rescales the error estimate for every pixel",
    )
    def test_one_cosmic_ray_does_not_rescale_the_whole_error_spectrum(
        self, noisy_spectrum
    ):
        clean = savitzky_golay(noisy_spectrum)[:, 1]
        hit = savitzky_golay(spike(noisy_spectrum, 300, 100.0))[:, 1]

        changed = np.abs(hit - clean) / clean > 0.1
        # Only pixels near the spike should move.
        assert changed.sum() < 0.05 * changed.size

    @pytest.mark.xfail(
        strict=True,
        reason="zero residual is floored at 1e-50 rather than at a physical "
        "noise level, so an interpolated-flat stretch gets a chi2 weight "
        "many orders of magnitude above every real pixel",
    )
    def test_flat_regions_do_not_get_infinite_statistical_weight(
        self, noisy_spectrum
    ):
        """A noiseless stretch must not dominate chi2.

        A detector gap patched by interpolation is exactly flat over its
        whole width, so this is a real input rather than a contrived one.
        """

        patched = noisy_spectrum.copy()
        patched[200:400, 1] = 1.0  # the interpolated-over gap

        err = savitzky_golay(patched)[:, 1]
        weight = 1.0 / err**2

        assert weight.max() / np.median(weight) < 1e6


class TestLinearError:
    """The ``"error_spectrum": "linear"`` option."""

    @pytest.mark.xfail(
        strict=True,
        reason="the residual `r` is assigned inside the per-bin loop and "
        "overwritten each pass, so every bin's sigma is computed against the "
        "LAST bin's straight-line fit",
    )
    def test_matches_a_per_bin_standard_deviation(self, noisy_spectrum):
        lam = noisy_spectrum[:, 0]
        flux = noisy_spectrum[:, 1]
        num = 10

        got = linear_error(noisy_spectrum)[:, 1]

        expected = []
        for chunk_lam, chunk_flux in zip(
            lam.reshape(-1, num), flux.reshape(-1, num)
        ):
            slope, intercept = np.polyfit(chunk_lam, chunk_flux, 1)
            residual = chunk_flux - (slope * chunk_lam + intercept)
            expected.append(np.std(residual, ddof=1))
        expected = np.repeat(expected, num)

        np.testing.assert_allclose(got, expected, rtol=0.05)

    def test_returns_one_error_per_input_pixel(self, noisy_spectrum):
        out = linear_error(noisy_spectrum)
        assert out.shape == (noisy_spectrum.shape[0], 2)


class TestBinning:
    @pytest.mark.xfail(
        strict=True,
        reason="the guard reads lam[15] - lam[16], which is negative for any "
        "ascending wavelength axis, so the 'already coarser than target' "
        "branch is unreachable -- and returns None if it were reached",
    )
    def test_spectrum_coarser_than_the_target_is_passed_through(self):
        """Asking for 10 A bins from an 82 A spectrum must not invent structure.

        The intended guard is "if the data are already coarser than the
        target, leave them alone", but it compares lam[15] - lam[16], which
        is negative whenever wavelength ascends. So the spectrum is pushed
        through binned_statistic onto a finer grid than it was measured on.
        """

        lam = np.linspace(4000.0, 8000.0, 50)  # ~82 A spacing
        spectrum = np.column_stack([lam, np.ones_like(lam)])

        out = bin_spectrum(spectrum, 10)

        assert out is not None
        np.testing.assert_allclose(np.asarray(out["lam_bin"]), lam, rtol=1e-9)

    @pytest.mark.xfail(
        strict=True,
        reason="bin_spectrum_bank indexes lam[15] and lam[16] unconditionally",
    )
    def test_short_spectrum_raises_something_meaningful(self):
        lam = np.linspace(4000.0, 4200.0, 10)
        spectrum = np.column_stack([lam, np.ones_like(lam)])

        with pytest.raises(ValueError):
            bin_spectrum_bank(spectrum, 10)

    def test_binning_preserves_the_wavelength_range(self, noisy_spectrum):
        out = bin_spectrum_bank(noisy_spectrum, 40)
        assert out[0, 0] >= noisy_spectrum[0, 0]
        assert out[-1, 0] <= noisy_spectrum[-1, 0]

    def test_binned_flux_is_median_normalised(self, noisy_spectrum):
        out = bin_spectrum_bank(noisy_spectrum, 40)
        np.testing.assert_allclose(np.nanmedian(out[:, 1]), 1.0, rtol=1e-9)


class TestTelluric:
    def test_masks_the_a_band(self, smooth_spectrum):
        out = remove_telluric(smooth_spectrum.copy())
        in_band = (out[:, 0] >= 7594) & (out[:, 0] <= 7680)
        assert np.isnan(out[in_band, 1]).all()
        assert np.isfinite(out[~in_band, 1]).all()

    @pytest.mark.xfail(
        strict=True,
        reason="remove_telluric writes the -10000 sentinel straight into the "
        "array it was handed, so the caller's spectrum is destroyed",
    )
    def test_does_not_mutate_its_argument(self, smooth_spectrum):
        original = smooth_spectrum.copy()
        remove_telluric(smooth_spectrum)
        np.testing.assert_array_equal(smooth_spectrum, original)

    @pytest.mark.xfail(
        strict=True,
        reason="-10000 is used as an in-band sentinel and then mapped to NaN "
        "everywhere it appears, so a genuine -10000 flux is silently masked",
    )
    def test_a_real_flux_of_minus_10000_survives(self, smooth_spectrum):
        spectrum = smooth_spectrum.copy()
        # Well away from the telluric band.
        spectrum[10, 1] = -10000.0

        out = remove_telluric(spectrum)

        assert not np.isnan(out[10, 1])


class TestNormalisation:
    @pytest.mark.xfail(
        strict=True,
        reason="flux is divided by its median with no guard, so a spectrum "
        "straddling zero is amplified by ~1e15",
    )
    def test_median_near_zero_does_not_explode(self, noisy_spectrum):
        flux = noisy_spectrum[:, 1]
        centred = flux - np.median(flux)  # median is now ~1e-17

        normalised = centred / np.nanmedian(centred)

        assert np.abs(normalised).max() < 1e6
