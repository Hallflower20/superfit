"""How the fit behaves when a spectrum contains extreme points.

Real spectra arrive with cosmic-ray hits, dead columns, negative sky
subtraction residuals and NaN gaps. Each test here drives one of those
through the code.

Every assertion in this file started life as a strict xfail recording a real
defect. They are guarantees now; if one starts failing, the corresponding
failure mode has come back.
"""

import warnings

import numpy as np
import pytest

from NGSF.error_routines import (
    apply_error_floor,
    linear_error,
    robust_scale,
    savitzky_golay,
)
from NGSF.Header_Binnings import bin_spectrum, bin_spectrum_bank, normalise_flux
from NGSF.SF_functions import remove_telluric


def spike(spectrum, index, height):
    """Return a copy of ``spectrum`` with one pixel driven to ``height``."""

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
        recovered = np.median(err) * robust_scale(flux)
        assert 0.01 < recovered < 0.04

    def test_a_single_nan_does_not_poison_the_whole_error_spectrum(
        self, noisy_spectrum
    ):
        """One masked pixel must not turn every uncertainty into NaN.

        It used to: the flux was normalised by .mean() and smoothed by
        savgol_filter, neither of which tolerates a NaN, so every chi2 in the
        run came out NaN.
        """

        contaminated = noisy_spectrum.copy()
        contaminated[300, 1] = np.nan

        err = savitzky_golay(contaminated)[:, 1]

        assert np.isfinite(err).all()
        assert (err > 0).all()

    def test_many_nan_gaps_still_give_a_usable_error(self, noisy_spectrum, rng):
        contaminated = noisy_spectrum.copy()
        contaminated[rng.choice(600, size=60, replace=False), 1] = np.nan

        err = savitzky_golay(contaminated)[:, 1]

        assert np.isfinite(err).all()
        clean = savitzky_golay(noisy_spectrum)[:, 1]
        # Same ballpark as the uncontaminated estimate.
        assert 0.5 < np.median(err) / np.median(clean) < 2.0

    def test_a_cosmic_ray_stays_local(self, noisy_spectrum):
        """A spike may inflate its own neighbourhood, not the whole spectrum.

        The old .mean() normalisation rescaled every pixel's error when a
        single hot pixel appeared. Contamination is now bounded by the width
        of the scatter window around the hit.
        """

        clean = savitzky_golay(noisy_spectrum)[:, 1]
        hit = savitzky_golay(spike(noisy_spectrum, 300, 100.0))[:, 1]

        changed = np.where(np.abs(hit - clean) / clean > 0.1)[0]

        assert changed.size > 0, "the spike should affect something"
        # Confined to the scatter window either side of pixel 300.
        assert changed.min() > 300 - 120
        assert changed.max() < 300 + 120

    def test_a_cosmic_ray_does_not_move_the_overall_error_scale(self, noisy_spectrum):
        clean = savitzky_golay(noisy_spectrum)[:, 1]
        hit = savitzky_golay(spike(noisy_spectrum, 300, 100.0))[:, 1]

        np.testing.assert_allclose(np.median(hit), np.median(clean), rtol=0.02)

    def test_flat_regions_do_not_get_infinite_statistical_weight(
        self, noisy_spectrum
    ):
        """A noiseless stretch must not dominate chi2.

        A detector gap patched by interpolation is exactly flat over its whole
        width. Its residual scatter is zero, and the old code floored that at
        1e-50 -- a chi2 weight of ~1e100 for pixels carrying no information.
        """

        patched = noisy_spectrum.copy()
        patched[200:400, 1] = 1.0

        err = savitzky_golay(patched)[:, 1]
        weight = 1.0 / err**2

        assert weight.max() / np.median(weight) < 1e6

    @pytest.mark.parametrize("true_sigma", [0.005, 0.02, 0.1])
    def test_recovers_a_known_noise_level(self, true_sigma):
        """Absolute calibration, not just self-consistency.

        Chi2 values are only meaningful if sigma is unbiased. Injecting known
        Gaussian noise on a smooth continuum should come back within a few
        percent. This is what says the estimator is honest rather than merely
        different from the previous one.
        """

        recovered = []
        for seed in range(8):
            rng = np.random.default_rng(seed)
            lam = np.linspace(4000.0, 8000.0, 1500)
            flux = 1.0 + 0.1 * np.sin(lam / 200.0)
            flux = flux + rng.normal(0.0, true_sigma, lam.size)
            spectrum = np.column_stack([lam, flux])

            err = savitzky_golay(spectrum)[:, 1]
            # savitzky_golay works on a normalised flux; undo that to compare.
            recovered.append(np.median(err) * robust_scale(flux))

        assert np.median(recovered) == pytest.approx(true_sigma, rel=0.1)

    def test_is_not_inflated_by_a_strong_emission_line(self):
        """A real line is signal, not noise, and must not raise sigma near it.

        The old estimator averaged residual**2 over a 100-pixel window, so
        every strong line inflated the error across its whole neighbourhood
        and quietly suppressed chi2 there.
        """

        rng = np.random.default_rng(11)
        lam = np.linspace(4000.0, 8000.0, 1500)
        noise = rng.normal(0.0, 0.02, lam.size)
        continuum = 1.0 + noise
        line = 5.0 * np.exp(-0.5 * ((lam - 6000.0) / 8.0) ** 2)

        without = savitzky_golay(np.column_stack([lam, continuum]))[:, 1]
        with_line = savitzky_golay(np.column_stack([lam, continuum + line]))[:, 1]

        near = np.abs(lam - 6000.0) < 200.0
        ratio = np.median(with_line[near]) / np.median(without[near])

        assert ratio < 2.0, "error near the line inflated by {:.1f}x".format(ratio)

    def test_rejects_a_spectrum_shorter_than_the_smoothing_window(self):
        short = np.column_stack([np.linspace(4000.0, 4100.0, 12), np.ones(12)])
        with pytest.raises(ValueError, match="smoothing window"):
            savitzky_golay(short)


class TestErrorFloor:
    def test_raises_zeros_to_a_fraction_of_the_typical_error(self):
        err = np.array([0.0, 0.1, 0.2, 0.1, 0.0])
        out = apply_error_floor(err, fraction=1e-2)
        assert (out > 0).all()
        assert out.max() == pytest.approx(0.2)

    def test_leaves_healthy_errors_alone(self):
        err = np.array([0.1, 0.2, 0.15])
        np.testing.assert_allclose(apply_error_floor(err), err)

    def test_all_zero_input_does_not_produce_zero_weights(self):
        out = apply_error_floor(np.zeros(10))
        assert np.isfinite(out).all()
        assert (out > 0).all()


class TestLinearError:
    """The ``"error_spectrum": "linear"`` option."""

    def test_matches_a_per_bin_standard_deviation(self, noisy_spectrum):
        """Each block's sigma comes from that block's own straight-line fit.

        It used to come from the *last* block's fit: the residual was assigned
        inside the loop and overwritten on every pass, making this error mode
        wrong by up to a factor of ~1000.
        """

        lam = noisy_spectrum[:, 0]
        flux = noisy_spectrum[:, 1]
        num = 10

        got = linear_error(noisy_spectrum, num=num)[:, 1]

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

    def test_handles_a_length_that_is_not_a_multiple_of_the_block(
        self, noisy_spectrum
    ):
        ragged = noisy_spectrum[:597]
        out = linear_error(ragged)
        assert out.shape == (597, 2)
        assert np.isfinite(out[:, 1]).all()

    def test_survives_nan_blocks(self, noisy_spectrum):
        contaminated = noisy_spectrum.copy()
        contaminated[100:130, 1] = np.nan

        out = linear_error(contaminated)

        assert np.isfinite(out[:, 1]).all()
        assert (out[:, 1] > 0).all()


class TestBinning:
    def test_spectrum_coarser_than_the_target_is_passed_through(self):
        """Asking for 10 A bins from an 82 A spectrum must not invent structure.

        The guard was written as lam[15] - lam[16], negative for any ascending
        wavelength axis, so the spectrum was pushed onto a grid finer than it
        was ever measured on.
        """

        lam = np.linspace(4000.0, 8000.0, 50)  # ~82 A spacing
        spectrum = np.column_stack([lam, np.ones_like(lam)])

        out = bin_spectrum(spectrum, 10)

        assert out is not None
        np.testing.assert_allclose(np.asarray(out["lam_bin"]), lam, rtol=1e-9)

    def test_short_spectrum_does_not_raise_indexerror(self):
        """lam[15] used to be indexed unconditionally."""

        lam = np.linspace(4000.0, 4200.0, 10)
        spectrum = np.column_stack([lam, np.ones_like(lam)])

        out = bin_spectrum_bank(spectrum, 10)

        assert out.shape == (10, 2)

    def test_two_pixel_spectrum_gives_a_clear_error(self):
        spectrum = np.column_stack([[4000.0], [1.0]])
        with pytest.raises(ValueError, match="at least two"):
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
        out = remove_telluric(smooth_spectrum)
        in_band = (out[:, 0] >= 7594) & (out[:, 0] <= 7680)
        assert np.isnan(out[in_band, 1]).all()
        assert np.isfinite(out[~in_band, 1]).all()

    def test_does_not_mutate_its_argument(self, smooth_spectrum):
        """It used to write a -10000 sentinel into the caller's own array."""

        original = smooth_spectrum.copy()
        remove_telluric(smooth_spectrum)
        np.testing.assert_array_equal(smooth_spectrum, original)

    def test_a_real_flux_of_minus_10000_survives(self, smooth_spectrum):
        """-10000 was a sentinel, so genuine -10000 flux was masked anywhere."""

        spectrum = smooth_spectrum.copy()
        spectrum[10, 1] = -10000.0  # well away from the telluric band

        out = remove_telluric(spectrum)

        assert out[10, 1] == -10000.0

    def test_keeps_a_shipped_error_column(self):
        """Three-column input used to come back with the error column dropped."""

        lam = np.linspace(7000.0, 8000.0, 200)
        spectrum = np.column_stack([lam, np.ones_like(lam), np.full(lam.size, 0.1)])

        out = remove_telluric(spectrum)

        assert out.shape == (200, 3)
        np.testing.assert_allclose(out[:, 2], 0.1)


class TestNormalisation:
    def test_normal_spectrum_gets_a_median_of_one(self, noisy_spectrum):
        out = normalise_flux(noisy_spectrum[:, 1])
        np.testing.assert_allclose(np.median(out), 1.0, rtol=1e-12)

    def test_median_near_zero_does_not_explode(self, noisy_spectrum):
        """Dividing by a ~1e-17 median used to amplify the spectrum by ~1e15."""

        centred = noisy_spectrum[:, 1] - np.median(noisy_spectrum[:, 1])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            normalised = normalise_flux(centred)

        assert np.abs(normalised).max() < 1e6

    def test_median_near_zero_warns(self, noisy_spectrum):
        centred = noisy_spectrum[:, 1] - np.median(noisy_spectrum[:, 1])

        with pytest.warns(RuntimeWarning, match="negligible"):
            normalise_flux(centred)

    def test_constant_spectrum_is_rejected(self):
        with pytest.raises(ValueError, match="constant"):
            normalise_flux(np.zeros(100))

    def test_ignores_nan(self, noisy_spectrum):
        flux = noisy_spectrum[:, 1].copy()
        flux[::10] = np.nan

        out = normalise_flux(flux)

        np.testing.assert_allclose(np.nanmedian(out), 1.0, rtol=1e-12)
