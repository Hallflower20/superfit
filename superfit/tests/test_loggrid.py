"""Log-wavelength grids and redshifting by array shift.

The whole speedup rests on one claim: shifting a template along a log grid
is the same operation as interpolating it onto a redshifted wavelength axis.
These tests check that claim directly, against ``np.interp``, because the
failure mode is quiet -- a half-bin error is a few hundred km/s, which does
not crash anything, it just makes every fit slightly wrong.

(A weight-swap bug of exactly that kind is why these exist.)
"""

import numpy as np
import pytest

from superfit.loggrid import (
    C_KM_S,
    LogGrid,
    RedshiftableTemplates,
    dlnlam_to_velocity,
    matched_velocity_resolution,
    velocity_to_dlnlam,
)


class TestConversions:
    def test_velocity_round_trip(self):
        assert dlnlam_to_velocity(velocity_to_dlnlam(500.0)) == pytest.approx(500.0)

    def test_one_bin_is_the_stated_velocity(self):
        """A step of dlnlam must correspond to c * dlnlam km/s."""

        dlnlam = velocity_to_dlnlam(300.0)
        grid = LogGrid(np.log(5000.0), dlnlam, 10)
        lam = grid.wavelength

        measured = C_KM_S * (lam[1] - lam[0]) / lam[0]
        assert measured == pytest.approx(300.0, rel=1e-3)

    @pytest.mark.parametrize("velocity", [-1.0, 0.0])
    def test_non_positive_velocity_is_rejected(self, velocity):
        with pytest.raises(ValueError, match="must be positive"):
            velocity_to_dlnlam(velocity)

    def test_matched_resolution_reproduces_the_linear_point_count(self):
        """A 10 A grid over the test spectrum's span is ~739 points."""

        lower, upper, resolution = 3053.4, 10447.7, 10
        n_linear = int((upper - lower) / resolution)

        velocity = matched_velocity_resolution(resolution, lower, upper)
        grid = LogGrid.spanning(lower, upper, velocity_to_dlnlam(velocity))

        assert abs(len(grid) - n_linear) <= 2

    def test_matched_resolution_is_a_few_hundred_km_s_for_optical(self):
        velocity = matched_velocity_resolution(10, 3000.0, 10000.0)
        assert 400 < velocity < 600


class TestLogGrid:
    def test_spacing_is_uniform_in_log(self):
        grid = LogGrid.spanning(3000.0, 10000.0, velocity_to_dlnlam(500.0))
        steps = np.diff(np.log(grid.wavelength))
        # 1e-9 of a 500 km/s bin is sub-mm/s; the residual spread is just the
        # exp/log round trip.
        np.testing.assert_allclose(steps, steps[0], rtol=1e-9)

    def test_spanning_covers_the_requested_range(self):
        grid = LogGrid.spanning(3000.0, 10000.0, velocity_to_dlnlam(500.0))
        lam = grid.wavelength
        assert lam[0] == pytest.approx(3000.0)
        assert lam[-1] >= 10000.0

    def test_shift_bins_is_zero_at_zero_redshift(self):
        grid = LogGrid.spanning(3000.0, 10000.0, velocity_to_dlnlam(500.0))
        assert grid.shift_bins(0.0) == 0.0

    def test_shift_bins_matches_the_definition(self):
        dlnlam = velocity_to_dlnlam(500.0)
        grid = LogGrid.spanning(3000.0, 10000.0, dlnlam)
        assert grid.shift_bins(0.127) == pytest.approx(np.log1p(0.127) / dlnlam)

    def test_extended_blueward_keeps_the_step_and_the_red_end(self):
        grid = LogGrid.spanning(3000.0, 10000.0, velocity_to_dlnlam(500.0))
        wider = grid.extended_blueward(50)

        assert wider.dlnlam == grid.dlnlam
        assert len(wider) == len(grid) + 50
        np.testing.assert_allclose(wider.wavelength[50:], grid.wavelength, rtol=1e-12)

    @pytest.mark.parametrize("bad", [dict(dlnlam=0), dict(dlnlam=-1e-3)])
    def test_rejects_a_non_positive_step(self, bad):
        with pytest.raises(ValueError, match="dlnlam must be positive"):
            LogGrid(np.log(3000.0), bad["dlnlam"], 10)


class TestRedshifting:
    """Shifting must agree with interpolating a redshifted axis."""

    @staticmethod
    def _bank(n_templates=4):
        """Templates whose structure is resolved by the grid.

        Deliberately smooth. Shifting on the log grid resamples twice --
        native to rest grid, then the shift -- so structure at the sampling
        scale gets averaged down, and comparing that against a single direct
        interpolation would be measuring the resampling, not the shift.
        ``test_unresolved_structure_is_smoothed`` covers that case on its own.
        """

        lam = np.linspace(2500.0, 12000.0, 4000)
        fluxes = [
            1.0 + 0.3 * np.sin(lam / (120.0 + 30 * i)) for i in range(n_templates)
        ]
        return [lam] * n_templates, fluxes

    @staticmethod
    def _observed_grid():
        return LogGrid.spanning(3200.0, 9500.0, velocity_to_dlnlam(400.0))

    @pytest.mark.parametrize("z", [0.0, 0.01, 0.0537, 0.1, 0.127, 0.3])
    def test_matches_direct_interpolation(self, z):
        """The reference: resample the template onto lam_obs / (1 + z)."""

        wavelengths, fluxes = self._bank()
        obs = self._observed_grid()
        bank = RedshiftableTemplates.from_templates(
            wavelengths, fluxes, obs, max_z=0.5
        )

        shifted = bank.at_redshift(z)

        rest_wanted = obs.wavelength / (1.0 + z)
        for i, (lam, flux) in enumerate(zip(wavelengths, fluxes)):
            expected = np.interp(rest_wanted, lam, flux, left=np.nan, right=np.nan)
            ok = np.isfinite(shifted[i]) & np.isfinite(expected)
            # Resolved structure survives the extra resampling at the 1e-3
            # level; a half-bin error would show up here as ~0.1.
            np.testing.assert_allclose(
                shifted[i][ok], expected[ok], rtol=5e-3, atol=5e-3
            )

    def test_unresolved_structure_is_smoothed(self):
        """Honest about the cost: shifting resamples twice.

        Structure finer than the grid -- pixel-to-pixel noise -- is averaged
        down rather than reproduced. Resolved features are not affected, so
        this bounds the effect rather than treating it as a defect.
        """

        obs = self._observed_grid()
        lam = np.linspace(2500.0, 12000.0, 4000)
        rng = np.random.default_rng(0)

        smooth = 1.0 + 0.3 * np.sin(lam / 100.0)
        noisy = 1.0 + 0.05 * rng.normal(size=lam.size)

        bank = RedshiftableTemplates.from_templates(
            [lam, lam], [smooth, noisy], obs, max_z=0.5
        )
        shifted = bank.at_redshift(0.127)
        rest_wanted = obs.wavelength / 1.127

        errors = []
        for row, flux in zip(shifted, (smooth, noisy)):
            direct = np.interp(rest_wanted, lam, flux, left=np.nan, right=np.nan)
            ok = np.isfinite(row) & np.isfinite(direct)
            errors.append(np.sqrt(np.mean((row[ok] - direct[ok]) ** 2)))

        smooth_rms, noisy_rms = errors
        assert smooth_rms < 1e-3, "resolved structure should survive"
        # Noise at the sampling scale is averaged, not reproduced.
        assert 0.5 < noisy_rms / 0.05 < 1.2

    def test_zero_redshift_is_the_identity(self):
        """z = 0 must return the templates untouched, not almost untouched."""

        wavelengths, fluxes = self._bank()
        obs = self._observed_grid()
        bank = RedshiftableTemplates.from_templates(
            wavelengths, fluxes, obs, max_z=0.5
        )

        at_rest = bank.at_redshift(0.0)
        directly = bank.values[:, bank.offset : bank.offset + len(obs)]

        np.testing.assert_array_equal(at_rest, directly)

    def test_an_exact_bin_shift_is_a_pure_slice(self):
        """A redshift worth a whole number of bins needs no blending at all."""

        obs = self._observed_grid()
        # Choose z so that ln(1+z) is exactly 7 bins.
        z = np.expm1(7 * obs.dlnlam)

        wavelengths, fluxes = self._bank()
        bank = RedshiftableTemplates.from_templates(
            wavelengths, fluxes, obs, max_z=0.5
        )

        shifted = bank.at_redshift(z)
        expected = bank.values[
            :, bank.offset - 7 : bank.offset - 7 + len(obs)
        ]

        np.testing.assert_allclose(shifted, expected, rtol=1e-12)

    def test_a_half_bin_shift_is_the_average_of_its_neighbours(self):
        """Pins the interpolation weights; swapping them passes every other test."""

        obs = self._observed_grid()
        z = np.expm1(3.5 * obs.dlnlam)

        wavelengths, fluxes = self._bank()
        bank = RedshiftableTemplates.from_templates(
            wavelengths, fluxes, obs, max_z=0.5
        )

        shifted = bank.at_redshift(z)
        n = len(obs)
        lo = bank.values[:, bank.offset - 4 : bank.offset - 4 + n]
        hi = bank.values[:, bank.offset - 3 : bank.offset - 3 + n]

        np.testing.assert_allclose(shifted, 0.5 * (lo + hi), rtol=1e-12)

    def test_shift_direction_is_towards_the_red(self):
        """A redshifted feature must move to longer wavelength, not shorter."""

        obs = LogGrid.spanning(4000.0, 8000.0, velocity_to_dlnlam(200.0))
        lam = np.linspace(3000.0, 12000.0, 20000)
        # A narrow line at 6000 A in the rest frame.
        flux = np.exp(-0.5 * ((lam - 6000.0) / 3.0) ** 2)

        bank = RedshiftableTemplates.from_templates([lam], [flux], obs, max_z=0.5)

        peak_rest = obs.wavelength[np.nanargmax(bank.at_redshift(0.0)[0])]
        peak_z = obs.wavelength[np.nanargmax(bank.at_redshift(0.1)[0])]

        assert peak_rest == pytest.approx(6000.0, rel=2e-3)
        assert peak_z == pytest.approx(6600.0, rel=2e-3)

    def test_every_template_shifts_by_the_same_amount(self):
        """The point of log space: one displacement for the whole bank."""

        obs = LogGrid.spanning(4000.0, 8000.0, velocity_to_dlnlam(200.0))
        lam = np.linspace(3000.0, 12000.0, 20000)
        fluxes = [
            np.exp(-0.5 * ((lam - centre) / 3.0) ** 2)
            for centre in (5000.0, 6000.0, 7000.0)
        ]

        bank = RedshiftableTemplates.from_templates(
            [lam] * 3, fluxes, obs, max_z=0.5
        )
        shifted = bank.at_redshift(0.1)

        for centre, row in zip((5000.0, 6000.0, 7000.0), shifted):
            assert obs.wavelength[np.nanargmax(row)] == pytest.approx(
                centre * 1.1, rel=3e-3
            )

    def test_uncovered_wavelengths_stay_nan(self):
        obs = LogGrid.spanning(3200.0, 9500.0, velocity_to_dlnlam(400.0))
        lam = np.linspace(5000.0, 6000.0, 500)

        bank = RedshiftableTemplates.from_templates(
            [lam], [np.ones_like(lam)], obs, max_z=0.5
        )
        shifted = bank.at_redshift(0.1)

        assert np.isnan(shifted[0, 0])
        assert np.isnan(shifted[0, -1])
        assert np.isfinite(shifted[0]).any()

    def test_beyond_the_padded_redshift_is_a_clear_error(self):
        wavelengths, fluxes = self._bank()
        obs = self._observed_grid()
        bank = RedshiftableTemplates.from_templates(
            wavelengths, fluxes, obs, max_z=0.1
        )

        assert bank.max_redshift == pytest.approx(0.1, abs=0.01)
        with pytest.raises(ValueError, match="outside this bank"):
            bank.at_redshift(0.9)

    def test_negative_max_z_is_rejected(self):
        wavelengths, fluxes = self._bank()
        with pytest.raises(ValueError, match="non-negative"):
            RedshiftableTemplates.from_templates(
                wavelengths, fluxes, self._observed_grid(), max_z=-0.1
            )


class TestBuildTasks:
    """Work is split so the redshift shift is done once per group."""

    @staticmethod
    def _tasks(n_z, n_ext, n_workers):
        from superfit.SF_functions import build_tasks

        return build_tasks(
            np.linspace(0.0, 0.2, n_z), np.linspace(-2.0, 2.0, n_ext), n_workers
        )

    def test_every_grid_point_appears_exactly_once(self):
        tasks = self._tasks(5, 21, 8)

        pairs = [(z, float(e)) for z, group in tasks for e in group]
        assert len(pairs) == 5 * 21
        assert len(set(pairs)) == 5 * 21

    def test_a_single_redshift_still_fills_the_pool(self):
        """One redshift must not collapse to one task and leave workers idle."""

        tasks = self._tasks(1, 21, 8)
        assert len(tasks) > 1

    def test_many_redshifts_group_extinction_together(self):
        """With plenty of redshifts, each task should carry several A_v values."""

        tasks = self._tasks(50, 21, 8)
        assert max(len(group) for _, group in tasks) > 1

    def test_each_task_is_one_redshift(self):
        for z, group in self._tasks(5, 21, 8):
            assert np.isscalar(z) or np.ndim(z) == 0
            assert len(group) >= 1


class TestRedshiftCostIsFlat:
    """The point of the change: trying more redshifts should not cost more
    per redshift than trying one, once the bank is resampled."""

    def test_shifting_does_not_depend_on_the_redshift(self):
        import time

        obs = LogGrid.spanning(3200.0, 9500.0, velocity_to_dlnlam(400.0))
        lam = np.linspace(2500.0, 12000.0, 3000)
        fluxes = [1.0 + 0.2 * np.sin(lam / (100.0 + i)) for i in range(200)]
        bank = RedshiftableTemplates.from_templates(
            [lam] * 200, fluxes, obs, max_z=0.5
        )

        def timed(z):
            bank.at_redshift(z)
            start = time.perf_counter()
            for _ in range(20):
                bank.at_redshift(z)
            return time.perf_counter() - start

        near, far = timed(0.001), timed(0.45)

        # Same slice-and-blend either way; allow a wide margin for a shared
        # machine, since the claim is "flat", not "identical to the microsecond".
        assert 0.2 < near / far < 5.0
