"""Estimating a per-pixel uncertainty for a spectrum that did not ship one.

Both estimators here have to survive real data: cosmic-ray hits, NaN gaps
where a detector column was masked, and stretches that were interpolated over
and so carry no scatter at all. The chi2 that consumes these weights each
pixel by 1 / sigma**2, so an underestimated sigma is not a small error -- one
pixel with sigma near zero silently decides the whole fit.
"""

import numpy as np
import scipy.signal as mf
from scipy.ndimage import median_filter, uniform_filter1d

# sigma = 1.4826 * MAD for Gaussian noise.
MAD_TO_SIGMA = 1.4826

# No pixel may claim an uncertainty smaller than this fraction of the typical
# one. Interpolated-flat regions have literally zero scatter; without a floor
# they get essentially infinite weight in the chi2.
MIN_ERROR_FRACTION = 1e-2


def robust_scale(values):
    """A median-based amplitude for ``values``, ignoring NaN.

    Used instead of the mean so that one cosmic ray cannot rescale the whole
    error spectrum.
    """

    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValueError("spectrum has no finite flux values")

    scale = np.median(np.abs(finite))
    if scale > 0:
        return scale

    # Flux straddling zero: fall back to its scatter.
    scale = MAD_TO_SIGMA * np.median(np.abs(finite - np.median(finite)))
    if scale > 0:
        return scale

    raise ValueError("spectrum flux is constant; cannot estimate an error")


def apply_error_floor(error, fraction=MIN_ERROR_FRACTION):
    """Raise implausibly small uncertainties to a fraction of the typical one."""

    error = np.asarray(error, dtype=float)
    positive = error[np.isfinite(error) & (error > 0)]

    if positive.size == 0:
        # Nothing to calibrate against; a flat unit error at least keeps the
        # chi2 finite and equally weighted.
        return np.ones_like(error)

    floor = fraction * np.median(positive)
    out = np.where(np.isfinite(error), error, np.nan)
    out = np.fmax(out, floor)
    return np.where(np.isfinite(out), out, floor)


def _interpolate_gaps(y):
    """Linearly fill NaN so a filter can run; returns (filled, was_valid)."""

    valid = np.isfinite(y)
    if valid.all():
        return y.astype(float), valid
    if not valid.any():
        raise ValueError("spectrum has no finite flux values")

    idx = np.arange(y.size)
    filled = y.astype(float).copy()
    filled[~valid] = np.interp(idx[~valid], idx[valid], y[valid])
    return filled, valid


def linear_error(spec_object, num=10):

    """
    Estimate sigma as the scatter of the flux about a straight line fitted
    within each consecutive block of ``num`` pixels.

    parameters
    ----------
    spec_object: (n, >=2) array of wavelength and flux
    num: block length in pixels

    returns
    -------
    (n, 2) array of wavelength and error.
    """

    lam = np.asarray(spec_object[:, 0], dtype=float)
    flux = np.asarray(spec_object[:, 1], dtype=float)

    n_full = (flux.size // num) * num
    remainder = flux.size - n_full

    if n_full == 0:
        raise ValueError(
            "spectrum has {} pixels, fewer than one {}-pixel block".format(
                flux.size, num
            )
        )

    lam_blocks = lam[:n_full].reshape(-1, num)
    flux_blocks = flux[:n_full].reshape(-1, num)

    sigma = np.full(lam_blocks.shape[0], np.nan)
    for i, (block_lam, block_flux) in enumerate(zip(lam_blocks, flux_blocks)):
        ok = np.isfinite(block_flux)
        if ok.sum() < 3:
            continue
        # Each block is measured against its OWN fit. The previous version
        # assigned the residual inside this loop and let it be overwritten, so
        # every block ended up scored against the last block's line.
        slope, intercept = np.polyfit(block_lam[ok], block_flux[ok], 1)
        residual = block_flux[ok] - (slope * block_lam[ok] + intercept)
        sigma[i] = np.std(residual, ddof=1)

    if np.isnan(sigma).all():
        raise ValueError("no block had enough finite pixels to estimate an error")

    # Blocks that were entirely masked inherit the typical scatter.
    sigma[np.isnan(sigma)] = np.nanmedian(sigma)

    error = np.repeat(sigma, num)
    if remainder:
        error = np.concatenate([error, np.full(remainder, error[-1])])

    return np.column_stack([lam, apply_error_floor(error)])


def savitzky_golay(spec, window=31, polyorder=3, scatter_window=100):

    """
    Estimate sigma as the local scatter of the flux about a Savitzky-Golay
    smoothing of itself.

    parameters
    ----------
    spec: (n, >=2) array of wavelength and flux
    window, polyorder: Savitzky-Golay smoothing parameters
    scatter_window: width in pixels over which the local scatter is measured

    returns
    -------
    (n, 2) array of wavelength and error.
    """

    x = np.asarray(spec[:, 0], dtype=float)
    y_raw = np.asarray(spec[:, 1], dtype=float)

    if y_raw.size < window:
        raise ValueError(
            "spectrum has {} pixels, fewer than the {}-pixel smoothing window".format(
                y_raw.size, window
            )
        )

    # Normalise by a median rather than a mean: the mean moves with a single
    # cosmic ray, which rescaled the error at every pixel, and it returns NaN
    # if any pixel is NaN, which turned the entire error spectrum into NaN.
    y = y_raw / robust_scale(y_raw)

    # savgol_filter cannot see through NaN, so bridge the gaps for smoothing
    # only and drop those pixels again when measuring the scatter.
    filled, valid = _interpolate_gaps(y)

    smooth = mf.savgol_filter(filled, window, polyorder, mode="nearest")
    residual = np.where(valid, y - smooth, np.nan)

    # A centred median of |residual| rather than a forward-looking mean of
    # residual**2: centred so the error is not offset by half a window, and a
    # median so a cosmic ray inflates its own neighbourhood rather than the
    # whole 100-pixel window it happens to fall in.
    abs_residual = np.where(valid, np.abs(residual), np.nan)
    filled_residual, _ = _interpolate_gaps(abs_residual)
    err_std = MAD_TO_SIGMA * median_filter(
        filled_residual, size=scatter_window, mode="nearest"
    )

    # A median of |residual| is zero across any interpolated-flat stretch, so
    # blend in the local mean to keep those regions from hitting exactly zero
    # before the floor is applied.
    mean_abs = uniform_filter1d(filled_residual, size=scatter_window, mode="nearest")
    err_std = np.where(err_std > 0, err_std, MAD_TO_SIGMA * mean_abs)

    return np.column_stack([x, apply_error_floor(err_std)])
