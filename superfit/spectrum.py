"""An observed spectrum, however it arrived.

The fitter used to take a path and re-read it from disk in four places,
which meant a spectrum already in memory -- from a survey pipeline, a
database, a simulation -- had to be written to a temporary file first, and
the binned version was written out and read straight back in. ``Spectrum``
is the single representation everything downstream works from; a file is
just one way to build one.
"""

import os

import numpy as np

from superfit.Header_Binnings import bin_spectrum_bank, kill_header, normalise_flux
from superfit.SF_functions import mask_host_lines, remove_telluric


class Spectrum:
    """Wavelength, flux, and optionally a flux error.

    Parameters
    ----------
    wavelength : array_like
        Wavelengths, in Angstroms, strictly increasing.
    flux : array_like
        Flux in any units; the fit is scale-free and normalises internally.
    error : array_like, optional
        Per-pixel flux uncertainty. Needed only for
        ``error_spectrum="included"``; otherwise it is estimated.
    name : str, optional
        Used to name output files and to fill the SPECTRUM column of the
        results table.
    """

    def __init__(self, wavelength, flux, error=None, name=None):
        wavelength = np.asarray(wavelength, dtype=float).ravel()
        flux = np.asarray(flux, dtype=float).ravel()

        if wavelength.size != flux.size:
            raise ValueError(
                "wavelength and flux have different lengths ({} vs {})".format(
                    wavelength.size, flux.size
                )
            )
        if wavelength.size == 0:
            raise ValueError("spectrum is empty")

        if error is not None:
            error = np.asarray(error, dtype=float).ravel()
            if error.size != flux.size:
                raise ValueError(
                    "error has length {} but flux has {}".format(
                        error.size, flux.size
                    )
                )

        if not np.isfinite(wavelength).all():
            raise ValueError("wavelength contains NaN or inf")

        order = np.argsort(wavelength)
        if not np.array_equal(order, np.arange(wavelength.size)):
            # Out-of-order input is a common export quirk rather than an
            # error; sort it and carry the other columns along.
            wavelength = wavelength[order]
            flux = flux[order]
            if error is not None:
                error = error[order]

        if np.any(np.diff(wavelength) == 0):
            raise ValueError("wavelength contains duplicate values")

        if not np.isfinite(flux).any():
            raise ValueError("flux is entirely NaN")

        self.wavelength = wavelength
        self.flux = flux
        self.error = error
        self.name = name or "spectrum"

    # -- constructors -----------------------------------------------------

    @classmethod
    def from_arrays(cls, wavelength, flux, error=None, name=None):
        """Explicit alias for the constructor, for symmetry with from_file."""

        return cls(wavelength, flux, error=error, name=name)

    @classmethod
    def from_file(cls, path, name=None):
        """Read a two- or three-column ascii spectrum, or a .csv.

        Header lines are stripped. A ``.csv`` is expected to carry
        ``wavelength``, ``flux`` and ``fluxerr`` columns.
        """

        data = kill_header(path)

        error = data[:, 2] if data.shape[1] > 2 else None
        stem = os.path.basename(path)
        stem = stem[: stem.rfind(".")] if "." in stem else stem

        return cls(data[:, 0], data[:, 1], error=error, name=name or stem)

    @classmethod
    def coerce(cls, value, name=None):
        """Accept a Spectrum, a path, or an (n, 2+) array and return a Spectrum."""

        if isinstance(value, cls):
            return value if name is None else value.renamed(name)

        if isinstance(value, (str, bytes, os.PathLike)):
            return cls.from_file(os.fspath(value), name=name)

        array = np.asarray(value, dtype=float)
        if array.ndim != 2 or array.shape[1] < 2:
            raise TypeError(
                "cannot interpret {!r} as a spectrum; pass a Spectrum, a file "
                "path, or an (n, 2) or (n, 3) array".format(type(value).__name__)
            )
        return cls(
            array[:, 0],
            array[:, 1],
            error=array[:, 2] if array.shape[1] > 2 else None,
            name=name,
        )

    # -- basics -----------------------------------------------------------

    def __len__(self):
        return self.wavelength.size

    def __repr__(self):
        return "Spectrum({!r}, {} points, {:.1f}-{:.1f} A, error={})".format(
            self.name,
            len(self),
            self.wavelength[0],
            self.wavelength[-1],
            "yes" if self.error is not None else "no",
        )

    def copy(self, **replace):
        """A copy with any of wavelength/flux/error/name replaced."""

        fields = {
            "wavelength": self.wavelength,
            "flux": self.flux,
            "error": self.error,
            "name": self.name,
        }
        fields.update(replace)
        return Spectrum(**fields)

    def renamed(self, name):
        return self.copy(name=name)

    def as_array(self):
        """(n, 2) or (n, 3) array, the layout the older functions expect."""

        if self.error is None:
            return np.column_stack([self.wavelength, self.flux])
        return np.column_stack([self.wavelength, self.flux, self.error])

    @property
    def median_spacing(self):
        return float(np.median(np.diff(self.wavelength)))

    @property
    def wavelength_range(self):
        return float(self.wavelength[0]), float(self.wavelength[-1])

    # -- preprocessing ----------------------------------------------------

    def normalised(self):
        """Flux divided by its median, guarded against a median near zero."""

        return self.copy(flux=normalise_flux(self.flux))

    def without_telluric(self):
        """The telluric A band blanked to NaN."""

        masked = remove_telluric(self.as_array())
        return self.copy(
            flux=masked[:, 1],
            error=masked[:, 2] if masked.shape[1] > 2 else None,
        )

    def without_host_lines(self, z):
        """Pixels inside host galaxy emission lines at redshift ``z`` removed."""

        kept = mask_host_lines(self.as_array(), float(np.ravel(z)[0]))
        return Spectrum(
            kept[:, 0],
            kept[:, 1],
            error=kept[:, 2] if kept.shape[1] > 2 else None,
            name=self.name,
        )

    def binned(self, resolution):
        """Median-binned onto a regular grid of width ``resolution`` Angstroms.

        The flux is reproduced bit for bit by ``bin_spectrum_bank``; an error
        column, if present, is carried through on the same bin edges as the
        root-mean-square of the errors falling in each bin, and scaled by the
        same median so flux and error stay on one scale.
        """

        flux_binned = bin_spectrum_bank(self.as_array()[:, :2], resolution)

        if self.error is None:
            return Spectrum(flux_binned[:, 0], flux_binned[:, 1], name=self.name)

        error_binned = self._bin_error(resolution)
        return Spectrum(
            flux_binned[:, 0],
            flux_binned[:, 1],
            error=error_binned,
            name=self.name,
        )

    def _bin_error(self, resolution):
        """Bin the error onto exactly the bins ``binned`` uses for the flux."""

        import math

        from scipy import stats

        lam = self.wavelength
        flux = self.flux

        # Mirror bin_spectrum_bank so the two cannot drift apart.
        n_bins = math.floor((lam[-1] - lam[0]) / resolution)
        bin_range = (lam.min(), lam.max())

        flux_bin, _, _ = stats.binned_statistic(
            lam, flux, statistic="median", range=bin_range, bins=n_bins
        )
        mean_square, _, _ = stats.binned_statistic(
            lam, self.error**2, statistic="mean", range=bin_range, bins=n_bins
        )

        keep = ~np.isnan(flux_bin)
        median_flux = np.nanmedian(flux_bin[keep])

        return np.sqrt(mean_square[keep]) / median_flux

    def interpolated_onto(self, grid):
        """Flux resampled onto ``grid``, NaN outside the observed range."""

        return np.interp(
            np.asarray(grid, dtype=float),
            self.wavelength,
            self.flux,
            left=np.nan,
            right=np.nan,
        )

    # -- checks -----------------------------------------------------------

    def check_long_enough(self, resolution, minimum_bins=35):
        """Raise if binning would leave too few points to fit meaningfully."""

        span = self.wavelength[-1] - self.wavelength[0]
        bin_width = 30 if self.median_spacing > 10 else 10
        n_bins = span / bin_width

        if n_bins < minimum_bins:
            raise ValueError(
                "Spectrum {!r} spans {:.0f} A, about {:.0f} bins of {} A. "
                "At least {} are needed to fit.".format(
                    self.name, span, n_bins, bin_width, minimum_bins
                )
            )
        return self
