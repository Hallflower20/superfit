"""Wavelength grids uniform in ln(lambda), and redshifting by shifting.

On a grid uniform in wavelength, redshifting a template rescales its
wavelength axis, so every template has to be interpolated onto the observed
grid again at every trial redshift -- roughly a thousand ``np.interp`` calls
per grid point.

On a grid uniform in ``ln(lambda)`` a redshift is a *translation*::

    ln(lambda_obs) = ln(lambda_rest) + ln(1 + z)

so it is the same displacement for every template and every wavelength.
Resample the bank onto a common rest-frame log grid once, and each trial
redshift becomes two array slices and a linear blend over the whole bank at
once, with no Python-level loop.

The grid is uniform in velocity rather than in Angstroms: a step of
``dlnlam`` corresponds to ``c * dlnlam`` km/s everywhere.
"""

import numpy as np

# Speed of light, km/s.
C_KM_S = 299792.458


def velocity_to_dlnlam(velocity_kms):
    """Convert a velocity resolution in km/s to a step in ln(lambda)."""

    if velocity_kms <= 0:
        raise ValueError(
            "velocity resolution must be positive, got {!r}".format(velocity_kms)
        )
    return float(velocity_kms) / C_KM_S


def dlnlam_to_velocity(dlnlam):
    """Convert a step in ln(lambda) to a velocity resolution in km/s."""

    return float(dlnlam) * C_KM_S


def matched_velocity_resolution(resolution_angstrom, lower, upper):
    """The velocity resolution giving as many bins as a linear grid would.

    A linear grid of ``resolution`` Angstroms across ``[lower, upper]`` holds
    ``(upper - lower) / resolution`` points; a log grid holds
    ``ln(upper / lower) / dlnlam``. Equating them keeps the sampling
    comparable, so an existing configuration written in Angstroms carries
    over to a log grid of the same size rather than silently changing how
    finely the fit is sampled.
    """

    if lower <= 0:
        raise ValueError("wavelength lower bound must be positive")
    if upper <= lower:
        raise ValueError("wavelength upper bound must exceed the lower bound")
    if resolution_angstrom <= 0:
        raise ValueError("resolution must be positive")

    dlnlam = np.log(upper / lower) * resolution_angstrom / (upper - lower)
    return dlnlam_to_velocity(dlnlam)


class LogGrid:
    """Wavelengths spaced uniformly in ln(lambda)."""

    def __init__(self, ln_start, dlnlam, n_points):
        if dlnlam <= 0:
            raise ValueError("dlnlam must be positive, got {!r}".format(dlnlam))
        if n_points < 2:
            raise ValueError("a grid needs at least two points")

        self.ln_start = float(ln_start)
        self.dlnlam = float(dlnlam)
        self.n_points = int(n_points)

    @classmethod
    def spanning(cls, lower, upper, dlnlam):
        """The smallest grid of step ``dlnlam`` covering ``[lower, upper]``."""

        if lower <= 0:
            raise ValueError("wavelength lower bound must be positive")
        if upper <= lower:
            raise ValueError("wavelength upper bound must exceed the lower bound")

        n = int(np.ceil(np.log(upper / lower) / dlnlam)) + 1
        return cls(np.log(lower), dlnlam, n)

    @property
    def ln_wavelength(self):
        return self.ln_start + self.dlnlam * np.arange(self.n_points)

    @property
    def wavelength(self):
        return np.exp(self.ln_wavelength)

    @property
    def velocity_resolution(self):
        return dlnlam_to_velocity(self.dlnlam)

    def __len__(self):
        return self.n_points

    def __repr__(self):
        lam = self.wavelength
        return "LogGrid({:.1f}-{:.1f} A, {} points, {:.0f} km/s)".format(
            lam[0], lam[-1], self.n_points, self.velocity_resolution
        )

    def shift_bins(self, z):
        """How many bins a redshift of ``z`` displaces a spectrum by."""

        return np.log1p(z) / self.dlnlam

    def extended_blueward(self, n_extra):
        """The same grid with ``n_extra`` bins added at the blue end."""

        return LogGrid(
            self.ln_start - n_extra * self.dlnlam,
            self.dlnlam,
            self.n_points + n_extra,
        )


class RedshiftableTemplates:
    """A bank of templates on a rest-frame log grid, shiftable in redshift.

    ``values`` is (n_templates, n_rest) sampled on ``rest_grid``, which must
    share its step with ``observed_grid`` and extend far enough blueward to
    cover the largest redshift asked for.
    """

    def __init__(self, values, rest_grid, observed_grid, offset):
        self.values = values
        self.rest_grid = rest_grid
        self.observed_grid = observed_grid
        # Index into the rest grid whose wavelength equals observed_grid[0].
        self.offset = int(offset)

    @classmethod
    def from_templates(cls, wavelengths, fluxes, observed_grid, max_z):
        """Resample templates onto a rest grid aligned with ``observed_grid``.

        Parameters
        ----------
        wavelengths, fluxes : sequences of 1-D arrays
            One pair per template; they need not share a sampling.
        observed_grid : LogGrid
            The grid the fit is evaluated on.
        max_z : float
            The largest redshift that will be requested. The rest grid is
            padded blueward by enough bins to cover it.
        """

        if max_z < 0:
            raise ValueError("max_z must be non-negative, got {!r}".format(max_z))

        # One extra bin so the linear blend can always reach its left neighbour.
        pad = int(np.ceil(observed_grid.shift_bins(max_z))) + 1
        rest_grid = observed_grid.extended_blueward(pad)
        rest_lam = rest_grid.wavelength

        values = np.empty((len(fluxes), len(rest_grid)), dtype=np.float64)
        for i, (lam, flux) in enumerate(zip(wavelengths, fluxes)):
            values[i] = np.interp(rest_lam, lam, flux, left=np.nan, right=np.nan)

        return cls(values, rest_grid, observed_grid, offset=pad)

    def __len__(self):
        return self.values.shape[0]

    def at_redshift(self, z, rows=None):
        """Templates redshifted to ``z``, sampled on the observed grid.

        Returns (n_templates, len(observed_grid)). The flux is *not* divided
        by (1 + z) and carries no extinction; both are cheap multiplies the
        caller applies, and keeping them out means this result can be reused
        across a whole extinction grid.

        ``rows`` is an optional ``slice`` selecting which templates to shift.
        A worker handling one block of a 15561-template bank has no use for
        the other 14500, and shifting them anyway is the largest single piece
        of duplicated memory traffic in a scan.
        """

        shift = self.observed_grid.shift_bins(z)
        whole = int(np.floor(shift))
        frac = shift - whole

        # Wanted rest index for observed pixel i is offset + i - shift, which
        # sits `frac` below (offset - whole + i). Blend that node with its
        # left neighbour.
        base = self.offset - whole
        n = len(self.observed_grid)

        if base - 1 < 0 or base + n > self.values.shape[1]:
            raise ValueError(
                "redshift {} needs rest-frame coverage outside this bank; it "
                "was built for z up to about {:.3f}".format(z, self.max_redshift)
            )

        block = self.values if rows is None else self.values[rows]

        left = block[:, base - 1 : base - 1 + n]
        right = block[:, base : base + n]

        return frac * left + (1.0 - frac) * right

    @property
    def max_redshift(self):
        """The largest redshift this bank was padded for."""

        return float(np.expm1((self.offset - 1) * self.rest_grid.dlnlam))
