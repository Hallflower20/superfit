"""Reading an observed spectrum off disk, in whatever form it arrived.

The reader used to be twenty lines that stripped comment characters and took
the first two whitespace-separated columns. That covers one export format.
Astronomers are handed FITS by every pipeline they use, wavelengths in nm or
microns by half of them, and inverse variance instead of an uncertainty by
the rest -- and a spectrum silently read a factor of ten off in wavelength
is a misclassification, not an error message.

This module handles:

* ascii, with the column layout given by position or by the names in a
  commented header line;
* csv, by column name;
* FITS binary tables, by column name, with names guessed from the usual
  spellings and units taken from TUNIT keywords;
* FITS images with a linear or log-linear wavelength WCS, which is how
  IRAF-descended pipelines write a one-dimensional spectrum;
* inverse variance, converted to an uncertainty;
* wavelengths in Angstroms, nm, microns, or log10(Angstrom).

Nothing here guesses silently where guessing could be wrong: a file whose
columns cannot be identified raises, and names the columns it did find.
"""

import os
import re

import numpy as np

# What a column name in a commented header line may look like.
_IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_.+-]*$")

# Multiply a wavelength in the named unit by this to get Angstroms.
WAVELENGTH_UNITS = {
    "": 1.0,
    "a": 1.0,
    "aa": 1.0,
    "ang": 1.0,
    "angstrom": 1.0,
    "angstroms": 1.0,
    "nm": 10.0,
    "nanometer": 10.0,
    "nanometers": 10.0,
    "um": 1.0e4,
    "µm": 1.0e4,
    "micron": 1.0e4,
    "microns": 1.0e4,
    "micrometer": 1.0e4,
    "mm": 1.0e7,
    "cm": 1.0e8,
    "m": 1.0e10,
}

# Units that mean "the stored number is log10 of a wavelength in Angstroms".
LOG_WAVELENGTH_UNITS = frozenset(["log10(aa)", "log10(angstrom)", "loglam", "log(aa)"])

# Column names seen in the wild, in the order they are preferred.
WAVELENGTH_NAMES = (
    "wavelength", "wave", "lam", "lambda", "wl", "awav", "spectral_axis", "loglam",
)
FLUX_NAMES = ("flux", "f_lambda", "flam", "fl", "spec", "intensity", "counts")
ERROR_NAMES = (
    "fluxerr", "flux_err", "error", "err", "sigma", "e_flux", "uncertainty",
    "noise", "stdev", "std",
)
IVAR_NAMES = ("ivar", "invvar", "inverse_variance", "flux_ivar")

FITS_SUFFIXES = (".fits", ".fit", ".fts", ".fz")


class SpectrumReadError(ValueError):
    """Raised when a file cannot be read as a spectrum."""


def wavelength_scale(unit):
    """Factor converting ``unit`` to Angstroms, or 'log' for a log axis."""

    if unit is None:
        return 1.0

    key = str(unit).strip().lower()
    if key in LOG_WAVELENGTH_UNITS:
        return "log"
    if key in WAVELENGTH_UNITS:
        return WAVELENGTH_UNITS[key]

    raise SpectrumReadError(
        "Unknown wavelength unit {!r}. Known units: {}.".format(
            unit, ", ".join(sorted(set(WAVELENGTH_UNITS) - {""}) + ["log10(AA)"])
        )
    )


def to_angstrom(wavelength, unit):
    """Convert a wavelength array to Angstroms."""

    scale = wavelength_scale(unit)
    wavelength = np.asarray(wavelength, dtype=float)
    if scale == "log":
        return 10.0**wavelength
    return wavelength * scale


def normalise_columns(columns):
    """Accept a sequence or a mapping and return a (wavelength, flux, error) dict."""

    if columns is None:
        return {}

    if isinstance(columns, dict):
        wanted = {}
        for key, value in columns.items():
            key = str(key).lower()
            if key in ("wavelength", "wave", "lam", "w"):
                wanted["wavelength"] = value
            elif key in ("flux", "f"):
                wanted["flux"] = value
            elif key in ("error", "err", "sigma", "e"):
                wanted["error"] = value
            elif key in ("ivar", "inverse_variance"):
                wanted["ivar"] = value
            else:
                raise SpectrumReadError(
                    "Unknown column role {!r}; expected wavelength, flux, "
                    "error or ivar.".format(key)
                )
        return wanted

    columns = list(columns)
    if not 2 <= len(columns) <= 3:
        raise SpectrumReadError(
            "Give two or three columns (wavelength, flux, and optionally "
            "error), got {}.".format(len(columns))
        )

    wanted = {"wavelength": columns[0], "flux": columns[1]}
    if len(columns) == 3:
        wanted["error"] = columns[2]
    return wanted


def _error_from_ivar(ivar):
    """Uncertainty from inverse variance, with non-positive entries dropped."""

    ivar = np.asarray(ivar, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        error = 1.0 / np.sqrt(ivar)
    return np.where(ivar > 0, error, np.nan)


# -- ascii -----------------------------------------------------------------


def _header_names(path, width):
    """Column names from a commented header line, if the file has one.

    The convention is a comment line whose tokens name the columns, as in
    ``# wavelength flux fluxerr``. Files carry prose comments too, and
    ``# instrument: a telescope`` is three tokens like the header above it,
    so a line counts only if every token also looks like a column name. Of
    those, the last is taken -- headers sit closest to the data. Consulted
    only when the caller asks for columns by name.
    """

    names = None
    with open(path, "r") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped[0] not in "#%@":
                break
            tokens = stripped.lstrip("#%@ ").split()
            if len(tokens) == width and all(_IDENTIFIER.match(t) for t in tokens):
                names = tokens
    return names


def _read_ascii(path, wanted):
    """Whitespace-separated columns, with comment and blank lines skipped."""

    rows = []
    with open(path, "r") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped[0] in "#%@" or stripped[0].isalpha():
                continue
            rows.append(stripped.split())

    if not rows:
        raise SpectrumReadError("{} has no data rows.".format(path))

    width = min(len(row) for row in rows)
    if width < 2:
        raise SpectrumReadError(
            "{} has only {} column(s); a spectrum needs wavelength and "
            "flux.".format(path, width)
        )

    try:
        data = np.array([[float(value) for value in row[:width]] for row in rows])
    except ValueError as exc:
        raise SpectrumReadError("{} has a non-numeric entry: {}".format(path, exc))

    names = None

    def index_of(role, default):
        nonlocal names
        requested = wanted.get(role)
        if requested is None:
            return default
        if isinstance(requested, (int, np.integer)):
            return int(requested)
        if names is None:
            names = [n.lower() for n in (_header_names(path, width) or [])]
        try:
            return names.index(str(requested).lower())
        except ValueError:
            raise SpectrumReadError(
                "{} has no column called {!r}. Its header names: {}. Columns "
                "can also be given by position, as 0, 1, 2.".format(
                    path, requested, ", ".join(names) or "none"
                )
            )

    lam_index = index_of("wavelength", 0)
    flux_index = index_of("flux", 1)
    error_index = index_of("error", 2 if width > 2 else None)
    ivar_index = index_of("ivar", None)

    for role, index in (("wavelength", lam_index), ("flux", flux_index)):
        if index >= width:
            raise SpectrumReadError(
                "{} asks for {} in column {}, but the file has {} "
                "columns.".format(path, role, index, width)
            )

    error = None
    if ivar_index is not None and ivar_index < width:
        error = _error_from_ivar(data[:, ivar_index])
    elif error_index is not None and error_index < width:
        error = data[:, error_index]

    return data[:, lam_index], data[:, flux_index], error, None


# -- csv -------------------------------------------------------------------


def _read_table(table, wanted, units, source):
    """Pull wavelength, flux and error out of a table with named columns."""

    lookup = {str(name).lower(): str(name) for name in table.colnames}

    def column(role, candidates, required):
        requested = wanted.get(role)
        if requested is not None:
            if isinstance(requested, (int, np.integer)):
                return table.colnames[int(requested)]
            key = str(requested).lower()
            if key not in lookup:
                raise SpectrumReadError(
                    "{} has no column called {!r}. It has: {}.".format(
                        source, requested, ", ".join(table.colnames)
                    )
                )
            return lookup[key]

        for candidate in candidates:
            if candidate in lookup:
                return lookup[candidate]

        if required:
            raise SpectrumReadError(
                "Could not find the {} column in {}. It has: {}. Name it "
                "explicitly with the column settings.".format(
                    role, source, ", ".join(table.colnames)
                )
            )
        return None

    lam_name = column("wavelength", WAVELENGTH_NAMES, True)
    flux_name = column("flux", FLUX_NAMES, True)
    ivar_name = column("ivar", IVAR_NAMES, False)
    error_name = None if ivar_name else column("error", ERROR_NAMES, False)

    wavelength = np.asarray(table[lam_name], dtype=float).ravel()
    flux = np.asarray(table[flux_name], dtype=float).ravel()

    error = None
    if ivar_name:
        error = _error_from_ivar(np.asarray(table[ivar_name], dtype=float).ravel())
    elif error_name:
        error = np.asarray(table[error_name], dtype=float).ravel()

    # A column named `loglam` carries log10(Angstrom) by convention, and
    # rarely says so in a unit keyword.
    unit = units.get(lam_name)
    if unit is None and lam_name.lower() == "loglam":
        unit = "log10(AA)"

    return wavelength, flux, error, unit


def _read_csv(path, wanted):
    from astropy.table import Table

    try:
        table = Table.read(path, format="csv")
    except Exception as exc:
        raise SpectrumReadError("Could not read {} as csv: {}".format(path, exc))

    units = {name: table[name].unit for name in table.colnames}
    return _read_table(table, wanted, units, path)


# -- FITS ------------------------------------------------------------------


def _wcs_wavelength(header, n_pixels):
    """The wavelength axis of a 1D image HDU, from its WCS keywords.

    Handles the linear case and the log-linear one that IRAF-descended
    pipelines write as ``CTYPE1 = 'AWAV-LOG'`` or ``DC-FLAG = 1``.
    """

    crval = header.get("CRVAL1")
    if crval is None:
        return None, None

    cdelt = header.get("CDELT1", header.get("CD1_1", header.get("CDELT", 1.0)))
    crpix = header.get("CRPIX1", 1.0)

    # FITS pixels are 1-based, and CRPIX names the pixel where CRVAL holds.
    axis = crval + (np.arange(n_pixels, dtype=float) + 1.0 - crpix) * float(cdelt)

    ctype = str(header.get("CTYPE1", "")).upper()
    if header.get("DC-FLAG", 0) or "LOG" in ctype:
        return 10.0**axis, "AA"

    return axis, header.get("CUNIT1")


def _read_fits(path, wanted, hdu_index):
    from astropy.io import fits

    with fits.open(path, memmap=False) as hdus:
        candidates = (
            [hdus[hdu_index]] if hdu_index is not None else list(hdus)
        )

        errors = []
        for hdu in candidates:
            if getattr(hdu, "data", None) is None:
                continue

            if hasattr(hdu, "columns"):
                from astropy.table import Table

                try:
                    table = Table(hdu.data)
                except Exception as exc:  # pragma: no cover - malformed table
                    errors.append(str(exc))
                    continue

                units = {}
                for i, name in enumerate(table.colnames, start=1):
                    units[name] = hdu.header.get("TUNIT{}".format(i))
                try:
                    return _read_table(table, wanted, units, path)
                except SpectrumReadError as exc:
                    errors.append(str(exc))
                    continue

            data = np.asarray(hdu.data)
            if data.ndim == 0:
                continue
            # Multispec files stack (flux, sky, error, ...) or several orders
            # along the first axis; the first row is the flux.
            flux = data if data.ndim == 1 else data.reshape(-1, data.shape[-1])[0]

            wavelength, unit = _wcs_wavelength(hdu.header, flux.size)
            if wavelength is None:
                errors.append(
                    "image HDU {!r} has no CRVAL1 to build a wavelength axis "
                    "from".format(hdu.name)
                )
                continue

            return wavelength, np.asarray(flux, dtype=float), None, unit

        raise SpectrumReadError(
            "No HDU in {} could be read as a spectrum.{}".format(
                path,
                "\n  " + "\n  ".join(errors) if errors else
                " It has no table or image data.",
            )
        )


# -- the entry point -------------------------------------------------------


def read_spectrum(path, columns=None, wavelength_unit=None, hdu=None):
    """Read a spectrum from ``path``.

    Parameters
    ----------
    path : str
        An ascii, csv or FITS file.
    columns : sequence or mapping, optional
        Which columns to use. A sequence names wavelength, flux and
        optionally error, either by name or by position. A mapping can also
        supply ``ivar`` instead of ``error``. Left out, the columns are
        identified by name for csv and FITS, and by position for ascii.
    wavelength_unit : str, optional
        The unit the wavelengths are in -- ``nm``, ``micron``,
        ``log10(AA)``, and so on. Left out, it is taken from the file's own
        unit keywords, and Angstroms assumed if it says nothing.
    hdu : int or str, optional
        Which FITS HDU to read. Left out, each is tried in turn.

    Returns
    -------
    (wavelength, flux, error)
        Wavelength in Angstroms; ``error`` is None when the file has none.
    """

    path = os.fspath(path)
    wanted = normalise_columns(columns)

    lowered = path.lower()
    if lowered.endswith(".gz"):
        lowered = lowered[:-3]

    if lowered.endswith(FITS_SUFFIXES):
        wavelength, flux, error, unit = _read_fits(path, wanted, hdu)
    elif lowered.endswith(".csv"):
        wavelength, flux, error, unit = _read_csv(path, wanted)
    else:
        wavelength, flux, error, unit = _read_ascii(path, wanted)

    # An explicit unit is the caller correcting the file, so it wins.
    wavelength = to_angstrom(wavelength, wavelength_unit or unit)

    return wavelength, np.asarray(flux, dtype=float), error
