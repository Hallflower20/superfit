import math
import warnings

import numpy as np
from astropy import table
from scipy import stats
from PyAstronomy import pyasl
from astropy.table import Table


def kill_header(file_name):

    """
    This function removes all entries beginning with '#' from a file with a header, keeping only the column

    data and saving it into a file.


    parameters
    ----------

    It takes one path (in the form of "/home/user/Dropbox/something") to pull and eliminate its header


    returns
    -------

    File without header.

    """

    if("csv" in file_name):
        file = Table.read(file_name, format='csv')

        lam_floats = file["wavelength"].data
        flux_floats = file["flux"].data
        fluxerr_floats = file["fluxerr"].data

        spectrum = np.array([lam_floats, flux_floats, fluxerr_floats]).T

        return spectrum

    lines = []

    file = open(file_name, "r")

    lines = file.readlines()

    lines = [i for i in lines if i]

    lines = [
        i
        for i in lines
        if not i[0].isalpha() and i[0] != "#" and i[0] != "%" and i[0] != "@"
    ]

    lines = [i for i in lines if i[0] != "\n"]

    lines = [s.strip("\n") for s in lines]  # remove empty lines

    lines = [s.replace("\n", "") for s in lines]  # replace with nothing


    columns = []

    for line in lines:
        ii = line.split()
        columns.append(ii)

    columns = np.array(columns)

    lam_floats = [float(i) for i in columns[:, 0]]
    flux_floats = [float(i) for i in columns[:, 1]]

    spectrum = np.array([lam_floats, flux_floats]).T

    return spectrum


def normalise_flux(flux):
    """Divide flux by its median, guarding the case where that median is ~0.

    A continuum-subtracted or sky-dominated spectrum can have a median
    indistinguishable from zero, and dividing by it amplified the spectrum by
    ~1e15. Fall back to the scatter when that happens.
    """

    flux = np.asarray(flux, dtype=float)
    finite = flux[np.isfinite(flux)]

    if finite.size == 0:
        raise ValueError("spectrum has no finite flux values")

    median = np.median(finite)
    mad = np.median(np.abs(finite - median))

    if np.abs(median) > 1e-3 * mad and median != 0:
        return flux / median

    warnings.warn(
        "median flux ({:.3g}) is negligible next to its scatter ({:.3g}); "
        "normalising by the scatter instead of the median".format(median, mad),
        RuntimeWarning,
        stacklevel=2,
    )

    if mad > 0:
        return flux / (1.4826 * mad)

    raise ValueError("spectrum flux is constant; cannot normalise it")


def median_spacing(lam):
    """Typical wavelength step, robust to a few irregular pixels."""

    lam = np.asarray(lam, dtype=float)
    if lam.size < 2:
        raise ValueError("spectrum needs at least two wavelength samples")
    return float(np.median(np.diff(lam)))


def bin_spectrum(spectrum, resolution):

    """

    Returns a median normalized flux, binned in a resolution given by the user.
    Parameters:
    -----------
    spectrum ‘array’: array of arrays containing the spectrum (with or without flux error).
    First array must be the wavelength, second the flux, (third the error on the flux).
    resolution ’int: the desired resolution, must match the units of wavelength in the spectrum file

    """

    lam = spectrum[:, 0]
    flux = spectrum[:, 1]

    if len(spectrum[0]) > 2:
        fluxerror = spectrum[:, 2]
    else:
        fluxerror = None

    # If the data are already coarser than the requested bin, rebinning would
    # invent structure that was never measured; hand them back as they are.
    # The old test read lam[15] - lam[16], which is negative for any ascending
    # wavelength axis, so this branch never ran -- and returned None when it
    # notionally did, because it had no return statement.
    if median_spacing(lam) >= resolution:
        passthrough = table.Table()
        passthrough["lam_bin"] = lam
        passthrough["bin_flux"] = normalise_flux(flux)
        if fluxerror is not None:
            passthrough["bin_fluxerror"] = fluxerror / np.nanmedian(flux)
        return passthrough
    else:
        number_of_bins = math.floor((lam[-1] - lam[0]) / resolution)
        flux_bin, bin_edge, index = stats.binned_statistic(
            lam,
            flux,
            statistic="median",
            range=(lam.min(), lam.max()),
            bins=number_of_bins,
        )
        bin_wavelength = [
            (bin_edge[i] + bin_edge[i + 1]) / 2 for i in range(len(bin_edge) - 1)
        ]

        # This is the condition I had to add to get rid of the NaNs, but I still don’t know why flux_bin has NaNs in the first place

        if fluxerror is not None:
            fluxerror_bin = []
            for i in range(len(bin_edge) - 1):
                error_squared = []
                for index, la in enumerate(lam):
                    if la >= bin_edge[i] and la < bin_edge[i + 1]:
                        error_squared.append(fluxerror[index] ** 2)
                    error = np.sqrt(np.sum(error_squared) / len(error_squared))
                fluxerror_bin.append(error)

        bin_wavelength = np.array(bin_wavelength)
        flux_bin = np.array(flux_bin)

        mask = [not np.isnan(x) for x in flux_bin]

        bin_wavelength = bin_wavelength[mask]
        flux_bin = flux_bin[mask]
        bin_spectra = table.Table()
        flux_bin = np.array(flux_bin)
        median_flux = np.nanmedian(flux_bin)
        flux_bin = flux_bin / median_flux
        bin_spectra["lam_bin"] = bin_wavelength
        bin_spectra["bin_flux"] = flux_bin
        if fluxerror is not None:
            fluxerror_bin = np.array(fluxerror_bin)[mask]
            fluxerror_bin = np.array(fluxerror_bin)
            fluxerror_bin = fluxerror_bin / median_flux
            bin_spectra["bin_fluxerror"] = fluxerror_bin
        bin_spectra = bin_spectra[bin_spectra["bin_flux"] != np.nan]

        return bin_spectra


def kill_header_and_bin(original, resolution=10, **kwargs):

    """

    Takes a path (in the form of "/home/user/Dropbox/something"), pulls a spectrum file that should consist of 2 or 3 columns (wavelength, flux and error)

    eliminates the header, bins the spectrum to a specific resolution and saves the file in the same directory from the original path, it adds "_binned.ascii"

    to the end of the file name


    ----------

    Outputs: astropy table with binned data, and saves file in the same path as the original with "_binned.ascii" in the name



    """

    saving_path = kwargs["save_bin"]

    noheader = kill_header(original)

    lam = [float(item) for item in noheader[:, 0]]
    flux = [float(item) for item in noheader[:, 1]]

    if len(noheader[0]) > 2:
        fluxerror = noheader[:, 2]
        spectrum = np.array([lam, flux, fluxerror]).T
    else:
        spectrum = np.array([lam, flux]).T

    bin_spec = bin_spectrum(spectrum, resolution)

    if np.min(np.diff(spectrum[:, 0])) > resolution:
        raise Exception(
            "The resolution you chose ({0} Ang) is less than a single bin ({1: .2f} ang). Decrease the resolution for this spectrum and try again".format(
                resolution, np.min(np.diff(spectrum[:, 0]))
            )
        )

    np.savetxt(saving_path, bin_spec, fmt="%s")

    return bin_spec, saving_path


def bin_spectrum_bank(spectrum, resolution):

    """

    Returns a median normalized flux, binned in a resolution given by the user. Modified to bin only
    the template bank since no error is involved.


    Parameters:
    -----------
    spectrum ‘array’: array of arrays containing the spectrum (with or without flux error).
    First array must be the wavelength, second the flux, (third the error on the flux).
    resolution ’int: the desired resolution, must match the units of wavelength in the spectrum file

    """
    # spectrum = np.loadtxt(spectrum)
    lam = spectrum[:, 0]
    flux = spectrum[:, 1]

    # See bin_spectrum: this guard used to be unreachable and indexed lam[15]
    # unconditionally, which raised IndexError on any spectrum under 17 pixels.
    if median_spacing(lam) >= resolution:
        return np.column_stack([lam, normalise_flux(flux)])

    else:
        number_of_bins = math.floor((lam[-1] - lam[0]) / resolution)
        flux_bin, bin_edge, index = stats.binned_statistic(
            lam,
            flux,
            statistic="median",
            range=(lam.min(), lam.max()),
            bins=number_of_bins,
        )
        bin_wavelength = [
            (bin_edge[i] + bin_edge[i + 1]) / 2 for i in range(len(bin_edge) - 1)
        ]

        bin_wavelength = np.array(bin_wavelength)
        flux_bin = np.array(flux_bin)

        mask = [not (np.isnan(x)) for x in flux_bin]

        bin_wavelength = bin_wavelength[mask]
        flux_bin = flux_bin[mask]

        bin_spectra = table.Table()
        flux_bin = np.array(flux_bin)
        median_flux = np.nanmedian(flux_bin)
        flux_bin = flux_bin / median_flux
        # bin_spectra['lam_bin'] = bin_wavelength
        # bin_spectra['bin_flux'] = flux_bin

        bin_spectra = np.array([bin_wavelength, flux_bin]).T

        return bin_spectra


# Nebular host-galaxy emission lines, rest frame, in air.
HOST_LINES = np.array(
    [
        6564.61,
        4862.69,
        3726.09,
        3729.88,
        5008.24,
        4960.30,
        6549.84,
        6585.23,
        6718.32,
        6732.71,
    ]
)

# airtovac2 is a slow PyAstronomy call and the lines never change, so convert
# once at import rather than on each of the ~1000 templates loaded per run.
HOST_LINES_VAC = pyasl.airtovac2(HOST_LINES)

# Half-width of the mask around each line: 400 km/s expressed as a redshift.
LINE_MASK_DISPERSION = 4e2 / 3e5


def host_line_ranges(z_obj=0.0):
    """(n_lines, 2) array of [low, high] wavelength bounds to mask at ``z_obj``."""

    centres = (1 + z_obj) * HOST_LINES_VAC
    return np.column_stack(
        [centres * (1 - LINE_MASK_DISPERSION), centres * (1 + LINE_MASK_DISPERSION)]
    )


def mask_host_lines(Data, z_obj=0.0):
    """Drop rows whose wavelength falls inside any host emission line.

    Vectorised over lines and pixels at once. The previous implementation
    looped over the ten lines and mapped a Python lambda across every pixel,
    which dominated the cost of loading the template bank.
    """

    ranges = host_line_ranges(z_obj)
    lam = Data[:, 0]

    inside = (lam[:, None] > ranges[None, :, 0]) & (lam[:, None] < ranges[None, :, 1])

    return Data[~inside.any(axis=1)]


def mask_lines_bank(Data, z_obj=0):
    """Mask host lines in a bank template. Bank objects are at rest, so z = 0."""

    return mask_host_lines(Data, z_obj)
