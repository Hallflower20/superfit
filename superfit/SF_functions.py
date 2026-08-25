import numpy as np
from scipy import interpolate
import extinction
from astropy import table
from astropy.table import Table
from astropy.io import ascii
import contextlib
import itertools
import os
from PyAstronomy import pyasl
import multiprocessing as mp
import threading
from tqdm import tqdm

from superfit.get_metadata import get_metadata
from superfit.error_routines import savitzky_golay, linear_error
from superfit.Header_Binnings import (
    bin_spectrum_bank,
    mask_host_lines,
    mask_lines_bank,
    normalise_flux,
)
from superfit.loggrid import RedshiftableTemplates
from superfit.packed import load_template
from superfit.paths import binning_dir

np.seterr(divide="ignore", invalid="ignore")



def redshifted_models(sn_bank, gal_bank, lam, z, extcon, R_v=3.1):
    """Supernova and galaxy models at one (redshift, extinction) point.

    Returns ``(sn, gal)`` shaped (1, n_sn, n_lam) and (n_gal, 1, n_lam), ready
    to broadcast into the chi2 grid.

    Both banks are dimmed by (1 + z). Only the supernova is reddened, and the
    extinction law is evaluated at the rest wavelength each observed pixel
    corresponds to, ``lam / (1 + z)``. That is a modelling choice, not an
    accident: reddening at the rest wavelength is host-galaxy dust. Galactic
    dust would act at the observed wavelength, and this single A_v absorbs
    both.

    Splitting the redshift from the extinction matters: the shift is shared
    by every template and by every extinction value, so a caller sweeping A_v
    at fixed z should hoist :func:`redshift_bank` out of its loop.
    """

    one_plus_z = 1.0 + z
    sn = redshift_bank(sn_bank, z)
    gal = redshift_bank(gal_bank, z)
    reddening = 10 ** (-0.4 * extcon * Alam(np.asarray(lam) / one_plus_z, R_v=R_v))

    return (sn * reddening)[np.newaxis, :, :], gal[:, np.newaxis, :]


def redshift_bank(bank, z):
    """Every template in ``bank`` at redshift ``z``, dimmed by (1 + z)."""

    return bank.at_redshift(z) / (1.0 + z)


# The telluric A band, in Angstroms.
TELLURIC_BAND = (7594.0, 7680.0)


def remove_telluric(spectrum):
    """Blank the telluric A band, returning a new array.

    The previous version wrote a -10000 sentinel into the caller's own array
    and then mapped every -10000 in the spectrum to NaN, so it both destroyed
    its input and masked any genuine flux that happened to equal -10000.
    Extra columns (a shipped error spectrum) are now carried through instead
    of being dropped.
    """

    out = np.array(spectrum, dtype=float, copy=True)

    in_band = (out[:, 0] >= TELLURIC_BAND[0]) & (out[:, 0] <= TELLURIC_BAND[1])
    out[in_band, 1] = np.nan

    return out


def Alam(lamin, A_v=1, R_v=3.1):

    """A_lambda in magnitudes from the CCM89 law, at ``A_v`` and ``R_v``.

    Magnitudes is what the name promises and what every caller assumes: each
    one forms the transmission itself, as ``10 ** (-0.4 * A_v * Alam(lam))``.

    This used to return ``extinction.apply(ccm89(...), ones)``, which is the
    *transmission* ``10 ** (-0.4 * A_lambda)`` rather than A_lambda. The
    callers then exponentiated a second time, computing
    ``10 ** (-0.4 * A_v * 10 ** (-0.4 * A_lambda))``. That curve removes more
    red light than blue -- the reddening ran backwards, so positive A_v made
    spectra bluer and negative A_v was doing the work extinction should.

    ``ccm89`` is linear in A_v, so the default ``A_v = 1`` gives the shape of
    the curve and a caller's own A_v scales it exactly. It is strict about
    dtype and wants a float64 array; ``np.asarray`` is the element-by-element
    ``[float(i) for i in lamin]`` conversion this used to do, without the
    Python loop.
    """

    return extinction.ccm89(np.asarray(lamin, dtype=np.float64), A_v, R_v)


def error_obj(kind, lam, object_to_fit):

    """
    Per-pixel uncertainty for the observation, resampled onto ``lam``.

    parameters
    ----------

    kind: "sg" (Savitzky-Golay), "linear", or "included" to use an error
        column supplied with the observation.
    lam: the wavelength grid the fit runs on.
    object_to_fit: a Spectrum, a path to one, or an (n, 2+) array.

    returns
    -------

    Error, one value per element of ``lam``, NaN outside the observation.

    """

    from superfit.spectrum import Spectrum

    spectrum = Spectrum.coerce(object_to_fit)

    # Median-normalise so the error is on the same scale as the flux the
    # chi2 sees. Already-normalised input is unaffected.
    object_spec = spectrum.normalised().as_array()

    if kind == "included":
        if object_spec.shape[1] < 3:
            raise ValueError(
                "error_spectrum='included' needs an uncertainty column, but "
                "{!r} has only wavelength and flux. Use 'sg' or 'linear', or "
                "supply error= when building the Spectrum.".format(spectrum.name)
            )
        error = np.column_stack([object_spec[:, 0], object_spec[:, 2]])

    elif kind == "linear":
        error = linear_error(object_spec)

    elif kind == "sg":
        error = savitzky_golay(object_spec)

    else:
        raise ValueError(
            "Unknown error_spectrum {!r}; expected 'sg', 'linear' or "
            "'included'.".format(kind)
        )

    object_err_interp = interpolate.interp1d(
        error[:, 0], error[:, 1], bounds_error=False, fill_value=np.nan
    )

    return object_err_interp(lam)


def solve_grid(sn, gal, int_obj, sigma, weighted=False):
    """Fit every (galaxy, supernova) pair at once and score each with chi2.

    For each pair this solves the same unweighted 2-parameter least squares
    problem as before -- observed ~= b * sn + d * gal, with negative
    amplitudes rejected -- and then evaluates the sigma-weighted chi2 of that
    solution, counting only pixels where the observation and both templates
    are defined.

    The direct form builds the full (n_gal, n_sn, n_lam) residual cube, and
    built it twice: once to count valid pixels and once for the chi2 itself.
    At the shipped bank size that is ten-odd 60 MB temporaries per grid point,
    which is what kept a redshift scan memory-bound rather than CPU-bound.

    Expanding the square

        chi2 = sum w*obj^2 - 2b sum w*obj*sn - 2d sum w*obj*gal
               + b^2 sum w*sn^2 + 2bd sum w*sn*gal + d^2 sum w*gal^2

    turns every term into a contraction over wavelength, i.e. one matrix
    product per term. Each sum has to run over the *same* pixel set -- the
    pixels valid for that particular pair -- which is why the masks are folded
    into the operands rather than applied afterwards. Peak memory drops from
    O(n_gal * n_sn * n_lam) to O(n_gal * n_sn).

    Parameters
    ----------
    sn : (1, n_sn, n_lam) array
    gal : (n_gal, 1, n_lam) array
    int_obj : (n_lam,) array
    sigma : (n_lam,) array
    weighted : bool
        When False (the default, and the historical behaviour) the amplitudes
        minimise the *unweighted* residual while the chi2 that ranks them is
        weighted by sigma -- so the reported chi2 is not the minimum of the
        reported model. When True the solve carries 1/sigma**2 too, which on
        the bundled spectrum lowers chi2 by a median 20% and can change which
        template wins. Exposed as the "weighted_solve" config option.

    Returns
    -------
    b, d, chi2, times : (n_gal, n_sn) arrays
        Supernova amplitude, galaxy amplitude, chi2, and the number of pixels
        each chi2 was accumulated over.
    """

    S = np.ascontiguousarray(sn[0], dtype=np.float64)  # (n_sn, n_lam)
    G = np.ascontiguousarray(gal[:, 0], dtype=np.float64)  # (n_gal, n_lam)
    obj = np.asarray(int_obj, dtype=np.float64)
    sig = np.asarray(sigma, dtype=np.float64)

    n_lam = S.shape[1]

    # Where each array actually carries a value.
    mS = np.isfinite(S)
    mG = np.isfinite(G)
    m_obj = np.isfinite(obj)

    S0 = np.where(mS, S, 0.0)
    G0 = np.where(mG, G, 0.0)
    obj0 = np.where(m_obj, obj, 0.0)

    # --- weighted contractions over wavelength ----------------------------
    # A pixel contributes only if the observation, the error and both
    # templates are defined there.
    valid_obs = m_obj & np.isfinite(sig)

    A = (mS & valid_obs).astype(np.float64)  # (n_sn, n_lam)
    B = (mG & valid_obs).astype(np.float64)  # (n_gal, n_lam)

    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where(valid_obs, 1.0 / sig**2, 0.0)

    times = B @ A.T

    SA = S0 * A
    GB = G0 * B

    t_oo = (B * (w * obj0 * obj0)) @ A.T
    t_os = B @ (SA * (w * obj0)).T
    t_og = (GB * (w * obj0)) @ A.T
    t_ss = B @ (SA * SA * w).T
    t_sg = (GB * w) @ SA.T
    t_gg = (GB * GB * w) @ A.T

    # --- amplitudes -------------------------------------------------------
    if weighted:
        # The normal equations for the chi2 that is actually reported: every
        # sum carries 1/sigma**2 and runs over the pair's own overlap. These
        # are the same six contractions the chi2 is built from.
        n_ss, n_gg, n_gs, n_so, n_go = t_ss, t_gg, t_sg, t_os, t_og
    else:
        # Historical behaviour: the amplitudes minimise the UNWEIGHTED
        # residual even though the chi2 they are scored by is weighted, so the
        # reported chi2 is not the minimum of the reported model. Kept as the
        # default so existing results reproduce; see the "weighted_solve"
        # option. Note sum(sn^2) runs over the supernova's own coverage rather
        # than its overlap with the galaxy, which is also preserved here.
        n_ss = (S0 * S0).sum(axis=1)[np.newaxis, :]          # (1, n_sn)
        n_gg = (G0 * G0).sum(axis=1)[:, np.newaxis]          # (n_gal, 1)
        n_gs = G0 @ S0.T                                      # (n_gal, n_sn)
        n_so = (S0 @ obj0)[np.newaxis, :]                     # (1, n_sn)
        n_go = (G0 @ obj0)[:, np.newaxis]                     # (n_gal, 1)

    with np.errstate(divide="ignore", invalid="ignore"):
        c = 1.0 / (n_ss * n_gg - n_gs**2)
        b = c * (n_gg * n_so - n_gs * n_go)
        d = c * (n_ss * n_go - n_gs * n_so)

    b = np.where(b < 0, np.nan, b)
    d = np.where(d < 0, np.nan, d)

    # --- chi2 -------------------------------------------------------------
    with np.errstate(invalid="ignore"):
        chi2 = (
            t_oo
            - 2.0 * b * t_os
            - 2.0 * d * t_og
            + b**2 * t_ss
            + 2.0 * b * d * t_sg
            + d**2 * t_gg
        )

    # Rounding can push an essentially perfect fit a hair below zero.
    chi2 = np.where(chi2 < 0, 0.0, chi2)

    # A rejected amplitude means the whole residual row was NaN, which the
    # direct form scored as zero valid pixels; the overlap cut then rejects it.
    rejected = ~np.isfinite(b) | ~np.isfinite(d)
    times = np.where(rejected, 0.0, times)
    chi2 = np.where(rejected, 0.0, chi2)

    return b, d, chi2, times


def core(
    int_obj,
    z,
    extcon,
    sn,
    gal,
    templates_sn_trunc,
    templates_gal_trunc,
    lam,
    iterations,
    sigma,
    **kwargs
):

    """

    Inputs:
    ------

    z - an array of redshifts

    extcon - array of values of A_v


    Outputs:
    --------


    Astropy table with the names for the best fit supernova and host galaxy,

    constants of proportionality for both the host galaxy and supernova templates,

    the value of chi2, the corresponding redshift and A_v.



    """

    kind = kwargs["kind"]
    original = kwargs["original"]
    minimum_overlap = kwargs["minimum_overlap"]

    # What goes in the SPECTRUM column. Falls back to the filename when the
    # observation came from disk and no name was given.
    name = kwargs.get("spectrum_name")
    if name is None:
        name = os.path.basename(original) if isinstance(original, str) else "spectrum"

    # sn and gal arrive already redshifted, extincted and resampled onto the
    # observed grid: on a log grid that is a shift shared by the whole bank,
    # so the caller does it once per redshift rather than once per grid point.

    b, d, chi2, times = solve_grid(
        sn, gal, int_obj, sigma, weighted=bool(kwargs.get("weighted_solve", False))
    )

    overlap = times / len(lam) > minimum_overlap

    # avoid short overlaps
    chi2[~overlap] = np.inf

    reduchi2 = chi2 / (times - 2) ** 2
    reduchi2 = np.where(reduchi2 == 0, 1e10, reduchi2)

    reduchi2_once = chi2 / (times - 2)
    reduchi2_once = np.where(reduchi2_once == 0, 1e10, reduchi2_once)

    # Flatten the matrix out and obtain indices corresponding values of proportionality constants
    reduchi2_1d = reduchi2.ravel()

    index = np.argsort(reduchi2_1d)

    # Column-wise accumulation, then ONE Table at the end. Building a
    # single-row astropy Table per result and stacking them cost more than
    # the chi2 it was reporting -- 11.6 ms against 9.7 ms per grid point on
    # the shipped bank, because each Table validates and converts 13 columns.

    redchi2 = []
    spectra = []
    galaxies = []
    supernovae = []
    const_sn = []
    const_gal = []
    phases = []
    bands = []
    frac_sn = []
    frac_gal = []
    chi2_dof = []
    chi2_dof2 = []

    for i in range(iterations):

        idx = np.unravel_index(index[i], reduchi2.shape)
        rchi2 = reduchi2[idx]

        redchi2.append(rchi2)

        supernova_file = templates_sn_trunc[idx[1]]
        host_galaxy_file = str(templates_gal_trunc[idx[0]])
        host_galaxy_file = host_galaxy_file[host_galaxy_file.rfind("/") + 1 :]

        bb = b[idx[0]][idx[1]]
        dd = d[idx[0]][idx[1]]

        # `sn` arrives already reddened -- it is the model the fitter scored.
        # This used to multiply by the extinction a second time, and at the
        # observed rather than the rest wavelength the model was reddened at,
        # so the reported flux split described a model nobody had fitted.
        sn_flux = sn[0, idx[1], :]
        gal_flux = gal[idx[0], 0, :]
        sn_contribution = bb * np.nanmean(sn_flux)
        gal_contribution = dd * np.nanmean(gal_flux)
        total = sn_contribution + gal_contribution

        ii = supernova_file.rfind(":")

        spectra.append(os.path.basename(name))
        galaxies.append(host_galaxy_file)
        supernovae.append(supernova_file)
        const_sn.append(bb)
        const_gal.append(dd)
        phases.append(supernova_file[ii + 1 : -1])
        bands.append(supernova_file[-1])
        frac_sn.append(sn_contribution / total)
        frac_gal.append(gal_contribution / total)
        chi2_dof.append(reduchi2_once[idx])
        chi2_dof2.append(reduchi2[idx])

    outputs = table.Table(
        [
            np.array(spectra, dtype="S200"),
            np.array(galaxies, dtype="S200"),
            np.array(supernovae, dtype="S200"),
            np.array(const_sn, dtype="f"),
            np.array(const_gal, dtype="f"),
            np.full(iterations, z, dtype="f"),
            np.full(iterations, extcon, dtype="f"),
            np.array(phases, dtype="S200"),
            np.array(bands, dtype="S200"),
            np.array(frac_sn, dtype="f"),
            np.array(frac_gal, dtype="f"),
            np.array(chi2_dof, dtype="f"),
            np.array(chi2_dof2, dtype="f"),
        ],
        names=(
            "SPECTRUM",
            "GALAXY",
            "SN",
            "CONST_SN",
            "CONST_GAL",
            "Z",
            "A_v",
            "Phase",
            "Band",
            "Frac(SN)",
            "Frac(gal)",
            "CHI2/dof",
            "CHI2/dof2",
        ),
    )

    return outputs, redchi2


def mask_gal_lines(Data, z_obj):
    """Mask host emission lines in the observed spectrum, at the object's z."""

    return mask_host_lines(Data, np.asarray(z_obj).reshape(-1)[0])

@contextlib.contextmanager
def _single_threaded_blas():
    """Pin BLAS to one thread for the duration of the block, if we can.

    chi2 is evaluated as several small matrix products. Every worker process
    would otherwise open its own BLAS thread pool sized to the whole machine,
    so the processes spend their time contending rather than computing.
    Degrades to a no-op when threadpoolctl is not installed.
    """

    try:
        from threadpoolctl import threadpool_limits
    except ImportError:
        yield
        return

    with threadpool_limits(limits=1):
        yield


# Read-only fit inputs shared with the worker processes.
#
# The previous version bound the whole template bank into a functools.partial
# and handed it to Pool.imap, so every (z, A_v) grid point pickled and shipped
# ~13 MB of templates to a worker and back. Instead the bank is published here
# before the pool is created; under fork the children inherit it for free.
_SHARED_STATE = None

# Held for as long as _SHARED_STATE is published; see all_parameter_space.
_FIT_LOCK = threading.Lock()


def _init_worker(state):
    """Fallback for start methods that do not inherit memory (spawn)."""

    global _SHARED_STATE
    _SHARED_STATE = state


def _fit_one_grid_point(args):
    """Fit the bank at one redshift, across a group of extinction values.

    The redshift shift is shared by the whole bank and does not depend on
    A_v, so it is done once here and reused across the group. That is what
    log binning buys: on a linear grid this step was a per-template
    interpolation repeated for every extinction value.
    """

    z, extinctions = args
    state = _SHARED_STATE

    lam = state["lam"]

    # Done once for the whole extinction group: the shift is shared.
    sn_at_z = redshift_bank(state["sn_bank"], z)[np.newaxis, :, :]
    gal_at_z = redshift_bank(state["gal_bank"], z)[:, np.newaxis, :]

    # Extinction acts at the template's rest wavelength, which for observed
    # pixel lam is lam / (1 + z). Evaluating the law there directly is exact
    # and avoids interpolating the extinction curve. Rest frame means this
    # models host-galaxy dust; see redshifted_models.
    alam_rest = Alam(lam / (1.0 + z), R_v=state["R_v"])

    results = []
    for extcon in extinctions:
        reddening = 10 ** (-0.4 * extcon * alam_rest)

        result, _ = core(
            state["int_obj"],
            z,
            extcon,
            sn_at_z * reddening,
            gal_at_z,
            state["sn_names"],
            state["gal_names"],
            lam,
            state["iterations"],
            state["sigma"],
            **state["kwargs"]
        )
        results.append(result)

    return results


# Past this many processes the fork and queue overhead costs more than the
# extra parallelism returns: a grid point is only tens of milliseconds of work.
# Measured on a 441-point scan over the shipped bank, 16-32 workers ran in
# ~2s while 244 took ~6s. Raise it with "n_cores" if your grid is much larger.
DEFAULT_MAX_WORKERS = 32


def build_tasks(redshift, extconstant, n_workers):
    """Split the (redshift, A_v) grid into units of work.

    A unit is one redshift plus a group of extinction values, because the
    redshift shift is shared across extinction and wants to be done once.
    Grouping too coarsely would leave workers idle when there is only one
    redshift, so the extinction axis is split into enough pieces to keep the
    pool busy and no more.
    """

    redshift = np.atleast_1d(redshift)
    extconstant = np.atleast_1d(extconstant)

    target_tasks = max(1, n_workers * 4)
    groups_per_z = int(np.ceil(target_tasks / len(redshift)))
    n_groups = max(1, min(len(extconstant), groups_per_z))

    extinction_groups = np.array_split(extconstant, n_groups)

    return [(float(z), group) for z in redshift for group in extinction_groups]


def resolve_worker_count(requested, n_tasks):
    """How many processes to actually start.

    Never more than there are grid points to evaluate, and never more than the
    CPUs this process is allowed to run on. `cores = 16` was hard-coded, which
    both oversubscribed a login node and ignored a large allocation.
    """

    try:
        available = len(os.sched_getaffinity(0))
    except AttributeError:
        available = mp.cpu_count()

    if requested is None or requested <= 0:
        requested = min(available, DEFAULT_MAX_WORKERS)
    else:
        requested = min(int(requested), available)

    return max(1, min(int(requested), n_tasks))

def all_parameter_space(
    int_obj,
    redshift,
    extconstant,
    templates_sn_trunc,
    templates_gal_trunc,
    lam,
    resolution,
    iterations,
    **kwargs
):

    """

    This function loops the core function of superfit over two user given arrays, one for redshift and one for

    the extinction constant, it then sorts all the chi2 values obtained and plots the curve that corresponds

    to the smallest one. This is not the recommended method to use, since it takes the longest time, it is

    rather a method to check results if there are any doubts with the two recommended methods.



    Parameters
    ----------

    Truncated SN and HG template libraries, extinction array and redshift array, lambda axis and **kwargs for the object path.



    Returns
    -------

    Astropy table with the best fit parameters: Host Galaxy and Supernova proportionality

    constants, redshift, extinction law constant and chi2 value, plots are optional.

    In this version for the fit the same SN can appear with two different redshifts (since it is a brute-force

    method in which we go over the whole parameter space we don't want to eliminate any results).





    For plotting: in order not to plot every single result the user can choose how many to plot, default

    set to the first three.


    """

    import time

    # Popped, not read: the rest of kwargs is forwarded to every worker, and
    # the whole Parameters object -- template lists, grids -- has no business
    # being pickled once per process.
    parameters = kwargs.pop("parameters")
    mask_galaxy_lines = parameters.mask_galaxy_lines
    metadata = get_metadata(parameters)

    print("superfit started")
    #print(len(templates_sn_trunc))
    start = time.time()

    # The full path to write, not a prefix to concatenate onto: see
    # superfit.output for why that distinction earned its own module.
    results_path = os.fspath(kwargs["results_path"])

    templates_sn_trunc_dict = {}
    templates_gal_trunc_dict = {}
    sn_spec_files = [str(x) for x in metadata.shorhand_dict.values()]
    path_dict = {}

    all_bank_files = [str(x) for x in metadata.dictionary_all_trunc_objects.values()]

    #print(len(all_bank_files))

    # Templates the bank lists but cannot supply. One bank's metadata names
    # eight spectra that are not on disk, and another has a hundred-odd whose
    # flux column carries the literal string "None" -- and a single one of
    # them used to end a 15000-template fit with a traceback from inside
    # np.loadtxt. A template that will not load is one template missing from
    # the comparison, not a reason to throw the other 15000 away, so they are
    # collected and reported together at the end.
    unreadable = []

    def read_sn(path, reader):
        try:
            return load_template(path, reader)
        except (OSError, ValueError) as exc:
            unreadable.append((path, "{}: {}".format(type(exc).__name__, exc)))
            return None

    if resolution == 10 or resolution == 30:

        for i in range(0, len(all_bank_files)):
            a = all_bank_files[i]

            full_name = a[a.find("sne") :]
            one_sn = read_sn(
                os.path.join(binning_dir(resolution), full_name), "loadtxt"
            )
            if one_sn is None:
                continue
            if mask_galaxy_lines:
                one_sn = mask_lines_bank(one_sn)

            idx = all_bank_files[i].rfind("/") + 1
            filename = all_bank_files[i][idx:]

            short_name = str(metadata.shorhand_dict[filename])

            path_dict[short_name] = all_bank_files[i]

            templates_sn_trunc_dict[short_name] = one_sn

    else:
        # Any other resolution: bin the original-resolution bank on the fly.
        for i in range(0, len(all_bank_files)):

            one_sn = read_sn(all_bank_files[i], "kill_header")
            if one_sn is None:
                continue
            if mask_galaxy_lines:
                one_sn = mask_lines_bank(one_sn)
            one_sn = bin_spectrum_bank(one_sn, resolution)

            idx = all_bank_files[i].rfind("/") + 1
            filename = all_bank_files[i][idx:]

            short_name = str(metadata.shorhand_dict[filename])

            path_dict[short_name] = all_bank_files[i]
            templates_sn_trunc_dict[short_name] = one_sn

    for i in range(0, len(templates_gal_trunc)):

        one_gal = load_template(templates_gal_trunc[i], "loadtxt")
        one_gal = bin_spectrum_bank(one_gal, resolution)
        templates_gal_trunc_dict[templates_gal_trunc[i]] = one_gal

    sn_spec_files = [x for x in path_dict.keys()]

    if unreadable:
        # Loud, and with examples: a bank quietly fitting against fewer
        # templates than it advertises is a result nobody can reproduce.
        print(
            "WARNING: {0} of {1} supernova templates could not be read and "
            "were left out of this fit. The bank lists them but cannot "
            "supply them; the classification below is against the remaining "
            "{2}.".format(
                len(unreadable), len(all_bank_files), len(sn_spec_files)
            )
        )
        for path, why in unreadable[:3]:
            print("    {0}: {1}".format(path, why))
        if len(unreadable) > 3:
            print("    ... and {0} more".format(len(unreadable) - 3))

    if not sn_spec_files:
        raise RuntimeError(
            "None of the {0} supernova templates this bank lists could be "
            "read, so there is nothing to fit against. The first failure was "
            "{1}".format(
                len(all_bank_files),
                unreadable[0][1] if unreadable else "not recorded",
            )
        )

    print(
        "Loaded {0} SN and {1} galaxy templates in {2: .1f}s".format(
            len(sn_spec_files), len(templates_gal_trunc_dict), time.time() - start
        )
    )

    # Resample the whole bank onto a rest-frame log grid aligned with the
    # observed one. From here a redshift is a shift, so this is the only time
    # anything is interpolated.
    observed_grid = kwargs["observed_grid"]
    max_z = float(np.max(redshift))

    bank_start = time.time()
    sn_bank = RedshiftableTemplates.from_templates(
        [templates_sn_trunc_dict[name][:, 0] for name in sn_spec_files],
        [templates_sn_trunc_dict[name][:, 1] for name in sn_spec_files],
        observed_grid,
        max_z,
    )
    gal_bank = RedshiftableTemplates.from_templates(
        [templates_gal_trunc_dict[name][:, 0] for name in templates_gal_trunc],
        [templates_gal_trunc_dict[name][:, 1] for name in templates_gal_trunc],
        observed_grid,
        max_z,
    )
    print(
        "Resampled onto a {0:.0f} km/s log grid ({1} bins) in {2: .1f}s".format(
            observed_grid.velocity_resolution, len(observed_grid),
            time.time() - bank_start,
        )
    )

    # The error spectrum depends only on the observed object, so derive it once
    # here rather than once per grid point inside core().
    sigma = error_obj(kwargs["kind"], lam, kwargs["original"])

    n_grid_points = len(redshift) * len(extconstant)
    n_workers = resolve_worker_count(kwargs.get("n_cores"), n_grid_points)
    params = build_tasks(redshift, extconstant, n_workers)

    fit_start = time.time()

    # _SHARED_STATE is how the bank reaches forked workers without being
    # pickled per task, which means it is one variable for the whole process.
    # Two threads calling this at once would publish over each other -- and
    # forking a pool from several threads at once deadlocks besides -- so
    # concurrent fits take turns here rather than corrupting each other. The
    # work inside is already spread across every core, so nothing is lost.
    global _SHARED_STATE

    with _FIT_LOCK:
        _SHARED_STATE = {
            "int_obj": int_obj,
            "sn_names": sn_spec_files,
            "gal_names": templates_gal_trunc,
            "sn_bank": sn_bank,
            "gal_bank": gal_bank,
            "lam": lam,
            "iterations": iterations,
            "sigma": sigma,
            "R_v": kwargs.get("R_v", 3.1),
            "kwargs": kwargs,
        }

        try:
            if n_workers == 1:
                # One grid point, or an explicit request for serial: skip the
                # pool entirely rather than paying to fork for a single task.
                results = [_fit_one_grid_point(p) for p in tqdm(params)]
            else:
                # Each worker's chi2 is a handful of small BLAS calls. Left to
                # themselves they would each spin up a full thread pool, so N
                # workers times N BLAS threads fight over the same cores.
                # Forked children inherit this limit.
                with _single_threaded_blas():
                    results = _run_pool(params, n_workers)
        finally:
            _SHARED_STATE = None

    print(
        "Fitted {0} grid points ({1} redshifts x {2} A_v) as {3} task(s) on "
        "{4} worker(s) in {5: .1f}s".format(
            n_grid_points, len(redshift), len(extconstant), len(params),
            n_workers, time.time() - fit_start,
        )
    )

    result = table.vstack([t for group in results for t in group])

    result.sort("CHI2/dof2")

    result = table.unique(result, keys="SN", keep="first")

    result.sort("CHI2/dof2")

    ascii.write(result, results_path, format="csv", fast_writer=False, overwrite=True)

    for message in grid_edge_warnings(result, redshift, extconstant):
        print("WARNING: " + message)

    end = time.time()
    print("Runtime: {0: .2f}s ".format(end - start))

    return


def grid_edge_warnings(result, redshift, extconstant, n_check=3):
    """Flag best fits that landed on the edge of the searched grid.

    A parameter pinned to the first or last value it was allowed to take
    usually means the real optimum lies outside the range, so the reported
    value is a boundary artefact rather than a measurement. The shipped
    A_v grid runs -2 to +2, and a top match at exactly -2 -- unphysical
    negative extinction, at the edge -- is easy to read straight past.

    Returns a list of human-readable warning strings.
    """

    messages = []
    if len(result) == 0:
        return messages

    grids = {"A_v": np.atleast_1d(extconstant), "Z": np.atleast_1d(redshift)}

    for column, grid in grids.items():
        if grid.size < 2 or column not in result.colnames:
            continue

        low, high = float(np.min(grid)), float(np.max(grid))
        values = np.asarray(result[column][:n_check], dtype=float)

        for rank, value in enumerate(values):
            if np.isclose(value, low):
                edge = "lower"
            elif np.isclose(value, high):
                edge = "upper"
            else:
                continue

            messages.append(
                "rank {0}: best-fit {1} = {2:g} sits on the {3} edge of the "
                "searched grid [{4:g}, {5:g}]; the true optimum is probably "
                "outside it, so widen the range before trusting this "
                "value".format(rank + 1, column, value, edge, low, high)
            )

    return messages


def _run_pool(params, n_workers):
    """Evaluate every grid point across ``n_workers`` processes."""

    ctx = mp.get_context("fork" if "fork" in mp.get_all_start_methods() else None)

    pool_kwargs = {}
    if ctx.get_start_method() != "fork":
        # Children will not inherit _SHARED_STATE; send it once per worker.
        pool_kwargs = {"initializer": _init_worker, "initargs": (_SHARED_STATE,)}

    # Hand out several grid points at a time so the queue round-trip is
    # amortised over real work.
    chunksize = max(1, len(params) // (n_workers * 4))

    with ctx.Pool(processes=n_workers, **pool_kwargs) as pool:
        return list(
            tqdm(
                pool.imap(_fit_one_grid_point, params, chunksize=chunksize),
                total=len(params),
            )
        )
