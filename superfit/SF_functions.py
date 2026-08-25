import numpy as np
from scipy import interpolate
import extinction
from astropy import table
from astropy.table import Table
from astropy.io import ascii
import contextlib
import itertools
import os
import time
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
from superfit import packed as packed_module
from superfit.packed import load_template
from superfit.paths import sne_dir

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


def solve_extinction_batch(sn, gal, int_obj, sigma, reddening, weighted=False):
    """Score every (extinction, galaxy, supernova) triple in one pass.

    :func:`solve_grid` does one extinction at a time, and doing that for a
    whole A_v grid repeats work that does not depend on A_v at all. Of the six
    weighted contractions the chi2 is built from, only three move with
    extinction -- ``t_os``, ``t_ss`` and ``t_sg``. ``times``, ``t_oo``,
    ``t_og`` and ``t_gg`` do not, because reddening is finite and positive
    everywhere, so it changes no template's coverage: the validity masks, the
    weights, and every contraction that does not touch the supernova flux are
    identical across the grid. The unweighted amplitude solve adds one more
    that does move (``n_gs``) and two that do not.

    So the A_v-independent halves are computed once, and each A_v-dependent
    contraction becomes a single matrix product with the extinction axis
    folded into its rows: stacking (n_av, n_gal, n_lam) into
    (n_av * n_gal, n_lam) turns twenty-one modest products into one large one.

    Parameters
    ----------
    sn : (n_sn, n_lam) array
        Templates at this redshift, *without* extinction applied.
    gal : (n_gal, n_lam) array
    int_obj, sigma : (n_lam,) arrays
    reddening : (n_av, n_lam) array
        The transmission for each A_v, ``10 ** (-0.4 * A_v * A_lambda)``.
    weighted : bool
        As in :func:`solve_grid`.

    Returns
    -------
    b, d, chi2, times : (n_av, n_gal, n_sn) arrays

    Agrees with :func:`solve_grid` applied one extinction at a time to about
    1e-14 relative, which is the reassociation of the same sums rather than a
    different calculation.
    """

    S = np.ascontiguousarray(sn, dtype=np.float64)
    G = np.ascontiguousarray(gal, dtype=np.float64)
    obj = np.asarray(int_obj, dtype=np.float64)
    sig = np.asarray(sigma, dtype=np.float64)
    red = np.ascontiguousarray(reddening, dtype=np.float64)

    n_av = red.shape[0]
    n_gal = G.shape[0]
    n_sn = S.shape[0]

    mS = np.isfinite(S)
    mG = np.isfinite(G)
    m_obj = np.isfinite(obj)

    S0 = np.where(mS, S, 0.0)
    G0 = np.where(mG, G, 0.0)
    obj0 = np.where(m_obj, obj, 0.0)

    valid_obs = m_obj & np.isfinite(sig)

    A = (mS & valid_obs).astype(np.float64)
    B = (mG & valid_obs).astype(np.float64)

    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where(valid_obs, 1.0 / sig**2, 0.0)

    SA = S0 * A
    GB = G0 * B

    # --- independent of A_v: once, then broadcast --------------------------
    times_2d = B @ A.T
    t_oo = ((B * (w * obj0 * obj0)) @ A.T)[np.newaxis, :, :]
    t_og = ((GB * (w * obj0)) @ A.T)[np.newaxis, :, :]
    t_gg = ((GB * GB * w) @ A.T)[np.newaxis, :, :]

    # --- moves with A_v: one stacked product each --------------------------
    def stacked(rows, scale):
        return (rows[np.newaxis, :, :] * scale[:, np.newaxis, :]).reshape(
            n_av * n_gal, -1
        )

    shape = (n_av, n_gal, n_sn)
    t_os = (stacked(B, red) @ (SA * (w * obj0)).T).reshape(shape)
    t_ss = (stacked(B, red * red) @ (SA * SA * w).T).reshape(shape)
    t_sg = (stacked(GB * w, red) @ SA.T).reshape(shape)

    if weighted:
        n_ss, n_gg, n_gs, n_so, n_go = t_ss, t_gg, t_sg, t_os, t_og
    else:
        # The historical unweighted solve; see solve_grid for why it is the
        # default. sum(sn^2) again runs over the supernova's own coverage.
        n_ss = ((red * red) @ (S0 * S0).T)[:, np.newaxis, :]
        n_gg = (G0 * G0).sum(axis=1)[np.newaxis, :, np.newaxis]
        n_gs = (stacked(G0, red) @ S0.T).reshape(shape)
        n_so = ((red * obj0) @ S0.T)[:, np.newaxis, :]
        n_go = (G0 @ obj0)[np.newaxis, :, np.newaxis]

    with np.errstate(divide="ignore", invalid="ignore"):
        c = 1.0 / (n_ss * n_gg - n_gs**2)
        b = c * (n_gg * n_so - n_gs * n_go)
        d = c * (n_ss * n_go - n_gs * n_so)

    b = np.where(b < 0, np.nan, b)
    d = np.where(d < 0, np.nan, d)

    with np.errstate(invalid="ignore"):
        chi2 = (
            t_oo
            - 2.0 * b * t_os
            - 2.0 * d * t_og
            + b**2 * t_ss
            + 2.0 * b * d * t_sg
            + d**2 * t_gg
        )

    chi2 = np.where(chi2 < 0, 0.0, chi2)

    rejected = ~np.isfinite(b) | ~np.isfinite(d)
    times = np.where(rejected, 0.0, np.broadcast_to(times_2d, shape))
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
        # basename, not a search for the last "/": a Windows bank path
        # separates with backslashes and would come through whole.
        host_galaxy_file = os.path.basename(str(templates_gal_trunc[idx[0]]))

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


# Supernova templates scored at once inside a worker. The batched kernel holds
# arrays of (n_av, n_gal, block), so the block is what bounds a worker's
# memory: 21 x 128 x 2048 float64 is 44 MB per term and there are a handful of
# terms live at once. Blocking is also what keeps the products large without
# letting the score arrays grow with the whole bank -- 21 x 128 x 15561 would
# be 335 MB per term, several times over.
DEFAULT_SN_BLOCK = 2048


def _top_k_per_extinction(reduchi2, k, sn_offset, n_sn_total):
    """Indices and values of the ``k`` smallest scores, per extinction.

    ``argpartition`` rather than a full sort: only ``k`` of these are ever
    looked at -- ten, by default -- and sorting the rest is the dominant cost
    once the bank is large. A full argsort of the 128 x 12316 scores the modern
    bank produces is 64 ms, which over a 21-point A_v grid is 1.3 s spent
    ordering candidates nobody reads; taking the top ten costs 6 ms.

    Ranks are returned as flat indices into the *whole* (n_gal, n_sn) grid, so
    a caller can merge blocks without knowing how the work was split. Ties are
    broken by that index, which makes the ordering independent of block size
    and of how the partition happened to land -- a full argsort left it to
    whatever quicksort did.
    """

    n_av, n_gal, n_block = reduchi2.shape
    flat = reduchi2.reshape(n_av, n_gal * n_block)
    take = min(k, flat.shape[1])

    # Column index within the block, mapped back to the full grid.
    local = np.argpartition(flat, take - 1, axis=1)[:, :take]
    values = np.take_along_axis(flat, local, axis=1)

    gal_index = local // n_block
    sn_index = local % n_block + sn_offset
    global_flat = gal_index * n_sn_total + sn_index

    # Sort the k candidates by (value, index). lexsort takes its last key as
    # the primary one and sorts along the last axis, so this is "by score,
    # ties by position". The tiebreak is what makes the ordering reproducible:
    # without it, equal scores would be ordered by wherever argpartition
    # happened to leave them, which depends on the block size.
    order = np.lexsort((global_flat, values))
    values = np.take_along_axis(values, order, axis=1)
    global_flat = np.take_along_axis(global_flat, order, axis=1)

    return global_flat, values


def _fit_one_block(args):
    """Score one block of supernova templates, at one redshift, for every A_v.

    Returns per-extinction candidate records rather than a table: numeric
    indices and numbers, which pickle cheaply and cost nothing to merge.
    Building an Astropy table here meant one per grid point, validated and
    type-converted thirteen columns at a time, and then thrown at vstack.
    """

    z, rows = args
    state = _SHARED_STATE

    lam = state["lam"]
    extinctions = state["extinctions"]
    n_sn_total = len(state["sn_names"])

    # Only this block is shifted. Shifting all 15561 templates in every task
    # was the largest duplicated cost in a scan.
    sn_at_z = state["sn_bank"].at_redshift(z, rows=rows) / (1.0 + z)
    gal_at_z = redshift_bank(state["gal_bank"], z)

    # Extinction acts at the template's rest wavelength, which for observed
    # pixel lam is lam / (1 + z). Evaluating the law there directly is exact
    # and avoids interpolating the extinction curve. Rest frame means this
    # models host-galaxy dust; see redshifted_models.
    alam_rest = Alam(lam / (1.0 + z), R_v=state["R_v"])
    reddening = 10.0 ** (-0.4 * np.asarray(extinctions)[:, np.newaxis] * alam_rest)

    b, d, chi2, times = solve_extinction_batch(
        sn_at_z,
        gal_at_z,
        state["int_obj"],
        state["sigma"],
        reddening,
        weighted=state["weighted"],
    )

    # Same scoring as core(): short overlaps are excluded, and a zero score
    # is a degenerate fit rather than a perfect one.
    overlap = times / len(lam) > state["minimum_overlap"]
    chi2 = np.where(overlap, chi2, np.inf)

    with np.errstate(divide="ignore", invalid="ignore"):
        dof = times - 2
        reduchi2 = chi2 / dof**2
        reduchi2_once = chi2 / dof

    reduchi2 = np.where(reduchi2 == 0, 1e10, reduchi2)
    reduchi2_once = np.where(reduchi2_once == 0, 1e10, reduchi2_once)

    ranked, values = _top_k_per_extinction(
        reduchi2, state["iterations"], rows.start, n_sn_total
    )

    # Mean flux of each candidate's model, for the Frac(SN) split. Taken from
    # the same reddened template the fitter scored, and computed as a
    # contraction rather than a nanmean per candidate: the numerator is
    # sum(reddening * flux) over the template's own coverage and the
    # denominator its pixel count, which is one small matrix product for the
    # whole block.
    mS = np.isfinite(sn_at_z)
    covered = mS.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        sn_mean = (reddening @ np.where(mS, sn_at_z, 0.0).T) / covered
    gal_mean = np.nanmean(gal_at_z, axis=1)

    records = []
    for a in range(len(extinctions)):
        gal_index = ranked[a] // n_sn_total
        sn_index = ranked[a] % n_sn_total
        local = sn_index - rows.start
        records.append(
            {
                "z": float(z),
                "extcon": float(extinctions[a]),
                "gal_index": gal_index,
                "sn_index": sn_index,
                "b": b[a][gal_index, local],
                "d": d[a][gal_index, local],
                "reduchi2": values[a],
                "reduchi2_once": reduchi2_once[a][gal_index, local],
                "sn_mean": sn_mean[a][local],
                "gal_mean": gal_mean[gal_index],
            }
        )

    return records


# Past this many processes the fork and queue overhead costs more than the
# extra parallelism returns: a grid point is only tens of milliseconds of work.
# Measured on a 441-point scan over the shipped bank, 16-32 workers ran in
# ~2s while 244 took ~6s. Raise it with "n_cores" if your grid is much larger.
DEFAULT_MAX_WORKERS = 32


def build_tasks(redshift, n_sn, n_workers, block=DEFAULT_SN_BLOCK):
    """Split the work into units of (one redshift, one block of supernovae).

    Not by extinction. The whole A_v grid is evaluated together now -- that is
    what makes the products large and lets the A_v-independent contractions be
    computed once -- so splitting it would undo the saving and, at an exact
    redshift, would also shift the same bank once per group.

    The block is shrunk when there is not enough other work to keep the pool
    busy: one redshift and one block would leave every worker but one idle.
    """

    redshift = np.atleast_1d(redshift)
    n_z = len(redshift)

    target_tasks = max(1, n_workers * 2)
    blocks_wanted = int(np.ceil(target_tasks / n_z))
    block = int(min(block, max(1, int(np.ceil(n_sn / max(1, blocks_wanted))))))

    bounds = list(range(0, n_sn, block)) + [n_sn]

    return [
        (float(z), slice(lo, hi))
        for z in redshift
        for lo, hi in zip(bounds[:-1], bounds[1:])
    ]


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

# Prepared banks this process has built, keyed by everything they depend on.
# One entry is the whole resampled bank -- 90 MB for the largest -- so this is
# deliberately a single slot rather than an unbounded cache: a batch fits the
# same bank over and over, and holding two of them costs more than rebuilding
# the one that was displaced.
_prepared = {}
_prepared_key = None


def forget_prepared_bank():
    """Drop the prepared bank this process is holding.

    For tests, and for a Session that has finished: the largest catalogue's
    resampled templates are around 90 MB.
    """

    global _prepared_key

    _prepared.clear()
    _prepared_key = None


def bank_state_token(bank_dir, revalidate):
    """A token that changes when the bank's contents change.

    Part of the prepared-bank key, and the reason a fit that reuses a prepared
    bank still notices an edited one. Without it the cache would answer from
    the first fit forever: the revalidation lives inside the loading, and a
    cache hit skips the loading, so the check would be skipped exactly when it
    was most needed.

    Revalidating costs the directory listing that validates each pack -- one,
    not one per template. ``revalidate=False`` returns a constant instead,
    which is a caller stating that it has pinned the bank; see Session.
    """

    if not revalidate:
        return "pinned"

    fingerprints = []
    try:
        directories = packed_module.packable_directories(bank_dir)
    except OSError:
        return "unknown"

    for relative in directories:
        pack = packed_module.open_pack(bank_dir, relative)
        fingerprints.append((relative, pack.fingerprint if pack else None))

    return tuple(fingerprints)


def prepared_bank_key(parameters, metadata, observed_grid, max_z, bank_state):
    """Everything a prepared bank depends on.

    The templates that go in depend on the bank -- its contents, by way of
    ``bank_state`` -- the binning, whether host lines are masked, and which
    objects the metadata scan selected. The resampling on top of that depends
    on the observed grid and on the largest redshift the grid has to reach.
    Nothing else, so two spectra sharing a wavelength range share the whole
    preparation, and two that do not share only the loading.
    """

    return (
        parameters.bank_dir,
        bank_state,
        parameters.resolution,
        bool(parameters.mask_galaxy_lines),
        parameters.metadata_key,
        len(metadata.dictionary_all_trunc_objects),
        tuple(str(x) for x in parameters.templates_gal_trunc),
        observed_grid.identity,
        round(float(max_z), 12),
    )


def prepare_bank(
    parameters,
    metadata,
    observed_grid,
    max_z,
    templates_gal_trunc,
    resolution,
    mask_galaxy_lines,
    revalidate_bank=True,
    quiet=False,
):
    """Load the bank and resample it onto a log grid, or reuse the last one.

    This is the expensive half of a fit that has nothing to do with the
    spectrum: on the largest bank it is a directory listing to validate the
    pack, 12316 template reads, host-line masking and the resample -- about
    2.6 s warm, and considerably more when the filesystem is cold. A batch of
    spectra on a common grid repeats all of it per spectrum unless it is kept,
    which is what Session exists to do.

    Returns ``(sn_bank, gal_bank, sn_names, unreadable)``.
    """

    global _prepared_key

    # Before the cache is consulted, not after: the bank's state is part of
    # the key, so a bank edited between two fits misses the cache rather than
    # being answered from it.
    packed_module.begin_fit(revalidate=revalidate_bank)
    key = prepared_bank_key(
        parameters,
        metadata,
        observed_grid,
        max_z,
        bank_state_token(parameters.bank_dir, revalidate_bank),
    )

    if _prepared_key == key and "bank" in _prepared:
        if not quiet:
            print("Reusing the prepared bank from the previous fit")
        return _prepared["bank"]

    load_start = time.time()

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

    # This fit's own bank. Everything below is resolved from it rather than
    # from the process-wide one, which another Superfit's construction can
    # have moved since this fit was set up.
    bank_dir = parameters.bank_dir
    sne_root = sne_dir(bank_dir=bank_dir)

    # 10 A and 30 A are pre-binned in the bank; anything else is binned from
    # the original resolution on the fly. Both walk the same list, so this is
    # one pass with the reader chosen up front rather than two copies of it.
    pre_binned = resolution in (10, 30)

    if pre_binned:
        binned_root = sne_dir(resolution, bank_dir=bank_dir)
        # relpath against the bank's own supernova root, not a search for the
        # substring "sne" in the whole path: a bank living under, say,
        # /home/snelling/ would have that search cut the path in the wrong
        # place, and on Windows the separators are backslashes.
        read_from = [
            os.path.join(binned_root, os.path.relpath(source, sne_root))
            for source in all_bank_files
        ]
        reader = "loadtxt"
    else:
        read_from = list(all_bank_files)
        reader = "kill_header"

    # One bulk read: the pack that covers these is resolved once for the whole
    # batch rather than re-resolved per template.
    source_of = dict(zip(read_from, all_bank_files))

    for path, one_sn, error in packed_module.load_templates(
        read_from, reader, bank_dir=bank_dir
    ):
        if error is not None:
            unreadable.append((path, "{}: {}".format(type(error).__name__, error)))
            continue

        if mask_galaxy_lines:
            one_sn = mask_lines_bank(one_sn)
        if not pre_binned:
            one_sn = bin_spectrum_bank(one_sn, resolution)

        source_path = source_of[path]
        short_name = str(metadata.shorhand_dict[os.path.basename(source_path)])

        path_dict[short_name] = source_path
        templates_sn_trunc_dict[short_name] = one_sn

    for i in range(0, len(templates_gal_trunc)):

        one_gal = load_template(
            templates_gal_trunc[i], "loadtxt", bank_dir=bank_dir
        )
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
            len(sn_spec_files), len(templates_gal_trunc_dict),
            time.time() - load_start,
        )
    )

    # Resample the whole bank onto a rest-frame log grid aligned with the
    # observed one. From here a redshift is a shift, so this is the only time
    # anything is interpolated.
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


    _prepared["bank"] = (sn_bank, gal_bank, sn_spec_files, unreadable)
    _prepared_key = key
    return _prepared["bank"]


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

    sn_bank, gal_bank, sn_spec_files, unreadable = prepare_bank(
        parameters,
        metadata,
        kwargs["observed_grid"],
        float(np.max(redshift)),
        templates_gal_trunc,
        resolution,
        mask_galaxy_lines,
        revalidate_bank=kwargs.get("revalidate_bank", True),
    )

    # The error spectrum depends only on the observed object, so derive it once
    # here rather than once per grid point inside core().
    sigma = error_obj(kwargs["kind"], lam, kwargs["original"])

    # What goes in the SPECTRUM column. Falls back to the filename when the
    # observation came from disk and no name was given.
    spectrum_name = kwargs.get("spectrum_name")
    if spectrum_name is None:
        original = kwargs["original"]
        spectrum_name = (
            os.path.basename(original) if isinstance(original, str) else "spectrum"
        )

    n_grid_points = len(redshift) * len(extconstant)
    n_workers = resolve_worker_count(kwargs.get("n_cores"), n_grid_points)
    params = build_tasks(redshift, len(sn_spec_files), n_workers)

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
            "extinctions": np.atleast_1d(extconstant).astype(float),
            "weighted": bool(kwargs.get("weighted_solve", False)),
            "minimum_overlap": kwargs["minimum_overlap"],
            "R_v": kwargs.get("R_v", 3.1),
        }

        try:
            if n_workers == 1:
                # One block, or an explicit request for serial: skip the pool
                # entirely rather than paying to fork for a single task.
                results = [_fit_one_block(p) for p in tqdm(params)]
            else:
                # Each worker's chi2 is a handful of BLAS calls. Left to
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

    result = assemble_results(
        merge_blocks(
            [record for group in results for record in group], iterations
        ),
        sn_spec_files,
        templates_gal_trunc,
        spectrum_name,
        iterations,
    )

    result.sort("CHI2/dof2")

    result = table.unique(result, keys="SN", keep="first")

    result.sort("CHI2/dof2")

    ascii.write(result, results_path, format="csv", fast_writer=False, overwrite=True)

    for message in grid_edge_warnings(result, redshift, extconstant):
        print("WARNING: " + message)

    print("Fit and write: {0: .2f}s ".format(time.time() - start))

    return


def merge_blocks(records, iterations):
    """Reduce each grid point's per-block candidates to its overall best.

    A worker sees one block of the bank, so its top ten are the best ten *in
    that block*. The grid point's answer is the best ten across all of them,
    which is this: group by (redshift, A_v), concatenate, and take the top ten
    again. Without it a run reports ten candidates per block rather than per
    grid point -- 394 rows where there should be 20, once duplicates are
    dropped.

    Ordering is by score and then by the position in the (galaxy, supernova)
    grid, the same total order each block used, so the result does not depend
    on how the bank was divided.
    """

    grouped = {}
    for record in records:
        grouped.setdefault((record["z"], record["extcon"]), []).append(record)

    merged = []
    for (z, extcon), group in grouped.items():
        # No shortcut for a single block. Each block is already trimmed to
        # `iterations`, so one block would come out right either way -- but a
        # function whose guarantee holds only when the work happened to be
        # split a particular way is a function whose guarantee does not hold.
        joined = {
            key: np.concatenate([np.atleast_1d(r[key]) for r in group])
            for key in (
                "gal_index", "sn_index", "b", "d",
                "reduchi2", "reduchi2_once", "sn_mean", "gal_mean",
            )
        }

        n_sn_total = int(joined["sn_index"].max()) + 1
        position = joined["gal_index"] * n_sn_total + joined["sn_index"]
        order = np.lexsort((position, joined["reduchi2"]))[:iterations]

        record = {"z": z, "extcon": extcon}
        record.update({key: value[order] for key, value in joined.items()})
        merged.append(record)

    return merged


def assemble_results(records, sn_names, gal_names, spectrum_name, iterations):
    """One results table, built once from every block's candidate records.

    Workers hand back numeric indices and numbers. Names, phases and bands are
    resolved here, and the table is constructed once for the whole run rather
    than once per grid point: each construction validates and type-converts
    thirteen columns, and at one per grid point that cost more than the chi2 it
    was reporting.
    """

    if not records:
        raise RuntimeError(
            "The fit produced no candidates at all. Every template was either "
            "rejected for too little overlap with the observation or had no "
            "usable amplitude solution."
        )

    z = np.concatenate([np.full(len(r["gal_index"]), r["z"]) for r in records])
    av = np.concatenate([np.full(len(r["gal_index"]), r["extcon"]) for r in records])
    gal_index = np.concatenate([r["gal_index"] for r in records])
    sn_index = np.concatenate([r["sn_index"] for r in records])
    b = np.concatenate([r["b"] for r in records])
    d = np.concatenate([r["d"] for r in records])
    reduchi2 = np.concatenate([r["reduchi2"] for r in records])
    reduchi2_once = np.concatenate([r["reduchi2_once"] for r in records])
    sn_mean = np.concatenate([r["sn_mean"] for r in records])
    gal_mean = np.concatenate([r["gal_mean"] for r in records])

    # A candidate whose score is the degenerate sentinel never belonged in the
    # table; a grid point with fewer than `iterations` real candidates used to
    # pad the table out with them.
    keep = np.isfinite(reduchi2) & (reduchi2 < 1e10)
    if not keep.any():
        keep = np.ones(len(reduchi2), dtype=bool)

    z, av = z[keep], av[keep]
    gal_index, sn_index = gal_index[keep], sn_index[keep]
    b, d = b[keep], d[keep]
    reduchi2, reduchi2_once = reduchi2[keep], reduchi2_once[keep]
    sn_mean, gal_mean = sn_mean[keep], gal_mean[keep]

    sn_lookup = np.asarray([str(x) for x in sn_names], dtype=object)
    gal_lookup = np.asarray(
        [os.path.basename(str(x)) for x in gal_names], dtype=object
    )

    supernovae = sn_lookup[sn_index]
    galaxies = gal_lookup[gal_index]

    # The shorthand carries "... phase-band : <phase><band>"; the phase is
    # everything after the colon bar the final character, which is the band.
    phases = [name[name.rfind(":") + 1 : -1] for name in supernovae]
    bands = [name[-1] for name in supernovae]

    sn_contribution = b * sn_mean
    gal_contribution = d * gal_mean
    with np.errstate(divide="ignore", invalid="ignore"):
        total = sn_contribution + gal_contribution
        frac_sn = sn_contribution / total
        frac_gal = gal_contribution / total

    return table.Table(
        [
            np.full(len(z), os.path.basename(str(spectrum_name)), dtype="S200"),
            np.array(galaxies, dtype="S200"),
            np.array(supernovae, dtype="S200"),
            b.astype("f"),
            d.astype("f"),
            z.astype("f"),
            av.astype("f"),
            np.array(phases, dtype="S200"),
            np.array(bands, dtype="S200"),
            frac_sn.astype("f"),
            frac_gal.astype("f"),
            reduchi2_once.astype("f"),
            reduchi2.astype("f"),
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
                pool.imap(_fit_one_block, params, chunksize=chunksize),
                total=len(params),
            )
        )
