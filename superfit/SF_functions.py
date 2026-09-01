import numpy as np
from scipy import interpolate
import extinction
from astropy import table
from astropy.io import ascii
import contextlib
import os
import multiprocessing as mp
import threading
from tqdm import tqdm

from superfit.get_metadata import get_metadata
from superfit.error_routines import savitzky_golay, linear_error
from superfit.Header_Binnings import (
    bin_spectrum_bank,
    mask_host_lines,
    mask_lines_bank,
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

    return GridSolver(sn, gal, int_obj, sigma, weighted=weighted).solve()


class GridSolver:
    """The chi2 grid at one redshift, solvable for many extinction values.

    :func:`solve_grid` on a reddened bank recomputes, at every A_v, work that
    A_v cannot change: reddening is a positive finite per-pixel multiplier, so
    it moves no validity mask, and half the contractions -- everything on the
    galaxy and observation side -- do not involve the supernova at all. This
    class computes those once, from the *unreddened* bank at one redshift, and
    :meth:`solve` folds a reddening vector into just the supernova-side terms.

    The per-A_v terms are evaluated by the same operations in the same order
    as :func:`solve_grid` always used, so the results are bit-identical to
    reddening the bank first; the regression suite holds either way. On the
    default 21-point A_v grid this removes roughly half the arithmetic from
    the hot loop.
    """

    def __init__(self, sn, gal, int_obj, sigma, weighted=False):
        S = np.ascontiguousarray(sn[0], dtype=np.float64)  # (n_sn, n_lam)
        G = np.ascontiguousarray(gal[:, 0], dtype=np.float64)  # (n_gal, n_lam)
        obj = np.asarray(int_obj, dtype=np.float64)
        sig = np.asarray(sigma, dtype=np.float64)

        self.weighted = bool(weighted)

        # Where each array actually carries a value. Reddening is finite and
        # positive, so these masks hold for every A_v that will be asked for.
        mS = np.isfinite(S)
        mG = np.isfinite(G)
        m_obj = np.isfinite(obj)

        self._S0 = np.where(mS, S, 0.0)
        G0 = np.where(mG, G, 0.0)
        obj0 = np.where(m_obj, obj, 0.0)

        # --- weighted contractions over wavelength ------------------------
        # A pixel contributes only if the observation, the error and both
        # templates are defined there.
        valid_obs = m_obj & np.isfinite(sig)

        A = (mS & valid_obs).astype(np.float64)  # (n_sn, n_lam)
        B = (mG & valid_obs).astype(np.float64)  # (n_gal, n_lam)

        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(valid_obs, 1.0 / sig**2, 0.0)

        self._times = B @ A.T

        self._SA0 = self._S0 * A  # supernova masked to the valid pixels
        GB = G0 * B

        # Everything the supernova does not appear in.
        self._t_oo = (B * (w * obj0 * obj0)) @ A.T
        self._t_og = (GB * (w * obj0)) @ A.T
        self._t_gg = (GB * GB * w) @ A.T

        # Operands for the supernova-side terms solve() rebuilds per A_v.
        self._B = B
        self._GBw = GB * w
        self._w = w
        self._wobj = w * obj0

        if not self.weighted:
            # Historical behaviour: the amplitudes minimise the UNWEIGHTED
            # residual even though the chi2 they are scored by is weighted, so
            # the reported chi2 is not the minimum of the reported model. Kept
            # as the default so existing results reproduce; see the
            # "weighted_solve" option. Note sum(sn^2) runs over the
            # supernova's own coverage rather than its overlap with the
            # galaxy, which is also preserved here.
            self._G0 = G0
            self._obj0 = obj0
            self._n_gg = (G0 * G0).sum(axis=1)[:, np.newaxis]  # (n_gal, 1)
            self._n_go = (G0 @ obj0)[:, np.newaxis]            # (n_gal, 1)

    def solve(self, reddening=None):
        """Amplitudes and chi2 with the supernovae dimmed by ``reddening``.

        ``reddening`` is a positive, finite transmission per observed pixel,
        or None for the bank exactly as it was given.

        Returns ``(b, d, chi2, times)`` as :func:`solve_grid` does.
        """

        if reddening is None:
            S0, SA = self._S0, self._SA0
        else:
            # Multiplying after the mask is exact: a masked entry is 0.0, and
            # 0.0 times a finite reddening is 0.0.
            S0 = self._S0 * reddening
            SA = self._SA0 * reddening

        t_os = self._B @ (SA * self._wobj).T
        t_ss = self._B @ (SA * SA * self._w).T
        t_sg = self._GBw @ SA.T

        # --- amplitudes ---------------------------------------------------
        if self.weighted:
            # The normal equations for the chi2 that is actually reported:
            # every sum carries 1/sigma**2 and runs over the pair's own
            # overlap. These are the same contractions the chi2 is built from.
            n_ss, n_gg, n_gs, n_so, n_go = (
                t_ss, self._t_gg, t_sg, t_os, self._t_og,
            )
        else:
            n_ss = (S0 * S0).sum(axis=1)[np.newaxis, :]  # (1, n_sn)
            n_gs = self._G0 @ S0.T                        # (n_gal, n_sn)
            n_so = (S0 @ self._obj0)[np.newaxis, :]       # (1, n_sn)
            n_gg = self._n_gg
            n_go = self._n_go

        with np.errstate(divide="ignore", invalid="ignore"):
            c = 1.0 / (n_ss * n_gg - n_gs**2)
            b = c * (n_gg * n_so - n_gs * n_go)
            d = c * (n_ss * n_go - n_gs * n_so)

        b = np.where(b < 0, np.nan, b)
        d = np.where(d < 0, np.nan, d)

        # --- chi2 ---------------------------------------------------------
        with np.errstate(invalid="ignore"):
            chi2 = (
                self._t_oo
                - 2.0 * b * t_os
                - 2.0 * d * self._t_og
                + b**2 * t_ss
                + 2.0 * b * d * t_sg
                + d**2 * self._t_gg
            )

        # Rounding can push an essentially perfect fit a hair below zero.
        chi2 = np.where(chi2 < 0, 0.0, chi2)

        # A rejected amplitude means the whole residual row was NaN, which the
        # direct form scored as zero valid pixels; the overlap cut then
        # rejects it.
        rejected = ~np.isfinite(b) | ~np.isfinite(d)
        times = np.where(rejected, 0.0, self._times)
        chi2 = np.where(rejected, 0.0, chi2)

        return b, d, chi2, times


class SingleTemplateSolver:
    """chi2 for ``observed ~= b * template``, every template at once.

    The stars and QSOs of the modern banks are fit ALONE -- no host galaxy
    underneath, and a QSO never offered as a host for a transient -- so their
    solve has one amplitude, not two. The conventions mirror
    :class:`GridSolver` exactly: the chi2 is weighted by 1/sigma**2 over the
    pixels where the observation, the error and the template are all defined;
    a negative amplitude is rejected; and under the legacy profile the
    amplitude minimises the *unweighted* residual (over the template's own
    coverage) while the reported chi2 is weighted, so category rows rank
    against supernova rows on the same footing.
    """

    def __init__(self, templates, int_obj, sigma, weighted=False):
        T = np.ascontiguousarray(templates, dtype=np.float64)  # (n_t, n_lam)
        obj = np.asarray(int_obj, dtype=np.float64)
        sig = np.asarray(sigma, dtype=np.float64)

        self.weighted = bool(weighted)

        mT = np.isfinite(T)
        m_obj = np.isfinite(obj)

        self._T0 = np.where(mT, T, 0.0)
        obj0 = np.where(m_obj, obj, 0.0)

        valid_obs = m_obj & np.isfinite(sig)
        A = (mT & valid_obs).astype(np.float64)

        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(valid_obs, 1.0 / sig**2, 0.0)

        self._times = A.sum(axis=1)
        self._TA = self._T0 * A
        self._t_oo = A @ (w * obj0 * obj0)
        self._w = w
        self._wobj = w * obj0
        self._obj0 = obj0

    def solve(self, reddening=None):
        """Amplitude and chi2 with the templates dimmed by ``reddening``.

        Returns ``(b, chi2, times)``, each (n_templates,).
        """

        if reddening is None:
            T0, TA = self._T0, self._TA
        else:
            T0 = self._T0 * reddening
            TA = self._TA * reddening

        t_os = TA @ self._wobj
        t_ss = (TA * TA) @ self._w

        if self.weighted:
            n_ss, n_so = t_ss, t_os
        else:
            n_ss = (T0 * T0).sum(axis=1)
            n_so = T0 @ self._obj0

        with np.errstate(divide="ignore", invalid="ignore"):
            b = n_so / n_ss

        b = np.where(b < 0, np.nan, b)

        with np.errstate(invalid="ignore"):
            chi2 = self._t_oo - 2.0 * b * t_os + b**2 * t_ss

        chi2 = np.where(chi2 < 0, 0.0, chi2)

        rejected = ~np.isfinite(b)
        times = np.where(rejected, 0.0, self._times)
        chi2 = np.where(rejected, 0.0, chi2)

        return b, chi2, times


def fit_standalone(
    names,
    arrays,
    int_obj,
    sigma,
    lam,
    observed_grid,
    redshifts,
    extinctions,
    R_v=3.1,
    weighted=False,
    minimum_overlap=0.7,
    spectrum_name="spectrum",
):
    """Fit a bank of host-less templates over a (z, A_v) grid.

    This is the fit behind the star and QSO categories: each template is
    matched to the observation on its own, over every redshift in
    ``redshifts`` (the caller passes ``[0.0]`` for stars, which are
    foreground) and every A_v in ``extinctions``, with the extinction law
    evaluated at the rest wavelength exactly as the supernova fit does.

    Returns an astropy Table with the same columns as :func:`core` -- one row
    per template at its best grid point, ``GALAXY`` fixed to ``"none"`` and
    the whole model flux attributed to the template -- or None when nothing
    survives the overlap cut. Rows are comparable with the supernova rows
    they will be ranked against: same chi2 conventions, same reduced-chi2
    denominators.
    """

    if not names:
        return None

    bank = RedshiftableTemplates.from_templates(
        [arrays[name][:, 0] for name in names],
        [arrays[name][:, 1] for name in names],
        observed_grid,
        float(np.max(redshifts)),
    )

    n = len(names)
    best_red2 = np.full(n, np.inf)
    best_red1 = np.full(n, np.inf)
    best_b = np.zeros(n)
    best_z = np.zeros(n)
    best_av = np.zeros(n)

    n_lam = len(lam)

    for z in np.atleast_1d(redshifts):
        z = float(z)
        solver = SingleTemplateSolver(
            redshift_bank(bank, z), int_obj, sigma, weighted=weighted
        )
        alam_rest = Alam(lam / (1.0 + z), R_v=R_v)

        for extcon in np.atleast_1d(extinctions):
            reddening = 10 ** (-0.4 * float(extcon) * alam_rest)
            b, chi2, times = solver.solve(reddening)

            overlap = times / n_lam > minimum_overlap
            chi2 = np.where(overlap, chi2, np.inf)

            with np.errstate(divide="ignore", invalid="ignore"):
                red2 = chi2 / (times - 2.0) ** 2
                red1 = chi2 / (times - 2.0)
            red2 = np.where(red2 == 0, 1e10, red2)
            red1 = np.where(red1 == 0, 1e10, red1)

            better = red2 < best_red2
            best_red2 = np.where(better, red2, best_red2)
            best_red1 = np.where(better, red1, best_red1)
            best_b = np.where(better, b, best_b)
            best_z = np.where(better, z, best_z)
            best_av = np.where(better, float(extcon), best_av)

    keep = np.isfinite(best_red2)
    if not keep.any():
        return None

    kept = [(i, names[i]) for i in np.flatnonzero(keep)]

    # Star names are metadata shorthands and carry a ": {phase}{band}" tail;
    # QSO names are plain and have neither.
    phases, bands = [], []
    for _i, name in kept:
        phase, band = split_phase_band(name) if ":" in name else ("u", "")
        phases.append(phase)
        bands.append(band)

    k = len(kept)
    return table.Table(
        [
            np.array([os.path.basename(spectrum_name)] * k, dtype="S200"),
            np.array(["none"] * k, dtype="S200"),
            np.array([name for _i, name in kept], dtype="S200"),
            np.array([best_b[i] for i, _n in kept], dtype="f"),
            np.zeros(k, dtype="f"),
            np.array([best_z[i] for i, _n in kept], dtype="f"),
            np.array([best_av[i] for i, _n in kept], dtype="f"),
            np.array(phases, dtype="S200"),
            np.array(bands, dtype="S200"),
            np.ones(k, dtype="f"),
            np.zeros(k, dtype="f"),
            np.array([best_red1[i] for i, _n in kept], dtype="f"),
            np.array([best_red2[i] for i, _n in kept], dtype="f"),
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
    reddening = kwargs.get("reddening")

    # What goes in the SPECTRUM column. Falls back to the filename when the
    # observation came from disk and no name was given.
    name = kwargs.get("spectrum_name")
    if name is None:
        name = os.path.basename(original) if isinstance(original, str) else "spectrum"

    # sn and gal arrive already redshifted, extincted and resampled onto the
    # observed grid: on a log grid that is a shift shared by the whole bank,
    # so the caller does it once per redshift rather than once per grid point.
    #
    # A caller sweeping A_v at one redshift passes solved= from a GridSolver,
    # which has already hoisted the extinction-independent half of the solve.
    solved = kwargs.get("solved")
    if solved is None:
        solved = solve_grid(
            sn, gal, int_obj, sigma,
            weighted=bool(kwargs.get("weighted_solve", False)),
        )
    b, d, chi2, times = solved

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

    # A bank smaller than `iterations` has fewer pairs than rows asked for;
    # indexing past the grid raised IndexError rather than returning them all.
    iterations = min(int(iterations), reduchi2_1d.size)

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

        # The model the fitter scored: either `sn` arrives already reddened,
        # or it arrives unreddened with the reddening alongside and only the
        # reported rows pay for the multiply. (This used to multiply an
        # already-reddened bank by the extinction a second time, and at the
        # observed rather than the rest wavelength the model was reddened at,
        # so the reported flux split described a model nobody had fitted.)
        sn_flux = sn[0, idx[1], :]
        if reddening is not None:
            sn_flux = sn_flux * reddening
        gal_flux = gal[idx[0], 0, :]
        sn_contribution = bb * np.nanmean(sn_flux)
        gal_contribution = dd * np.nanmean(gal_flux)
        total = sn_contribution + gal_contribution

        phase, band = split_phase_band(supernova_file)

        spectra.append(os.path.basename(name))
        galaxies.append(host_galaxy_file)
        supernovae.append(supernova_file)
        const_sn.append(bb)
        const_gal.append(dd)
        phases.append(phase)
        bands.append(band)
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


def split_phase_band(shorthand):
    """The (phase, band) a template's shorthand name ends with.

    The shorthand get_metadata builds ends ``": {phase}{band}"``, where phase
    is a number or ``"u"`` for unknown and band is a filter letter that can
    be absent. Splitting blindly at the last character used to hand the last
    digit of the phase to the Band column whenever the band was empty --
    every object the bank's phase table does not list.
    """

    tail = str(shorthand)[str(shorthand).rfind(":") + 1 :].strip()
    if tail and tail != "u" and not tail[-1].isdigit():
        return tail[:-1], tail[-1]
    return tail, ""


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

    # Everything A_v cannot change -- masks, weights, the galaxy-side
    # contractions -- is computed once here and shared across the group.
    solver = GridSolver(
        sn_at_z,
        gal_at_z,
        state["int_obj"],
        state["sigma"],
        weighted=bool(state["kwargs"].get("weighted_solve", False)),
    )

    results = []
    for extcon in extinctions:
        reddening = 10 ** (-0.4 * extcon * alam_rest)

        result, _ = core(
            state["int_obj"],
            z,
            extcon,
            sn_at_z,
            gal_at_z,
            state["sn_names"],
            state["gal_names"],
            lam,
            state["iterations"],
            state["sigma"],
            solved=solver.solve(reddening),
            reddening=reddening,
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
    redshift shift -- and now the whole extinction-independent half of the
    solve, see GridSolver -- is shared across extinction and wants to be done
    once. Grouping too coarsely would leave workers idle when there is only
    one redshift, so the extinction axis is split into enough pieces to keep
    the pool busy and no more; splitting beyond that repeats the shared work,
    which is why a single worker never splits at all.
    """

    redshift = np.atleast_1d(redshift)
    extconstant = np.atleast_1d(extconstant)

    if n_workers <= 1:
        # Serial: any split of the extinction axis is pure repetition.
        groups_per_z = 1
    else:
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

def _load_metadata_templates(
    metadata, resolution, mask_galaxy_lines, bank_dir, unreadable
):
    """Load every template a metadata scan names, prepared as the fit uses them.

    10 A and 30 A come from the pre-binned bank; anything else is read at the
    original resolution and binned on the fly. Host lines are masked when the
    fit masks them. The relative path is taken against the bank's own
    supernova root, not a search for the substring "sne" in the whole path: a
    bank living under, say, /home/snelling/ would have that search cut the
    path in the wrong place, and on Windows the separators are backslashes.

    Returns ``(names, arrays)``: shorthand names in bank order, and the
    prepared (n, 2) array for each. Templates that will not load are appended
    to ``unreadable`` rather than raised, so one bad file costs one template,
    not the fit.
    """

    sne_root = sne_dir(bank_dir=bank_dir)
    pre_binned = resolution in (10, 30)
    if pre_binned:
        binned_root = sne_dir(resolution, bank_dir=bank_dir)

    names = []
    arrays = {}

    for source_path in (
        str(x) for x in metadata.dictionary_all_trunc_objects.values()
    ):
        try:
            if pre_binned:
                relative = os.path.relpath(source_path, sne_root)
                one = load_template(
                    os.path.join(binned_root, relative), "loadtxt",
                    bank_dir=bank_dir,
                )
            else:
                one = load_template(source_path, "kill_header", bank_dir=bank_dir)
        except (OSError, ValueError) as exc:
            unreadable.append(
                (source_path, "{}: {}".format(type(exc).__name__, exc))
            )
            continue

        if mask_galaxy_lines:
            one = mask_lines_bank(one)
        if not pre_binned:
            one = bin_spectrum_bank(one, resolution)

        short_name = str(metadata.shorhand_dict[os.path.basename(source_path)])
        names.append(short_name)
        arrays[short_name] = one

    return names, arrays


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

    templates_gal_trunc_dict = {}

    all_bank_files = [str(x) for x in metadata.dictionary_all_trunc_objects.values()]

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

    # Revalidate any pack this process already has open, once, against the
    # bank on disk. A long-lived process fitting many spectra must not keep
    # reading a pack built before the bank changed under it.
    packed_module.begin_fit()

    sn_spec_files, templates_sn_trunc_dict = _load_metadata_templates(
        metadata, resolution, mask_galaxy_lines, bank_dir, unreadable
    )

    for i in range(0, len(templates_gal_trunc)):

        one_gal = load_template(
            templates_gal_trunc[i], "loadtxt", bank_dir=bank_dir
        )
        one_gal = bin_spectrum_bank(one_gal, resolution)
        templates_gal_trunc_dict[templates_gal_trunc[i]] = one_gal

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

    # --- the star and QSO categories ---------------------------------------
    # Standalone single-template fits, ranked into the same table as the
    # SN + host rows. Banks of a few dozen templates, solved vectorised in
    # this process: no pool is worth forking for them. The legacy bank has
    # neither, so both resolve off there and nothing below runs.
    category_tables = _fit_categories(
        parameters,
        metadata_kwargs=dict(
            resolution=resolution,
            mask_galaxy_lines=mask_galaxy_lines,
            bank_dir=bank_dir,
        ),
        int_obj=int_obj,
        sigma=sigma,
        lam=lam,
        observed_grid=observed_grid,
        redshift=redshift,
        extconstant=extconstant,
        R_v=kwargs.get("R_v", 3.1),
        weighted=bool(kwargs.get("weighted_solve", False)),
        minimum_overlap=kwargs["minimum_overlap"],
        spectrum_name=kwargs.get("spectrum_name") or "spectrum",
    )

    result = table.vstack(
        [t for group in results for t in group] + category_tables
    )

    result.sort("CHI2/dof2")

    result = table.unique(result, keys="SN", keep="first")

    result.sort("CHI2/dof2")

    ascii.write(result, results_path, format="csv", fast_writer=False, overwrite=True)

    for message in grid_edge_warnings(result, redshift, extconstant):
        print("WARNING: " + message)

    end = time.time()
    print("Runtime: {0: .2f}s ".format(end - start))

    return


def _warn_unreadable_category(unreadable, category):
    if not unreadable:
        return
    print(
        "WARNING: {0} {1} template(s) could not be read and were left out "
        "of this fit:".format(len(unreadable), category)
    )
    for path, why in unreadable[:3]:
        print("    {0}: {1}".format(path, why))
    if len(unreadable) > 3:
        print("    ... and {0} more".format(len(unreadable) - 3))


def _fit_categories(
    parameters,
    metadata_kwargs,
    int_obj,
    sigma,
    lam,
    observed_grid,
    redshift,
    extconstant,
    R_v,
    weighted,
    minimum_overlap,
    spectrum_name,
):
    """Fit the star and QSO categories, returning their result tables.

    Only the categories this fit's bank supplies and this fit's settings ask
    for; see Parameters.fit_stars / fit_qsos. Each is a standalone fit --
    :func:`fit_standalone` -- so no supernova is ever placed on top of a
    QSO, and no star or QSO gets a host underneath it.
    """

    import time

    resolution = metadata_kwargs["resolution"]
    mask_galaxy_lines = metadata_kwargs["mask_galaxy_lines"]
    bank_dir = metadata_kwargs["bank_dir"]

    tables = []

    if parameters.fit_stars and parameters.star_types:
        from superfit.get_metadata import metadata_for_types

        start = time.time()
        unreadable = []
        star_metadata = metadata_for_types(parameters, parameters.star_types)
        names, arrays = _load_metadata_templates(
            star_metadata, resolution, mask_galaxy_lines, bank_dir, unreadable
        )
        _warn_unreadable_category(unreadable, "star")

        # Stars are foreground, so their redshift is pinned to zero. The A_v
        # grid still applies: at z = 0 the law acts at the observed
        # wavelength, which is exactly Galactic dust toward a star.
        result = fit_standalone(
            names, arrays, int_obj, sigma, lam, observed_grid,
            [0.0], extconstant,
            R_v=R_v, weighted=weighted, minimum_overlap=minimum_overlap,
            spectrum_name=spectrum_name,
        )
        if result is not None:
            tables.append(result)
        print(
            "Fitted {0} star templates at z = 0 in {1: .1f}s".format(
                len(names), time.time() - start
            )
        )

    if parameters.fit_qsos and parameters.templates_qso:
        start = time.time()
        unreadable = []
        names = []
        arrays = {}

        for path in parameters.templates_qso:
            try:
                one = load_template(path, "loadtxt", bank_dir=bank_dir)
            except (OSError, ValueError) as exc:
                unreadable.append(
                    (path, "{}: {}".format(type(exc).__name__, exc))
                )
                continue
            if mask_galaxy_lines:
                one = mask_lines_bank(one)
            one = bin_spectrum_bank(one, resolution)

            name = "QSO/" + os.path.basename(str(path))
            names.append(name)
            arrays[name] = one

        _warn_unreadable_category(unreadable, "QSO")

        # QSOs are extragalactic: the full redshift grid, like a supernova --
        # just never with a host underneath and never under a transient.
        result = fit_standalone(
            names, arrays, int_obj, sigma, lam, observed_grid,
            redshift, extconstant,
            R_v=R_v, weighted=weighted, minimum_overlap=minimum_overlap,
            spectrum_name=spectrum_name,
        )
        if result is not None:
            tables.append(result)
        print(
            "Fitted {0} QSO templates over {1} redshift(s) in {2: .1f}s".format(
                len(names), len(np.atleast_1d(redshift)), time.time() - start
            )
        )

    return tables


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
