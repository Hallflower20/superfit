# Methods

What the fit actually does, and the parts of it worth understanding before
quoting a result.

## How a fit is put together

`Superfit.run()` prepares the observation and then hands it to
`all_parameter_space`, which evaluates every (redshift, A_v) grid point
across a process pool. Within a grid point, the fit solves every
(galaxy, supernova) pair at once as a set of matrix products and scores each
with chi2 — that is where nearly all the time goes. A `GridSolver` is built
once per redshift and computes everything A_v cannot change (the validity
masks, the weights, every galaxy-side contraction); each A_v then only pays
for the supernova-side terms, which is bit-identical to solving the reddened
bank from scratch and roughly three times faster over the default A_v grid.

```
Spectrum  ──▶  mask host lines / telluric  ──▶  normalise  ──▶  int_obj
    │                                                             │
    └──▶  bin to `resolution`  ──▶  error model  ──▶  sigma       │
                                                        │         │
                          all_parameter_space  ◀────────┴─────────┘
                                    │
              for each z:  GridSolver  ──▶  for each A_v:  solve
                                    │             ──▶  top `iterations`
                          rank, dedupe by SN, write results.csv
```

The quantities involved:

- `int_obj` — the observation, masked, normalised and interpolated onto the
  fitting grid.
- `redshift` — the redshift or redshifts to optimise over.
- `extconstant` — the A_v values to optimise over.
- `templates_sn_trunc`, `templates_gal_trunc` — the supernova and host
  galaxy templates in play, after `temp_sn_tr` / `temp_gal_tr` narrowing.
- `lam` — the wavelength grid the fit runs on.
- `sigma` — the per-pixel uncertainty, from the chosen error model.
- `minimum_overlap` — the least fraction of `lam` a template must cover to be
  scored at all.

The results table keeps at most `iterations` templates per grid point, ranks
everything by reduced chi2, and deduplicates by supernova template keeping
the best-scoring appearance of each. The same supernova can therefore not
appear twice at two different redshifts; the best one wins.

## Why the fit runs in log wavelength

The fitting grid is uniform in `ln(lambda)`, not in Angstroms, because a
redshift is then a *translation*:

    ln(lambda_obs) = ln(lambda_rest) + ln(1 + z)

The displacement is the same for every template and every wavelength, so the
bank is resampled once onto a rest-frame log grid and each trial redshift is
two array slices and a linear blend across the whole bank at once. On a
linear grid every template had to be interpolated onto the observed axis
again at every trial redshift — about a thousand `np.interp` calls per grid
point.

A consequence worth knowing: the grid is uniform in velocity, so its width in
Angstroms grows with wavelength. At 500 km/s a bin is 8 A at 5000 A and 15 A
at 9000 A, where a linear 10 A grid was 10 A everywhere. The blue end is
sampled more finely and the red end less.

Shifting also resamples twice — once onto the rest grid, once for the shift —
so structure finer than a bin is averaged rather than reproduced. For
resolved features the two agree to about 1 part in 10^3; for noise at the
sampling scale it is a smoothing.

## Extinction

A_v is the CCM89 extinction in magnitudes at V, applied to the supernova
template at its **rest** wavelength. That is a modelling choice, not an
accident: reddening at the rest wavelength is host-galaxy dust. Any Milky
Way component is absorbed into the same single term.

For an observed pixel at `lambda`, the law is evaluated at
`lambda / (1 + z)` — exactly the rest wavelength that pixel corresponds to —
so the curve itself is never interpolated.

Only the supernova is reddened; the galaxy template is not. Both are dimmed
by `(1 + z)`.

Negative A_v is not physical extinction. It is kept in the default grid as
slack for a template that is redder than the object, and so that a best fit
at A_v = 0 is an interior point of the searched range rather than an edge —
see the warnings below.

`R_v` shapes the curve that A_v scales; 3.1 is the diffuse Milky Way
average and the usual default.

## A note on the weighting

In the `legacy` profile — the default — the supernova and galaxy amplitudes
are solved by minimising the **unweighted** residual, while the chi2 that
ranks the templates is weighted by the error spectrum. The reported chi2 is
therefore not the minimum of the reported model.

The `modern` profile (`weighted_solve`) carries `1/sigma^2` into the solve as
well, which on the bundled test spectrum lowers chi2 by a median 20% and
reorders the tail of the results. The top match is unchanged there, but it
need not be in general.

`legacy` is the default so that existing results reproduce. If you are
starting fresh, `modern` is the statistically consistent choice:

```bash
superfit fit spectrum.flm --z 0.127 --profile modern
```

Whichever you use, the run's `config.json` records it.

## Warnings you should not ignore

The run prints a warning when a best-fit `A_v` or `Z` lands on the first or
last value of the grid it was allowed to search:

    WARNING: rank 1: best-fit A_v = -2 sits on the lower edge of the searched
    grid [-2, 2]; the true optimum is probably outside it, so widen the range
    before trusting this value

A parameter pinned to the edge of its range is a boundary artefact, not a
measurement. Widen `--av-range` (or the redshift range) and refit.

The top three matches are checked, so a warning about rank 3 does not
invalidate rank 1 — but it does say the grid is too narrow somewhere.

## The template bank

The bank contains major subclasses — calcium-rich supernovae, type II
flashers, TDEs, SLSN-I and II, among others — in separate folders, for more
accurate classification. The default binning is 10 A.

superfit is only as good as the template bank it uses. `superfit bank status`
prints the checksum of the one you have; quote it when you publish.

The bank has two top-level folders:

- `original_resolution` — the raw spectra, wavelengths in the observed
  frame. Also holds the wiserep metadata files for each object (name,
  redshift, observation date), which the fit reads to build phases and
  display names.
- `binnings` — the same spectra binned and redshift-corrected, at 10 A and
  30 A. A fit at either of those resolutions uses these directly; any other
  resolution bins `original_resolution` on the fly, which is slower.

Host galaxy lines can be masked in the templates as well as in the object,
which is what `mask_galaxy_lines` does. The lines involved are listed in
`mask_lines_bank` in `superfit/Header_Binnings.py`.

## Reading a spectrum

Wavelengths are converted to Angstroms on the way in, from whatever unit the
file declares. Inverse variance is converted to an uncertainty, with
non-positive entries becoming NaN rather than infinity — a pixel with zero
weight has no measured uncertainty, and pretending otherwise puts an
infinity into the chi2.

A `Spectrum` sorts descending wavelengths rather than complaining, and
refuses input that cannot be a spectrum — mismatched lengths, duplicate
wavelengths, all-NaN flux — naming the problem. It also refuses a spectrum
too short to fit meaningfully, rather than returning a fit based on thirty
pixels.

Normalisation divides the flux by its median, with a guard: a
continuum-subtracted or sky-dominated spectrum can have a median
indistinguishable from zero, and dividing by it amplified the spectrum by a
factor of 10^15. When that happens the scatter is used instead, and a warning
is issued.
