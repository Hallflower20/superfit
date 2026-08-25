# Configuration reference

Every setting has a default, so you only have to state what you care about —
in practice the spectrum and the redshift. Anything you do not set takes its
value from `superfit.config.DEFAULT_CONFIG`.

The *complete* effective configuration is written to `config.json` in the
run directory, so a one-flag command is as reproducible as a full JSON file.

## How to pass settings

Settings can come as keyword arguments, a dict, a JSON string, a path to a
JSON file, or command-line flags:

```python
from superfit import Superfit

Superfit(spectrum, z=0.127, resolution=30)
Superfit(spectrum, config={"resolution": 30, "minimum_overlap": 0.5})
Superfit(spectrum, config="parameters.json", z=0.127)   # file, one override
```

Precedence, lowest to highest: defaults, then the profile, then `config`,
then keyword arguments.

Some shorthand is accepted, because `z=0.127` is what you mean and
`{"use_exact_z": True, "z_exact": 0.127}` is what the fitter needs:

| Shorthand | Real key |
| --- | --- |
| `z` | `z_exact` (and sets `use_exact_z`) |
| `output_dir` | `saving_results_path` |
| `n_plots` | `how_many_plots` |
| `error_model` | `error_spectrum` |
| `cores` | `n_cores` |
| `resolution_angstrom` | `resolution` |

Unknown keys are rejected rather than ignored, and combinations that cannot
work — masking host lines while scanning redshift, a negative resolution, an
inverted A_v range — are caught up front rather than failing inside a worker
process partway through a run.

### Booleans

On/off settings accept `true`/`false` and the legacy `0`/`1` (and `yes`/`no`,
`on`/`off`) interchangeably, and are normalised to booleans at load time.
Files written years ago and files written today both work.

### Comments

JSON has no comment syntax, so `superfit` ignores any key beginning with an
underscore. That is how `superfit config create` writes a parameter file
that explains itself:

```json
{
    "_resolution": "Binning resolution in Angstroms.",
    "resolution": 10
}
```

## Profiles

A profile is a named set of scientific defaults, applied over the built-in
defaults and under everything else.

| Profile | `weighted_solve` | For |
| --- | --- | --- |
| `legacy` (default) | off | Reproducing published superfit results. |
| `modern` | on | New work; the statistically consistent choice. |

A setting given directly wins over the profile it was applied over, and the
recorded profile then reads `custom` — `{"profile": "legacy",
"weighted_solve": true}` is not a legacy run, and `config.json` has to
describe the settings that came out. A `custom` config reloads unchanged.

```bash
superfit fit spectrum.flm --z 0.1 --profile modern
```

```python
Superfit(spectrum, z=0.1, profile="modern")
```

The profile in force is recorded in the run's `config.json`, so a result says
which science produced it. See [methods](methods.md#a-note-on-the-weighting)
for what the difference actually is.

---

## What to fit

`object_to_fit` — path to the spectrum. `None` when the spectrum was passed
directly as arrays or as a `Spectrum`, in which case the spectrum's name is
recorded here instead.

### Reading the spectrum

These apply only when `object_to_fit` is a path. All default to `None`,
meaning "work it out from the file".

`spectrum_columns` — which columns hold wavelength, flux and error. A list of
two or three entries, by name (csv, FITS) or by position (ascii); or a
mapping, which may also name `ivar` instead of `error`:

```python
Superfit("spec.fits", spectrum_columns=["WAVE", "FLUX", "IVAR"])
Superfit("wide.txt", spectrum_columns=[1, 2, 4])
Superfit("spec.fits", spectrum_columns={"wavelength": "lam", "ivar": "weight"})
```

Left unset, columns are matched against the usual spellings:

| Role | Names recognised |
| --- | --- |
| wavelength | `wavelength`, `wave`, `lam`, `lambda`, `wl`, `awav`, `spectral_axis`, `loglam` |
| flux | `flux`, `f_lambda`, `flam`, `fl`, `spec`, `intensity`, `counts` |
| error | `fluxerr`, `flux_err`, `error`, `err`, `sigma`, `e_flux`, `uncertainty`, `noise`, `stdev`, `std` |
| inverse variance | `ivar`, `invvar`, `inverse_variance`, `flux_ivar` |

For ascii files, columns are taken by position — 0, 1, 2 — unless the file
has a commented header line naming them, as in `# wavelength flux fluxerr`,
in which case names work too. Every data row must have the same number of
columns; a ragged file is refused, naming the line, rather than silently
narrowed to its shortest row.

Inverse variance is converted to an uncertainty; entries that are zero or
negative become NaN rather than infinity. Naming an error column explicitly
takes precedence over an `ivar` column that happens to be in the same file.

`wavelength_unit` — `AA`, `nm`, `micron`, `mm`, `cm`, `m`, or `log10(AA)`.
Left unset, taken from the file's `TUNIT`/`CUNIT1` keyword, and Angstroms
assumed if the file says nothing. A column named `loglam` is understood as
log10(Angstrom) even without a unit. An explicit value wins over the file, so
this is also how you correct a mislabelled one.

`spectrum_hdu` — which FITS HDU to read, by index or by name. Left unset,
each is tried in turn: binary tables by column name, then one-dimensional
image HDUs by their wavelength WCS (`CRVAL1` plus `CDELT1` or `CD1_1`, with
`CRPIX1`, linear or log-linear).

An image HDU with no dispersion keyword is refused rather than assumed to be
1 A per pixel, and a *two*-dimensional image is only read when you name its
HDU — for a multispec stack the first row is the flux, but for a long-slit
or drizzled frame it is sky, and that is not a guess worth making for you.

## Redshift

`use_exact_z` — `true` to fit at a single redshift, `false` to scan a range.

`z_exact` — the redshift to fit at, when `use_exact_z` is true.

`z_range_begin`, `z_range_end`, `z_int` — the scan: first redshift, last
redshift, and step.

Masking host galaxy lines needs a single redshift to place them at, so
`mask_galaxy_lines` and a scan are rejected together. Set an exact redshift,
or turn the masking off.

## The wavelength grid

`resolution` — the nominal resolution of the fit in Angstroms, default 10.
The fit itself runs on a grid uniform in `ln(lambda)`
(see [methods](methods.md#why-the-fit-runs-in-log-wavelength)), so this
chooses that grid's step unless `velocity_resolution` overrides it. Values
of 10 and 30 use the pre-binned bank directly; anything else bins the
original-resolution bank on the fly, which is slower.

`velocity_resolution` — the step of the fitting grid, in km/s. Left unset, it
is derived from `resolution` so the log grid holds as many bins as a linear
grid of that many Angstroms would have over the same span — about 500 km/s
for the usual 10 A across the optical.

`lower_lam`, `upper_lam` — the wavelength range to fit over. If they are
equal, the range is taken from the object: its own span, padded by 300 A at
each end.

## Preprocessing

`error_spectrum` — how the per-pixel uncertainty is obtained.

| Value | Meaning |
| --- | --- |
| `sg` (default) | Savitzky-Golay estimate from the spectrum itself. Recommended. |
| `linear` | Linear fit every ten points. |
| `included` | The observation's own uncertainty column. |

`included` fails clearly if the spectrum has no uncertainty — supply one with
`error=` when building the `Spectrum`, or name the column with
`spectrum_columns`.

`mask_galaxy_lines` — mask host galaxy emission lines, in both the templates
and the object. Needs a single redshift.

`mask_telluric` — mask the telluric A band, 7594–7680 A in the observer's
frame.

`minimum_overlap` — the least fraction of the fitting grid a template must
cover to be scored at all, between 0 and 1. Templates covering less are
given infinite chi2. Keep it near the default of 0.7.

## Templates

`bank` — which template bank to fit against, by name: `legacy`,
`modern-curated`, `modern`, or anything registered with `superfit bank
install <name> --from <directory>`. Empty (the default) means "whatever this
machine resolves to", which is what every fit did when there was only one
bank. See [docs/cli.md](cli.md#superfit-bank-list).

Which bank produced a classification is part of the result, so the name is
written into the run's `config.json`. Naming a bank that is not installed
stops the fit rather than falling back to a different one.

`temp_sn_tr`, `temp_gal_tr` — which supernova subtypes and host galaxy types
to consider. Both default to everything in the bank, which is the
recommendation; narrowing them is for asking a specific question, not for
speed.

`epoch_low`, `epoch_high` — restrict templates to a phase window, in days
from maximum light. Equal values (the default, both 0) mean no restriction.

## Extinction

`Alam_low`, `Alam_high`, `Alam_interval` — the A_v grid to search, in
magnitudes at V. Default -2 to 2 in steps of 0.2.

`R_v` — total-to-selective extinction ratio for the CCM89 law. Default 3.1,
the diffuse Milky Way average. Lower it for denser sightlines.

A_v is CCM89 extinction applied to the supernova template at its **rest**
wavelength, so it models host-galaxy dust; any Milky Way component is
absorbed into the same single term. Negative A_v is not physical extinction —
it is kept in the default grid as slack for a template redder than the
object, and so that a best fit at zero extinction is an interior point rather
than a grid edge. See [methods](methods.md#extinction).

## Fitting

`iterations` — how many templates each grid point contributes to the result
table. Default 10. Raise it if the final list is shorter than you want; the
table is deduplicated by SN, so many grid points agreeing on the same
template yields few rows.

`n_cores` — worker processes for the (redshift, A_v) grid. `0` (the default)
means "use the CPUs this process is allowed on", capped at 32 — past that,
process startup costs more than the extra parallelism returns. The pool is
never larger than the number of grid points.

`weighted_solve` — carry the `1/sigma^2` weights into the amplitude solve as
well as the chi2. Off in the `legacy` profile, on in `modern`. See
[methods](methods.md#a-note-on-the-weighting).

`profile` — `legacy` or `modern`; see above.

## Output

`saving_results_path` — where run directories are created. Empty (the
default) means the working directory. Each fit gets its own subdirectory
under this, named after the spectrum.

`overwrite` — replace a run directory that already holds a `results.csv`.
Off by default: losing a fit to a re-run aimed at the wrong place is not
recoverable.

`how_many_plots` — plot this many of the best fits. **Default 0** — no plots
unless you ask.

`show_plot` — display the plots as well as saving them. Default off. Only
does anything with an interactive matplotlib backend; `superfit doctor`
says whether yours is one.

`show_plot_png` — save plots as PNG rather than PDF. Default off.

## The result

`run()` returns a `FitResult`:

```python
result = Superfit(spectrum, z=0.127, n_plots=2).run()

result.results        # the ranked pandas DataFrame, best match first
result.best           # its first row, as a Series
result.directory      # the run directory, as a pathlib.Path
result.artifacts      # every file written, as Paths
result.plots          # just the plots
result.artifact("results.csv")

len(result)           # these fall through to the DataFrame
result["CHI2/dof"]
result.head()
```

Columns of the table:

| Column | Meaning |
| --- | --- |
| `SPECTRUM` | The object's name. |
| `GALAXY` | Host galaxy template. |
| `SN` | Supernova template, as `type/object/instrument phase-band`. |
| `Z` | Redshift of this match. |
| `A_v` | Extinction of this match. |
| `CONST_SN`, `CONST_GAL` | Fitted amplitudes of the two components. |
| `Frac(SN)`, `Frac(gal)` | Fraction of the model flux from each; they sum to 1. |
| `CHI2/dof`, `CHI2/dof2` | Reduced chi2. The table is ranked by `CHI2/dof2`. |
