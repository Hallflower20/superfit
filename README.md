# superfit

Spectral classification of supernovae. `superfit` fits an observed spectrum
against a bank of supernova and host-galaxy templates over a grid of
redshift and extinction, and ranks the matches by chi2.

Formerly published as **NGSF** (Next Generation SuperFit). `import NGSF`
still works and warns; new code should use `superfit`.

## Install

```bash
pip install superfit
superfit bank install
```

`bank install` downloads the 74 MB template bank, checks it against a known
checksum, and unpacks it into your platform's per-user data directory, which
`superfit` searches automatically. Nothing else needs setting up.

Check it worked:

```bash
superfit doctor
```

If anything is missing — a dependency, the bank, a writable output
directory — `doctor` says which, and exits non-zero.

## Quick start

One spectrum, one redshift:

```bash
superfit fit spectrum.flm --z 0.127
```

```
Best match: Ic/1994I/KAST phase-band : -2.57B
  z = 0.127   A_v = 1.2   chi2/dof = 7.244   SN contributes 87%

Written to ./spectrum
  config.json
  binned.txt
  results.csv
```

Or, if you don't know the redshift, search for it:

```bash
superfit fit spectrum.flm --scan-z 0.0 0.3 --plots 3
```

From Python, on arrays you already have — nothing has to touch disk:

```python
from superfit import Superfit

result = Superfit(wavelength=lam, flux=flux, error=err, z=0.127).run()

print(result.best["SN"], result.best["CHI2/dof"])
print(result[["SN", "GALAXY", "A_v", "CHI2/dof"]].head())
```

`result` is a `FitResult`: the ranked table (a pandas DataFrame, best match
first, also at `result.results`), the directory the run was written to, and
the list of files in it. It indexes and iterates like the DataFrame it
wraps.

From a file instead, in ascii, csv or FITS:

```python
from superfit import Spectrum, Superfit

spectrum = Spectrum.from_file("SN2021urb.flm")
result = Superfit(spectrum, z=0.127).run()
```

A `Spectrum` is reusable, so trying a few redshifts costs one read:

```python
for z in (0.12, 0.127, 0.13):
    result = Superfit(spectrum, z=z, output_dir=f"z{z}").run()
    print(z, result.best["CHI2/dof"])
```

Each fit owns its own settings, so building several and running them in any
order — or concurrently — does what it looks like it does.

## What you get

Each fit writes a directory named after the spectrum:

```
<output_dir>/<spectrum name>/
    results.csv      the ranked table, best match first
    config.json      every effective setting, enough to reproduce the run
    binned.txt       the binned observation
    bestfit_1.pdf    the best fits, if you asked for plots
    bestfit_2.pdf
```

A directory that already holds a result is not overwritten unless you pass
`--overwrite`.

![Output](ZTF18abokyfk_20180925_P60_v1_10_0.png)

The plot shows the input object in red and the combined SN + host galaxy
template in green. The legend gives the SN type, the host type, and the
percentage of the flux the SN template contributes; the redshift is above
the plot.

## Documentation

- **[Command line](docs/cli.md)** — every command and flag, with examples.
- **[Configuration](docs/configuration.md)** — every setting, what it means,
  and the Python API for passing them.
- **[Methods](docs/methods.md)** — what the fit actually does: the log
  wavelength grid, the extinction model, the weighting, the template bank,
  and the warnings you should not ignore.
- **[Development](docs/development.md)** — building an environment, running
  the tests, cutting a release.

## Requirements

Python 3.9+, with numpy, scipy, astropy, pandas, matplotlib, extinction,
PyAstronomy and tqdm. `pip install superfit` handles these.

Optional: `pip install superfit[speed]` adds `threadpoolctl`, which stops
worker processes contending over BLAS threads. Without it the fit is
correct, just slower.

Tested on both an older stack (numpy 1.24, scipy 1.8, astropy 5.2, pandas
1.5) and a current one (numpy 2.5, scipy 1.18, astropy 8.0, pandas 3.0),
producing byte-identical results on either.

## Citing

superfit is only as good as the template bank it uses. If you publish a
classification, say which bank version you used — `superfit bank status`
prints its checksum.
