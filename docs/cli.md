# The `superfit` command

Four commands: `fit`, `bank`, `doctor`, `config`. Everything has a default,
so the common case is one line.

```bash
superfit --help
superfit fit --help
```

---

## `superfit fit`

Fit one spectrum against the template bank.

```bash
superfit fit spectrum.flm --z 0.127
```

That is the whole thing. Everything below is optional.

### Redshift

| Flag | Meaning |
| --- | --- |
| `--z Z` | Fit at exactly this redshift. The usual case. |
| `--scan-z BEGIN END` | Search this range instead of fixing z. |
| `--z-step STEP` | Step for the scan. Default 0.01. |

`--z` and `--scan-z` are mutually exclusive, and `--z-step` only means
something with `--scan-z`; either mistake is refused rather than ignored.

`--scan-z` turns host-line masking off, and says so. The mask places galaxy
emission lines at the object's redshift, and a scan has no single one. You
can still mask explicitly at a fixed redshift with `--z`.

```bash
superfit fit spectrum.flm --scan-z 0.0 0.3 --z-step 0.005
```

### Reading the spectrum

Ascii, csv and FITS are all read. Columns and units are worked out from the
file where possible, and can be stated when they cannot.

| Flag | Meaning |
| --- | --- |
| `--columns W,F[,E]` | Which columns hold wavelength, flux and error. Names for csv and FITS, positions for ascii. |
| `--wavelength-unit UNIT` | `AA`, `nm`, `micron`, `mm`, `log10(AA)`. Default: whatever the file says, else Angstroms. |
| `--hdu N` | Which FITS HDU to read, by index or name. Default: try each. |
| `--name NAME` | Name for the spectrum. Default: the filename. |

```bash
superfit fit spec.fits --z 0.1 --columns WAVE,FLUX,IVAR
superfit fit spec.fits --z 0.1 --wavelength-unit nm --hdu SPECTRUM
superfit fit wide.txt  --z 0.1 --columns 1,2,4
```

`ivar` named as the third column is understood as inverse variance and
converted to an uncertainty. See [configuration](configuration.md#reading-the-spectrum)
for the full list of names that are recognised automatically.

### The fit

| Flag | Meaning |
| --- | --- |
| `--resolution A` | Binning resolution in Angstroms. Default 10. |
| `--error-model {sg,linear,included}` | How the error spectrum is obtained. Default `sg`. |
| `--profile {legacy,modern}` | Named set of scientific defaults. Default `legacy`. |
| `--weighted-solve` | Shorthand for `--profile modern`. |
| `--rv R_V` | Total-to-selective extinction ratio. Default 3.1. |
| `--av-range LOW HIGH` | Extinction range to search. Default -2 to 2. |
| `--epochs LOW HIGH` | Restrict templates to a phase window, in days. Equal values mean no restriction. |
| `--no-mask-galaxy-lines` | Do not mask host emission lines. |
| `--no-mask-telluric` | Do not mask the telluric A band. |
| `--n-cores N` | Worker processes. 0 (default) means auto. |

`--error-model included` uses the observation's own uncertainty column, and
fails clearly if the file does not have one.

`legacy` and `modern` differ in one thing: whether the amplitude solve is
weighted consistently with the chi2. `legacy` is the default so that
published results reproduce. See
[methods](methods.md#a-note-on-the-weighting).

### Output

| Flag | Meaning |
| --- | --- |
| `--output DIR`, `-o DIR` | Where to create the run directory. Default: here. |
| `--overwrite` | Replace a run directory that already holds results. |
| `--plots N` | Plot the N best fits. Default 0. |
| `--show` | Display the plots as well as saving them. Implies `--plots 1`. |
| `--png` | Save plots as PNG rather than PDF. |
| `--quiet` | Print only the best match. |

A fit writes `<output_dir>/<spectrum name>/` and refuses to overwrite one
that already holds a `results.csv`. This is deliberate: a re-run aimed at
the wrong place used to destroy the earlier answer silently.

```bash
superfit fit spectrum.flm --z 0.127 --plots 3 --output results/
superfit fit spectrum.flm --z 0.127 --overwrite
```

### Using a parameter file

For runs that need to be reproducible, or that set more than a handful of
things, keep the settings in JSON:

```bash
superfit fit spectrum.flm --config parameters.json
superfit fit --config parameters.json          # the file names the spectrum
```

Command-line flags override the file. `superfit parameters.json` — the form
that predates the subcommands — still works and means the second line above.

The complete effective configuration is written to `config.json` in the run
directory whatever route you took, so a one-flag command is as reproducible
as a full JSON file.

---

## `superfit bank`

### `superfit bank install`

Download, verify and unpack the template bank.

```bash
superfit bank install
```

It goes into the platform's per-user data directory
(`~/.local/share/superfit/bank` on Linux, `~/Library/Application
Support/superfit/bank` on macOS, `%LOCALAPPDATA%\superfit\bank` on Windows),
which `superfit` searches automatically. No environment variable needed.

| Flag | Meaning |
| --- | --- |
| `--dir PATH` | Install somewhere else. |
| `--archive PATH` | Use a zip already on disk instead of downloading. |
| `--url URL` | Download from somewhere else. |
| `--overwrite` | Replace an existing installation. |
| `--no-verify` | Skip the checksum check. |
| `--quiet` | No progress bars. |
| `--no-pack` | Skip the packing step below. |

The archive is checked against a checksum pinned in the package, so a
truncated or substituted download is caught immediately rather than showing
up as a strange fit result later. The unpack is staged and moved into place,
so an interrupted install cannot leave a half-populated bank behind.

On a machine with no direct internet access, download the archive elsewhere
and point at it:

```bash
superfit bank install --archive /path/to/supyfit_bank.zip
```

### `superfit bank pack`

Concatenate each template directory into a single `.npy`, so a fit opens two
files instead of a thousand.

```bash
superfit bank pack
```

`bank install` does this for you. Run it by hand for a bank installed some
other way — unzipped, copied from a colleague, or pointed at with
`$SUPERFIT_BANK_DIR` — or after editing templates.

| Flag | Meaning |
| --- | --- |
| `--dir PATH` | Pack this bank instead of the one a fit would use. |
| `--quiet` | Do not list what was packed. |

A fit reads roughly a thousand template files. Parsing them is cheap;
*opening* them is not, because on a parallel filesystem every open is a round
trip to a metadata server. Measured on Perlmutter's CFS, the first fit of the
day spent 10.9 s getting the 10 Å bank off disk, of which 0.67 s was parsing.
From the pack the same templates arrive in about 40 ms.

The pack goes next to the bank when that directory is writable, so a shared
bank is packed once for everyone, and in the per-user data directory when it
is not. It stores the raw float64 each template parses to and nothing derived
from it, so results are bit-identical either way, and no fit setting — the
epoch window, the template selection, `mask_galaxy_lines` — can make it stale.
Only the bank can, and a pack whose directory listing no longer matches is
ignored in favour of the text. Packing is entirely optional: without it,
fits read the text files exactly as before.

### `superfit bank status`

Say which bank a fit would use, where it came from, and whether it is whole.

```bash
$ superfit bank status
Template bank: /home/you/.local/share/superfit/bank
  located via: found on the search path
  contents:    35 supernova types, 12 galaxy templates
  packed:      all 8 directories, in /home/you/.local/share/superfit/bank/packed
  installed:   2026-08-24T04:31:07+00:00
  sha256:      42689295b35b77568e9f831925344eea42ebe9c3bdfe61b027df4e8b2c367ce8
```

Exits non-zero if there is no bank, or if the one it found is incomplete.

Where superfit looks, in order:

1. `$SUPERFIT_BANK_DIR` (or the older `$NGSF_BANK_DIR`)
2. `./bank`
3. the repository root, for a source checkout
4. next to the installed package
5. the per-user data directory `bank install` writes to

A bank sitting next to the spectra you are working on therefore wins over
the installed one, which is usually what you want when comparing banks.

---

## `superfit doctor`

Check the whole installation and say what is wrong with it.

```bash
$ superfit doctor
superfit doctor

superfit 0.1.0 from /home/you/envs/superfit/lib/python3.13/site-packages/superfit
python 3.13.12

[ok  ] python version
[ok  ] numpy  -- 2.5.2
...
[ok  ] template bank  -- /home/you/.local/share/superfit/bank (35 SN types)
[ok  ] install directory  -- /home/you/.local/share/superfit/bank
[ok  ] matplotlib backend  -- Agg
[note] Agg cannot display plots; --show will do nothing, but saved plots are unaffected

Everything checks out, with 1 note(s).
```

Checks the Python version, every dependency and its version, the optional
ones, the template bank, whether the install directory is writable, and
whether the matplotlib backend can actually display anything. Exits non-zero
if any check fails; notes do not fail it.

Run this first when something is not working.

---

## `superfit config`

### `superfit config create`

Write a parameter file that explains itself.

```bash
superfit config create parameters.json
```

JSON has no comment syntax, so each setting is preceded by an
underscore-prefixed sibling key holding its explanation; `superfit` ignores
any key starting with `_`.

```json
{
    "_z_exact": "The object's redshift.",
    "z_exact": 0.0,
    "_resolution": "Binning resolution in Angstroms. 10 and 30 use the pre-binned bank.",
    "resolution": 10
}
```

| Flag | Meaning |
| --- | --- |
| `--full` | Include every setting, not just the common ones. |
| `--profile P` | Start from a named profile. |
| `--force` | Overwrite an existing file. |

### `superfit config show`

Print the configuration a fit would run with, after defaults, profile and
file are merged.

```bash
superfit config show
superfit config show parameters.json
superfit config show --profile modern
```

Useful for answering "what is this actually going to do" without running it.

---

## Exit codes

| Code | Meaning |
| --- | --- |
| 0 | Success. |
| 1 | A real error: bad setting, missing file, refused overwrite, missing bank. |
| 2 | Nothing to do — a command or action was not given, and help was printed. |
| 130 | Interrupted. |
