# Development

## Building an environment from source

A self-contained virtualenv, which keeps the dependencies out of your user
site-packages (`~/.local`) where they can collide with other projects:

```bash
python -m venv /path/to/envs/superfit        # NOT --system-site-packages
source /path/to/envs/superfit/bin/activate
pip install -e ".[dev]"
```

Building *without* `--system-site-packages` is deliberate. A venv that
inherits the system packages also picks up `~/.local`, and a stale package
there shadows the system one — the usual symptom is an import error from
inside matplotlib that has nothing to do with your code.

Then get the template bank:

```bash
superfit bank install
```

In a source checkout the bank is also found at `./bank`, so unzipping the
archive in the repository root works too, and takes precedence.

### On NERSC (Perlmutter)

An environment built this way lives at

```
/global/cfs/cdirs/desicollab/users/xhall/envs/superfit
```

Activate and go, from any working directory:

```bash
source /global/cfs/cdirs/desicollab/users/xhall/envs/superfit/bin/activate
superfit fit /path/to/spectrum.flm --z 0.127
```

It is an editable install of the checkout at
`/global/cfs/cdirs/desicollab/users/xhall/GitHub/superfit`, so the template
bank is found automatically and `git pull` takes effect without
reinstalling. Moving or renaming that checkout breaks the link; re-run
`pip install -e .` if you do.

## Tests

```bash
pytest                      # everything, ~35s
pytest -m "not endtoend"    # skip the tests that need the template bank
```

The suite:

| File | What it covers |
| --- | --- |
| `test_core_kernel.py` | The fitting linear algebra, pinned against an independent loop-based reference so optimisation work can be checked against something other than itself. |
| `test_robustness.py` | Extreme inputs — cosmic rays, NaN gaps, interpolated-flat stretches, spectra straddling zero — through the error and binning routines. |
| `test_loggrid.py` | The log wavelength grid and the redshift-as-translation machinery. |
| `test_spectrum.py`, `test_io.py` | Reading spectra: ascii, csv, FITS tables and images, units, inverse variance. |
| `test_config.py` | Defaults, aliases, profiles, booleans, and the validation that rejects unusable combinations up front. |
| `test_cli.py` | What each command-line flag turns into. A flag that maps to the wrong setting is a fit that answers a different question. |
| `test_output.py` | Run directories, the overwrite guard, and the `FitResult` object. |
| `test_bank.py` | Installing the bank: checksums, staging, hostile archives. Builds its own tiny archive; never touches the network. |
| `test_api.py` | Fits the same data as a file, as arrays, as a `Spectrum`, as a raw array and as FITS, and asserts they all agree. Also that two fits in one process cannot disturb each other. |
| `test_regression.py` | The real pipeline over the real bank, against a committed golden result. |

`test_regression.py` separates "which template won and in what order" from
"what were the numbers", so a failure says which kind of change happened.

If you change the science on purpose, regenerate the golden file in the same
commit and say why:

```bash
python superfit/tests/test_regression.py --regenerate
```

## Things to be careful about

**The template bank must never end up in the sdist.** It is a 74 MB download,
not source. CI greps the built distribution and fails if it is there.

**Settings belong to a fit, not to the process.** `Parameters` is built by
and owned by the `Superfit` that needs it, and is frozen once built. This was
not always true, and the consequence was that building a second fit silently
retuned the first — see the commit history. Anything that reaches for
process-global configuration is reintroducing that bug.

**Output paths are `pathlib.Path`s, assembled by `superfit.output`.** Not
string prefixes concatenated onto filenames. That is what
`--out results` writing `resultsSN2021urb.csv` into the parent directory
looked like.

**Close your figures.** Anything that draws with pyplot must close what it
drew, or a batch of a few hundred fits accumulates every figure it has
made.

## Building a release

```bash
python -m build
python -m twine check dist/*
```

Check the sdist does not contain the bank:

```bash
tar tzf dist/superfit-*.tar.gz | grep -i bank/ && echo "BAD" || echo "ok"
```

## Updating the pinned bank checksum

`superfit.bank.BANK_SHA256` pins the sha256 of the published archive, so a
truncated or substituted download is caught at install time. If the
published bank is updated, verify the new archive and update the constant in
the same commit:

```bash
python -c "
from superfit import bank
import tempfile, os
p = os.path.join(tempfile.mkdtemp(), 'b.zip')
bank.download(bank.BANK_URL, p)
print(bank.sha256_of(p))
"
```

Note that the server hosting the archive answers 403 to urllib's default
User-Agent, which is why `bank.download` sets one. A plain
`urllib.request.urlopen(url)` will appear to be a network failure.
