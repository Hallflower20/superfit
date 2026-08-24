"""Where one fit's output goes, and what it produced.

Output used to be a string prefix that four call sites concatenated file
names onto. That has three failure modes and all of them bite in practice:
``--out results`` (no trailing slash) wrote ``resultsSN2021urb.csv`` into the
parent directory; a directory that did not exist yet made the fit run to
completion and then fail to save it; and a second fit of the same spectrum
overwrote the first without a word.

A fit now writes into its own directory, created up front, with fixed names
inside it::

    <output_dir>/<spectrum name>/
        results.csv        the ranked table
        config.json        every effective setting, for reproducing the run
        binned.txt         the binned observation
        bestfit_1.pdf      the best fits, if plotting was asked for
        bestfit_2.pdf

Fixed names make a run directory self-describing and scriptable. Refusing to
overwrite one, unless overwriting was asked for, makes a result you already
have safe from a re-run you meant to point somewhere else.
"""

import os
import re
from pathlib import Path

RESULTS_CSV = "results.csv"
USED_CONFIG_JSON = "config.json"
BINNED_TXT = "binned.txt"
PLOT_STEM = "bestfit_"

# Files a previous run of this fit would have left. Used only to clear a
# directory the caller explicitly asked to overwrite -- never anything wider,
# and never the directory itself.
_OUR_FILES = (RESULTS_CSV, USED_CONFIG_JSON, BINNED_TXT)
_OUR_PLOTS = re.compile(r"^{}\d+\.(pdf|png)$".format(re.escape(PLOT_STEM)))

# Anything that would let a spectrum name escape its parent directory or
# confuse a shell. Spectrum names come from filenames and FITS headers, so
# they are not always tame.
_UNSAFE = re.compile(r"[^A-Za-z0-9._+-]+")


class OutputExistsError(FileExistsError):
    """Raised rather than overwrite results that are already on disk."""


def safe_name(name):
    """A directory-safe version of a spectrum name."""

    cleaned = _UNSAFE.sub("_", str(name).strip()).strip("._")
    return cleaned or "spectrum"


class RunDirectory:
    """One fit's output directory: created, checked, and named.

    Parameters
    ----------
    base : str, Path, or None
        Where run directories are created. ``None`` or empty means the
        current working directory.
    name : str
        The spectrum's name; becomes the directory name.
    overwrite : bool
        Replace results already in the directory. Without it, a directory
        that already holds a results table is an error rather than a
        silent overwrite.
    """

    def __init__(self, base, name, overwrite=False):
        base = Path(base).expanduser() if base else Path.cwd()
        self.name = safe_name(name)
        self.path = base.expanduser() / self.name
        self._prepare(bool(overwrite))

    def _prepare(self, overwrite):
        results = self.path / RESULTS_CSV

        if results.exists() and not overwrite:
            raise OutputExistsError(
                "{} already holds a fit of {!r}. Pass overwrite=True (or "
                "--overwrite) to replace it, or point --output somewhere "
                "else.".format(self.path, self.name)
            )

        try:
            self.path.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise OutputExistsError(
                "Could not create the output directory {}: {}".format(self.path, exc)
            )

        if overwrite:
            self._clear()

    def _clear(self):
        """Remove what a previous run of this fit left, and nothing else.

        Plots are numbered by rank, so a re-run that asks for fewer of them
        would otherwise leave the tail of the old run behind, mixed in with
        the new one and indistinguishable from it.
        """

        for entry in self.path.iterdir():
            if entry.is_file() and (
                entry.name in _OUR_FILES or _OUR_PLOTS.match(entry.name)
            ):
                entry.unlink()

    # -- the files themselves ---------------------------------------------

    @property
    def results_csv(self):
        return self.path / RESULTS_CSV

    @property
    def used_config_json(self):
        return self.path / USED_CONFIG_JSON

    @property
    def binned_txt(self):
        return self.path / BINNED_TXT

    def plot(self, rank, png=False):
        """Path for the plot of the ``rank``-th best fit, counting from 1."""

        return self.path / "{}{}.{}".format(PLOT_STEM, rank, "png" if png else "pdf")

    def __fspath__(self):
        return os.fspath(self.path)

    def __str__(self):
        return str(self.path)

    def __repr__(self):
        return "RunDirectory({!r})".format(str(self.path))


class FitResult:
    """The ranked table, plus where everything a fit produced was written.

    Indexing, iteration and attribute access fall through to the underlying
    ``pandas.DataFrame``, so ``result["CHI2/dof"]`` and ``len(result)`` mean
    what they used to when ``run()`` returned the frame itself.
    """

    def __init__(self, results, directory, artifacts=()):
        self.results = results
        self.directory = Path(str(directory))
        #: Every file this fit wrote, in the order it wrote them.
        self.artifacts = [Path(str(p)) for p in artifacts]

    @property
    def best(self):
        """The top-ranked match, as a pandas Series."""

        return self.results.iloc[0]

    @property
    def plots(self):
        return [p for p in self.artifacts if _OUR_PLOTS.match(p.name)]

    def artifact(self, name):
        """The written file called ``name``, or None if it was not written."""

        for path in self.artifacts:
            if path.name == name:
                return path
        return None

    # -- behave like the DataFrame it wraps --------------------------------

    def __getitem__(self, key):
        return self.results[key]

    def __len__(self):
        return len(self.results)

    def __iter__(self):
        return iter(self.results)

    def __contains__(self, key):
        return key in self.results

    def __getattr__(self, name):
        # Only reached for names this class does not define, so `results`
        # itself can never recurse -- it is set in __init__.
        try:
            return getattr(self.__dict__["results"], name)
        except KeyError:  # pragma: no cover - during unpickling
            raise AttributeError(name)

    def __repr__(self):
        if len(self) == 0:
            return "<FitResult: no surviving templates, in {}>".format(self.directory)
        best = self.best
        return "<FitResult: {} at z={:g}, chi2/dof={:.3g}, {} rows, in {}>".format(
            best["SN"], best["Z"], best["CHI2/dof"], len(self), self.directory
        )
