"""Filesystem locations for the superfit template bank and metadata.

Historically these paths were hard-coded absolute strings pointing at one
developer's home directory, which made the package unrunnable anywhere else.
Everything is now resolved relative to the installed package, the current
working directory, or an explicit ``SUPERFIT_BANK_DIR`` override.

Resolution is deferred until something actually asks for a path. ``import
superfit`` must succeed on a machine with no template bank -- otherwise a
fresh ``pip install superfit`` cannot even be imported to read its own
documentation.
"""

import os
import sys

# Environment overrides, newest name first. NGSF_BANK_DIR is the name this
# used to have and is still honoured.
BANK_DIR_ENV_VARS = ("SUPERFIT_BANK_DIR", "NGSF_BANK_DIR")

# Names a bank, rather than giving its path: "legacy", "modern-curated".
# Lower precedence than SUPERFIT_BANK_DIR, which is a specific directory and
# so a more specific instruction.
BANK_NAME_ENV_VAR = "SUPERFIT_BANK"

BANK_URL = "https://www.wiserep.org/sites/default/files/supyfit_bank.zip"

# Directory containing this file, i.e. the installed superfit package.
PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

# One level up: the repository root for a source checkout.
REPO_DIR = os.path.dirname(PACKAGE_DIR)

# MJD-of-maximum table ships inside the package, so it needs no searching.
MJD_MAX_BRIGHTNESS_CSV = os.path.join(PACKAGE_DIR, "mjd_of_maximum_brightness.csv")

_cached_bank_dir = None


def _env_override():
    """An explicitly configured bank directory, if any.

    A variable that is set but points nowhere is an error rather than a
    reason to fall back: someone who names a bank means that bank, and
    quietly fitting against a different one is worse than stopping.
    """

    for var in BANK_DIR_ENV_VARS:
        override = os.environ.get(var)
        if not override:
            continue
        if not os.path.isdir(override):
            raise FileNotFoundError(
                "{} is set to {!r}, which is not a directory.".format(var, override)
            )
        return override
    return None


def _candidate_bank_dirs():
    """Places a template bank may live, in order of precedence."""

    yield os.path.join(os.getcwd(), "bank")
    yield os.path.join(REPO_DIR, "bank")
    yield os.path.join(PACKAGE_DIR, "bank")
    # Where `superfit bank install` puts it. Last, so a bank sitting next to
    # the spectra you are working on still wins -- but present, so a plain
    # `pip install superfit && superfit bank install` needs no further setup.
    yield _user_data_bank_dir()


def _user_data_bank_dir():
    """The per-user data directory ``superfit bank install`` writes to.

    Duplicated from superfit.bank rather than imported: this module is the
    one every other module depends on, and it must not start importing
    things back.
    """

    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~\\AppData\\Local")
    elif sys.platform == "darwin":
        base = os.path.expanduser("~/Library/Application Support")
    else:
        base = os.environ.get("XDG_DATA_HOME") or os.path.expanduser("~/.local/share")
    return os.path.join(base, "superfit", "bank")


def bank_dir_for_name(name):
    """The directory holding the bank called ``name``.

    Raises rather than falling back: someone who names a bank means that
    bank, and quietly fitting against a different one is worse than stopping.
    """

    from superfit.bank import UnknownBank, resolve_bank

    directory = resolve_bank(name)
    if directory is None:
        raise FileNotFoundError(
            "The bank named {!r} is not installed on this machine.\n\n"
            "Install it with:\n\n    superfit bank install {}\n\n"
            "or, for a bank that already exists on disk:\n\n"
            "    superfit bank install {} --from <directory>\n\n"
            "`superfit bank list` shows what is available.".format(name, name, name)
        )
    return directory


def find_bank_dir(use_cache=True, name=None):
    """Return the template bank directory, or raise with a useful message.

    ``name`` selects a bank by name and overrides everything else -- it is
    the fit saying which bank it wants, which beats any ambient setting.

    ``use_cache=False`` re-resolves without touching the cache. It used to
    overwrite it, so merely *inspecting* the bank -- ``superfit bank status``,
    ``superfit doctor`` -- silently discarded a directory the caller had
    chosen with :func:`set_bank_dir`, and the next fit ran against a
    different bank than the one just reported.
    """

    if name:
        return os.path.abspath(bank_dir_for_name(name))

    if use_cache and _cached_bank_dir is not None:
        return _cached_bank_dir

    def remember(path):
        global _cached_bank_dir
        path = os.path.abspath(path)
        if use_cache:
            _cached_bank_dir = path
        return path

    override = _env_override()
    if override is not None:
        return remember(override)

    named = os.environ.get(BANK_NAME_ENV_VAR)
    if named:
        return remember(bank_dir_for_name(named))

    tried = []
    for candidate in _candidate_bank_dirs():
        tried.append(candidate)
        if os.path.isdir(candidate):
            return remember(candidate)

    raise FileNotFoundError(
        "Could not locate the superfit template bank. Looked in:\n  "
        + "\n  ".join(tried)
        + "\n\nInstall it with:\n\n    superfit bank install\n\n"
        "or download it from "
        + BANK_URL
        + " and point SUPERFIT_BANK_DIR at it.\n\n"
        "`superfit bank list` shows the banks superfit knows by name."
    )


# The phase table's first column names the object; the rest describe its
# maximum. A table without the name column cannot be looked up in, whatever
# else it contains.
PHASE_TABLE_NAME_COLUMN = "Name"


def phase_table_is_usable(path):
    """Whether a phase table has the ``Name`` column callers key on.

    Returns ``(ok, reason)``. Only the header is read.

    This exists because a table can be present, non-empty and still
    unusable: dropping the name column leaves a file whose every row parses
    and whose every key is a maximum-light MJD, so nothing raises and no
    object ever matches. Silent, and it degrades classification rather than
    stopping it, which is the worst way for this to fail.
    """

    try:
        with open(path) as handle:
            header = handle.readline()
    except OSError as exc:
        return False, "cannot be read ({})".format(exc)

    if not header.strip():
        return False, "is empty"

    columns = [c.strip().strip('"') for c in header.split(",")]
    if columns[0] != PHASE_TABLE_NAME_COLUMN:
        return False, (
            "has no {!r} column -- its header is {!r}. Every lookup is by "
            "object name, so a table without that column silently matches "
            "nothing".format(PHASE_TABLE_NAME_COLUMN, header.strip())
        )

    return True, None


def mjd_max_brightness_csv(bank_dir=None, warn=True):
    """The phase table to use: the bank's own if it has one, else the package's.

    A bank that carries its own table knows the objects in it, which the copy
    shipped with the package cannot -- that one lists the 189 objects the
    legacy bank was built from, and a bank of 15561 templates from elsewhere
    needs its own. This is how such a table gets used, without a caller
    assigning over module globals to arrange it.

    A bank-local table that is not usable is refused rather than used: the
    package copy giving too few phases is a visible shortfall, while a
    malformed table giving none at all looks exactly like success.
    """

    if bank_dir is None:
        try:
            bank_dir = find_bank_dir()
        except FileNotFoundError:
            return MJD_MAX_BRIGHTNESS_CSV

    local = os.path.join(bank_dir, "mjd_of_maximum_brightness.csv")
    if not os.path.isfile(local):
        return MJD_MAX_BRIGHTNESS_CSV

    ok, reason = phase_table_is_usable(local)
    if ok:
        return local

    if warn:
        import warnings

        warnings.warn(
            "Ignoring the phase table shipped with this bank: {} {}. Falling "
            "back to the copy inside the package, which only knows the "
            "objects the legacy bank was built from -- so templates from "
            "this bank will report an unknown phase. Fix the table in the "
            "bank to get phases back.".format(local, reason),
            RuntimeWarning,
            stacklevel=2,
        )

    return MJD_MAX_BRIGHTNESS_CSV


def set_bank_dir(path):
    """Use ``path`` as the template bank for the rest of this process."""

    global _cached_bank_dir

    if not os.path.isdir(path):
        raise FileNotFoundError("No such template bank directory: {}".format(path))
    _cached_bank_dir = os.path.abspath(path)
    return _cached_bank_dir


def bank_is_available():
    """True when a template bank can be located, without raising if not."""

    try:
        find_bank_dir()
    except FileNotFoundError:
        return False
    return True


# Every helper below takes an optional ``bank_dir``, and a fit always passes
# its own. Falling back to the process-wide search is for interactive use and
# for the command line, where there is only one bank in play.
#
# A fit must not read that global. Two Superfits with different banks exist
# at once often enough -- a batch loop, a notebook, a comparison between two
# banks -- and a bank pinned process-wide is one that the *next* constructor
# moves out from under the fit already holding it. That fit then reads one
# bank's supernovae against another's galaxies and reports the name it was
# told, which is worse than either bank alone. The fit lock does not help:
# the bank is chosen during construction, outside it.


def original_resolution_dir(bank_dir=None):
    return os.path.join(bank_dir or find_bank_dir(), "original_resolution")


def binnings_dir(bank_dir=None):
    return os.path.join(bank_dir or find_bank_dir(), "binnings")


def binning_dir(resolution, bank_dir=None):
    """Directory holding the bank pre-binned to ``resolution`` Angstroms."""

    return os.path.join(binnings_dir(bank_dir), "{}A".format(resolution))


def sne_dir(resolution=None, bank_dir=None):
    """Supernova template directory, binned or at original resolution."""

    if resolution is None:
        return os.path.join(original_resolution_dir(bank_dir), "sne")
    return os.path.join(binning_dir(resolution, bank_dir), "sne")


def gal_dir(resolution=None, bank_dir=None):
    """Host galaxy template directory, binned or at original resolution."""

    if resolution is None:
        return os.path.join(original_resolution_dir(bank_dir), "gal")
    return os.path.join(binning_dir(resolution, bank_dir), "gal")


# Module-level constants stay readable as attributes (paths.BANK_DIR), but
# resolve on access rather than at import.
_LAZY_ATTRS = {
    "BANK_DIR": find_bank_dir,
    "ORIGINAL_RESOLUTION_DIR": original_resolution_dir,
    "BINNINGS_DIR": binnings_dir,
}


def __getattr__(name):
    if name in _LAZY_ATTRS:
        return _LAZY_ATTRS[name]()
    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name))


def __dir__():
    return sorted(list(globals()) + list(_LAZY_ATTRS))
