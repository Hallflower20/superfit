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

# Environment overrides, newest name first. NGSF_BANK_DIR is the name this
# used to have and is still honoured.
BANK_DIR_ENV_VARS = ("SUPERFIT_BANK_DIR", "NGSF_BANK_DIR")

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


def find_bank_dir(use_cache=True):
    """Return the template bank directory, or raise with a useful message."""

    global _cached_bank_dir

    if use_cache and _cached_bank_dir is not None:
        return _cached_bank_dir

    override = _env_override()
    if override is not None:
        _cached_bank_dir = os.path.abspath(override)
        return _cached_bank_dir

    tried = []
    for candidate in _candidate_bank_dirs():
        tried.append(candidate)
        if os.path.isdir(candidate):
            _cached_bank_dir = os.path.abspath(candidate)
            return _cached_bank_dir

    raise FileNotFoundError(
        "Could not locate the superfit template bank. Looked in:\n  "
        + "\n  ".join(tried)
        + "\n\nDownload it from "
        + BANK_URL
        + "\nand unzip it into the working directory, or point "
        "SUPERFIT_BANK_DIR at it."
    )


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


def original_resolution_dir():
    return os.path.join(find_bank_dir(), "original_resolution")


def binnings_dir():
    return os.path.join(find_bank_dir(), "binnings")


def binning_dir(resolution):
    """Directory holding the bank pre-binned to ``resolution`` Angstroms."""

    return os.path.join(binnings_dir(), "{}A".format(resolution))


def sne_dir(resolution=None):
    """Supernova template directory, binned or at original resolution."""

    if resolution is None:
        return os.path.join(original_resolution_dir(), "sne")
    return os.path.join(binning_dir(resolution), "sne")


def gal_dir(resolution=None):
    """Host galaxy template directory, binned or at original resolution."""

    if resolution is None:
        return os.path.join(original_resolution_dir(), "gal")
    return os.path.join(binning_dir(resolution), "gal")


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
