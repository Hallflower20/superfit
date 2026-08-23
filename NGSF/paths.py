"""Filesystem locations for the NGSF template bank and metadata.

Historically these paths were hard-coded absolute strings pointing at one
developer's home directory, which made the package unrunnable anywhere else.
Everything is now resolved relative to the installed package, the current
working directory, or an explicit ``NGSF_BANK_DIR`` override.
"""

import os

# Directory containing this file, i.e. the installed NGSF package.
PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

# One level up: the repository root for a source checkout.
REPO_DIR = os.path.dirname(PACKAGE_DIR)


def _candidate_bank_dirs():
    """Places a template bank may live, in order of precedence."""

    override = os.environ.get("NGSF_BANK_DIR")
    if override:
        yield override

    yield os.path.join(os.getcwd(), "bank")
    yield os.path.join(REPO_DIR, "bank")
    yield os.path.join(PACKAGE_DIR, "bank")


def find_bank_dir():
    """Return the template bank directory, or raise with a useful message."""

    tried = []
    for candidate in _candidate_bank_dirs():
        tried.append(candidate)
        if os.path.isdir(candidate):
            return os.path.abspath(candidate)

    raise FileNotFoundError(
        "Could not locate the NGSF template bank. Looked in:\n  "
        + "\n  ".join(tried)
        + "\n\nDownload it from "
        "https://www.wiserep.org/sites/default/files/supyfit_bank.zip "
        "and unzip it next to this package, or set NGSF_BANK_DIR to point at it."
    )


BANK_DIR = find_bank_dir()

ORIGINAL_RESOLUTION_DIR = os.path.join(BANK_DIR, "original_resolution")
BINNINGS_DIR = os.path.join(BANK_DIR, "binnings")

# MJD-of-maximum table ships inside the package.
MJD_MAX_BRIGHTNESS_CSV = os.path.join(PACKAGE_DIR, "mjd_of_maximum_brightness.csv")


def binning_dir(resolution):
    """Directory holding the bank pre-binned to ``resolution`` Angstroms."""

    return os.path.join(BINNINGS_DIR, "{}A".format(resolution))


def sne_dir(resolution=None):
    """Supernova template directory, binned or at original resolution."""

    if resolution is None:
        return os.path.join(ORIGINAL_RESOLUTION_DIR, "sne")
    return os.path.join(binning_dir(resolution), "sne")


def gal_dir(resolution=None):
    """Host galaxy template directory, binned or at original resolution."""

    if resolution is None:
        return os.path.join(ORIGINAL_RESOLUTION_DIR, "gal")
    return os.path.join(binning_dir(resolution), "gal")
