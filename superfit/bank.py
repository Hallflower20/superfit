"""Getting the template bank onto the machine.

The bank is a 74 MB archive of supernova and host-galaxy spectra that is not
shipped with the package. Installing it used to be a paragraph of README:
curl the zip, unzip it, work out where to put it, export an environment
variable. Every one of those steps is a place to stop, and the failure at
the end -- "Could not locate the superfit template bank" -- arrives long
after the step that caused it.

``superfit bank install`` does the whole thing: download, check, unpack into
a standard per-user data directory that :mod:`superfit.paths` already knows
to look in, and record what was installed so ``superfit bank status`` can
say whether it is the bank it should be.
"""

import hashlib
import json
import os
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

from superfit import paths

BANK_URL = paths.BANK_URL

# sha256 of the archive published at BANK_URL. Recorded so that a truncated
# or tampered download is caught here rather than as a strange fit result
# three weeks later. Pass expected_sha256=None to install without checking.
BANK_SHA256 = "42689295b35b77568e9f831925344eea42ebe9c3bdfe61b027df4e8b2c367ce8"

# The server behind BANK_URL answers 403 to urllib's default User-Agent and
# 200 to a browser's. Not a policy worth arguing with; just send one.
USER_AGENT = "superfit/{} (+https://github.com/Hallflower20/superfit)".format(
    getattr(sys.modules.get("superfit"), "__version__", "0")
)

MANIFEST_NAME = "INSTALLED.json"

# What a usable bank has in it.
REQUIRED_SUBDIRS = ("original_resolution", "binnings")

# Archive members that are packaging debris rather than templates.
_JUNK = ("__MACOSX/", ".DS_Store")


class BankError(RuntimeError):
    """Raised when the bank cannot be installed or is not what it claims."""


def user_data_dir():
    """The per-user data directory this platform expects applications to use."""

    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~\\AppData\\Local")
    elif sys.platform == "darwin":
        base = os.path.expanduser("~/Library/Application Support")
    else:
        base = os.environ.get("XDG_DATA_HOME") or os.path.expanduser("~/.local/share")
    return Path(base) / "superfit"


def default_install_dir():
    """Where ``superfit bank install`` puts the bank if not told otherwise."""

    return user_data_dir() / "bank"


# Banks superfit knows by name.
#
# There is more than one now, and they differ scientifically rather than
# incidentally: the legacy bank is what published superfit results were
# produced against, and the modern ones are a larger, provenance-tracked
# rebuild with DESI DR1 host galaxies. Which one a fit used is part of the
# result, so it is a named choice recorded in the run's config rather than
# an environment variable someone set months ago.
#
# `url` is None for a bank that is not published yet: those are installed
# from a directory with `superfit bank install <name> --from PATH`.
KNOWN_BANKS = {
    "legacy": {
        "summary": "The published superfit bank",
        "description": (
            "1007 supernova templates and 12 host galaxies. What superfit "
            "has always shipped, and what published classifications were "
            "produced against."
        ),
        "url": BANK_URL,
        "sha256": BANK_SHA256,
    },
    "modern-curated": {
        "summary": "Modern empirical bank, balanced",
        "description": (
            "2474 supernova templates and 128 DESI DR1 host galaxies, "
            "curated for balance across type and phase. The recommended "
            "modern default."
        ),
        "url": None,
        "sha256": None,
    },
    "modern": {
        "summary": "Modern empirical bank, full",
        "description": (
            "15561 supernova templates and 128 DESI DR1 host galaxies: "
            "every template that passed quality cuts, unbalanced. Slower, "
            "and biased towards whichever types were observed most."
        ),
        "url": None,
        "sha256": None,
    },
}

# What `superfit bank install` and a fit with no bank named will use.
DEFAULT_BANK = "legacy"

# Where the name -> directory registry lives. A bank installed from a
# directory is recorded here rather than copied: the modern banks are 0.4 and
# 1.3 GB, and on a shared filesystem one copy per user is a waste of the
# quota it comes out of.
REGISTRY_NAME = "banks.json"


class UnknownBank(BankError):
    """Raised when a bank is asked for by a name superfit does not know."""


def registry_path():
    """The file recording which named banks are installed where."""

    return user_data_dir() / REGISTRY_NAME


def read_registry():
    """Name -> directory, for every bank registered on this machine."""

    try:
        data = json.loads(registry_path().read_text())
    except (OSError, ValueError):
        return {}

    if not isinstance(data, dict):
        return {}
    return {str(k): str(v) for k, v in data.get("banks", {}).items()}


def write_registry(banks):
    """Record the name -> directory map, atomically."""

    path = registry_path()
    path.parent.mkdir(parents=True, exist_ok=True)

    scratch = path.with_name(path.name + ".tmp")
    scratch.write_text(json.dumps({"banks": dict(banks)}, indent=2, sort_keys=True))
    scratch.replace(path)
    return path


def register(name, directory):
    """Point ``name`` at ``directory`` for future fits."""

    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise BankError("No such bank directory: {}".format(directory))

    banks = read_registry()
    banks[name] = str(directory)
    write_registry(banks)
    return directory


def unregister(name):
    """Forget ``name``. Does not delete anything on disk."""

    banks = read_registry()
    removed = banks.pop(name, None)
    if removed is not None:
        write_registry(banks)
    return removed


def bank_install_dir(name):
    """Where ``superfit bank install <name>`` unpacks a downloaded bank.

    The legacy bank keeps the unnamed location it has always used, so an
    existing installation is still found after upgrading.
    """

    if name == DEFAULT_BANK:
        return default_install_dir()
    return user_data_dir() / "banks" / name


def resolve_bank(name):
    """The directory for a named bank, or None if it is not installed.

    The registry wins over the install directory, so a bank registered from
    a directory on a shared filesystem is used in preference to a stale copy
    someone unpacked earlier.
    """

    if name not in KNOWN_BANKS:
        registered = read_registry().get(name)
        if registered is None:
            raise UnknownBank(
                "No bank named {!r}. Known names: {}. Any other name has to "
                "be registered first with `superfit bank install <name> "
                "--from <directory>`.".format(
                    name, ", ".join(sorted(KNOWN_BANKS))
                )
            )
        return registered if os.path.isdir(registered) else None

    registered = read_registry().get(name)
    if registered and os.path.isdir(registered):
        return registered

    installed = bank_install_dir(name)
    return str(installed) if installed.is_dir() else None


def installed_banks():
    """Every bank this machine can use, as ``{name: directory}``."""

    found = {}

    for name in KNOWN_BANKS:
        directory = None
        try:
            directory = resolve_bank(name)
        except BankError:
            pass
        if directory:
            found[name] = directory

    for name, directory in read_registry().items():
        if os.path.isdir(directory):
            found[name] = directory

    return found

def sha256_of(path, block=1 << 20):
    """Hex sha256 of a file, read in chunks so a 74 MB zip is not held in RAM."""

    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(block), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url, destination, quiet=False):
    """Fetch ``url`` to ``destination``, showing progress."""

    from urllib.error import HTTPError, URLError
    from urllib.request import Request, urlopen

    request = Request(url, headers={"User-Agent": USER_AGENT})

    try:
        response = urlopen(request)
    except HTTPError as exc:
        raise BankError(
            "The template bank download returned HTTP {} for {}. The archive "
            "may have moved; see the project README, or download it by hand "
            "and pass --archive.".format(exc.code, url)
        )
    except URLError as exc:
        raise BankError(
            "Could not reach {}: {}. If this machine has no direct internet "
            "access, download the archive elsewhere and pass "
            "--archive.".format(url, exc.reason)
        )

    total = int(response.headers.get("Content-Length") or 0)

    with response, open(destination, "wb") as out:
        if quiet:
            shutil.copyfileobj(response, out)
        else:
            from tqdm import tqdm

            with tqdm(
                total=total or None, unit="B", unit_scale=True, desc="template bank"
            ) as bar:
                for chunk in iter(lambda: response.read(1 << 20), b""):
                    out.write(chunk)
                    bar.update(len(chunk))

    size = os.path.getsize(destination)
    if total and size != total:
        raise BankError(
            "Download is {} bytes but the server said {}. The transfer was "
            "cut short; try again.".format(size, total)
        )
    return destination


def _members_to_extract(archive):
    """Template members of the archive, with any leading ``bank/`` stripped.

    Yields ``(member, relative_path)``. Refuses absolute paths and ``..``
    segments so a hostile archive cannot write outside the target directory.
    """

    for member in archive.infolist():
        name = member.filename
        if any(part in name for part in _JUNK) or name.endswith("/"):
            continue

        relative = name[5:] if name.startswith("bank/") else name
        if not relative:
            continue

        parts = Path(relative).parts
        if os.path.isabs(relative) or ".." in parts:
            raise BankError(
                "Refusing to unpack {!r}: the archive tries to write outside "
                "the target directory.".format(name)
            )
        yield member, Path(*parts)


def unpack(archive_path, destination, quiet=False):
    """Unpack the bank archive into ``destination``."""

    destination = Path(destination)

    try:
        archive = zipfile.ZipFile(archive_path)
    except zipfile.BadZipFile:
        raise BankError(
            "{} is not a readable zip archive. If it was downloaded through a "
            "browser or a proxy, it may be an error page rather than the "
            "bank.".format(archive_path)
        )

    with archive:
        members = list(_members_to_extract(archive))
        if not members:
            raise BankError("{} contains no template files.".format(archive_path))

        progress = members
        if not quiet:
            from tqdm import tqdm

            progress = tqdm(members, desc="unpacking", unit="file")

        for member, relative in progress:
            target = destination / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(member) as source, open(target, "wb") as out:
                shutil.copyfileobj(source, out)

    missing = [d for d in REQUIRED_SUBDIRS if not (destination / d).is_dir()]
    if missing:
        raise BankError(
            "Unpacked {} but it has no {} directory. This does not look like "
            "the superfit template bank.".format(
                archive_path, " or ".join(missing)
            )
        )
    return destination


def write_manifest(destination, url, checksum, source):
    """Record what was installed, so ``bank status`` can describe it."""

    import datetime

    manifest = {
        "url": url,
        "sha256": checksum,
        "source": str(source),
        "installed_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "superfit_version": getattr(
            sys.modules.get("superfit"), "__version__", "unknown"
        ),
    }
    path = Path(destination) / MANIFEST_NAME
    path.write_text(json.dumps(manifest, indent=2))
    return manifest


def read_manifest(directory):
    """The install manifest for a bank directory, or None if there is none."""

    path = Path(directory) / MANIFEST_NAME
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def install(
    destination=None,
    url=BANK_URL,
    archive=None,
    overwrite=False,
    expected_sha256=BANK_SHA256,
    quiet=False,
):
    """Download, check and unpack the template bank.

    Parameters
    ----------
    destination : path, optional
        Where to install. Defaults to the per-user data directory, which
        :mod:`superfit.paths` searches, so no environment variable is needed.
    archive : path, optional
        A zip already on disk; skips the download. For machines with no
        direct internet access, and for re-installing without re-fetching.
    overwrite : bool
        Replace an installation that is already there.
    expected_sha256 : str or None
        Checksum the archive must have. None skips the check.

    Returns
    -------
    dict
        The manifest that was written.
    """

    destination = Path(destination) if destination else default_install_dir()

    if destination.is_dir() and any(destination.iterdir()) and not overwrite:
        raise BankError(
            "{} already exists. Pass --overwrite to replace it, or "
            "`superfit bank status` to see what is there.".format(destination)
        )

    scratch = None
    try:
        if archive is None:
            scratch = tempfile.mkdtemp(prefix="superfit-bank-")
            archive = Path(scratch) / "supyfit_bank.zip"
            if not quiet:
                print("Downloading the template bank from {}".format(url))
            download(url, archive, quiet=quiet)
            source = url
        else:
            archive = Path(archive).expanduser()
            if not archive.is_file():
                raise BankError("No such archive: {}".format(archive))
            source = str(archive)

        checksum = sha256_of(archive)
        if expected_sha256 and checksum != expected_sha256:
            raise BankError(
                "Checksum mismatch for {}:\n  expected {}\n  got      {}\n"
                "The archive is not the one this version of superfit knows "
                "about. Re-run to retry the download, or pass "
                "--no-verify if you mean to install a different "
                "bank.".format(archive, expected_sha256, checksum)
            )

        # Unpack beside the destination and move into place, so an
        # interrupted install cannot leave a half-populated bank that
        # `find_bank_dir` would then happily fit against.
        destination.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(
            tempfile.mkdtemp(prefix=".{}-".format(destination.name), dir=destination.parent)
        )
        try:
            unpack(archive, staging, quiet=quiet)
            write_manifest(staging, url, checksum, source)

            if destination.exists():
                shutil.rmtree(destination)
            staging.replace(destination)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
    finally:
        if scratch:
            shutil.rmtree(scratch, ignore_errors=True)

    paths.set_bank_dir(destination)
    return read_manifest(destination)


def install_from_directory(name, source, copy=False, quiet=False):
    """Make an existing bank directory available under ``name``.

    By default this records where the bank is rather than copying it. The
    modern banks are 0.4 and 1.3 GB; on a shared filesystem, one copy per
    user is a waste of the quota it comes out of, and a bank that everyone
    reads from one place is also a bank everyone packs once. ``copy=True``
    materialises it under the per-user data directory for a bank that lives
    somewhere temporary.

    Returns the directory the name now resolves to.
    """

    source = Path(source).expanduser()
    if not source.is_dir():
        raise BankError("No such directory: {}".format(source))

    missing = [d for d in REQUIRED_SUBDIRS if not (source / d).is_dir()]
    if missing:
        raise BankError(
            "{} has no {} directory, so it does not look like a superfit "
            "template bank.".format(source, " or ".join(missing))
        )

    if not copy:
        directory = register(name, source)
        if not quiet:
            print("Registered {!r} -> {}".format(name, directory))
        return directory

    destination = bank_install_dir(name)
    destination.parent.mkdir(parents=True, exist_ok=True)

    staging = Path(
        tempfile.mkdtemp(prefix=".{}-".format(destination.name), dir=destination.parent)
    )
    try:
        if not quiet:
            print("Copying {} -> {}".format(source, destination))
        # dirs_exist_ok: staging is already there, made by mkdtemp.
        shutil.copytree(source, staging, dirs_exist_ok=True, symlinks=True)
        if destination.exists():
            shutil.rmtree(destination)
        staging.replace(destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    register(name, destination)
    return destination


def status(directory=None):
    """Describe the template bank this machine would use.

    Returns a dict rather than printing, so the CLI can format it and
    ``superfit doctor`` can judge it.
    """

    report = {
        "found": False,
        "path": None,
        "source": None,
        "complete": False,
        "missing": list(REQUIRED_SUBDIRS),
        "n_sn_types": 0,
        "n_galaxy_templates": 0,
        "manifest": None,
        "install_dir": str(default_install_dir()),
        "packed": None,
        "error": None,
    }

    if directory is not None:
        found = Path(directory)
        report["source"] = "given on the command line"
        if not found.is_dir():
            report["error"] = "{} is not a directory".format(found)
            return report
    else:
        for var in paths.BANK_DIR_ENV_VARS:
            if os.environ.get(var):
                report["source"] = "${}".format(var)
                break
        try:
            found = Path(paths.find_bank_dir(use_cache=False))
        except FileNotFoundError as exc:
            report["error"] = str(exc)
            return report
        if report["source"] is None:
            report["source"] = "found on the search path"

    report["found"] = True
    report["path"] = str(found)
    report["missing"] = [d for d in REQUIRED_SUBDIRS if not (found / d).is_dir()]
    report["complete"] = not report["missing"]
    report["manifest"] = read_manifest(found)

    sne = found / "original_resolution" / "sne"
    gal = found / "binnings" / "10A" / "gal"
    # Dotfiles skipped: the published archive carries .DS_Store files, and
    # counting those as templates makes the report disagree with the fit.
    if sne.is_dir():
        report["n_sn_types"] = sum(
            1 for p in sne.iterdir() if p.is_dir() and not p.name.startswith(".")
        )
    if gal.is_dir():
        report["n_galaxy_templates"] = sum(
            1 for p in gal.iterdir() if p.is_file() and not p.name.startswith(".")
        )

    report["packed"] = pack_status(found)

    return report


def pack_status(directory):
    """Which of a bank's template directories have a usable pack.

    A stale pack counts as unpacked: it is what a fit would do with it.
    Reported so that "why is every run ten seconds slow" has an answer in
    ``superfit bank status`` rather than needing a profiler.
    """

    from superfit import packed

    directory = str(directory)
    report = {"root": None, "packed": [], "unpacked": []}

    for relative in packed.packable_directories(directory):
        if packed.open_pack(directory, relative) is None:
            report["unpacked"].append(relative)
        else:
            report["packed"].append(relative)

    if report["packed"]:
        for root in packed.pack_roots(directory):
            if os.path.isdir(root):
                report["root"] = root
                break

    return report
