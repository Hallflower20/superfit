"""The template bank packed into one array per directory.

A fit reads roughly a thousand template files, each a few tens of kilobytes
of text. Parsing them is not the expensive part -- ``np.loadtxt`` gets through
the whole supernova bank in 0.7 s. *Opening* them is. On a parallel filesystem
every open is a round trip to a metadata server, and a thousand of those cost
several seconds before any arithmetic happens. Measured on Perlmutter's CFS,
reading the 30 A bank: 1022 files in 7.0 s cold, 0.17 s once the page cache is
warm. The 10 A bank inside a real fit cost 10.9 s cold, of which parsing was
0.67 s and host-line masking 0.07 s.

So this module trades a thousand files for two: every template in a directory
is concatenated into one ``.npy``, and a JSON index records where each one
starts. Reading it back is one memory-mapped open and a slice.

What gets packed is the *raw* array each template parses to, and nothing
derived from it. Host-line masking, binning and the log-grid resample still
happen at fit time, where they cost 0.07 s, 0.18 s and 0.02 s respectively.
That is deliberate. It keeps the pack independent of every fit setting -- the
epoch window, the selected types, ``mask_galaxy_lines`` -- so the only thing
its validity depends on is the bank itself, and no configuration change can
leave a fit reading a pack that does not match what it asked for.

The values stored are the float64 the text parsed to, so a fit against the
pack gives bit-identical results to a fit against the text. The pack is an
optimisation and never a requirement: nothing here raises into a fit, and a
missing, stale, or unreadable pack simply falls back to reading the text.
"""

import hashlib
import json
import os
import shutil
import tempfile

import numpy as np

# Directory name the pack is written into, both inside a writable bank and
# inside the per-user data directory.
PACK_DIRNAME = "packed"

# Bumped when the on-disk layout changes in a way an older reader would get
# wrong. A pack whose format does not match is ignored, not repaired.
FORMAT_VERSION = 1

INDEX_SUFFIX = ".index.json"

# Files that live among the templates but are not templates. params.py
# filters the same names out of the template lists it builds; the two lists
# have to agree, or the pack would either miss templates the fit wants or
# spend its time failing to parse photometry PDFs.
NON_TEMPLATE_NAMES = ("wiserep_spectra.csv", "info", "photometry", "photometry.pdf")


class PackError(RuntimeError):
    """Raised when a pack cannot be built. Never raised when reading one."""


def is_template(name):
    """True when a filename under a template directory is a template."""

    if name.startswith("."):
        return False
    return not any(junk in name for junk in NON_TEMPLATE_NAMES)


def _read_with_loadtxt(path):
    return np.loadtxt(path)


def _read_with_kill_header(path):
    from superfit.Header_Binnings import kill_header

    return kill_header(path)


# How each template directory is read at fit time. The name is recorded in the
# index and checked on lookup: a caller that reads a directory some other way
# falls back to the text rather than being handed an array built by a reader
# that might not agree with it. ``original_resolution/sne`` carries headers
# and needs kill_header; everything else the fit reads with np.loadtxt.
READERS = {
    "loadtxt": _read_with_loadtxt,
    "kill_header": _read_with_kill_header,
}


def reader_name_for(relative_directory):
    """Which reader the fit uses for a directory, relative to the bank root."""

    parts = relative_directory.replace("\\", "/").split("/")
    if parts[0] == "original_resolution" and parts[-1] == "sne":
        return "kill_header"
    return "loadtxt"


def packable_directories(bank_dir):
    """Template directories in ``bank_dir``, as paths relative to it.

    Both layouts, and every binning that happens to be installed, so a bank
    with a 20 A binning gets one too without this list being edited.
    """

    found = []

    for kind in ("sne", "gal"):
        if os.path.isdir(os.path.join(bank_dir, "original_resolution", kind)):
            found.append("original_resolution/" + kind)

    binnings = os.path.join(bank_dir, "binnings")
    if os.path.isdir(binnings):
        for entry in sorted(os.listdir(binnings)):
            for kind in ("sne", "gal"):
                if os.path.isdir(os.path.join(binnings, entry, kind)):
                    found.append("binnings/{}/{}".format(entry, kind))

    return found


def list_templates(directory):
    """``(relative path, size)`` for every template under ``directory``, sorted.

    ``os.scandir`` rather than ``os.walk``: walking the whole bank with a
    ``stat`` per file takes 1.3 s, and with scandir 0.19 s, because scandir
    keeps the directory entry the readdir already returned.
    """

    found = []
    stack = [directory]

    while stack:
        current = stack.pop()
        with os.scandir(current) as entries:
            for entry in entries:
                if entry.name.startswith("."):
                    continue
                if entry.is_dir(follow_symlinks=False):
                    stack.append(entry.path)
                elif is_template(entry.name):
                    relative = os.path.relpath(entry.path, directory)
                    found.append((relative.replace(os.sep, "/"), entry.stat().st_size))

    found.sort()
    return found


def fingerprint(listing):
    """A digest of a directory listing, used to notice an edited bank.

    Names and sizes, not contents: hashing 30 MB of templates would cost more
    than the reads the pack exists to avoid, while a template that changed
    without changing length is not a thing that happens to a published bank
    by accident. Cheap enough to check on every fit -- the listing it is
    computed from takes about 40 ms for the supernova bank.
    """

    digest = hashlib.sha256()
    digest.update("{}\n".format(FORMAT_VERSION).encode())
    for relative, size in listing:
        digest.update("{}\0{}\n".format(relative, size).encode())
    return digest.hexdigest()


def _bank_token(bank_dir):
    """A short stable name for a bank directory, for the per-user pack path."""

    return hashlib.sha256(os.path.abspath(bank_dir).encode()).hexdigest()[:16]


def _user_pack_root(bank_dir):
    from superfit.bank import user_data_dir

    return os.path.join(str(user_data_dir()), PACK_DIRNAME, _bank_token(bank_dir))


def pack_roots(bank_dir):
    """Where a pack for ``bank_dir`` may live, in the order they are tried.

    Beside the bank first, because a bank shared between users is worth
    packing once; then the per-user data directory, which is where the pack
    goes when the bank itself is read-only.
    """

    return [os.path.join(bank_dir, PACK_DIRNAME), _user_pack_root(bank_dir)]


def writable_pack_root(bank_dir):
    """Where a pack for ``bank_dir`` should be written."""

    beside = os.path.join(bank_dir, PACK_DIRNAME)
    if os.access(bank_dir, os.W_OK):
        return beside
    return _user_pack_root(bank_dir)


def _index_path(root, relative_directory):
    return os.path.join(root, relative_directory + INDEX_SUFFIX)


def _array_path(root, relative_directory):
    return os.path.join(root, relative_directory + ".npy")


def pack_directory(bank_dir, relative_directory, root=None, quiet=True):
    """Pack one template directory. Returns the index that was written.

    Written to a scratch name and renamed into place, so a fit reading the
    pack concurrently sees either the old one or the new one and never half
    of either.
    """

    directory = os.path.join(bank_dir, relative_directory)
    root = root or writable_pack_root(bank_dir)
    reader_name = reader_name_for(relative_directory)
    reader = READERS[reader_name]

    listing = list_templates(directory)
    if not listing:
        raise PackError("No templates found in {}".format(directory))

    names = []
    arrays = []
    skipped = []

    for relative, _size in listing:
        try:
            array = np.asarray(reader(os.path.join(directory, relative)), dtype=float)
        except Exception as exc:  # noqa: BLE001 -- any parse failure is the same answer
            skipped.append((relative, "{}: {}".format(type(exc).__name__, exc)))
            continue

        if array.ndim != 2 or array.shape[0] == 0 or array.shape[1] < 2:
            skipped.append((relative, "not a 2-column spectrum"))
            continue

        names.append(relative)
        arrays.append(array)

    if not arrays:
        raise PackError("No readable templates in {}".format(directory))

    # Ragged widths are padded with NaN and the true width recorded, so a
    # template carrying an error column keeps it and its neighbours are not
    # silently widened. Every template in the published bank is two columns,
    # so in practice nothing is padded.
    widths = [int(a.shape[1]) for a in arrays]
    width = max(widths)

    lengths = [int(a.shape[0]) for a in arrays]
    offsets = np.concatenate([[0], np.cumsum(lengths)]).astype(np.int64)

    packed = np.full((int(offsets[-1]), width), np.nan, dtype=np.float64)
    for array, start, stop in zip(arrays, offsets[:-1], offsets[1:]):
        packed[start:stop, : array.shape[1]] = array

    index = {
        "format": FORMAT_VERSION,
        "directory": relative_directory.replace(os.sep, "/"),
        "reader": reader_name,
        "fingerprint": fingerprint(listing),
        "n_templates": len(names),
        "width": width,
        "names": names,
        "offsets": [int(x) for x in offsets],
        "widths": widths,
        "skipped": [name for name, _why in skipped],
    }

    array_path = _array_path(root, relative_directory)
    index_path = _index_path(root, relative_directory)
    os.makedirs(os.path.dirname(array_path), exist_ok=True)

    scratch = tempfile.mkdtemp(prefix=".packing-", dir=os.path.dirname(array_path))
    try:
        staged_array = os.path.join(scratch, "templates.npy")
        staged_index = os.path.join(scratch, "templates.json")
        np.save(staged_array, packed)
        with open(staged_index, "w") as handle:
            json.dump(index, handle)

        # The array goes first: a reader checks for the index and would
        # otherwise find one pointing at an array that is not there yet.
        os.replace(staged_array, array_path)
        os.replace(staged_index, index_path)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    if not quiet:
        print(
            "  {}: {} templates, {:.1f} MB{}".format(
                relative_directory,
                len(names),
                os.path.getsize(array_path) / 1e6,
                "" if not skipped else ", {} unreadable".format(len(skipped)),
            )
        )

    return index


def pack(bank_dir=None, root=None, quiet=True):
    """Pack every template directory in a bank. Returns the indexes written."""

    if bank_dir is None:
        from superfit.paths import find_bank_dir

        bank_dir = find_bank_dir()

    directories = packable_directories(bank_dir)
    if not directories:
        raise PackError(
            "{} has no template directories to pack. Is it a superfit "
            "bank?".format(bank_dir)
        )

    root = root or writable_pack_root(bank_dir)
    if not quiet:
        print("Packing {} into {}".format(bank_dir, root))

    return [
        pack_directory(bank_dir, relative, root=root, quiet=quiet)
        for relative in directories
    ]


class _Pack:
    """One packed directory, opened for reading."""

    def __init__(self, array_path, index):
        self._array_path = array_path
        self._array = None
        self.reader = index["reader"]
        self.fingerprint = index["fingerprint"]
        self._span = {
            name: (int(start), int(stop), int(w))
            for name, start, stop, w in zip(
                index["names"], index["offsets"][:-1], index["offsets"][1:],
                index["widths"],
            )
        }

    def get(self, relative):
        span = self._span.get(relative)
        if span is None:
            return None

        if self._array is None:
            # Memory-mapped, so opening the 10 A supernova pack touches the
            # header and nothing else; the pages come in as templates are
            # sliced out of it.
            self._array = np.load(self._array_path, mmap_mode="r")

        start, stop, width = span
        # A copy, not a view onto the map: callers mask and bin these in
        # place, and a fit must not be able to write through to the pack --
        # nor to keep the whole mapping alive by holding one template.
        return np.array(self._array[start:stop, :width], dtype=np.float64)


def open_pack(bank_dir, relative_directory, verify=True):
    """Open the pack for one directory, or return None if there is not a usable one.

    ``verify`` re-lists the source directory and compares the fingerprint, so
    a bank edited since it was packed falls back to the text instead of being
    fitted against a stale copy.
    """

    for root in pack_roots(bank_dir):
        index_path = _index_path(root, relative_directory)
        array_path = _array_path(root, relative_directory)

        try:
            with open(index_path) as handle:
                index = json.load(handle)
        except (OSError, ValueError):
            continue

        if index.get("format") != FORMAT_VERSION or not os.path.isfile(array_path):
            continue

        if verify:
            source = os.path.join(bank_dir, relative_directory)
            try:
                current = fingerprint(list_templates(source))
            except OSError:
                continue
            if current != index.get("fingerprint"):
                continue

        try:
            return _Pack(array_path, index)
        except (KeyError, TypeError, ValueError):
            continue

    return None


# Packs opened by this process, keyed by (bank directory, relative directory).
# A fit asks for the same directory a thousand times, and the fingerprint
# check -- cheap, but 40 ms -- should happen once per fit, not once per
# template. None is cached too: a bank with no pack must not be re-checked
# a thousand times either. Populated before the worker pool forks, so the
# children inherit it rather than each re-verifying the same bank.
_open_packs = {}

# packable_directories() lists the binnings directory, which is one more
# round trip to the filesystem than this lookup can afford to make per
# template.
_packable = {}


def _pack_for(path):
    """The pack covering ``path``, and the path relative to its directory."""

    from superfit.paths import find_bank_dir

    try:
        bank_dir = find_bank_dir()
    except FileNotFoundError:
        return None, None

    absolute = os.path.abspath(path)
    relative_to_bank = os.path.relpath(absolute, bank_dir)
    if relative_to_bank.startswith(os.pardir):
        return None, None

    if bank_dir not in _packable:
        _packable[bank_dir] = frozenset(packable_directories(bank_dir))
    packable = _packable[bank_dir]

    parts = relative_to_bank.replace(os.sep, "/").split("/")
    for depth in (3, 2):
        if len(parts) <= depth:
            continue
        relative_directory = "/".join(parts[:depth])
        if relative_directory not in packable:
            continue

        key = (bank_dir, relative_directory)
        if key not in _open_packs:
            _open_packs[key] = open_pack(bank_dir, relative_directory)

        pack_for_dir = _open_packs[key]
        if pack_for_dir is not None:
            return pack_for_dir, "/".join(parts[depth:])

    return None, None


def forget_open_packs():
    """Drop every pack and directory listing this process has cached.

    For tests, and for anything that writes a pack in a process that has
    already looked for one.
    """

    _open_packs.clear()
    _packable.clear()


def load_template(path, reader="loadtxt"):
    """One bank template, from the pack when there is a usable one.

    ``reader`` names how the caller would read the text, and a pack built by
    a different reader is not used: the point of the pack is to hand back
    exactly what the caller would have parsed, and two readers do not always
    agree about a file. Falls back to reading the text for any reason at all
    -- no pack, a stale one, a template added since it was built.
    """

    pack_for_dir, relative = _pack_for(path)

    if pack_for_dir is not None and pack_for_dir.reader == reader:
        array = pack_for_dir.get(relative)
        if array is not None:
            return array

    return READERS[reader](path)
