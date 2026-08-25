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
#
# 2: the array name carries a generation and the index names it, so
#    publishing is one atomic rename; the index also records rows/dtype so it
#    can be checked against the array; and the fingerprint covers mtime as
#    well as size. A version-1 pack is not stale, it is a different format,
#    and saying so is what makes `bank status` report "not packed" rather
#    than implying the bank was edited.
FORMAT_VERSION = 2

INDEX_SUFFIX = ".index.json"

# Files that live among the templates but are not templates. params.py
# filters the same names out of the template lists it builds; the two lists
# have to agree, or the pack would either miss templates the fit wants or
# spend its time failing to parse photometry PDFs.
NON_TEMPLATE_NAMES = ("wiserep_spectra.csv", "info", "photometry", "photometry.pdf")

# Processes used to parse templates, when the caller does not say. Deliberately
# modest: packing is a setup step that often runs on a shared login node, and
# quietly taking 200 cores there to save a minute is not a trade to make on
# the user's behalf. `superfit bank pack --jobs N` raises it.
DEFAULT_PACK_JOBS = 8


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
    """``(relative path, size, mtime_ns)`` per template under ``directory``, sorted.

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
                    stat = entry.stat()
                    found.append(
                        (
                            relative.replace(os.sep, "/"),
                            stat.st_size,
                            stat.st_mtime_ns,
                        )
                    )

    found.sort()
    return found


def fingerprint(listing):
    """A digest of a directory listing, used to notice an edited bank.

    Names, sizes and modification times -- not contents. Hashing a gigabyte
    of templates would cost more than the reads the pack exists to avoid, and
    this has to be cheap enough to check at every fit boundary. Size alone was
    not enough: editing a flux value in place usually leaves the byte count
    exactly as it was, and the stale pack then stayed "valid" while the text
    said something else. mtime_ns closes that, so the remaining hole is an
    edit that preserves both length and timestamp, which takes deliberate
    effort rather than an accident.

    For a bank whose reproducibility has to be provable rather than probable,
    record its checksum: a pack is only ever as trustworthy as the cheapest
    check that validates it.
    """

    digest = hashlib.sha256()
    digest.update("{}\n".format(FORMAT_VERSION).encode())
    for relative, size, mtime_ns in listing:
        digest.update("{}\0{}\0{}\n".format(relative, size, mtime_ns).encode())
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


def _array_path(root, relative_directory, generation):
    """Where one *generation* of a packed array lives.

    The array name carries the generation, and the index names the array. That
    makes the index the single thing that has to change atomically, so a
    reader either sees the whole old pack or the whole new one.

    Replacing a fixed ``sne.npy`` and then its index was two renames, and a
    reader arriving between them saw the new array described by the old index
    -- offsets pointing into different data, which is not a failure that
    announces itself. It reads as a fit against subtly wrong templates.
    """

    parent, name = os.path.split(relative_directory)
    return os.path.join(root, parent, "{}.{}.npy".format(name, generation))


def _generation_of(index):
    """The array name an index refers to, or None for one that names none."""

    array = index.get("array")
    if not isinstance(array, str) or not array or "/" in array or "\\" in array:
        return None
    return array


def _sweep_old_generations(root, relative_directory, keep):
    """Delete arrays for this directory that no index refers to any more.

    On POSIX a reader that has the old array open, or mapped, keeps reading it
    after the unlink, so this cannot pull the ground out from under a fit in
    progress. Windows refuses the unlink instead, which is why failure here is
    ignored: a leftover array costs disk, and nothing else.
    """

    parent, name = os.path.split(relative_directory)
    directory = os.path.join(root, parent)
    prefix = name + "."

    try:
        entries = os.listdir(directory)
    except OSError:
        return

    for entry in entries:
        if entry == keep or not entry.startswith(prefix) or not entry.endswith(".npy"):
            continue
        try:
            os.remove(os.path.join(directory, entry))
        except OSError:
            pass


def _parse_one(args):
    """Read one template, and digest its bytes. Module level, so it pickles.

    The digest is of the file exactly as it sits on disk, so `bank verify` can
    later prove the pack still describes the text rather than merely agreeing
    with its size and timestamp.
    """

    directory, relative, reader_name = args
    path = os.path.join(directory, relative)

    try:
        with open(path, "rb") as handle:
            digest = hashlib.sha256(handle.read()).hexdigest()
    except OSError as exc:
        return relative, None, None, "{}: {}".format(type(exc).__name__, exc)

    try:
        array = np.asarray(READERS[reader_name](path), dtype=float)
    except Exception as exc:  # noqa: BLE001 -- any parse failure is the same answer
        return relative, None, None, "{}: {}".format(type(exc).__name__, exc)

    if array.ndim != 2 or array.shape[0] == 0 or array.shape[1] < 2:
        return relative, None, None, "not a 2-column spectrum"

    return relative, array, digest, None


def resolve_pack_jobs(requested):
    """How many processes to parse with.

    Capped low by default. Packing is a setup step that runs on whatever
    machine the user happens to be on -- often a shared login node -- and
    taking every core there to save a minute is not a trade worth making
    silently. ``--jobs`` raises it.
    """

    if requested is not None and requested > 0:
        return int(requested)

    try:
        available = len(os.sched_getaffinity(0))
    except AttributeError:
        available = os.cpu_count() or 1

    return max(1, min(available, DEFAULT_PACK_JOBS))


def _parse_all(directory, listing, reader_name, jobs, quiet):
    """Read every template in ``listing``, in order, across ``jobs`` processes.

    Parsing dominates packing -- kill_header runs a Python ``float()`` per
    value and manages 27 files a second, so the 15,561 original-resolution
    templates of the largest bank are ten minutes of pure parsing on one
    core. It is also embarrassingly parallel, so it is not spent on one core
    unless asked.
    """

    work = [(directory, entry[0], reader_name) for entry in listing]

    if jobs > 1 and len(work) > jobs:
        import multiprocessing as mp

        ctx = mp.get_context("fork" if "fork" in mp.get_all_start_methods() else None)
        chunksize = max(1, min(64, len(work) // (jobs * 8) or 1))
        with ctx.Pool(processes=jobs) as pool:
            results = pool.imap(_parse_one, work, chunksize=chunksize)
            yield from _with_progress(results, len(work), directory, quiet)
    else:
        yield from _with_progress(
            (_parse_one(item) for item in work), len(work), directory, quiet
        )


def _with_progress(results, total, directory, quiet):
    if quiet:
        yield from results
        return

    from tqdm import tqdm

    with tqdm(results, total=total, unit="tmpl", desc=os.path.basename(directory)) as bar:
        yield from bar


def pack_directory(bank_dir, relative_directory, root=None, quiet=True, jobs=None):
    """Pack one template directory. Returns the index that was written.

    Written to a scratch name and renamed into place, so a fit reading the
    pack concurrently sees either the old one or the new one and never half
    of either.
    """

    directory = os.path.join(bank_dir, relative_directory)
    root = root or writable_pack_root(bank_dir)
    reader_name = reader_name_for(relative_directory)

    listing = list_templates(directory)
    if not listing:
        raise PackError("No templates found in {}".format(directory))

    names = []
    arrays = []
    skipped = []
    content = hashlib.sha256()

    for relative, array, digest, why in _parse_all(
        directory, listing, reader_name, resolve_pack_jobs(jobs), quiet
    ):
        if array is None:
            skipped.append((relative, why))
        else:
            names.append(relative)
            arrays.append(array)
            # In listing order, which _parse_all preserves, so the digest is
            # reproducible rather than dependent on which worker finished first.
            content.update("{}\0{}\n".format(relative, digest).encode())

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
    # Each template is dropped as it is copied in. The largest bank's
    # original-resolution directory is 334 MB packed, and holding the list of
    # pieces alongside the finished array would make that 670 MB for no
    # reason.
    for i, (start, stop) in enumerate(zip(offsets[:-1], offsets[1:])):
        array = arrays[i]
        packed[start:stop, : array.shape[1]] = array
        arrays[i] = None

    digest = fingerprint(listing)

    # The generation names the array file, so publishing a new pack never
    # writes over the array an existing reader is holding. Derived from the
    # fingerprint, so packing the same unchanged directory twice is idempotent.
    generation = digest[:16]

    index = {
        "format": FORMAT_VERSION,
        "directory": relative_directory.replace(os.sep, "/"),
        "reader": reader_name,
        "fingerprint": digest,
        # sha256 over the bytes of every template that went in. Not checked on
        # the hot path -- that would mean re-reading the whole bank, which is
        # the cost the pack exists to remove -- but `superfit bank verify`
        # recomputes it, which is how a pack's agreement with the text becomes
        # provable rather than probable.
        "content": content.hexdigest(),
        "array": os.path.basename(_array_path(root, relative_directory, generation)),
        "n_templates": len(names),
        "width": width,
        "rows": int(offsets[-1]),
        "dtype": "float64",
        "names": names,
        "offsets": [int(x) for x in offsets],
        "widths": widths,
        "skipped": [name for name, _why in skipped],
    }

    array_path = _array_path(root, relative_directory, generation)
    index_path = _index_path(root, relative_directory)
    os.makedirs(os.path.dirname(array_path), exist_ok=True)

    scratch = tempfile.mkdtemp(prefix=".packing-", dir=os.path.dirname(array_path))
    try:
        staged_array = os.path.join(scratch, "templates.npy")
        staged_index = os.path.join(scratch, "templates.json")
        np.save(staged_array, packed)
        with open(staged_index, "w") as handle:
            json.dump(index, handle)

        # The array lands under a name nothing refers to yet, so it is
        # complete and fsynced before any index points at it. Replacing the
        # index is then the one step that publishes the new generation, and it
        # is a single atomic rename.
        os.replace(staged_array, array_path)
        os.replace(staged_index, index_path)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    _sweep_old_generations(root, relative_directory, os.path.basename(array_path))

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


def pack(bank_dir=None, root=None, quiet=True, jobs=None, only=None):
    """Pack every template directory in a bank. Returns the indexes written.

    ``only`` restricts the work to the directories named in it, which is how
    a caller packs the binnings it will actually fit at without paying for
    the original-resolution bank it will not.
    """

    if bank_dir is None:
        from superfit.paths import find_bank_dir

        bank_dir = find_bank_dir()

    directories = packable_directories(bank_dir)
    if only is not None:
        wanted = set(only)
        directories = [d for d in directories if d in wanted]

    if not directories:
        raise PackError(
            "{} has no template directories to pack. Is it a superfit "
            "bank?".format(bank_dir)
        )

    root = root or writable_pack_root(bank_dir)
    jobs = resolve_pack_jobs(jobs)

    if not quiet:
        print(
            "Packing {} into {} ({} process{})".format(
                bank_dir, root, jobs, "" if jobs == 1 else "es"
            )
        )

    return [
        pack_directory(bank_dir, relative, root=root, quiet=quiet, jobs=jobs)
        for relative in directories
    ]


def verify(bank_dir=None, quiet=True, jobs=None):
    """Re-read a bank's text and prove each pack still describes it.

    The per-fit check compares names, sizes and modification times, which is
    cheap enough to run at every fit boundary and catches every edit that
    changes any of the three. What it cannot catch is an edit that lands
    inside the filesystem's timestamp granularity of the pack being built --
    /tmp resolves mtime to about 8 ms, so two writes in one tick share a
    timestamp.

    This closes that by rehashing every template. It costs a full read of the
    bank, which is exactly what the pack exists to avoid, so it is a command
    someone runs rather than something a fit does.

    Returns ``[(relative_directory, ok, detail), ...]``.
    """

    if bank_dir is None:
        from superfit.paths import find_bank_dir

        bank_dir = find_bank_dir()

    results = []

    for relative in packable_directories(bank_dir):
        index = None
        for root in pack_roots(bank_dir):
            try:
                with open(_index_path(root, relative)) as handle:
                    index = json.load(handle)
                break
            except (OSError, ValueError):
                continue

        if index is None:
            results.append((relative, None, "not packed"))
            continue

        recorded = index.get("content")
        if not recorded:
            results.append(
                (relative, None, "packed before content digests were recorded")
            )
            continue

        directory = os.path.join(bank_dir, relative)
        reader_name = index.get("reader", reader_name_for(relative))
        listing = list_templates(directory)

        content = hashlib.sha256()
        for name, _array, digest, why in _parse_all(
            directory, listing, reader_name, resolve_pack_jobs(jobs), quiet
        ):
            if why is None:
                content.update("{}\0{}\n".format(name, digest).encode())

        ok = content.hexdigest() == recorded
        results.append(
            (
                relative,
                ok,
                "matches the text" if ok else "DOES NOT match the text -- repack",
            )
        )

    return results


class PackUnusable(Exception):
    """Internal: this pack cannot be trusted, so read the text instead."""


def _check_index(index):
    """The spans an index describes, or raise PackUnusable.

    Everything the index claims is checked against everything else it claims,
    before any of it is used to slice an array. An index that disagrees with
    itself -- offsets that run backwards, a name count that does not match the
    offset count, a width wider than the array -- would otherwise hand back a
    slice of some other template, which is the one failure mode a fit cannot
    detect for itself.
    """

    try:
        names = index["names"]
        offsets = index["offsets"]
        widths = index["widths"]
        width = int(index["width"])
        rows = int(index["rows"])
    except (KeyError, TypeError, ValueError) as exc:
        raise PackUnusable("index is missing fields: {}".format(exc))

    if not isinstance(names, list) or not isinstance(offsets, list):
        raise PackUnusable("index names/offsets are not lists")
    if len(offsets) != len(names) + 1 or len(widths) != len(names):
        raise PackUnusable(
            "index describes {} names, {} offsets and {} widths".format(
                len(names), len(offsets), len(widths)
            )
        )
    if not names:
        raise PackUnusable("index describes no templates")
    if offsets[0] != 0 or offsets[-1] != rows:
        raise PackUnusable("index offsets do not span 0..rows")
    if any(b < a for a, b in zip(offsets, offsets[1:])):
        raise PackUnusable("index offsets are not monotonic")
    if any(not 0 < w <= width for w in widths):
        raise PackUnusable("index has a width outside 1..{}".format(width))

    return {
        name: (int(start), int(stop), int(w))
        for name, start, stop, w in zip(names, offsets[:-1], offsets[1:], widths)
    }


def _open_array(array_path, index):
    """Memory-map a packed array, checking it is the one the index describes.

    A truncated or corrupt array is the case that used to escape: the index
    parsed, so the pack looked usable, and the failure surfaced from inside
    np.load or from slicing -- in the middle of a fit, as a traceback, with
    the text sitting right there unread. Checked here so the caller can fall
    back instead.
    """

    try:
        array = np.load(array_path, mmap_mode="r")
    except Exception as exc:  # noqa: BLE001 -- a bad header raises many things
        raise PackUnusable("array will not open: {}".format(exc))

    expected = (int(index["rows"]), int(index["width"]))
    if array.shape != expected:
        raise PackUnusable(
            "array is {} but the index describes {}".format(array.shape, expected)
        )
    if array.dtype != np.float64:
        raise PackUnusable("array is {}, not float64".format(array.dtype))

    # np.load trusts the header's shape; a file cut short after it was written
    # then faults on access rather than on open, and a SIGBUS from a mapped
    # page is not something a fit can catch. Compare the size on disk instead.
    try:
        on_disk = os.path.getsize(array_path)
    except OSError as exc:
        raise PackUnusable("array cannot be stat'd: {}".format(exc))
    if on_disk < array.nbytes:
        raise PackUnusable(
            "array is truncated: {} bytes on disk, {} needed".format(
                on_disk, array.nbytes
            )
        )

    return array


class _Pack:
    """One packed directory, opened for reading."""

    def __init__(self, array_path, index, source_directory):
        self.reader = index["reader"]
        self.fingerprint = index["fingerprint"]
        self.source_directory = source_directory
        self._array_path = array_path
        self._index = index
        self._span = _check_index(index)
        self._array = None
        # Cleared at each fit boundary; see stale().
        self._revalidate = False

    def _mapped(self):
        if self._array is None:
            self._array = _open_array(self._array_path, self._index)
        return self._array

    def get(self, relative):
        """One template, or None if this pack does not hold it.

        Raises PackUnusable rather than anything else, so a caller has exactly
        one thing to catch to fall back to the text.
        """

        span = self._span.get(relative)
        if span is None:
            return None

        start, stop, width = span
        array = self._mapped()

        try:
            # A copy, not a view onto the map: callers mask and bin these in
            # place, and a fit must not be able to write through to the pack
            # -- nor keep the whole mapping alive by holding one template.
            return np.array(array[start:stop, :width], dtype=np.float64)
        except Exception as exc:  # noqa: BLE001 -- a bad map raises many things
            raise PackUnusable("reading {}: {}".format(relative, exc))

    def stale(self):
        """Whether the source directory has changed since this was opened.

        Checked at most once per fit boundary, not per template: re-listing
        the largest bank's supernova directory is about 1.8 s, which is worth
        paying once against 77 s of reading text, and not worth paying 15000
        times.
        """

        if not self._revalidate:
            return False
        self._revalidate = False

        try:
            current = fingerprint(list_templates(self.source_directory))
        except OSError:
            return True
        return current != self.fingerprint


def open_pack(bank_dir, relative_directory, verify=True):
    """Open the pack for one directory, or return None if there is not a usable one.

    ``verify`` re-lists the source directory and compares the fingerprint, so
    a bank edited since it was packed falls back to the text instead of being
    fitted against a stale copy.

    Returns None for every way a pack can be unusable -- absent, stale,
    corrupt, self-contradictory, truncated. None of them is an error worth
    raising: the text is still there.
    """

    source = os.path.join(bank_dir, relative_directory)

    for root in pack_roots(bank_dir):
        index_path = _index_path(root, relative_directory)

        try:
            with open(index_path) as handle:
                index = json.load(handle)
        except (OSError, ValueError):
            continue

        if not isinstance(index, dict) or index.get("format") != FORMAT_VERSION:
            continue

        array_name = _generation_of(index)
        if array_name is None:
            continue

        parent = os.path.dirname(_index_path(root, relative_directory))
        array_path = os.path.join(parent, array_name)
        if not os.path.isfile(array_path):
            continue

        if verify:
            try:
                current = fingerprint(list_templates(source))
            except OSError:
                continue
            if current != index.get("fingerprint"):
                continue

        try:
            pack = _Pack(array_path, index, source)
            # Map it now rather than on first use, so a corrupt array is one
            # more reason to fall back here instead of a surprise later.
            pack._mapped()
            return pack
        except (PackUnusable, KeyError, TypeError, ValueError):
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


def _pack_for(path, bank_dir=None):
    """The pack covering ``path``, and the path relative to its directory.

    ``bank_dir`` is the bank the caller is fitting against. A fit always
    passes its own: resolving the process-wide bank here would let a second
    Superfit's construction decide which pack this fit reads.
    """

    if bank_dir is None:
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

        # Marked for revalidation by begin_fit(). A bank edited between two
        # fits in one process must not be read out of the pack built before
        # the edit; the cache used to hold the first fit's answer for the life
        # of the process, so only a caller that knew to clear it saw the edit.
        if pack_for_dir is not None and pack_for_dir.stale():
            pack_for_dir = open_pack(bank_dir, relative_directory)
            _open_packs[key] = pack_for_dir

        if pack_for_dir is not None:
            return pack_for_dir, "/".join(parts[depth:])

    return None, None


def begin_fit():
    """Mark every open pack for one revalidation against its source.

    Called at the start of a fit. The check itself is deferred to the next
    lookup and happens at most once per pack per fit, so a fit that reads
    15000 templates pays for one directory listing rather than 15000.
    """

    for pack in _open_packs.values():
        if pack is not None:
            pack._revalidate = True

    # A pack that was absent last time may have been built since, and the set
    # of packable directories can grow when a binning is added.
    for key, pack in list(_open_packs.items()):
        if pack is None:
            del _open_packs[key]
    _packable.clear()


def forget_open_packs():
    """Drop every pack and directory listing this process has cached.

    For tests, and for anything that writes a pack in a process that has
    already looked for one.
    """

    _open_packs.clear()
    _packable.clear()


def load_template(path, reader="loadtxt", bank_dir=None):
    """One bank template, from the pack when there is a usable one.

    ``reader`` names how the caller would read the text, and a pack built by
    a different reader is not used: the point of the pack is to hand back
    exactly what the caller would have parsed, and two readers do not always
    agree about a file.

    ``bank_dir`` is the bank to look in, and a fit always passes its own.

    Falls back to reading the text for every reason a pack can fail -- absent,
    stale, corrupt, truncated, a template added since it was built. The pack
    is an optimisation, so a broken one costs time and never an answer.
    """

    pack_for_dir, relative = _pack_for(path, bank_dir)

    if pack_for_dir is not None and pack_for_dir.reader == reader:
        try:
            array = pack_for_dir.get(relative)
        except PackUnusable:
            # Poison it for the rest of this fit rather than retrying a broken
            # map 15000 times, then read the text.
            for key, value in list(_open_packs.items()):
                if value is pack_for_dir:
                    _open_packs[key] = None
            array = None
        if array is not None:
            return array

    return READERS[reader](path)
