"""Packing the template bank into one array per directory.

The pack exists to make a fit faster, so the thing worth testing is not that
it is fast but that it is *invisible*: a fit reading it must get exactly the
array it would have parsed from the text, and every way the pack can be wrong
-- absent, stale, half-written, built by a reader that disagrees -- has to end
in reading the text rather than in a wrong answer.

Everything here builds a small bank in a tmp_path, so none of it needs the
real 74 MB one.
"""

import json
import os
import time

import numpy as np
import pytest

from superfit import packed, paths


TEMPLATE_A = "4000.0 1.5\n4010.0 1.6\n4020.0 1.7\n"
TEMPLATE_B = "5000.0 0.5\n5010.0 0.25\n"
GALAXY = "4000.0 2.0\n4010.0 2.5\n"

SN_A = "binnings/10A/sne/Ia-norm/SN2011fe/a.ascii"
SN_B = "binnings/10A/sne/Ic/SN1994I/b.ascii"


def array_path(bank_dir, relative_directory):
    """The packed array a directory's index currently names.

    Looked up through the index rather than assumed, because the array name
    carries a generation: that is what makes publishing a new pack a single
    atomic rename of the index.
    """

    index_path = bank_dir / "packed" / (relative_directory + ".index.json")
    index = json.loads(index_path.read_text())
    return index_path.parent / index["array"]


@pytest.fixture
def bank_dir(tmp_path, monkeypatch):
    """A bank with both layouts, pointed at by SUPERFIT_BANK_DIR."""

    root = tmp_path / "bank"

    files = {
        SN_A: TEMPLATE_A,
        SN_B: TEMPLATE_B,
        "binnings/10A/gal/E": GALAXY,
        "original_resolution/sne/Ia-norm/SN2011fe/a.ascii": "# a header\n" + TEMPLATE_A,
        "original_resolution/gal/E": GALAXY,
        # Not templates. params.py filters these out of the fit's own template
        # list, and the pack has to agree with it.
        "original_resolution/sne/Ia-norm/SN2011fe/wiserep_spectra.csv": "junk\n",
        "original_resolution/sne/Ia-norm/SN2011fe/info": "-\n",
    }

    for relative, content in files.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)

    monkeypatch.setenv("SUPERFIT_BANK_DIR", str(root))
    monkeypatch.setattr(paths, "_cached_bank_dir", None)
    packed.forget_open_packs()
    yield root
    packed.forget_open_packs()


class TestWhatGetsPacked:
    def test_every_template_directory_in_the_bank(self, bank_dir):
        assert packed.packable_directories(str(bank_dir)) == [
            "original_resolution/sne",
            "original_resolution/gal",
            "binnings/10A/sne",
            "binnings/10A/gal",
        ]

    def test_a_binning_that_happens_to_be_installed(self, bank_dir):
        """The list is discovered, not hard-coded: a 20 A bank gets one too."""

        extra = bank_dir / "binnings" / "20A" / "sne" / "Ic" / "SN1994I"
        extra.mkdir(parents=True)
        (extra / "b.ascii").write_text(TEMPLATE_B)

        assert "binnings/20A/sne" in packed.packable_directories(str(bank_dir))

    def test_non_templates_are_left_out(self, bank_dir):
        index = packed.pack_directory(str(bank_dir), "original_resolution/sne")

        assert index["names"] == ["Ia-norm/SN2011fe/a.ascii"]

    def test_the_reader_matches_how_the_fit_reads_that_directory(self):
        """Only the headered supernova directory needs kill_header."""

        assert packed.reader_name_for("original_resolution/sne") == "kill_header"
        assert packed.reader_name_for("original_resolution/gal") == "loadtxt"
        assert packed.reader_name_for("binnings/10A/sne") == "loadtxt"

    def test_an_empty_directory_is_an_error_not_an_empty_pack(self, bank_dir):
        empty = bank_dir / "binnings" / "40A" / "sne"
        empty.mkdir(parents=True)

        with pytest.raises(packed.PackError):
            packed.pack_directory(str(bank_dir), "binnings/40A/sne")


class TestReadingBack:
    def test_a_packed_template_is_what_the_text_would_have_parsed(self, bank_dir):
        expected = {
            relative: np.loadtxt(bank_dir / relative) for relative in (SN_A, SN_B)
        }

        packed.pack(str(bank_dir))
        packed.forget_open_packs()

        for relative, wanted in expected.items():
            got = packed.load_template(str(bank_dir / relative), "loadtxt")
            # Bit-for-bit, not close: the whole point is that a fit against
            # the pack gives the same answer as a fit against the text.
            assert np.array_equal(got, wanted)

    def test_a_header_is_stripped_the_same_way(self, bank_dir):
        from superfit.Header_Binnings import kill_header

        path = bank_dir / "original_resolution/sne/Ia-norm/SN2011fe/a.ascii"
        expected = kill_header(str(path))

        packed.pack(str(bank_dir))
        packed.forget_open_packs()

        assert np.array_equal(
            packed.load_template(str(path), "kill_header"), expected
        )

    def test_the_caller_cannot_write_through_to_the_pack(self, bank_dir):
        """A fit masks host lines in place; that must not reach the array."""

        packed.pack(str(bank_dir))
        packed.forget_open_packs()

        first = packed.load_template(str(bank_dir / SN_A), "loadtxt")
        first[0, 1] = np.nan

        second = packed.load_template(str(bank_dir / SN_A), "loadtxt")
        assert second[0, 1] == 1.5


class TestFallingBackToTheText:
    def test_with_no_pack_at_all(self, bank_dir):
        expected = np.loadtxt(bank_dir / SN_A)

        assert np.array_equal(
            packed.load_template(str(bank_dir / SN_A), "loadtxt"), expected
        )

    def test_when_the_bank_has_been_edited_since_it_was_packed(self, bank_dir):
        packed.pack(str(bank_dir))
        packed.forget_open_packs()

        # A template whose length changed: the pack now describes a bank that
        # is not there any more, and using it would fit the old spectrum.
        (bank_dir / SN_A).write_text(TEMPLATE_A + "4030.0 1.8\n")

        got = packed.load_template(str(bank_dir / SN_A), "loadtxt")
        assert got.shape == (4, 2)
        assert np.array_equal(got, np.loadtxt(bank_dir / SN_A))

    def test_when_a_template_was_added_after_packing(self, bank_dir):
        packed.pack(str(bank_dir))
        packed.forget_open_packs()

        added = bank_dir / "binnings/10A/sne/Ic/SN1994I/c.ascii"
        added.write_text(TEMPLATE_B)

        assert np.array_equal(
            packed.load_template(str(added), "loadtxt"), np.loadtxt(added)
        )

    def test_when_the_reader_disagrees(self, bank_dir):
        """A pack built by one reader is not handed to a caller using another."""

        packed.pack(str(bank_dir))
        packed.forget_open_packs()

        # The 10 A bank is packed with loadtxt. Asking for it as kill_header
        # must read the text, not hand back the loadtxt array.
        got = packed.load_template(str(bank_dir / SN_A), "kill_header")
        from superfit.Header_Binnings import kill_header

        assert np.array_equal(got, kill_header(str(bank_dir / SN_A)))

    def test_when_the_format_version_moves_on(self, bank_dir):
        packed.pack(str(bank_dir))
        index_path = bank_dir / "packed" / "binnings/10A/sne.index.json"
        index = json.loads(index_path.read_text())
        index["format"] = packed.FORMAT_VERSION + 1
        index_path.write_text(json.dumps(index))
        packed.forget_open_packs()

        assert packed.open_pack(str(bank_dir), "binnings/10A/sne") is None

    def test_when_the_index_is_corrupt(self, bank_dir):
        packed.pack(str(bank_dir))
        (bank_dir / "packed" / "binnings/10A/sne.index.json").write_text("{not json")
        packed.forget_open_packs()

        assert packed.open_pack(str(bank_dir), "binnings/10A/sne") is None
        assert np.array_equal(
            packed.load_template(str(bank_dir / SN_A), "loadtxt"),
            np.loadtxt(bank_dir / SN_A),
        )

    def test_when_the_array_is_missing(self, bank_dir):
        packed.pack(str(bank_dir))
        array_path(bank_dir, "binnings/10A/sne").unlink()
        packed.forget_open_packs()

        assert packed.open_pack(str(bank_dir), "binnings/10A/sne") is None

    def test_when_the_array_is_truncated(self, bank_dir):
        """The case that used to escape: index fine, array cut short."""

        packed.pack(str(bank_dir))
        path = array_path(bank_dir, "binnings/10A/sne")
        data = path.read_bytes()
        path.write_bytes(data[: len(data) - 16])
        packed.forget_open_packs()

        assert packed.open_pack(str(bank_dir), "binnings/10A/sne") is None
        assert np.array_equal(
            packed.load_template(str(bank_dir / SN_A), "loadtxt"),
            np.loadtxt(bank_dir / SN_A),
        )

    def test_when_the_array_header_is_garbage(self, bank_dir):
        packed.pack(str(bank_dir))
        array_path(bank_dir, "binnings/10A/sne").write_bytes(b"not an npy file")
        packed.forget_open_packs()

        assert packed.open_pack(str(bank_dir), "binnings/10A/sne") is None
        assert np.array_equal(
            packed.load_template(str(bank_dir / SN_A), "loadtxt"),
            np.loadtxt(bank_dir / SN_A),
        )

    def test_when_the_index_disagrees_with_the_array(self, bank_dir):
        """Offsets that describe a different array would slice the wrong rows."""

        packed.pack(str(bank_dir))
        index_path = bank_dir / "packed" / "binnings/10A/sne.index.json"
        index = json.loads(index_path.read_text())
        index["rows"] = index["rows"] + 100
        index["offsets"][-1] = index["rows"]
        index_path.write_text(json.dumps(index))
        packed.forget_open_packs()

        assert packed.open_pack(str(bank_dir), "binnings/10A/sne") is None

    @pytest.mark.parametrize(
        "damage",
        [
            pytest.param(lambda i: i.update(offsets=i["offsets"][:-1]), id="short"),
            pytest.param(lambda i: i.update(offsets=[0, 9, 3]), id="backwards"),
            pytest.param(lambda i: i.update(widths=[99, 99]), id="too-wide"),
            pytest.param(lambda i: i.update(names=[]), id="no-names"),
        ],
    )
    def test_an_index_that_disagrees_with_itself(self, bank_dir, damage):
        packed.pack(str(bank_dir))
        index_path = bank_dir / "packed" / "binnings/10A/sne.index.json"
        index = json.loads(index_path.read_text())
        damage(index)
        index_path.write_text(json.dumps(index))
        packed.forget_open_packs()

        assert packed.open_pack(str(bank_dir), "binnings/10A/sne") is None

    def test_for_a_path_outside_the_bank(self, bank_dir, tmp_path):
        """load_template is given bank paths, but must not guess at others."""

        stray = tmp_path / "elsewhere.ascii"
        stray.write_text(TEMPLATE_B)

        assert np.array_equal(
            packed.load_template(str(stray), "loadtxt"), np.loadtxt(stray)
        )


class TestPublication:
    def test_the_array_name_carries_a_generation(self, bank_dir):
        packed.pack(str(bank_dir))

        name = array_path(bank_dir, "binnings/10A/sne").name
        assert name.startswith("sne.") and name.endswith(".npy")
        assert name != "sne.npy"

    def test_repacking_unchanged_content_is_idempotent(self, bank_dir):
        packed.pack(str(bank_dir))
        first = array_path(bank_dir, "binnings/10A/sne").name

        packed.pack(str(bank_dir))

        assert array_path(bank_dir, "binnings/10A/sne").name == first

    def test_repacking_changed_content_publishes_a_new_generation(self, bank_dir):
        packed.pack(str(bank_dir))
        first = array_path(bank_dir, "binnings/10A/sne")

        (bank_dir / SN_A).write_text(TEMPLATE_A + "4030.0 1.8\n")
        packed.pack(str(bank_dir))
        second = array_path(bank_dir, "binnings/10A/sne")

        assert second.name != first.name
        assert second.is_file()

    def test_a_reader_holding_the_old_array_is_not_overwritten(self, bank_dir):
        """The old generation stays readable; only the index switches over."""

        packed.pack(str(bank_dir))
        old = array_path(bank_dir, "binnings/10A/sne")
        held = np.load(old, mmap_mode="r")
        first_row = np.array(held[0])

        (bank_dir / SN_A).write_text(TEMPLATE_A + "4030.0 1.8\n")
        packed.pack(str(bank_dir))

        # Deleted from the directory, but the open mapping still reads what it
        # always read -- no writer scribbled a different array over it.
        assert np.array_equal(np.array(held[0]), first_row)

    def test_old_generations_are_swept_up(self, bank_dir):
        packed.pack(str(bank_dir))
        (bank_dir / SN_A).write_text(TEMPLATE_A + "4030.0 1.8\n")
        packed.pack(str(bank_dir))

        directory = bank_dir / "packed" / "binnings" / "10A"
        arrays = sorted(p.name for p in directory.glob("sne.*.npy"))
        assert len(arrays) == 1


class TestRevalidationAtFitBoundaries:
    """A pack opened once must not be trusted forever.

    These deliberately do not call forget_open_packs(): the earlier tests did,
    which meant they proved the fingerprint works while saying nothing about
    whether a long-lived process ever checks it again.
    """

    def test_an_edit_between_fits_is_noticed(self, bank_dir):
        packed.pack(str(bank_dir))
        packed.forget_open_packs()

        # Fit one: reads the pack.
        packed.begin_fit()
        assert np.array_equal(
            packed.load_template(str(bank_dir / SN_A), "loadtxt"),
            np.loadtxt(bank_dir / SN_A),
        )

        # The bank changes underneath, without the pack being rebuilt.
        (bank_dir / SN_A).write_text(TEMPLATE_A + "4030.0 1.8\n")

        # Fit two: must see the text, not the pack built before the edit.
        packed.begin_fit()
        got = packed.load_template(str(bank_dir / SN_A), "loadtxt")
        assert got.shape == (4, 2)
        assert np.array_equal(got, np.loadtxt(bank_dir / SN_A))

    def test_an_edit_that_keeps_the_byte_count_is_still_noticed(self, bank_dir):
        """Size alone missed this: editing a value in place rarely resizes."""

        packed.pack(str(bank_dir))
        packed.forget_open_packs()
        packed.begin_fit()
        packed.load_template(str(bank_dir / SN_A), "loadtxt")

        # Past the filesystem's mtime granularity, which is what the check
        # rests on. /tmp resolves to about 8 ms, so an edit in the same tick
        # as the pack build shares its timestamp and is genuinely invisible
        # here -- that is the hole `superfit bank verify` exists to close.
        time.sleep(0.05)

        original = (bank_dir / SN_A).read_text()
        edited = original.replace("1.5", "9.9")
        assert len(edited) == len(original)
        (bank_dir / SN_A).write_text(edited)

        packed.begin_fit()
        got = packed.load_template(str(bank_dir / SN_A), "loadtxt")
        assert got[0, 1] == 9.9


class TestVerify:
    """The provable check, for when probable is not enough."""

    def test_a_fresh_pack_matches_the_text(self, bank_dir):
        packed.pack(str(bank_dir))

        results = packed.verify(str(bank_dir))

        assert results
        assert all(ok for _relative, ok, _detail in results)

    def test_an_edit_the_cheap_check_could_miss_is_still_caught(self, bank_dir):
        """Same length, same timestamp -- only the bytes give it away."""

        packed.pack(str(bank_dir))

        path = bank_dir / SN_A
        stat = path.stat()
        original = path.read_text()
        edited = original.replace("1.5", "9.9")
        assert len(edited) == len(original)
        path.write_text(edited)
        # Put the timestamp back, so size and mtime both agree with the pack.
        os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))

        assert packed.fingerprint(
            packed.list_templates(str(bank_dir / "binnings/10A/sne"))
        ) == json.loads(
            (bank_dir / "packed" / "binnings/10A/sne.index.json").read_text()
        )["fingerprint"], "the cheap check should be fooled here"

        results = dict(
            (relative, ok) for relative, ok, _detail in packed.verify(str(bank_dir))
        )
        assert results["binnings/10A/sne"] is False

    def test_an_unpacked_directory_is_reported_not_failed(self, bank_dir):
        results = dict(
            (relative, ok) for relative, ok, _detail in packed.verify(str(bank_dir))
        )

        assert set(results.values()) == {None}

    def test_within_one_fit_the_pack_is_listed_once(self, bank_dir, monkeypatch):
        """Revalidation is per fit, not per template."""

        packed.pack(str(bank_dir))
        packed.forget_open_packs()

        calls = []
        real = packed.list_templates
        monkeypatch.setattr(
            packed,
            "list_templates",
            lambda d: (calls.append(d), real(d))[1],
        )

        packed.begin_fit()
        for _ in range(5):
            packed.load_template(str(bank_dir / SN_A), "loadtxt")
            packed.load_template(str(bank_dir / SN_B), "loadtxt")

        # One listing to open the pack; none of the ten lookups adds another.
        assert len(calls) == 1

    def test_a_pack_built_between_fits_starts_being_used(self, bank_dir):
        packed.forget_open_packs()

        packed.begin_fit()
        assert packed.load_template(str(bank_dir / SN_A), "loadtxt") is not None

        packed.pack(str(bank_dir))

        packed.begin_fit()
        assert packed.open_pack(str(bank_dir), "binnings/10A/sne") is not None
        assert np.array_equal(
            packed.load_template(str(bank_dir / SN_A), "loadtxt"),
            np.loadtxt(bank_dir / SN_A),
        )


class TestWhereThePackGoes:
    def test_beside_the_bank_when_it_is_writable(self, bank_dir):
        assert packed.writable_pack_root(str(bank_dir)) == str(bank_dir / "packed")

    def test_in_the_user_data_directory_when_the_bank_is_read_only(
        self, bank_dir, monkeypatch
    ):
        monkeypatch.setattr(packed.os, "access", lambda path, mode: False)

        root = packed.writable_pack_root(str(bank_dir))
        assert root != str(bank_dir / "packed")
        assert "packed" in root

    def test_both_are_searched_when_reading(self, bank_dir):
        roots = packed.pack_roots(str(bank_dir))

        assert roots[0] == str(bank_dir / "packed")
        assert len(roots) == 2


class TestFingerprint:
    def test_it_notices_a_changed_size(self, bank_dir):
        directory = str(bank_dir / "binnings" / "10A" / "sne")
        before = packed.fingerprint(packed.list_templates(directory))

        (bank_dir / SN_A).write_text(TEMPLATE_A + "4030.0 1.8\n")

        assert packed.fingerprint(packed.list_templates(directory)) != before

    def test_it_notices_a_removed_template(self, bank_dir):
        directory = str(bank_dir / "binnings" / "10A" / "sne")
        before = packed.fingerprint(packed.list_templates(directory))

        (bank_dir / SN_B).unlink()

        assert packed.fingerprint(packed.list_templates(directory)) != before

    def test_it_does_not_depend_on_the_order_the_filesystem_lists_them(
        self, bank_dir
    ):
        """A digest that moved with readdir order would invalidate itself."""

        directory = str(bank_dir / "binnings" / "10A" / "sne")
        listing = packed.list_templates(directory)

        assert listing == sorted(listing)
        assert packed.fingerprint(listing) == packed.fingerprint(
            packed.list_templates(directory)
        )


class TestBankStatus:
    def test_it_says_when_nothing_is_packed(self, bank_dir):
        from superfit import bank

        report = bank.pack_status(str(bank_dir))

        assert report["packed"] == []
        assert len(report["unpacked"]) == 4

    def test_it_says_when_everything_is(self, bank_dir):
        from superfit import bank

        packed.pack(str(bank_dir))
        packed.forget_open_packs()

        report = bank.pack_status(str(bank_dir))

        assert report["unpacked"] == []
        assert len(report["packed"]) == 4
        assert report["root"] == str(bank_dir / "packed")

    def test_a_stale_pack_counts_as_unpacked(self, bank_dir):
        from superfit import bank

        packed.pack(str(bank_dir))
        (bank_dir / SN_A).write_text(TEMPLATE_A + "4030.0 1.8\n")
        packed.forget_open_packs()

        report = bank.pack_status(str(bank_dir))

        assert "binnings/10A/sne" in report["unpacked"]
