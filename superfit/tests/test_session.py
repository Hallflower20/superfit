"""Reusing a prepared bank across fits, and what that must not hide.

Preparing the bank -- validating the pack, reading every template, masking,
resampling -- is the expensive half of a fit and has nothing to do with the
spectrum. Keeping it is what makes a batch fast. The risk is that keeping it
also keeps a stale answer, so most of what is tested here is when the cache
must *miss*.
"""

import numpy as np
import pytest

from superfit import SF_functions
from superfit.session import Session


class FakeGrid:
    def __init__(self, identity=(1.0, 2.0, 3)):
        self.identity = identity


class FakeParameters:
    def __init__(self, **kw):
        self.bank_dir = kw.get("bank_dir", "/bank")
        self.resolution = kw.get("resolution", 10)
        self.mask_galaxy_lines = kw.get("mask_galaxy_lines", True)
        self.metadata_key = kw.get("metadata_key", ("/bank", "table", ("Ia",), 0, 0))
        self.templates_gal_trunc = kw.get("templates_gal_trunc", ["/bank/gal/E"])


class FakeMetadata:
    def __init__(self, n=5):
        self.dictionary_all_trunc_objects = {i: i for i in range(n)}


def key(**kw):
    grid = kw.pop("grid", FakeGrid())
    max_z = kw.pop("max_z", 0.1)
    state = kw.pop("bank_state", "pinned")
    n_objects = kw.pop("n_objects", 5)
    return SF_functions.prepared_bank_key(
        FakeParameters(**kw), FakeMetadata(n_objects), grid, max_z, state
    )


class TestPreparedBankKey:
    def test_the_same_inputs_give_the_same_key(self):
        assert key() == key()

    @pytest.mark.parametrize(
        "change",
        [
            {"bank_dir": "/other"},
            {"resolution": 30},
            {"mask_galaxy_lines": False},
            {"metadata_key": ("/bank", "table", ("Ic",), 0, 0)},
            {"templates_gal_trunc": ["/bank/gal/E", "/bank/gal/Sa"]},
            {"n_objects": 6},
            {"grid": FakeGrid((9.0, 2.0, 3))},
            {"max_z": 0.2},
            {"bank_state": (("binnings/10A/sne", "abc"),)},
        ],
        ids=[
            "bank", "resolution", "masking", "metadata", "galaxies",
            "object-count", "grid", "max-z", "bank-contents",
        ],
    )
    def test_anything_the_preparation_depends_on_changes_the_key(self, change):
        assert key(**change) != key()

    def test_the_grid_is_part_of_it(self):
        """Templates resampled onto one grid mean nothing on another."""

        assert key(grid=FakeGrid((1.0, 2.0, 4))) != key(grid=FakeGrid((1.0, 2.0, 3)))


class TestBankStateToken:
    def test_pinning_returns_a_constant(self, tmp_path):
        a = SF_functions.bank_state_token(str(tmp_path), revalidate=False)
        b = SF_functions.bank_state_token(str(tmp_path), revalidate=False)

        assert a == b == "pinned"

    def test_an_unpacked_bank_still_produces_a_token(self, tmp_path):
        for relative in ("binnings/10A/sne/Ia-norm/x", "binnings/10A/gal"):
            (tmp_path / relative).mkdir(parents=True)
        (tmp_path / "binnings/10A/sne/Ia-norm/x/a.ascii").write_text(
            "4000 1.0\n4010 1.1\n"
        )

        token = SF_functions.bank_state_token(str(tmp_path), revalidate=True)

        assert token != "pinned"

    def test_it_follows_the_packs(self, tmp_path, monkeypatch):
        """A repacked bank must not be answered from a prepared bank."""

        from superfit import packed, paths

        for relative in ("binnings/10A/sne/Ia-norm/x", "binnings/10A/gal"):
            (tmp_path / relative).mkdir(parents=True)
        template = tmp_path / "binnings/10A/sne/Ia-norm/x/a.ascii"
        template.write_text("4000 1.0\n4010 1.1\n")
        (tmp_path / "binnings/10A/gal/E").write_text("4000 2.0\n4010 2.1\n")
        (tmp_path / "original_resolution/gal").mkdir(parents=True)
        (tmp_path / "original_resolution/gal/E").write_text("4000 2.0\n4010 2.1\n")
        (tmp_path / "original_resolution/sne/Ia-norm/x").mkdir(parents=True)
        (tmp_path / "original_resolution/sne/Ia-norm/x/a.ascii").write_text(
            "4000 1.0\n4010 1.1\n"
        )

        monkeypatch.setenv("SUPERFIT_BANK_DIR", str(tmp_path))
        monkeypatch.setattr(paths, "_cached_bank_dir", None)
        packed.forget_open_packs()

        packed.pack(str(tmp_path))
        packed.forget_open_packs()
        before = SF_functions.bank_state_token(str(tmp_path), revalidate=True)

        import time

        time.sleep(0.05)
        template.write_text("4000 9.9\n4010 9.9\n")
        packed.forget_open_packs()
        after = SF_functions.bank_state_token(str(tmp_path), revalidate=True)

        assert before != after
        packed.forget_open_packs()


class TestForgetPreparedBank:
    def test_it_empties_the_slot(self):
        SF_functions._prepared["bank"] = ("sn", "gal", [], [])
        SF_functions._prepared_key = ("something",)

        SF_functions.forget_prepared_bank()

        assert SF_functions._prepared == {}
        assert SF_functions._prepared_key is None


class TestSessionSettings:
    def test_settings_are_shared_across_fits(self):
        session = Session(bank="legacy", resolution=30)

        assert session.settings == {"bank": "legacy", "resolution": 30}

    def test_settings_can_be_changed(self):
        session = Session(resolution=10).set(resolution=30, z=0.1)

        assert session.settings["resolution"] == 30
        assert session.settings["z"] == 0.1

    def test_the_returned_settings_are_a_copy(self):
        session = Session(resolution=10)

        session.settings["resolution"] = 999

        assert session.settings["resolution"] == 10

    def test_a_closed_session_refuses_to_fit(self):
        session = Session()
        session.close()

        with pytest.raises(RuntimeError, match="closed"):
            session.fit("spectrum.flm")

    def test_it_works_as_a_context_manager(self):
        with Session(resolution=10) as session:
            assert session.settings["resolution"] == 10

        with pytest.raises(RuntimeError):
            session.fit("x.flm")

    def test_revalidate_asks_for_one_more_check(self):
        session = Session()
        session._validated = True

        session.revalidate()

        assert session._validated is False

    def test_the_repr_says_how_many_fits_and_which_bank(self):
        session = Session(bank="modern")

        assert "modern" in repr(session)
        assert "0 fit" in repr(session)
