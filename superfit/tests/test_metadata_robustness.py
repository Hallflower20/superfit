"""What the bank scan does when a bank is not shaped like the legacy one.

The legacy bank is 188 objects, every one of them listed in the phase table
and every one with a populated wiserep CSV, so none of these cases could
arise. The modern banks are 1479 and 10069 objects built from several
surveys, and all of them arise there: 51 objects whose CSV is a header with
no rows, and thousands the phase table does not cover.

Each of these used to end the whole fit with an exception from inside a
dictionary lookup. None of them is an error -- they all mean "this object has
less metadata than the legacy bank had", which the code already had an answer
for ('u' for an unknown phase). These tests pin that answer down.
"""

import numpy as np
import pytest

from superfit.get_metadata import Metadata, _as_float


WISEREP_HEADER = (
    "Obj. ID,IAU name,Obj. Type,Redshift,JD,Obs-date,Instrument,Ascii file,"
    "WL Medium,Obj. RA,Obj. DEC\n"
)
WISEREP_ROW = (
    "1,SN2011fe,Ia-norm,0.001,2455815.5,2011-09-01,KAST,a.ascii,Air,210.77,54.27\n"
)

PHASE_TABLE = "Name,mjd_peak,band_peak,isupperlimit\nSN2011fe,55814.0,B,0\n"


class FakeParameters:
    """Only the three attributes the scan reads."""

    def __init__(self, epoch_low=0, epoch_high=0, temp_sn_tr=("Ia-norm",)):
        self.epoch_low = epoch_low
        self.epoch_high = epoch_high
        self.temp_sn_tr = list(temp_sn_tr)


@pytest.fixture
def bank_with(tmp_path, monkeypatch):
    """Build a one-object bank, and point the scan at it."""

    def build(wiserep, phase_table=PHASE_TABLE, obj="SN2011fe"):
        objdir = tmp_path / "original_resolution" / "sne" / "Ia-norm" / obj
        objdir.mkdir(parents=True)
        (objdir / "a.ascii").write_text("4000 1.0\n")
        (objdir / "wiserep_spectra.csv").write_text(wiserep)

        table = tmp_path / "mjd_of_maximum_brightness.csv"
        table.write_text(phase_table)

        import superfit.get_metadata as gm

        monkeypatch.setattr(
            gm, "sne_dir", lambda *a, **k: str(tmp_path / "original_resolution" / "sne")
        )
        monkeypatch.setattr(gm, "mjd_max_brightness_csv", lambda *a, **k: str(table))
        return tmp_path

    return build


class TestAsFloat:
    @pytest.mark.parametrize("value", [None, "", "  ", "nope", float("nan"), np.nan])
    def test_things_that_are_not_a_number(self, value):
        assert _as_float(value) is None

    @pytest.mark.parametrize(
        "value,expected", [("55814.0", 55814.0), (3, 3.0), (-1, -1.0), ("-1", -1.0)]
    )
    def test_things_that_are(self, value, expected):
        assert _as_float(value) == expected


class TestEmptyWiserep:
    def test_a_header_with_no_rows_does_not_end_the_fit(self, bank_with):
        """One bank has 51 of these; .iloc[0] on them raised IndexError."""

        bank_with(WISEREP_HEADER)

        metadata = Metadata(FakeParameters())

        assert metadata.dictionary_all_trunc_objects == {}

    def test_the_object_is_recorded_as_having_no_metadata(self, bank_with):
        root = bank_with(WISEREP_HEADER)

        metadata = Metadata(FakeParameters())

        assert any("SN2011fe" in p for p in metadata.no_wiserep)

    def test_a_populated_one_still_works(self, bank_with):
        bank_with(WISEREP_HEADER + WISEREP_ROW)

        metadata = Metadata(FakeParameters())

        assert list(metadata.dictionary_all_trunc_objects) == ["a.ascii"]


class TestMissingPhase:
    def test_an_object_the_table_does_not_list_gets_an_unknown_phase(self, bank_with):
        """Thousands of objects in the modern banks; this used to KeyError."""

        bank_with(
            WISEREP_HEADER + WISEREP_ROW,
            phase_table="Name,mjd_peak,band_peak,isupperlimit\nsomeone-else,1.0,B,0\n",
        )

        metadata = Metadata(FakeParameters())

        assert "phase-band : u" in metadata.shorhand_dict["a.ascii"]

    def test_the_minus_one_sentinel_still_means_unknown(self, bank_with):
        bank_with(
            WISEREP_HEADER + WISEREP_ROW,
            phase_table="Name,mjd_peak,band_peak,isupperlimit\nSN2011fe,-1,B,0\n",
        )

        metadata = Metadata(FakeParameters())

        assert "phase-band : u" in metadata.shorhand_dict["a.ascii"]

    def test_a_known_maximum_still_gives_a_phase(self, bank_with):
        bank_with(WISEREP_HEADER + WISEREP_ROW)

        metadata = Metadata(FakeParameters())

        # JD 2455815.5 against MJD 55814.0 -> 2455814.5, so one day after max.
        assert "phase-band : 1.0B" in metadata.shorhand_dict["a.ascii"]

    def test_a_spectrum_with_no_epoch_of_its_own_is_unknown_not_nan(self, bank_with):
        """695 rows in one bank have a blank JD; nan propagated silently."""

        bank_with(WISEREP_HEADER + WISEREP_ROW.replace("2455815.5", ""))

        metadata = Metadata(FakeParameters())

        assert "phase-band : u" in metadata.shorhand_dict["a.ascii"]

    def test_an_unknown_phase_is_excluded_by_an_epoch_window(self, bank_with):
        """'u' is not a number; it cannot satisfy a window and must not try."""

        bank_with(
            WISEREP_HEADER + WISEREP_ROW,
            phase_table="Name,mjd_peak,band_peak,isupperlimit\nsomeone-else,1.0,B,0\n",
        )

        metadata = Metadata(FakeParameters(epoch_low=-5, epoch_high=5))

        assert metadata.dictionary_all_trunc_objects == {}

    def test_a_missing_band_leaves_the_label_readable(self, bank_with):
        bank_with(
            WISEREP_HEADER + WISEREP_ROW,
            phase_table="Name,mjd_peak,band_peak,isupperlimit\nother,1.0,B,0\n",
        )

        metadata = Metadata(FakeParameters())

        label = metadata.shorhand_dict["a.ascii"]
        assert label.endswith("phase-band : u")
