"""Choosing between template banks by name.

There is more than one bank now, and they differ scientifically: which one
produced a classification is part of the result. So a bank is named, the name
is recorded in the run's config, and naming one that is not there stops rather
than quietly fitting against whichever bank the machine happened to find.

Everything here builds small banks in a tmp_path; none of it needs a real one.
"""

import json

import pytest

from superfit import bank, paths


TEMPLATE = "4000.0 1.5\n4010.0 1.6\n"
PHASE_TABLE = "Name,mjd_peak,band_peak,isupperlimit\n2011fe,55814.0,B,0\n"

# The shape both modern banks ship: the Name column dropped, so every key is
# a maximum-light MJD and no object ever matches.
HEADLESS_PHASE_TABLE = "mjd_peak,band_peak,isupperlimit\n55814.0,B,0\n"


def make_bank(root, phase_table=None):
    """The smallest directory tree that counts as a template bank."""

    for relative in (
        "binnings/10A/sne/Ia-norm/SN2011fe/a.ascii",
        "binnings/10A/gal/E",
        "original_resolution/sne/Ia-norm/SN2011fe/a.ascii",
        "original_resolution/gal/E",
    ):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(TEMPLATE)

    if phase_table is not None:
        (root / "mjd_of_maximum_brightness.csv").write_text(phase_table)

    return root


@pytest.fixture
def home(tmp_path, monkeypatch):
    """An isolated per-user data directory, so the registry is not the real one."""

    monkeypatch.setattr(bank.sys, "platform", "linux")
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "share"))
    monkeypatch.delenv("SUPERFIT_BANK_DIR", raising=False)
    monkeypatch.delenv("NGSF_BANK_DIR", raising=False)
    monkeypatch.delenv("SUPERFIT_BANK", raising=False)
    monkeypatch.setattr(paths, "_cached_bank_dir", None)
    return tmp_path


class TestRegistry:
    def test_a_bank_starts_out_not_installed(self, home):
        assert bank.resolve_bank("modern") is None
        assert bank.installed_banks() == {}

    def test_registering_points_a_name_at_a_directory(self, home):
        directory = make_bank(home / "somewhere")

        bank.register("modern", directory)

        assert bank.resolve_bank("modern") == str(directory.resolve())
        assert bank.installed_banks()["modern"] == str(directory.resolve())

    def test_it_survives_being_written_and_read_again(self, home):
        make_bank(home / "b")
        bank.register("modern", home / "b")

        assert json.loads(bank.registry_path().read_text())["banks"]["modern"]

    def test_forgetting_a_name_leaves_the_bank_on_disk(self, home):
        directory = make_bank(home / "b")
        bank.register("modern", directory)

        bank.unregister("modern")

        assert bank.resolve_bank("modern") is None
        assert (directory / "binnings").is_dir()

    def test_an_unknown_name_is_an_error_not_a_shrug(self, home):
        with pytest.raises(bank.UnknownBank):
            bank.resolve_bank("no-such-bank")

    def test_a_registered_name_need_not_be_one_superfit_ships(self, home):
        directory = make_bank(home / "mine")
        bank.register("mine", directory)

        assert bank.resolve_bank("mine") == str(directory.resolve())

    def test_a_registry_pointing_at_a_deleted_bank_reads_as_absent(self, home):
        directory = make_bank(home / "b")
        bank.register("modern", directory)

        import shutil

        shutil.rmtree(directory)

        assert bank.resolve_bank("modern") is None


class TestInstallFromDirectory:
    def test_it_records_the_path_rather_than_copying(self, home):
        source = make_bank(home / "shared")

        directory = bank.install_from_directory("modern", source, quiet=True)

        assert directory == source.resolve()
        assert not (bank.user_data_dir() / "banks" / "modern").exists()

    def test_copy_materialises_it_under_the_data_directory(self, home):
        source = make_bank(home / "shared")

        directory = bank.install_from_directory(
            "modern", source, copy=True, quiet=True
        )

        assert directory == bank.bank_install_dir("modern")
        assert (directory / "binnings" / "10A" / "gal" / "E").read_text() == TEMPLATE

    def test_a_directory_that_is_not_a_bank_is_refused(self, home):
        (home / "junk").mkdir()

        with pytest.raises(bank.BankError, match="does not look like"):
            bank.install_from_directory("modern", home / "junk", quiet=True)

    def test_a_missing_directory_is_refused(self, home):
        with pytest.raises(bank.BankError, match="No such directory"):
            bank.install_from_directory("modern", home / "nope", quiet=True)


class TestSelection:
    def test_a_name_beats_the_search_path(self, home, monkeypatch):
        named = make_bank(home / "named")
        make_bank(home / "cwd_bank")
        bank.register("modern", named)
        monkeypatch.setenv("SUPERFIT_BANK_DIR", str(home / "cwd_bank"))

        assert paths.find_bank_dir(name="modern") == str(named.resolve())

    def test_the_environment_name_is_honoured(self, home, monkeypatch):
        named = make_bank(home / "named")
        bank.register("modern", named)
        monkeypatch.setenv("SUPERFIT_BANK", "modern")

        assert paths.find_bank_dir(use_cache=False) == str(named.resolve())

    def test_an_explicit_directory_beats_an_environment_name(self, home, monkeypatch):
        named = make_bank(home / "named")
        explicit = make_bank(home / "explicit")
        bank.register("modern", named)
        monkeypatch.setenv("SUPERFIT_BANK", "modern")
        monkeypatch.setenv("SUPERFIT_BANK_DIR", str(explicit))

        assert paths.find_bank_dir(use_cache=False) == str(explicit)

    def test_naming_a_bank_that_is_not_installed_stops(self, home):
        with pytest.raises(FileNotFoundError, match="not installed"):
            paths.find_bank_dir(name="modern")


class TestBankIsPerFitNotPerProcess:
    """Constructing a second fit must not move the first one's bank.

    `a = Superfit(bank="legacy"); b = Superfit(bank="modern"); a.run()` used to
    fit `a` against modern's supernovae and legacy's galaxies, and label the
    result "legacy". Bank selection happens during construction, so the fit
    lock does not help; the only fix is for nothing about a fit to be read
    from process-wide state.
    """

    def test_resolving_a_named_bank_does_not_move_the_process_bank(
        self, home, monkeypatch
    ):
        make_bank(home / "searched")
        named = make_bank(home / "named")
        bank.register("modern", named)
        monkeypatch.setenv("SUPERFIT_BANK_DIR", str(home / "searched"))

        before = paths.find_bank_dir(use_cache=False)
        paths.find_bank_dir(name="modern")

        assert paths.find_bank_dir(use_cache=False) == before

    def test_two_parameters_keep_their_own_banks(self, home, monkeypatch):
        from superfit.params import Parameters

        legacy = make_bank(home / "legacy_bank", phase_table=PHASE_TABLE)
        modern = make_bank(home / "modern_bank", phase_table=PHASE_TABLE)
        bank.register("legacy-test", legacy)
        bank.register("modern-test", modern)

        base = dict(
            object_to_fit="", lower_lam=4000, upper_lam=5000, resolution=10,
            mask_galaxy_lines=False,
        )

        a = Parameters(dict(base, bank="legacy-test"))
        b = Parameters(dict(base, bank="modern-test"))

        # a's identity must be unchanged by b having been built afterwards.
        assert a.bank_dir == str(legacy.resolve())
        assert b.bank_dir == str(modern.resolve())
        assert a.phase_table.startswith(a.bank_dir)
        assert b.phase_table.startswith(b.bank_dir)

    def test_the_metadata_cache_key_follows_the_fit_not_the_process(
        self, home, monkeypatch
    ):
        """Two banks must not share a cached scan."""

        from superfit.params import Parameters

        legacy = make_bank(home / "legacy_bank")
        modern = make_bank(home / "modern_bank")
        bank.register("legacy-test", legacy)
        bank.register("modern-test", modern)

        base = dict(
            object_to_fit="", lower_lam=4000, upper_lam=5000, resolution=10,
            mask_galaxy_lines=False,
        )

        a = Parameters(dict(base, bank="legacy-test"))
        b = Parameters(dict(base, bank="modern-test"))

        assert a.metadata_key != b.metadata_key
        assert a.metadata_key[0] == a.bank_dir

    def test_template_lists_come_from_the_named_bank(self, home):
        from superfit.params import Parameters

        legacy = make_bank(home / "legacy_bank")
        modern = make_bank(home / "modern_bank")
        # Give modern a second galaxy so the two are distinguishable.
        (modern / "binnings" / "10A" / "gal" / "Sa").write_text(TEMPLATE)
        bank.register("legacy-test", legacy)
        bank.register("modern-test", modern)

        base = dict(
            object_to_fit="", lower_lam=4000, upper_lam=5000, resolution=10,
            mask_galaxy_lines=False, temp_gal_tr=["E", "Sa"],
        )

        a = Parameters(dict(base, bank="legacy-test"))
        b = Parameters(dict(base, bank="modern-test"))

        assert len(a.templates_gal_trunc) == 1
        assert len(b.templates_gal_trunc) == 2
        assert all(str(legacy) in str(p) for p in a.templates_gal_trunc)


class TestPhaseTable:
    def test_a_banks_own_table_is_preferred(self, home):
        directory = make_bank(home / "b", phase_table=PHASE_TABLE)

        assert paths.mjd_max_brightness_csv(str(directory)) == str(
            directory / "mjd_of_maximum_brightness.csv"
        )

    def test_a_bank_without_one_falls_back_to_the_package_copy(self, home):
        directory = make_bank(home / "b")

        assert (
            paths.mjd_max_brightness_csv(str(directory))
            == paths.MJD_MAX_BRIGHTNESS_CSV
        )

    def test_a_table_with_no_name_column_is_refused_loudly(self, home):
        """Every lookup is by name; without that column nothing ever matches."""

        directory = make_bank(home / "b", phase_table=HEADLESS_PHASE_TABLE)

        with pytest.warns(RuntimeWarning, match="Name"):
            table = paths.mjd_max_brightness_csv(str(directory))

        assert table == paths.MJD_MAX_BRIGHTNESS_CSV

    def test_an_empty_table_is_refused(self, home):
        directory = make_bank(home / "b", phase_table="")

        with pytest.warns(RuntimeWarning):
            assert (
                paths.mjd_max_brightness_csv(str(directory))
                == paths.MJD_MAX_BRIGHTNESS_CSV
            )

    def test_the_check_reports_why(self, home):
        directory = make_bank(home / "b", phase_table=HEADLESS_PHASE_TABLE)
        path = str(directory / "mjd_of_maximum_brightness.csv")

        ok, reason = paths.phase_table_is_usable(path)

        assert not ok
        assert "Name" in reason

    def test_a_good_table_passes(self, home):
        directory = make_bank(home / "b", phase_table=PHASE_TABLE)
        path = str(directory / "mjd_of_maximum_brightness.csv")

        assert paths.phase_table_is_usable(path) == (True, None)
