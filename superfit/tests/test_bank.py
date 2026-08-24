"""Installing the template bank.

Everything here works on a small archive built in the test, so none of it
touches the network or the real 74 MB download. What is being checked is the
part that can go wrong quietly: a truncated or substituted archive, a
half-finished install left behind by an interruption, and an archive that
tries to write outside the directory it was pointed at.
"""

import io
import json
import zipfile

import pytest

from superfit import bank


def make_archive(path, entries, prefix="bank/"):
    """A zip shaped like the published one: everything under ``bank/``."""

    with zipfile.ZipFile(path, "w") as archive:
        for name, content in entries.items():
            archive.writestr(prefix + name, content)
    return path


@pytest.fixture
def tiny_bank(tmp_path):
    """The smallest archive that counts as a template bank."""

    return make_archive(
        tmp_path / "bank.zip",
        {
            "original_resolution/sne/Ia-norm/SN2011fe/spectrum.dat": "4000 1.0\n",
            "original_resolution/gal/E": "4000 1.0\n",
            "binnings/10A/sne/Ia-norm/SN2011fe/spectrum.dat": "4000 1.0\n",
            "binnings/10A/gal/E": "4000 1.0\n",
            ".DS_Store": "junk",
        },
    )


class TestInstallDirectory:
    def test_it_is_under_the_platform_data_directory(self, monkeypatch):
        monkeypatch.setattr(bank.sys, "platform", "linux")
        monkeypatch.setenv("XDG_DATA_HOME", "/data")

        assert bank.default_install_dir().as_posix() == "/data/superfit/bank"

    def test_paths_looks_there(self, monkeypatch, tmp_path):
        """An install with no environment variable still has to be found."""

        from superfit import paths

        monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path))
        monkeypatch.setattr(paths.sys, "platform", "linux")

        assert paths._user_data_bank_dir() == str(tmp_path / "superfit" / "bank")
        assert paths._user_data_bank_dir() in list(paths._candidate_bank_dirs())


class TestInstall:
    def test_it_unpacks_a_usable_bank(self, tiny_bank, tmp_path):
        destination = tmp_path / "installed"

        bank.install(destination, archive=tiny_bank, expected_sha256=None, quiet=True)

        assert (destination / "original_resolution" / "gal" / "E").is_file()
        assert (destination / "binnings" / "10A" / "gal" / "E").is_file()

    def test_the_leading_bank_directory_is_stripped(self, tiny_bank, tmp_path):
        """The archive holds `bank/binnings/...`; a nested bank/bank/ is wrong."""

        destination = tmp_path / "installed"
        bank.install(destination, archive=tiny_bank, expected_sha256=None, quiet=True)

        assert not (destination / "bank").exists()

    def test_packaging_debris_is_left_out(self, tmp_path):
        archive = make_archive(
            tmp_path / "b.zip",
            {
                "original_resolution/gal/E": "x",
                "binnings/10A/gal/E": "x",
                ".DS_Store": "junk",
            },
        )
        with zipfile.ZipFile(archive, "a") as handle:
            handle.writestr("__MACOSX/._bank", "junk")

        destination = tmp_path / "installed"
        bank.install(destination, archive=archive, expected_sha256=None, quiet=True)

        assert not (destination / "__MACOSX").exists()
        assert not (destination / ".DS_Store").exists()

    def test_a_manifest_records_what_was_installed(self, tiny_bank, tmp_path):
        destination = tmp_path / "installed"

        manifest = bank.install(
            destination, archive=tiny_bank, expected_sha256=None, quiet=True
        )

        assert manifest["sha256"] == bank.sha256_of(tiny_bank)
        assert manifest["source"] == str(tiny_bank)
        on_disk = json.loads((destination / bank.MANIFEST_NAME).read_text())
        assert on_disk == manifest

    def test_the_checksum_is_checked(self, tiny_bank, tmp_path):
        with pytest.raises(bank.BankError, match="Checksum mismatch"):
            bank.install(
                tmp_path / "installed",
                archive=tiny_bank,
                expected_sha256="0" * 64,
                quiet=True,
            )

    def test_a_failed_checksum_installs_nothing(self, tiny_bank, tmp_path):
        destination = tmp_path / "installed"

        with pytest.raises(bank.BankError):
            bank.install(
                destination, archive=tiny_bank, expected_sha256="0" * 64, quiet=True
            )

        assert not destination.exists()

    def test_an_existing_installation_is_not_replaced_silently(
        self, tiny_bank, tmp_path
    ):
        destination = tmp_path / "installed"
        destination.mkdir()
        (destination / "keep").write_text("mine")

        with pytest.raises(bank.BankError, match="already exists"):
            bank.install(
                destination, archive=tiny_bank, expected_sha256=None, quiet=True
            )

        assert (destination / "keep").read_text() == "mine"

    def test_overwrite_replaces_it(self, tiny_bank, tmp_path):
        destination = tmp_path / "installed"
        destination.mkdir()
        (destination / "stale") .write_text("old")

        bank.install(
            destination,
            archive=tiny_bank,
            expected_sha256=None,
            overwrite=True,
            quiet=True,
        )

        assert not (destination / "stale").exists()
        assert (destination / "binnings").is_dir()

    def test_an_archive_that_is_not_a_bank_is_refused(self, tmp_path):
        archive = make_archive(tmp_path / "b.zip", {"readme.txt": "hello"})

        with pytest.raises(bank.BankError, match="does not look like"):
            bank.install(
                tmp_path / "installed",
                archive=archive,
                expected_sha256=None,
                quiet=True,
            )

    def test_a_partial_install_is_not_left_behind(self, tmp_path):
        """An interrupted unpack must not leave a bank the fitter would use."""

        archive = make_archive(tmp_path / "b.zip", {"readme.txt": "hello"})
        destination = tmp_path / "installed"

        with pytest.raises(bank.BankError):
            bank.install(
                destination, archive=archive, expected_sha256=None, quiet=True
            )

        assert not destination.exists()
        assert list(tmp_path.glob(".installed-*")) == []

    def test_a_missing_archive_says_so(self, tmp_path):
        with pytest.raises(bank.BankError, match="No such archive"):
            bank.install(tmp_path / "x", archive=tmp_path / "nope.zip", quiet=True)

    def test_something_that_is_not_a_zip_says_so(self, tmp_path):
        not_a_zip = tmp_path / "b.zip"
        not_a_zip.write_text("<html>404 Not Found</html>")

        with pytest.raises(bank.BankError, match="not a readable zip"):
            bank.install(
                tmp_path / "installed",
                archive=not_a_zip,
                expected_sha256=None,
                quiet=True,
            )


class TestHostileArchives:
    @pytest.mark.parametrize("name", ["../escaped", "bank/../../escaped"])
    def test_paths_that_escape_the_target_are_refused(self, tmp_path, name):
        archive = tmp_path / "b.zip"
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr("bank/binnings/10A/gal/E", "x")
            handle.writestr(name, "gotcha")

        with pytest.raises(bank.BankError, match="outside the target directory"):
            bank.install(
                tmp_path / "installed",
                archive=archive,
                expected_sha256=None,
                quiet=True,
            )


class TestChecksum:
    def test_it_matches_hashlib(self, tmp_path):
        import hashlib

        path = tmp_path / "f"
        path.write_bytes(b"superfit" * 1000)

        assert bank.sha256_of(path) == hashlib.sha256(b"superfit" * 1000).hexdigest()

    def test_it_reads_in_chunks(self, tmp_path):
        """A 74 MB archive must not be held in memory to be hashed."""

        path = tmp_path / "f"
        path.write_bytes(b"x" * 5000)

        assert bank.sha256_of(path, block=64) == bank.sha256_of(path)


class TestStatus:
    def test_it_describes_an_installed_bank(self, tiny_bank, tmp_path):
        destination = tmp_path / "installed"
        bank.install(destination, archive=tiny_bank, expected_sha256=None, quiet=True)

        report = bank.status(destination)

        assert report["found"] and report["complete"]
        assert report["missing"] == []
        assert report["n_sn_types"] == 1
        assert report["n_galaxy_templates"] == 1
        assert report["manifest"]["sha256"]

    def test_it_names_what_is_missing(self, tmp_path):
        half = tmp_path / "half"
        (half / "binnings").mkdir(parents=True)

        report = bank.status(half)

        assert report["found"]
        assert not report["complete"]
        assert report["missing"] == ["original_resolution"]

    def test_a_directory_that_is_not_there(self, tmp_path):
        report = bank.status(tmp_path / "nope")

        assert not report["found"]
        assert "not a directory" in report["error"]

    def test_it_says_where_the_bank_would_go(self):
        assert bank.status(None)["install_dir"].endswith("bank")


class TestDownloadErrors:
    def test_an_http_error_explains_the_way_out(self, tmp_path, monkeypatch):
        from urllib.error import HTTPError

        def boom(*args, **kwargs):
            raise HTTPError("url", 404, "Not Found", {}, io.BytesIO(b""))

        monkeypatch.setattr("urllib.request.urlopen", boom)

        with pytest.raises(bank.BankError, match="--archive"):
            bank.download("http://example.invalid/b.zip", tmp_path / "b.zip")

    def test_an_unreachable_host_explains_the_way_out(self, tmp_path, monkeypatch):
        from urllib.error import URLError

        def boom(*args, **kwargs):
            raise URLError("no route to host")

        monkeypatch.setattr("urllib.request.urlopen", boom)

        with pytest.raises(bank.BankError, match="--archive"):
            bank.download("http://example.invalid/b.zip", tmp_path / "b.zip")

    def test_a_user_agent_is_sent(self, tmp_path, monkeypatch):
        """The server behind the real URL answers 403 to urllib's default."""

        seen = {}

        class Response(io.BytesIO):
            headers = {"Content-Length": "2"}

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

        def fake_urlopen(request, *args, **kwargs):
            seen["user_agent"] = request.get_header("User-agent")
            return Response(b"ok")

        monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)

        bank.download("http://example.invalid/b.zip", tmp_path / "b.zip", quiet=True)

        assert "superfit" in seen["user_agent"]
