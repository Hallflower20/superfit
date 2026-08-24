"""Run directories and the object a fit returns.

These are the parts that decide where a result lands and whether an existing
one survives, so they are tested without the template bank -- they have to
hold for every fit, and nobody should have to download 74 MB to find out.
"""

import numpy as np
import pandas as pd
import pytest

from superfit.output import (
    BINNED_TXT,
    RESULTS_CSV,
    USED_CONFIG_JSON,
    FitResult,
    OutputExistsError,
    RunDirectory,
    safe_name,
)


@pytest.fixture
def table():
    return pd.DataFrame(
        {
            "SN": ["Ia-norm/SN2011fe/x", "II/SN2013ej/y"],
            "Z": [0.127, 0.127],
            "CHI2/dof": [1.2, 3.4],
        }
    )


class TestSafeName:
    @pytest.mark.parametrize(
        "given,expected",
        [
            ("SN2021urb", "SN2021urb"),
            ("SN2021urb_2021-08-06_Keck1.flm", "SN2021urb_2021-08-06_Keck1.flm"),
            ("a spectrum", "a_spectrum"),
            ("../../etc/passwd", "etc_passwd"),
            ("/absolute/path", "absolute_path"),
            ("", "spectrum"),
            ("...", "spectrum"),
        ],
    )
    def test_names_stay_inside_one_directory(self, given, expected):
        assert safe_name(given) == expected


class TestRunDirectory:
    def test_the_directory_is_created(self, tmp_path):
        run = RunDirectory(tmp_path / "does" / "not" / "exist", "SN2021urb")

        assert run.path.is_dir()
        assert run.path == tmp_path / "does" / "not" / "exist" / "SN2021urb"

    def test_each_fit_gets_its_own_folder(self, tmp_path):
        a = RunDirectory(tmp_path, "SN2021urb")
        b = RunDirectory(tmp_path, "SN2011fe")

        assert a.path != b.path
        assert a.path.parent == b.path.parent == tmp_path

    def test_no_base_means_the_working_directory(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)

        assert RunDirectory("", "SN2021urb").path == tmp_path / "SN2021urb"

    def test_a_bare_name_is_not_glued_onto_the_parent(self, tmp_path):
        """`--out results` used to write `resultsSN2021urb.csv` one level up."""

        run = RunDirectory(tmp_path / "results", "SN2021urb")

        assert run.results_csv == tmp_path / "results" / "SN2021urb" / RESULTS_CSV

    def test_existing_results_are_refused(self, tmp_path):
        first = RunDirectory(tmp_path, "SN2021urb")
        first.results_csv.write_text("old, precious\n")

        with pytest.raises(OutputExistsError, match="already holds a fit"):
            RunDirectory(tmp_path, "SN2021urb")

        assert first.results_csv.read_text() == "old, precious\n"

    def test_the_refusal_names_the_way_out(self, tmp_path):
        RunDirectory(tmp_path, "SN2021urb").results_csv.write_text("x")

        with pytest.raises(OutputExistsError, match="--overwrite"):
            RunDirectory(tmp_path, "SN2021urb")

    def test_overwrite_replaces_the_run(self, tmp_path):
        first = RunDirectory(tmp_path, "SN2021urb")
        first.results_csv.write_text("old\n")
        first.binned_txt.write_text("old\n")

        second = RunDirectory(tmp_path, "SN2021urb", overwrite=True)

        assert second.path == first.path
        assert not second.results_csv.exists()
        assert not second.binned_txt.exists()

    def test_overwrite_clears_stale_plots(self, tmp_path):
        """A re-run asking for fewer plots must not leave the old tail behind."""

        run = RunDirectory(tmp_path, "SN2021urb")
        for rank in (1, 2, 3):
            run.plot(rank).write_text("old")

        RunDirectory(tmp_path, "SN2021urb", overwrite=True)

        assert not list(run.path.glob("bestfit_*"))

    def test_overwrite_leaves_everything_else_alone(self, tmp_path):
        run = RunDirectory(tmp_path, "SN2021urb")
        notes = run.path / "notes.txt"
        notes.write_text("mine")
        (run.path / "subdir").mkdir()

        RunDirectory(tmp_path, "SN2021urb", overwrite=True)

        assert notes.read_text() == "mine"
        assert (run.path / "subdir").is_dir()

    def test_a_directory_with_no_results_is_not_a_collision(self, tmp_path):
        """An empty or half-written directory is not a result worth guarding."""

        (tmp_path / "SN2021urb").mkdir()

        RunDirectory(tmp_path, "SN2021urb")

    def test_plot_paths_are_ranked_from_one(self, tmp_path):
        run = RunDirectory(tmp_path, "SN2021urb")

        assert run.plot(1).name == "bestfit_1.pdf"
        assert run.plot(2, png=True).name == "bestfit_2.png"

    def test_the_file_names_are_fixed(self, tmp_path):
        run = RunDirectory(tmp_path, "SN2021urb")

        assert run.results_csv.name == RESULTS_CSV
        assert run.used_config_json.name == USED_CONFIG_JSON
        assert run.binned_txt.name == BINNED_TXT


class TestFitResult:
    def test_it_indexes_like_the_frame_it_wraps(self, tmp_path, table):
        result = FitResult(table, tmp_path)

        assert len(result) == 2
        assert list(result["Z"]) == [0.127, 0.127]
        assert list(result.columns) == ["SN", "Z", "CHI2/dof"]
        assert "SN" in result

    def test_best_is_the_top_row(self, tmp_path, table):
        assert FitResult(table, tmp_path).best["SN"] == "Ia-norm/SN2011fe/x"

    def test_it_lists_what_was_written(self, tmp_path, table):
        artifacts = [tmp_path / RESULTS_CSV, tmp_path / "bestfit_1.pdf"]
        result = FitResult(table, tmp_path, artifacts)

        assert [p.name for p in result.artifacts] == [RESULTS_CSV, "bestfit_1.pdf"]
        assert [p.name for p in result.plots] == ["bestfit_1.pdf"]
        assert result.artifact(RESULTS_CSV) == tmp_path / RESULTS_CSV
        assert result.artifact("nothing.txt") is None

    def test_the_frame_is_still_reachable(self, tmp_path, table):
        result = FitResult(table, tmp_path)

        assert result.results is table
        np.testing.assert_allclose(result.results["CHI2/dof"], [1.2, 3.4])

    def test_repr_says_what_won(self, tmp_path, table):
        assert "SN2011fe" in repr(FitResult(table, tmp_path))

    def test_repr_survives_an_empty_result(self, tmp_path, table):
        assert "no surviving" in repr(FitResult(table.iloc[:0], tmp_path))
