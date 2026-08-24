"""The command line: what each flag turns into, and what still works.

None of this needs the template bank. The point is the translation layer --
a flag that quietly maps to the wrong setting is a fit that answers a
different question, and that is exactly the kind of bug that survives a
green end-to-end suite.
"""

import json

import pytest

from superfit.cli import build_parser, fit_overrides, rewrite_legacy
from superfit.config import ConfigError, load_config


def parse(*argv):
    return build_parser().parse_args(list(argv))


def overrides_for(*argv):
    """The configuration a `superfit fit ...` command line asks for."""

    args = parse("fit", "spectrum.flm", *argv)
    overrides, _ = fit_overrides(args)
    return load_config(**overrides)


class TestTheOneLinerWorks:
    """`superfit fit spectrum.flm --z 0.127` is the whole point."""

    def test_a_spectrum_and_a_redshift_is_enough(self):
        args = parse("fit", "spectrum.flm", "--z", "0.127")

        assert args.spectrum == "spectrum.flm"
        assert args.z == 0.127
        assert args.config is None

    def test_it_becomes_an_exact_redshift_fit(self):
        config = overrides_for("--z", "0.127")

        assert config["z_exact"] == 0.127
        assert config["use_exact_z"] is True

    def test_everything_else_takes_its_default(self):
        from superfit.config import DEFAULT_CONFIG

        config = overrides_for("--z", "0.127")

        for key in ("resolution", "error_spectrum", "minimum_overlap", "R_v"):
            assert config[key] == DEFAULT_CONFIG[key]


class TestRedshift:
    def test_scan_z_sets_the_range(self):
        config = overrides_for("--scan-z", "0.0", "0.3")

        assert config["use_exact_z"] is False
        assert config["z_range_begin"] == 0.0
        assert config["z_range_end"] == 0.3

    def test_z_step_sets_the_step(self):
        assert overrides_for("--scan-z", "0", "0.3", "--z-step", "0.005")["z_int"] == (
            0.005
        )

    def test_scan_z_turns_off_host_line_masking(self):
        """Otherwise every scan is rejected by the validator, which helps nobody."""

        args = parse("fit", "s.flm", "--scan-z", "0", "0.3")
        overrides, notes = fit_overrides(args)

        assert overrides["mask_galaxy_lines"] is False
        assert any("masking" in note for note in notes)

    def test_the_note_is_not_printed_when_the_user_said_it_first(self):
        args = parse("fit", "s.flm", "--scan-z", "0", "0.3", "--no-mask-galaxy-lines")
        _, notes = fit_overrides(args)

        assert notes == []

    def test_z_and_scan_z_together_is_refused(self):
        args = parse("fit", "s.flm", "--z", "0.1", "--scan-z", "0", "0.3")

        with pytest.raises(ValueError, match="not both"):
            fit_overrides(args)

    def test_z_step_without_a_scan_is_refused(self):
        args = parse("fit", "s.flm", "--z", "0.1", "--z-step", "0.01")

        with pytest.raises(ValueError, match="only means something with --scan-z"):
            fit_overrides(args)


class TestFlagsMapToSettings:
    @pytest.mark.parametrize(
        "argv,key,value",
        [
            (["--resolution", "30"], "resolution", 30),
            (["--error-model", "linear"], "error_spectrum", "linear"),
            (["--error-model", "included"], "error_spectrum", "included"),
            (["--plots", "5"], "how_many_plots", 5),
            (["--n-cores", "4"], "n_cores", 4),
            (["--rv", "2.5"], "R_v", 2.5),
            (["--output", "results"], "saving_results_path", "results"),
            (["--overwrite"], "overwrite", True),
            (["--png"], "show_plot_png", True),
            (["--no-mask-telluric"], "mask_telluric", False),
            (["--no-mask-galaxy-lines"], "mask_galaxy_lines", False),
        ],
    )
    def test_flag(self, argv, key, value):
        assert overrides_for("--z", "0.1", *argv)[key] == value

    def test_av_range(self):
        config = overrides_for("--z", "0.1", "--av-range", "-1", "3")

        assert (config["Alam_low"], config["Alam_high"]) == (-1.0, 3.0)

    def test_epochs(self):
        config = overrides_for("--z", "0.1", "--epochs", "-10", "20")

        assert (config["epoch_low"], config["epoch_high"]) == (-10.0, 20.0)

    def test_show_implies_something_to_show(self):
        """Plotting is off by default, so --show alone would show nothing."""

        assert overrides_for("--z", "0.1", "--show")["how_many_plots"] >= 1

    def test_show_respects_an_explicit_count(self):
        config = overrides_for("--z", "0.1", "--show", "--plots", "4")

        assert config["how_many_plots"] == 4

    def test_output_directory_is_not_glued_to_the_filename(self):
        """`--out results` used to need a trailing slash or it wrote elsewhere."""

        assert overrides_for("--z", "0.1", "-o", "results")[
            "saving_results_path"
        ] == "results"


class TestProfiles:
    def test_legacy_is_the_default(self):
        config = load_config()

        assert config["profile"] == "legacy"
        assert config["weighted_solve"] is False

    def test_modern_weights_the_solve(self):
        config = load_config(profile="modern")

        assert config["weighted_solve"] is True

    def test_weighted_solve_is_shorthand_for_modern(self):
        args = parse("fit", "s.flm", "--z", "0.1", "--weighted-solve")

        assert args.weighted_solve is True

    def test_a_profile_named_in_a_file_is_honoured(self):
        assert load_config({"profile": "modern"})["weighted_solve"] is True

    def test_an_unknown_profile_is_named_in_the_error(self):
        with pytest.raises(ConfigError, match="beta"):
            load_config(profile="beta")

    def test_the_profile_is_recorded(self):
        """So a run directory says which science it was run with."""

        assert load_config(profile="modern")["profile"] == "modern"

    def test_a_setting_that_contradicts_the_profile_is_not_recorded_as_it(self):
        """config.json has to describe the settings that came out."""

        assert load_config({"weighted_solve": True})["profile"] == "custom"
        assert load_config(profile="modern", weighted_solve=False)["profile"] == (
            "custom"
        )

    def test_a_custom_config_still_reloads(self):
        """Which is the whole point of writing config.json."""

        recorded = load_config({"weighted_solve": True})

        assert load_config(recorded) == recorded


class TestBackwardsCompatibility:
    def test_a_bare_json_path_still_means_a_config(self):
        assert rewrite_legacy(["parameters.json"]) == [
            "fit", "--config", "parameters.json",
        ]

    def test_a_bare_json_string_still_means_a_config(self):
        assert rewrite_legacy(['{"resolution": 10}']) == [
            "fit", "--config", '{"resolution": 10}',
        ]

    def test_old_flags_are_carried_along(self):
        assert rewrite_legacy(["params.json", "--out", "here"]) == [
            "fit", "--config", "params.json", "--out", "here",
        ]

    def test_subcommands_are_left_alone(self):
        assert rewrite_legacy(["fit", "s.flm"]) == ["fit", "s.flm"]
        assert rewrite_legacy(["bank", "status"]) == ["bank", "status"]

    def test_options_are_left_alone(self):
        assert rewrite_legacy(["--help"]) == ["--help"]
        assert rewrite_legacy([]) == []

    @pytest.mark.parametrize("argv", [["--object", "s.flm"], ["--out", "here"]])
    def test_the_previous_flag_names_still_parse(self, argv):
        parse("fit", *argv)

    def test_object_still_names_the_spectrum(self):
        args = parse("fit", "--object", "s.flm", "--z", "0.1")

        assert (args.spectrum or args.spectrum_legacy) == "s.flm"

    def test_out_still_names_the_output_directory(self):
        args = parse("fit", "s.flm", "--out", "here")
        overrides, _ = fit_overrides(args)

        assert overrides["output_dir"] == "here"

    def test_no_plots_still_silences_plotting(self):
        assert overrides_for("--z", "0.1", "--no-plots")["how_many_plots"] == 0


class TestSubcommandWiring:
    @pytest.mark.parametrize(
        "argv,handler",
        [
            (["fit", "s.flm"], "run_fit"),
            (["bank", "install"], "run_bank_install"),
            (["bank", "status"], "run_bank_status"),
            (["doctor"], "run_doctor"),
            (["config", "create"], "run_config_create"),
            (["config", "show"], "run_config_show"),
        ],
    )
    def test_each_command_reaches_its_handler(self, argv, handler):
        assert parse(*argv).handler.__name__ == handler

    def test_no_command_prints_help_rather_than_crashing(self, capsys):
        from superfit.cli import main

        assert main([]) == 2
        assert "superfit fit spectrum.flm" in capsys.readouterr().out

    def test_an_unknown_command_is_rejected(self):
        with pytest.raises(SystemExit):
            parse("frobnicate")


class TestConfigCreate:
    def run(self, *argv):
        from superfit.cli import main

        return main(["config", "create", *argv])

    def test_it_writes_a_file_that_loads(self, tmp_path, capsys):
        path = tmp_path / "parameters.json"

        assert self.run(str(path)) == 0
        capsys.readouterr()

        written = json.loads(path.read_text())
        load_config(written)

    def test_the_comments_are_ignored_when_it_loads(self, tmp_path, capsys):
        path = tmp_path / "parameters.json"
        self.run(str(path))
        capsys.readouterr()

        written = json.loads(path.read_text())
        assert any(key.startswith("_") for key in written), "nothing explains itself"

        config = load_config(written)
        assert not any(key.startswith("_") for key in config)

    def test_every_setting_is_explained(self, tmp_path, capsys):
        path = tmp_path / "parameters.json"
        self.run(str(path))
        capsys.readouterr()

        written = json.loads(path.read_text())
        for key in written:
            if key.startswith("_"):
                continue
            assert "_" + key in written, "{} has no explanation".format(key)

    def test_full_covers_every_setting(self, tmp_path, capsys):
        from superfit.config import DEFAULT_CONFIG

        path = tmp_path / "parameters.json"
        self.run(str(path), "--full")
        capsys.readouterr()

        written = json.loads(path.read_text())
        real = {k for k in written if not k.startswith("_")}
        assert real == set(DEFAULT_CONFIG)

    def test_an_existing_file_is_not_clobbered(self, tmp_path, capsys):
        path = tmp_path / "parameters.json"
        path.write_text("mine")

        assert self.run(str(path)) == 1
        assert path.read_text() == "mine"
        assert "--force" in capsys.readouterr().err

    def test_force_overwrites(self, tmp_path, capsys):
        path = tmp_path / "parameters.json"
        path.write_text("mine")

        assert self.run(str(path), "--force") == 0
        assert json.loads(path.read_text())

    def test_the_profile_is_carried_into_the_file(self, tmp_path, capsys):
        path = tmp_path / "parameters.json"
        self.run(str(path), "--profile", "modern")
        capsys.readouterr()

        assert json.loads(path.read_text())["profile"] == "modern"


class TestErrorsAreReported:
    def test_a_bad_setting_exits_nonzero_with_a_message(self, capsys):
        from superfit.cli import main

        assert main(["fit", "s.flm", "--z", "0.1", "--resolution", "0"]) == 1
        assert "resolution must be positive" in capsys.readouterr().err

    def test_no_spectrum_and_no_config_is_a_clear_error(self, capsys):
        from superfit.cli import main

        assert main(["fit"]) == 1
        assert "No spectrum to fit" in capsys.readouterr().err
