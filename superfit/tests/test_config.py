"""Configuration loading, defaulting and validation.

Every key has a default, so the point of these tests is that the defaults
are complete and that a configuration which cannot work is rejected at load
time with a message naming the problem -- rather than thirty seconds later,
inside a worker process, as a KeyError or an IndexError.
"""

import json

import pytest

from superfit.config import (
    DEFAULT_CONFIG,
    ConfigError,
    load_config,
    read_config,
    resolve_aliases,
)


class TestDefaults:
    def test_an_empty_config_is_complete(self):
        config = load_config()
        assert set(config) == set(DEFAULT_CONFIG)

    def test_defaults_are_not_shared_between_calls(self):
        """A mutable default leaking between fits would be a nasty bug."""

        first = load_config()
        first["temp_sn_tr"].append("bogus")

        assert "bogus" not in load_config()["temp_sn_tr"]
        assert "bogus" not in DEFAULT_CONFIG["temp_sn_tr"]

    def test_overrides_win_over_defaults(self):
        assert load_config(resolution=30)["resolution"] == 30

    def test_keyword_overrides_win_over_the_source(self):
        assert load_config({"resolution": 30}, resolution=10)["resolution"] == 10

    def test_none_valued_overrides_are_ignored(self):
        """argparse hands us None for every flag the user did not pass."""

        assert load_config(resolution=None)["resolution"] == (
            DEFAULT_CONFIG["resolution"]
        )


class TestSources:
    def test_dict(self):
        assert load_config({"resolution": 25})["resolution"] == 25

    def test_json_string(self):
        assert load_config('{"resolution": 25}')["resolution"] == 25

    def test_json_file(self, tmp_path):
        path = tmp_path / "params.json"
        path.write_text(json.dumps({"resolution": 25}))
        assert load_config(str(path))["resolution"] == 25

    def test_unreadable_source_is_reported_with_its_value(self):
        with pytest.raises(ConfigError, match="not_a_config"):
            load_config("not_a_config")

    def test_read_config_applies_no_defaults(self):
        assert read_config({"resolution": 25}) == {"resolution": 25}


class TestAliases:
    def test_z_maps_to_exact_redshift(self):
        config = load_config(z=0.127)
        assert config["z_exact"] == 0.127
        assert config["use_exact_z"] == 1

    def test_a_scan_range_implies_not_exact(self):
        config = load_config(z_range_begin=0.0, z_range_end=0.2, mask_galaxy_lines=0)
        assert config["use_exact_z"] == 0

    def test_explicit_use_exact_z_is_respected(self):
        config = load_config(z=0.1, use_exact_z=1)
        assert config["use_exact_z"] == 1

    @pytest.mark.parametrize(
        "alias,target,value",
        [
            ("output_dir", "saving_results_path", "out/"),
            ("n_plots", "how_many_plots", 3),
            ("error_model", "error_spectrum", "linear"),
            ("cores", "n_cores", 8),
        ],
    )
    def test_shorthand_names(self, alias, target, value):
        assert load_config(**{alias: value})[target] == value

    def test_resolve_aliases_leaves_real_keys_alone(self):
        assert resolve_aliases({"resolution": 10}) == {"resolution": 10}


class TestValidation:
    def test_unknown_key_is_rejected_and_named(self):
        with pytest.raises(ConfigError, match="resolutoin"):
            load_config(resolutoin=10)

    def test_unknown_key_error_lists_the_valid_ones(self):
        with pytest.raises(ConfigError, match="minimum_overlap"):
            load_config(nonsense=1)

    @pytest.mark.parametrize(
        "kwargs,match",
        [
            (dict(resolution=0), "resolution must be positive"),
            (dict(resolution=-5), "resolution must be positive"),
            (dict(minimum_overlap=1.5), "between 0 and 1"),
            (dict(minimum_overlap=-0.1), "between 0 and 1"),
            (dict(error_spectrum="magic"), "must be 'sg'"),
            (dict(Alam_interval=0), "Alam_interval must be positive"),
            (dict(Alam_low=2, Alam_high=-2), "below Alam_low"),
            (dict(temp_sn_tr=[]), "nothing to fit"),
            (dict(temp_gal_tr=[]), "nothing to fit"),
        ],
    )
    def test_unusable_values_are_rejected(self, kwargs, match):
        with pytest.raises(ConfigError, match=match):
            load_config(**kwargs)

    def test_masking_host_lines_during_a_scan_is_rejected(self):
        """The lines are placed at the object's z, so there must be only one."""

        with pytest.raises(ConfigError, match="mask_galaxy_lines needs a single"):
            load_config(use_exact_z=0, mask_galaxy_lines=1)

    def test_scan_with_masking_off_is_fine(self):
        load_config(use_exact_z=0, mask_galaxy_lines=0)

    @pytest.mark.parametrize(
        "kwargs,match",
        [
            (dict(use_exact_z=0, mask_galaxy_lines=0, z_int=0), "z_int must be"),
            (
                dict(
                    use_exact_z=0,
                    mask_galaxy_lines=0,
                    z_range_begin=0.5,
                    z_range_end=0.1,
                ),
                "below z_range_begin",
            ),
        ],
    )
    def test_scan_ranges_are_checked(self, kwargs, match):
        with pytest.raises(ConfigError, match=match):
            load_config(**kwargs)

    def test_the_shipped_parameters_file_is_valid(self, base_parameters):
        """The example config must survive its own validator."""

        load_config(base_parameters)
