"""Configuration defaults and loading for superfit.

Every setting has a default here, so a caller only has to state what they
actually care about -- in practice the spectrum and the redshift. The
effective configuration (defaults plus overrides) is what gets recorded in
``*_used.json`` alongside the results, so a run stays reproducible even
though the input was three keys long.
"""

import copy
import json

# Every supernova subtype in the shipped bank.
ALL_SN_TYPES = [
    "IIb-flash",
    "computed",
    "Ia 02es-like",
    "Ia-02cx like",
    "TDE He",
    "Ca-Ia",
    "Ia-CSM-(ambigious)",
    "II",
    "super_chandra",
    "SLSN-II",
    "IIn",
    "FBOT",
    "Ibn",
    "SLSN-IIn",
    "Ia 91T-like",
    "IIb",
    "TDE H",
    "SN - Imposter",
    "II-flash",
    "ILRT",
    "Ia 99aa-like",
    "Ic",
    "SLSN-I",
    "Ia-pec",
    "Ib",
    "Ia-CSM",
    "Ia-norm",
    "SLSN-Ib",
    "TDE H+He",
    "Ia 91bg-like",
    "Ca-Ib",
    "Ia-rapid",
    "Ic-BL",
    "Ic-pec",
    "SLSN-IIb",
]

# Every host galaxy template in the shipped bank.
ALL_GALAXY_TYPES = [
    "E",
    "S0",
    "Sa",
    "Sb",
    "SB1",
    "SB2",
    "SB3",
    "SB4",
    "SB5",
    "SB6",
    "Sc",
]


DEFAULT_CONFIG = {
    # --- what to fit -----------------------------------------------------
    # None means "the spectrum was supplied directly, not as a path".
    "object_to_fit": None,
    # --- redshift --------------------------------------------------------
    "use_exact_z": 1,
    "z_exact": 0.0,
    "z_range_begin": 0.0,
    "z_range_end": 0.1,
    "z_int": 0.01,
    # --- extinction grid -------------------------------------------------
    "Alam_low": -2.0,
    "Alam_high": 2.0,
    "Alam_interval": 0.2,
    # --- templates -------------------------------------------------------
    "temp_sn_tr": ALL_SN_TYPES,
    "temp_gal_tr": ALL_GALAXY_TYPES,
    "resolution": 10,
    "epoch_low": 0,
    "epoch_high": 0,
    # --- wavelength range; equal values mean "match the object" ----------
    "lower_lam": 0,
    "upper_lam": 0,
    # --- preprocessing ---------------------------------------------------
    "error_spectrum": "sg",
    "mask_galaxy_lines": 1,
    "mask_telluric": 1,
    "minimum_overlap": 0.7,
    # --- fitting ---------------------------------------------------------
    "iterations": 10,
    "n_cores": 0,
    "weighted_solve": 0,
    # --- output ----------------------------------------------------------
    "saving_results_path": "",
    "show_plot": 0,
    "show_plot_png": False,
    "how_many_plots": 0,
}


# Shorthand accepted by Superfit(...) and mapped onto real config keys.
# These exist because "z=0.127" is what someone means, and
# {"use_exact_z": 1, "z_exact": 0.127} is what the fitter needs.
ALIASES = {
    "z": "z_exact",
    "resolution_angstrom": "resolution",
    "output_dir": "saving_results_path",
    "n_plots": "how_many_plots",
    "error_model": "error_spectrum",
    "cores": "n_cores",
}


class ConfigError(ValueError):
    """Raised for a configuration that cannot be used as given."""


def parse_json_string(text):
    """Return the parsed object, or False when ``text`` is not JSON."""

    try:
        return json.loads(text)
    except (ValueError, TypeError):
        return False


def parse_json_file(path):
    """Return the parsed file contents, or False when unreadable."""

    try:
        with open(path, "r") as handle:
            return json.load(handle)
    except Exception:
        return False


def read_config(source):
    """Read a raw configuration from a dict, a JSON string, or a file path.

    No defaults are applied; see :func:`load_config`.
    """

    if source is None:
        return {}

    if isinstance(source, dict):
        return dict(source)

    parsed = parse_json_string(source)
    if isinstance(parsed, dict):
        return parsed

    parsed = parse_json_file(source)
    if isinstance(parsed, dict):
        return parsed

    raise ConfigError(
        "Could not read superfit parameters from {!r}. Expected a dict, a "
        "JSON string, or a path to a JSON file.".format(source)
    )


def resolve_aliases(overrides):
    """Expand shorthand keys, and turn ``z=`` into an exact-redshift request."""

    resolved = {}
    for key, value in overrides.items():
        if value is None:
            continue
        resolved[ALIASES.get(key, key)] = value

    # Naming a single redshift implies wanting exactly that redshift, unless
    # the caller also asked for a scan.
    if "z" in overrides and "use_exact_z" not in overrides:
        resolved.setdefault("use_exact_z", 1)

    # Conversely, asking for a scan range implies not wanting an exact z.
    scan_keys = {"z_range_begin", "z_range_end", "z_int"}
    if scan_keys & set(resolved) and "use_exact_z" not in overrides and "z" not in overrides:
        resolved.setdefault("use_exact_z", 0)

    return resolved


def load_config(source=None, **overrides):
    """Build a complete configuration.

    Parameters
    ----------
    source : dict, str, or None
        A parameter dict, a JSON string, or a path to a JSON file.
    **overrides
        Individual settings, taking precedence over ``source``. Shorthand
        names in :data:`ALIASES` are accepted.

    Returns
    -------
    dict
        Every key in :data:`DEFAULT_CONFIG`, with the caller's values merged
        over the defaults.
    """

    config = copy.deepcopy(DEFAULT_CONFIG)
    config.update(read_config(source))
    config.update(resolve_aliases(overrides))

    unknown = set(config) - set(DEFAULT_CONFIG)
    if unknown:
        raise ConfigError(
            "Unknown configuration key(s): {}. Valid keys are: {}".format(
                ", ".join(sorted(repr(k) for k in unknown)),
                ", ".join(sorted(DEFAULT_CONFIG)),
            )
        )

    validate_config(config)
    return config


def validate_config(config):
    """Check a merged configuration, raising ConfigError on anything unusable."""

    if config["resolution"] <= 0:
        raise ConfigError(
            "resolution must be positive, got {!r}".format(config["resolution"])
        )

    if not 0 <= config["minimum_overlap"] <= 1:
        raise ConfigError(
            "minimum_overlap is a fraction between 0 and 1, got {!r}".format(
                config["minimum_overlap"]
            )
        )

    if config["error_spectrum"] not in ("sg", "linear", "included"):
        raise ConfigError(
            "error_spectrum must be 'sg', 'linear' or 'included', got {!r}".format(
                config["error_spectrum"]
            )
        )

    if config["Alam_interval"] <= 0:
        raise ConfigError(
            "Alam_interval must be positive, got {!r}".format(config["Alam_interval"])
        )

    if config["Alam_high"] < config["Alam_low"]:
        raise ConfigError(
            "Alam_high ({}) is below Alam_low ({})".format(
                config["Alam_high"], config["Alam_low"]
            )
        )

    if not config["use_exact_z"]:
        if config["z_int"] <= 0:
            raise ConfigError(
                "z_int must be positive for a redshift scan, got {!r}".format(
                    config["z_int"]
                )
            )
        if config["z_range_end"] < config["z_range_begin"]:
            raise ConfigError(
                "z_range_end ({}) is below z_range_begin ({})".format(
                    config["z_range_end"], config["z_range_begin"]
                )
            )

    # Masking host lines shifts them to the object's redshift, which is only
    # defined if there is a single redshift to shift to.
    if config["mask_galaxy_lines"] and not config["use_exact_z"]:
        raise ConfigError(
            "mask_galaxy_lines needs a single redshift, but a scan was "
            "requested. Either set an exact z, or pass mask_galaxy_lines=0."
        )

    for name in ("temp_sn_tr", "temp_gal_tr"):
        if not config[name]:
            raise ConfigError("{} is empty; there would be nothing to fit".format(name))

    return config
