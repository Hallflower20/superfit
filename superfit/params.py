import glob
import numpy as np
import os
import sys
import json
from astropy.table import Table
from superfit import config as config_module
from superfit.auxiliary import select_templates
from superfit.paths import gal_dir, sne_dir


# Configuration reading and defaulting live in superfit.config; these names
# are kept here because older code imports them from this module.
parseJsonString = config_module.parse_json_string
parseJsonFile = config_module.parse_json_file
load_config = config_module.load_config


# The active configuration. Importing superfit must never require command-line
# arguments -- that made the package impossible to use from a notebook, from
# another script, or from a test. A configuration given on argv is still
# picked up automatically so `python run.py parameters.json` keeps working.
data = None
_cached_parameters = None
# The spectrum the active configuration refers to, when it was supplied as
# arrays rather than as a path for Parameters to go and read.
_active_spectrum = None


def set_config(source, spectrum=None):
    """Install ``source`` as the active configuration and return it."""

    global data, _cached_parameters, _active_spectrum
    data = load_config(source)
    _cached_parameters = None
    _active_spectrum = spectrum
    return data


def get_parameters():
    """Return the active :class:`Parameters`, building it at most once.

    ``Parameters`` globs the template bank and reads the object spectrum, so
    the previous habit of constructing it at import time in four different
    modules did that work four times per run.
    """

    global _cached_parameters
    if _cached_parameters is None:
        if data is None:
            raise RuntimeError(
                "No superfit parameters loaded. Call superfit.params.set_config(...) "
                "with a dict or a path to a JSON file, pass one on the command "
                "line, or use the Superfit(config=...) argument."
            )
        _cached_parameters = Parameters(data, spectrum=_active_spectrum)
    return _cached_parameters


class _LazyParameters:
    """Attribute proxy so module-level ``parameters.foo`` resolves on first use.

    Existing call sites read ``parameters.resolution`` and friends at call
    time; this keeps them working without requiring a configuration to exist
    at import time.
    """

    def __getattr__(self, name):
        return getattr(get_parameters(), name)

    def __repr__(self):
        if data is None or _cached_parameters is None:
            return "<superfit parameters: not loaded>"
        return "<superfit parameters for {}>".format(_cached_parameters.object_to_fit)


parameters = _LazyParameters()


if len(sys.argv) > 1:
    try:
        set_config(sys.argv[1])
    except ValueError:
        # argv[1] is not a superfit config -- e.g. pytest arguments. Stay unloaded.
        pass


class Parameters:
    def __init__(self, data, spectrum=None):
        """Derived quantities for one fit.

        Parameters
        ----------
        data : dict
            A configuration; missing keys are filled from DEFAULT_CONFIG.
        spectrum : Spectrum, optional
            The observation being fitted. Only needed when the wavelength
            grid is to be derived from it (``lower_lam == upper_lam``); it
            saves reading the file again, and it is the only way to fit a
            spectrum that never existed as a file.
        """

        data = config_module.load_config(data)

        # Keep the merged config so it can be echoed back into *_used.json.
        self.config = data

        self.object_to_fit = data["object_to_fit"]
        self.save_results_path = data["saving_results_path"]

        self.use_exact_z = data["use_exact_z"]
        self.z_exact = data["z_exact"]
        self.z_range_begin = data["z_range_begin"]
        self.z_range_end = data["z_range_end"]
        self.z_int = data["z_int"]

        if self.use_exact_z:
            self.redshift = np.array([self.z_exact])
        else:
            z_num = int((self.z_range_end - self.z_range_begin) / self.z_int) + 1
            self.redshift = np.linspace(self.z_range_begin, self.z_range_end, z_num)

        self.mask_galaxy_lines = data["mask_galaxy_lines"]
        self.mask_telluric = data["mask_telluric"]

        if self.mask_galaxy_lines == 1 and len(self.redshift) != 1:
            raise Exception(
                "Make sure to pick an exact value for z in order to mask the host lines accordingly!"
            )

        # Epochs
        self.epoch_high = data["epoch_high"]
        self.epoch_low = data["epoch_low"]

        # Chose minimum overlap
        self.minimum_overlap = data["minimum_overlap"]

        # Number of steps for A_v (do not change)
        self.Alam_high = data["Alam_high"]
        self.Alam_low = data["Alam_low"]
        self.Alam_interval = data["Alam_interval"]

        alam_num = int((self.Alam_high - self.Alam_low) / self.Alam_interval) + 1
        self.extconstant = np.linspace(self.Alam_low, self.Alam_high, alam_num)

        # Library to look at
        self.temp_gal_tr = data["temp_gal_tr"]
        self.temp_sn_tr = data["temp_sn_tr"]

        self.resolution = data["resolution"]
        self.upper = data["upper_lam"]
        self.lower = data["lower_lam"]

        if self.upper == self.lower:

            wavelength = self._observed_wavelength(spectrum)

            # Historical: the grid starts from the SECOND sample, not the
            # first. Preserved deliberately -- it shifts the start by one
            # pixel, which the 300 A padding swamps, but changing it would
            # move lam and so every chi2 in the output.
            self.lower = wavelength[1] - 300
            self.upper = wavelength[-1] + 300

            interval = int((self.upper - self.lower) / self.resolution)
            self.lam = np.linspace(self.lower, self.upper, interval)

        else:

            self.upper = data["upper_lam"]
            self.lower = data["lower_lam"]
            interval = int((self.upper - self.lower) / self.resolution)
            self.lam = np.linspace(self.lower, self.upper, interval)

        # Kind of error spectrum ('SG', 'linear' or 'included')
        self.kind = data["error_spectrum"]

        # Show plot?
        self.show = data["show_plot"]

        # Allow png output.
        if "show_plot_png" in data:
            self.show_plot_png = data["show_plot_png"]
        else:
            self.show_plot_png = False

        # How many results to plot?
        self.n = data["how_many_plots"]

        self.iterations = data.get("iterations", 10)

        # Worker processes for the (z, A_v) grid. 0 or absent means "use every
        # CPU this process is allowed on"; the pool is capped at the number of
        # grid points regardless.
        self.n_cores = data.get("n_cores", 0)

        # Solve for the SN/galaxy amplitudes using the same 1/sigma**2 weights
        # the chi2 uses. Off by default: it is the statistically consistent
        # choice, but it shifts every chi2 and can change the best-fit
        # template, so switching it on is a deliberate act.
        self.weighted_solve = data.get("weighted_solve", 0)

        # Template library

        if self.resolution == 10 or self.resolution == 30:
            templates_gal = glob.glob(os.path.join(gal_dir(self.resolution), "*"))
            templates_gal = [
                x for x in templates_gal if "CVS" not in x and "README" not in x
            ]
            templates_gal = np.array(templates_gal)

            templates_sn = glob.glob(
                os.path.join(sne_dir(self.resolution), "**", "**", "*")
            )

            templates_sn = [
                x
                for x in templates_sn
                if "wiserep_spectra.csv" not in x
                and "info" not in x
                and "photometry" not in x
                and "photometry.pdf" not in x
            ]
            templates_sn = np.array(templates_sn)

        else:
            templates_gal = glob.glob(os.path.join(gal_dir(), "*"))
            templates_gal = [
                x for x in templates_gal if "CVS" not in x and "README" not in x
            ]
            templates_gal = np.array(templates_gal)

            templates_sn = glob.glob(os.path.join(sne_dir(), "**", "**", "*"))
            templates_sn = [
                x
                for x in templates_sn
                if "wiserep_spectra.csv" not in x
                and "info" not in x
                and "photometry" not in x
                and "photometry.pdf" not in x
            ]
            templates_sn = np.array(templates_sn)

        self.templates_sn_trunc = select_templates(templates_sn, self.temp_sn_tr)
        self.templates_gal_trunc = select_templates(templates_gal, self.temp_gal_tr)

    def _observed_wavelength(self, spectrum):
        """The observed wavelength axis, from the Spectrum or from the file."""

        if spectrum is not None:
            return np.asarray(spectrum.wavelength, dtype=float)

        if not self.object_to_fit:
            raise config_module.ConfigError(
                "The wavelength grid is derived from the observation "
                "(lower_lam == upper_lam), but no spectrum was supplied and "
                "'object_to_fit' is unset. Pass a spectrum, or set "
                "lower_lam/upper_lam explicitly."
            )

        from superfit.spectrum import Spectrum

        return Spectrum.from_file(self.object_to_fit).wavelength
