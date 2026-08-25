"""Derived quantities for one fit.

A :class:`Parameters` belongs to the fit that built it. It used to be
process-global state, installed by ``set_config()`` and read back through a
module-level proxy, which meant constructing a second ``Superfit`` silently
retuned the first one -- so ``a = Superfit(s1, z=0.1); b = Superfit(s2,
z=0.2); a.run()`` fitted s1 at z=0.2. Two fits in one process, and any kind
of concurrency, were quietly wrong. Nothing here is global any more.
"""

import glob
import numpy as np
import os
from superfit import config as config_module
from superfit.auxiliary import select_templates
from superfit.loggrid import (
    LogGrid,
    matched_velocity_resolution,
    velocity_to_dlnlam,
)
from superfit.paths import gal_dir, sne_dir


# Configuration reading and defaulting live in superfit.config; these names
# are kept here because older code imports them from this module.
parseJsonString = config_module.parse_json_string
parseJsonFile = config_module.parse_json_file
load_config = config_module.load_config


class Parameters:
    """Everything one fit needs, derived once from one configuration.

    Instances are frozen after construction: the fit reads them from several
    modules and from worker processes, and a setting that can change halfway
    through is a setting that cannot be trusted in the results.
    """

    _frozen = False

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
        self.overwrite = data["overwrite"]

        # Resolved before anything below reads the bank -- the template lists
        # at the end of this method glob it, and get_metadata walks it.
        #
        # A named bank is pinned process-wide rather than carried on this
        # object, because the path helpers every reader goes through
        # (sne_dir, gal_dir) take no parameters. That matches the constraint
        # the fit already has: superfit serialises concurrent fits around
        # _FIT_LOCK, so two banks are not in play at once. Naming no bank
        # changes nothing and leaves the search path alone.
        self.bank = data.get("bank") or ""
        if self.bank:
            from superfit.paths import bank_dir_for_name, set_bank_dir

            self.bank_dir = set_bank_dir(bank_dir_for_name(self.bank))
        else:
            from superfit.paths import find_bank_dir

            self.bank_dir = find_bank_dir()

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

        if self.mask_galaxy_lines and len(self.redshift) != 1:
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

        # Shape of the CCM89 curve the A_v grid scales. Applied at the
        # template's rest wavelength, so A_v is host-galaxy dust.
        self.R_v = data["R_v"]

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
            # pixel, which the 300 A padding swamps.
            self.lower = wavelength[1] - 300
            self.upper = wavelength[-1] + 300

        else:

            self.upper = data["upper_lam"]
            self.lower = data["lower_lam"]

        # The fit runs on a grid uniform in ln(lambda), so that a redshift is
        # a shift along the axis rather than a rescaling of it. See
        # superfit.loggrid.
        if self.lower <= 0:
            raise config_module.ConfigError(
                "The fitting grid starts at {:.1f} A, which a logarithmic "
                "grid cannot represent. Set lower_lam explicitly to a "
                "positive wavelength.".format(self.lower)
            )

        self.velocity_resolution = data["velocity_resolution"]
        if self.velocity_resolution is None:
            self.velocity_resolution = matched_velocity_resolution(
                self.resolution, self.lower, self.upper
            )

        self.observed_grid = LogGrid.spanning(
            self.lower, self.upper, velocity_to_dlnlam(self.velocity_resolution)
        )
        self.lam = self.observed_grid.wavelength

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

        self._frozen = True

    def __setattr__(self, name, value):
        if self._frozen:
            raise AttributeError(
                "Parameters is immutable once built; {!r} cannot be changed. "
                "Build a new Superfit with the settings you want.".format(name)
            )
        object.__setattr__(self, name, value)

    def __repr__(self):
        return "<Parameters for {!r}: {} SN and {} galaxy templates>".format(
            self.object_to_fit,
            len(self.templates_sn_trunc),
            len(self.templates_gal_trunc),
        )

    @property
    def metadata_key(self):
        """What a bank metadata scan actually depends on.

        Used to cache the scan across fits: it walks every object directory
        and parses ~190 CSVs, and two fits that differ only in redshift want
        the same answer.

        The bank directory is part of the key as well as the settings --
        the scan reads the bank, so pointing at a different one has to
        invalidate it.
        """

        from superfit.paths import find_bank_dir

        return (
            find_bank_dir(),
            tuple(self.temp_sn_tr),
            self.epoch_low,
            self.epoch_high,
        )

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
