"""Derived quantities for one fit.

A :class:`Parameters` belongs to the fit that built it. It used to be
process-global state, installed by ``set_config()`` and read back through a
module-level proxy, which meant constructing a second ``Superfit`` silently
retuned the first one -- so ``a = Superfit(s1, z=0.1); b = Superfit(s2,
z=0.2); a.run()`` fitted s1 at z=0.2. Two fits in one process, and any kind
of concurrency, were quietly wrong. Nothing here is global any more.
"""

import functools
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
from superfit import paths as paths_module
from superfit.paths import gal_dir, sne_dir


# Configuration reading and defaulting live in superfit.config; these names
# are kept here because older code imports them from this module.
parseJsonString = config_module.parse_json_string
parseJsonFile = config_module.parse_json_file
load_config = config_module.load_config


def is_qso_name(path):
    """Whether a galaxy-directory template is a QSO, by its name.

    The modern banks ship their QSO templates in the galaxy directory as
    ``DESI_QSO_<targetid>``. A QSO is never a host: offering one under a
    transient would fit supernovae on top of quasars, which is exactly the
    combination the QSO category exists to replace. So anything the galaxy
    glob finds with a QSO name is pulled out of the host list unconditionally
    and fit on its own instead.
    """

    return "QSO" in os.path.basename(str(path)).upper()


def is_star_type(name):
    """Whether a supernova-type directory holds stars, by its name.

    The modern banks ship stellar templates as type directories named
    ``star-A`` .. ``star-WD`` alongside the supernova types. They are a
    category of their own: fit alone, with no host, at redshift zero.
    """

    return str(name).lower().startswith("star")


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
        # Resolved onto *this object* and nowhere else. An earlier version
        # called set_bank_dir() here, which made the choice process-wide, so
        # constructing a second Superfit moved the bank out from under the
        # first: `a = Superfit(bank="legacy"); b = Superfit(bank="modern");
        # a.run()` fitted a against modern's supernovae and legacy's galaxies
        # and labelled the result "legacy". Every bank-dependent path is now
        # derived from self.bank_dir, and the process-wide search is only for
        # a caller that named no bank.
        self.bank = data.get("bank") or ""
        if self.bank:
            self.bank_dir = paths_module.bank_dir_for_name(self.bank)
        else:
            self.bank_dir = paths_module.find_bank_dir()

        # The phase table belongs to the bank, so it is settled here too
        # rather than re-resolved from the global by whoever reads it.
        self.phase_table = paths_module.mjd_max_brightness_csv(self.bank_dir)

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

        # Globbed out of self.bank_dir, never out of the process-wide bank.
        binned = self.resolution if self.resolution in (10, 30) else None
        self._binned_resolution = binned

        templates_gal = glob.glob(
            os.path.join(gal_dir(binned, self.bank_dir), "*")
        )
        templates_gal = [
            x for x in templates_gal if "CVS" not in x and "README" not in x
        ]

        # QSOs live in the galaxy directory but are not hosts; see
        # is_qso_name. Split them out BEFORE the type selection, so no
        # temp_gal_tr setting can put a transient on top of a quasar.
        self.templates_qso = sorted(x for x in templates_gal if is_qso_name(x))
        templates_gal = np.array([x for x in templates_gal if not is_qso_name(x)])

        self.templates_gal_trunc = select_templates(templates_gal, self.temp_gal_tr)

        if len(self.templates_gal_trunc) == 0:
            raise config_module.ConfigError(
                "No galaxy template in {} matches temp_gal_tr={!r}. Every "
                "transient is fit as SN + host, so an empty host list fits "
                "nothing. (QSO templates do not count: they are never "
                "hosts.)".format(
                    gal_dir(binned, self.bank_dir), list(self.temp_gal_tr)
                )
            )

        # Star types the bank supplies: type directories named star-*.
        try:
            sne_root = sne_dir(bank_dir=self.bank_dir)
            self.star_types = sorted(
                entry
                for entry in os.listdir(sne_root)
                if is_star_type(entry)
                and os.path.isdir(os.path.join(sne_root, entry))
            )
        except OSError:
            self.star_types = []

        # "auto" means "when the bank supplies them", which is what keeps the
        # legacy bank exactly as it was: it has neither, so both resolve off.
        self.fit_stars = self._resolve_category(
            data["fit_stars"], "fit_stars", self.star_types,
            "star templates (sne/star-* type directories)",
        )
        self.fit_qsos = self._resolve_category(
            data["fit_qsos"], "fit_qsos", self.templates_qso,
            "QSO templates (gal/*QSO* files)",
        )

        self._frozen = True

    def _resolve_category(self, setting, key, supply, what):
        """Turn a true/false/"auto" category setting into a plain bool."""

        if setting == "auto":
            return bool(supply)
        if setting and not supply:
            raise config_module.ConfigError(
                "{}=true, but the bank at {} has no {}. The legacy bank does "
                "not carry them; fit against one of the modern banks, or "
                "leave {} on 'auto'.".format(key, self.bank_dir, what, key)
            )
        return bool(setting)

    @functools.cached_property
    def templates_sn_trunc(self):
        """Supernova template paths matching the selected types.

        The fit itself does not read this -- it resolves its supernova list
        from the bank metadata, which matches types exactly -- so the glob of
        the whole supernova tree it takes to build (a thousand-odd directory
        entries, several seconds on a cold parallel filesystem) is deferred
        until something actually asks.
        """

        templates_sn = glob.glob(
            os.path.join(
                sne_dir(self._binned_resolution, self.bank_dir), "**", "**", "*"
            )
        )
        templates_sn = [
            x
            for x in templates_sn
            if "wiserep_spectra.csv" not in x
            and "info" not in x
            and "photometry" not in x
            and "photometry.pdf" not in x
        ]

        return select_templates(np.array(templates_sn), self.temp_sn_tr)

    def __setattr__(self, name, value):
        if self._frozen:
            raise AttributeError(
                "Parameters is immutable once built; {!r} cannot be changed. "
                "Build a new Superfit with the settings you want.".format(name)
            )
        object.__setattr__(self, name, value)

    def __repr__(self):
        # Deliberately does not touch templates_sn_trunc: a repr should not
        # cost a walk of the template bank.
        return "<Parameters for {!r}: {} SN types and {} galaxy templates>".format(
            self.object_to_fit,
            len(self.temp_sn_tr),
            len(self.templates_gal_trunc),
        )

    @property
    def metadata_key(self):
        """What a bank metadata scan actually depends on.

        Used to cache the scan across fits: it walks every object directory
        and parses ~190 CSVs, and two fits that differ only in redshift want
        the same answer.

        This fit's own bank directory and phase table are part of the key as
        well as the settings: the scan reads both, so two fits against
        different banks must not share an answer. It used to read the
        process-wide bank here, which meant the key described whichever bank
        was resolved last rather than the one this fit is about to use.
        """

        return (
            self.bank_dir,
            self.phase_table,
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
