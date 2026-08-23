import os
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.ndimage import gaussian_filter1d
from astropy.table import Table


from superfit.SF_functions import Alam, all_parameter_space, remove_telluric, mask_gal_lines
from superfit.config import ConfigError, load_config
from superfit.Header_Binnings import kill_header, kill_header_and_bin, normalise_flux
from superfit.error_routines import linear_error, savitzky_golay
from superfit.get_metadata import get_metadata
from superfit.params import parameters, get_parameters, set_config
from superfit.paths import gal_dir, sne_dir
from superfit.spectrum import Spectrum



def _looks_like_config(value):
    """True when a first positional argument is a configuration, not a spectrum.

    A dict is unambiguous. A string is a config if it parses as JSON or ends
    in .json; anything else is treated as a path to a spectrum.
    """

    if isinstance(value, dict):
        return True
    if isinstance(value, str):
        if value.strip().startswith("{"):
            return True
        return value.lower().endswith(".json")
    return False


class Superfit:
    """One fit of one observed spectrum against the template bank.

    The observation can be supplied in whatever form it is already in::

        Superfit(wavelength=lam, flux=flux, error=err, z=0.127)
        Superfit(Spectrum.from_file("SN2021urb.flm"), z=0.127)
        Superfit("SN2021urb.flm", z=0.127)
        Superfit(numpy_array_of_shape_n_by_3, z=0.127)

    and the configuration as a dict, a JSON string, a path to a JSON file,
    or simply as keyword arguments::

        Superfit(spectrum, config="parameters.json")
        Superfit(spectrum, config={"resolution": 30}, z=0.127)
        Superfit("parameters.json")     # legacy: config only, spectrum named in it

    Anything not specified takes its value from
    :data:`superfit.config.DEFAULT_CONFIG`.
    """

    def __init__(
        self,
        spectrum=None,
        config=None,
        wavelength=None,
        flux=None,
        error=None,
        name=None,
        **overrides
    ):
        # `Superfit(config_dict)` and `Superfit("parameters.json")` predate
        # the spectrum argument and still mean the config.
        if config is None and _looks_like_config(spectrum):
            config, spectrum = spectrum, None

        spectrum = self._resolve_spectrum(
            spectrum, wavelength=wavelength, flux=flux, error=error, name=name
        )

        merged = load_config(config, **overrides)
        if spectrum is None:
            if not merged["object_to_fit"]:
                raise ConfigError(
                    "No spectrum to fit. Pass a Spectrum, a file path, or "
                    "wavelength= and flux= arrays, or set 'object_to_fit' in "
                    "the configuration."
                )
            spectrum = Spectrum.from_file(merged["object_to_fit"])
        elif not merged["object_to_fit"]:
            # Recorded in *_used.json so the run stays self-describing.
            merged["object_to_fit"] = spectrum.name

        set_config(merged, spectrum=spectrum)

        self.observation = spectrum.check_long_enough(parameters.resolution)

        self.name = spectrum.name
        self.name_no_extension = spectrum.name
        self.original_path_name = merged["object_to_fit"]

        self.lamda = spectrum.wavelength
        self.flux = spectrum.flux
        self.spectrum = spectrum.as_array()

        prefix = parameters.save_results_path
        self.binned_name = prefix + self.name_no_extension + "_binned.txt"
        self.results_name = prefix + self.name_no_extension
        self.results_path = self.results_name + ".csv"

        # What the chi2 is measured against: masked, normalised, resampled.
        masked = spectrum
        if parameters.mask_galaxy_lines:
            masked = masked.without_host_lines(parameters.redshift)
        if parameters.mask_telluric:
            masked = masked.without_telluric()

        self.masked = masked.normalised()
        self.int_obj = self.masked.interpolated_onto(parameters.lam)

        # The error spectrum is measured on the binned, *unmasked* observation.
        self.binned = spectrum.binned(parameters.resolution)

        self.metadata = get_metadata()

        self._write_used_config()

    @staticmethod
    def _resolve_spectrum(spectrum, wavelength, flux, error, name):
        """Turn whichever of the input forms was used into a Spectrum."""

        if wavelength is not None or flux is not None:
            if spectrum is not None:
                raise TypeError(
                    "Pass either a spectrum or wavelength=/flux= arrays, not both."
                )
            if wavelength is None or flux is None:
                raise TypeError(
                    "wavelength= and flux= must be given together."
                )
            return Spectrum(wavelength, flux, error=error, name=name)

        if spectrum is None:
            return None

        return Spectrum.coerce(spectrum, name=name)

    def _write_used_config(self):
        """Record the effective configuration next to the results."""

        used_json = parameters.save_results_path + "{}_used.json".format(
            self.name_no_extension
        )
        try:
            with open(used_json, "w") as handle:
                json.dump(parameters.config, handle, indent=2, default=str)
        except OSError as exc:
            # Not being able to write the audit file should not lose the fit.
            print("WARNING: could not write {}: {}".format(used_json, exc))

    def plot(self):

        plt.figure(figsize=(7 * np.sqrt(2), 7))
        plt.title(str(self.name), fontsize=17)
        plt.ylabel("Flux", fontsize=16)
        plt.xlabel("Lamda", fontsize=16)
        plt.plot(self.lamda, self.flux, "k")

    def mask_telluric(self):

        masked_spectrum = remove_telluric(self.spectrum)
        plt.figure(figsize=(7 * np.sqrt(2), 7))
        plt.title("Masked Telluric for " + str(self.name), fontsize=17)
        plt.ylabel("Flux", fontsize=16)
        plt.xlabel("Lamda", fontsize=16)
        plt.plot(masked_spectrum[:, 0], masked_spectrum[:, 1], "k")

    def mask_galaxy_lines(self):

        if parameters.use_exact_z != 1:
            raise Exception(
                "Make sure to pick an exact value for z in order to mask the host lines accordingly!"
            )
        else:
            Data = mask_gal_lines(self.spectrum, z_obj=parameters.redshift[0])
            plt.figure(figsize=(7 * np.sqrt(2), 7))
            plt.ylabel("Flux arbitrary", fontsize=14)
            plt.xlabel("Lamda", fontsize=14)
            plt.title(
                "Galaxy lines masked at z=" + str(parameters.redshift[0]),
                fontsize=15,
                fontweight="bold",
            )
            plt.plot(
                self.lamda, self.flux / np.median(self.flux), "r", label=str(self.name)
            )
            plt.plot(
                Data[:, 0],
                Data[:, 1] / np.median(Data[:, 1]),
                "b",
                label="Masked object",
            )
            plt.legend(framealpha=1, frameon=True, fontsize=12)
            # plt.savefig(str(self.name) + '_masked.pdf' )

    def sg_error(self):

        Data = self.spectrum
        error = savitzky_golay(self.spectrum)[:, 1]

        plt.figure(figsize=(7 * np.sqrt(2), 7))
        plt.ylabel("Flux arbitrary", fontsize=14)
        plt.xlabel("Lamda", fontsize=14)
        plt.title("Savitzky-Golay error estimation", fontsize=15, fontweight="bold")
        plt.fill_between(
            Data[:, 0],
            Data[:, 1] / np.median(Data[:, 1]) - error,
            Data[:, 1] / np.median(Data[:, 1]) + error,
            color="#FF4500",
            label="error",
        )
        plt.plot(
            Data[:, 0], Data[:, 1] / np.median(Data[:, 1]), "k", label=str(self.name)
        )
        plt.legend(framealpha=1, frameon=True, fontsize=12)
        # plt.savefig(str(self.name) + '_sg.pdf' )

    def linear_error(self):

        Data = self.spectrum
        error = linear_error(self.spectrum)[:, 1]

        plt.figure(figsize=(7 * np.sqrt(2), 7))
        plt.ylabel("Flux arbitrary", fontsize=14)
        plt.xlabel("Lamda", fontsize=14)
        plt.title("Linear error estimation", fontsize=15, fontweight="bold")
        plt.fill_between(
            Data[:, 0],
            (Data[:, 1] - error) / np.median(Data[:, 1]),
            (Data[:, 1] + error) / np.median(Data[:, 1]),
            color="#03AC13",
            label="error",
        )
        plt.plot(
            Data[:, 0], Data[:, 1] / np.median(Data[:, 1]), "k", label=str(self.name)
        )
        plt.legend(framealpha=1, frameon=True, fontsize=12)

    def mask_gal_lines_and_telluric(self):

        Data_masked = mask_gal_lines(self.name, z_obj=parameters.redshift)
        masked_spectrum = remove_telluric(Data_masked)

        plt.figure(figsize=(7 * np.sqrt(2), 7))
        plt.title("Masked Telluric and Galaxy lines", fontsize=17)
        plt.ylabel("Flux", fontsize=16)
        plt.xlabel("Lamda", fontsize=16)
        plt.plot(
            masked_spectrum[:, 0], masked_spectrum[:, 1], "k", label=str(self.name)
        )
        plt.legend(framealpha=1, frameon=True, fontsize=12)

    def run(self, save_binned=True):
        """Fit the spectrum and return the ranked results.

        Parameters
        ----------
        save_binned : bool
            Also write the binned observation next to the results. It is an
            intermediate, kept because it is useful to inspect; the fit no
            longer reads it back.

        Returns
        -------
        pandas.DataFrame
            One row per surviving template, best match first.
        """

        print(
            "Running optimization for spectrum: {0} with resolution = {1} Å".format(
                self.name_no_extension, parameters.resolution
            )
        )

        if save_binned:
            self._write_binned()

        all_parameter_space(
            self.int_obj,
            parameters.redshift,
            parameters.extconstant,
            parameters.templates_sn_trunc,
            parameters.templates_gal_trunc,
            parameters.lam,
            parameters.resolution,
            parameters.iterations,
            kind=parameters.kind,
            original=self.binned,
            spectrum_name=self.name,
            save=self.results_name,
            show=parameters.show,
            minimum_overlap=parameters.minimum_overlap,
            n_cores=parameters.n_cores,
            observed_grid=parameters.observed_grid,
            weighted_solve=parameters.weighted_solve,
        )

        self.results = pd.read_csv(self.results_path)
        self._plot_best_fits()
        return self.results

    def _write_binned(self):
        """Save the binned observation as two columns of text."""

        try:
            np.savetxt(
                self.binned_name,
                np.column_stack([self.binned.wavelength, self.binned.flux]),
                fmt="%s",
            )
        except OSError as exc:
            print("WARNING: could not write {}: {}".format(self.binned_name, exc))

    def superfit(self):
        """Backwards-compatible alias for :meth:`run`."""

        return self.run()

    def _plot_best_fits(self):

        result_number = 0

        if parameters.n > len(self.results):

            result_number = result_number + len(self.results)

        elif len(self.results) >= parameters.n:

            result_number = result_number + parameters.n

        for j in range(result_number):

            row = self.results.iloc[j]

            hg_name = row["GALAXY"]
            short_name = row["SN"]
            bb = row["CONST_SN"]
            dd = row["CONST_GAL"]
            z = row["Z"]
            extmag = row["A_v"]
            sn_cont = row["Frac(SN)"]

            # Get all names from the dictionary
            full_names = [str(x) for x in self.metadata.shorhand_dict.keys()]
            short_names = [str(x) for x in self.metadata.shorhand_dict.values()]

            # print(full_names)

            for i in range(0, len(short_names)):
                if str(short_names[i]) == str(short_name):
                    sn_best_fullname = full_names[i]
                    sn_short_name = short_names[i]
                    idx = sn_short_name.rfind("/")
                    subtype = sn_short_name[:idx]

            int_obj = self.int_obj

            sn_name = os.path.join(sne_dir(10), subtype, sn_best_fullname)
            hg_name = os.path.join(gal_dir(10), hg_name)

            # print(sn_name)

            nova = kill_header(sn_name)
            nova[:, 1] = nova[:, 1] / np.nanmedian(nova[:, 1])

            host = np.loadtxt(hg_name)
            host[:, 1] = host[:, 1] / np.nanmedian(host[:, 1])

            # Interpolate supernova and host galaxy
            # redshifted_nova   =  nova[:,0]*(z+1)
            # extinct_nova      =  nova[:,1]*10**(-0.4*extmag * Alam(nova[:,0]))/(1+z)

            # reshifted_host    =  host[:,0]*(z+1)
            # reshifted_hostf   =  host[:,1]/(z+1)

            redshifted_nova = nova[:, 0] * (z + 1)
            extinct_nova = (
                nova[:, 1] * 10 ** (-0.4 * extmag * Alam(nova[:, 0])) / (z + 1)
            )

            reshifted_host = host[:, 0] * (z + 1)
            reshifted_hostf = host[:, 1] / (z + 1)

            nova_int = interpolate.interp1d(
                redshifted_nova, extinct_nova, bounds_error=False, fill_value="nan"
            )
            host_int = interpolate.interp1d(
                reshifted_host, reshifted_hostf, bounds_error=False, fill_value="nan"
            )
            host_nova = bb * nova_int(parameters.lam) + dd * host_int(parameters.lam)

            sn_type = short_name[: short_name.find("/")]
            hg_name = hg_name[hg_name.rfind("/") + 1 :]
            subclass = short_name[short_name.find("/") + 1 : short_name.rfind("/")]
            phase = str(short_name[short_name.rfind(":") + 1 : -1])

            plt.figure(figsize=(8 * np.sqrt(2), 8))
            plt.plot(parameters.lam, int_obj, "r", label="Input object: " + self.name)
            plt.plot(
                parameters.lam,
                host_nova,
                "g",
                label="SN: "
                + sn_type
                + " - "
                + subclass
                + " - Phase: "
                + phase
                + "\nHost: "
                + str(hg_name)
                + "\nSN contrib: {0: .1f}%".format(100 * sn_cont),
            )
            plt.legend(framealpha=1, frameon=True, fontsize=12)
            plt.ylabel("Flux arbitrary", fontsize=14)
            plt.xlabel("Lamda", fontsize=14)
            plt.title("Best fit for z = " + str(z), fontsize=15, fontweight="bold")

            if parameters.show_plot_png:
                plt.savefig(self.results_name + "_" + str(j) + ".png")
            else:
                plt.savefig(self.results_name + "_" + str(j) + ".pdf")

            if parameters.show == 1:
                plt.show()

    # NOTE: there is no results() method. `run()` assigns self.results, which
    # would shadow any method of that name anyway; the old one was
    # unreachable and returned itself.

    def any_result(self, j):

        row = self.results.iloc[j]

        hg_name = row["GALAXY"]
        short_name = row["SN"]
        bb = row["CONST_SN"]
        dd = row["CONST_GAL"]
        z = row["Z"]
        extmag = row["A_v"]
        sn_cont = row["Frac(SN)"]

        # Get all names from the dictionary
        full_names = [str(x) for x in self.metadata.shorhand_dict.keys()]
        short_names = [str(x) for x in self.metadata.shorhand_dict.values()]

        for i in range(0, len(short_names)):
            if str(short_names[i]) == str(short_name):
                sn_best_fullname = full_names[i]
                sn_short_name = short_names[i]
                idx = sn_short_name.rfind("/")
                subtype = sn_short_name[:idx]

        int_obj = self.int_obj

        sn_name = os.path.join(sne_dir(10), subtype, sn_best_fullname)
        hg_name = os.path.join(gal_dir(10), hg_name)

        nova = kill_header(sn_name)
        nova[:, 1] = nova[:, 1] / np.nanmedian(nova[:, 1])

        host = np.loadtxt(hg_name)
        host[:, 1] = host[:, 1] / np.nanmedian(host[:, 1])

        # Interpolate supernova and host galaxy
        redshifted_nova = nova[:, 0] * (z + 1)
        extinct_nova = nova[:, 1] * 10 ** (-0.4 * extmag * Alam(nova[:, 0])) / (z + 1)

        reshifted_host = host[:, 0] * (z + 1)
        reshifted_hostf = host[:, 1] / (z + 1)

        nova_int = interpolate.interp1d(
            redshifted_nova, extinct_nova, bounds_error=False, fill_value="nan"
        )
        host_int = interpolate.interp1d(
            reshifted_host, reshifted_hostf, bounds_error=False, fill_value="nan"
        )
        host_nova = bb * nova_int(parameters.lam) + dd * host_int(parameters.lam)

        sn_type = short_name[: short_name.find("/")]
        hg_name = hg_name[hg_name.rfind("/") + 1 :]
        subclass = short_name[short_name.find("/") + 1 : short_name.rfind("/")]
        phase = str(short_name[short_name.rfind(":") + 1 : -1])
        plt.figure(figsize=(8 * np.sqrt(2), 8))
        plt.plot(parameters.lam, int_obj, "r", label="Input object: " + self.name)
        plt.plot(
            parameters.lam,
            host_nova,
            "g",
            label="SN: "
            + sn_type
            + " - "
            + subclass
            + " - Phase: "
            + phase
            + "\nHost: "
            + str(hg_name)
            + "\nSN contrib: {0: .1f}%".format(100 * sn_cont),
        )
        plt.legend(framealpha=1, frameon=True, fontsize=12)
        plt.ylabel("Flux arbitrary", fontsize=14)
        plt.xlabel("Lamda", fontsize=14)
        plt.title("Best fit for z = " + str(z), fontsize=15, fontweight="bold")

        if parameters.show_plot_png:
            plt.savefig(self.results_name + "_" + str(j) + ".png")
        else:
            plt.savefig(self.results_name + "_" + str(j) + ".pdf")

        if parameters.show == 1:
            plt.show()

    def convolution(self):

        obj_res = 100

        obj_med = np.median(self.lamda)
        width = obj_med / obj_res
        sig = width / (2 * np.sqrt(2 * np.log(2)))

        filtered = gaussian_filter1d(self.flux, sig)
        plt.plot(self.lamda, filtered)
