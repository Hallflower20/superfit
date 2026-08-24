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
from superfit.output import FitResult, RunDirectory
from superfit.params import Parameters
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

    Each instance owns its settings, as ``self.parameters``. Building a
    second Superfit does not disturb the first, so several fits can be set up
    and run in any order, or concurrently.
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

        # Settings first, spectrum second. Merging and validating is pure and
        # instant; reading a spectrum is I/O. A misspelled setting should be
        # reported as a misspelled setting, not hidden behind whatever the
        # filesystem happens to say about a path that was also wrong.
        merged = load_config(config, **overrides)

        spectrum = self._resolve_spectrum(
            spectrum, wavelength=wavelength, flux=flux, error=error, name=name
        )

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

        # This fit's own settings. Nothing about them is global, so a second
        # Superfit built after this one cannot change what this one does.
        self.parameters = parameters = Parameters(merged, spectrum=spectrum)

        self.observation = spectrum.check_long_enough(parameters.resolution)

        self.name = spectrum.name
        self.name_no_extension = spectrum.name
        self.original_path_name = merged["object_to_fit"]

        self.lamda = spectrum.wavelength
        self.flux = spectrum.flux
        self.spectrum = spectrum.as_array()

        # Created now rather than at save time: a directory that cannot be
        # written, or that already holds a fit of this spectrum, is worth
        # hearing about before the fit runs, not thirty seconds after.
        self.output = RunDirectory(
            parameters.save_results_path, self.name, overwrite=parameters.overwrite
        )
        self.results_path = self.output.results_csv
        self._artifacts = []

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

        self.metadata = get_metadata(parameters)

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

        used_json = self.output.used_config_json
        try:
            with open(used_json, "w") as handle:
                json.dump(self.parameters.config, handle, indent=2, default=str)
        except OSError as exc:
            # Not being able to write the audit file should not lose the fit.
            print("WARNING: could not write {}: {}".format(used_json, exc))
        else:
            self._artifacts.append(used_json)

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

        parameters = self.parameters
        if not parameters.use_exact_z:
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

        Data_masked = mask_gal_lines(self.name, z_obj=self.parameters.redshift)
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
        FitResult
            The ranked table -- one row per surviving template, best match
            first -- together with the directory the run was written to and
            every file in it. It indexes and iterates like the DataFrame it
            wraps, which is also available as ``.results``.
        """

        parameters = self.parameters

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
            results_path=self.results_path,
            show=parameters.show,
            minimum_overlap=parameters.minimum_overlap,
            n_cores=parameters.n_cores,
            observed_grid=parameters.observed_grid,
            weighted_solve=parameters.weighted_solve,
            R_v=parameters.R_v,
            parameters=parameters,
        )

        self.results = pd.read_csv(self.results_path)
        self._artifacts.append(self.results_path)

        self._plot_best_fits()

        self.result = FitResult(self.results, self.output, self._artifacts)
        return self.result

    def _write_binned(self):
        """Save the binned observation as two columns of text."""

        binned_txt = self.output.binned_txt
        try:
            np.savetxt(
                binned_txt,
                np.column_stack([self.binned.wavelength, self.binned.flux]),
                fmt="%s",
            )
        except OSError as exc:
            print("WARNING: could not write {}: {}".format(binned_txt, exc))
        else:
            self._artifacts.append(binned_txt)

    def superfit(self):
        """Backwards-compatible alias for :meth:`run`."""

        return self.run()

    def _plot_best_fits(self):
        """Plot the top ``how_many_plots`` matches into the run directory."""

        for j in range(min(self.parameters.n, len(self.results))):
            self._artifacts.append(self.plot_rank(j))

    # NOTE: there is no results() method. `run()` assigns self.results, which
    # would shadow any method of that name anyway; the old one was
    # unreachable and returned itself.

    def plot_rank(self, j):
        """Plot the ``j``-th ranked match against the observation.

        ``j`` counts from 0, as the results table does; the file is named
        for the rank, counting from 1. Returns the path it wrote.
        """

        parameters = self.parameters
        row = self.results.iloc[j]

        short_name = row["SN"]
        bb = row["CONST_SN"]
        dd = row["CONST_GAL"]
        z = row["Z"]
        extmag = row["A_v"]
        sn_cont = row["Frac(SN)"]

        # The results table carries the shorthand name; the file it came from
        # is the key of the same entry. Inverting the dict once beats the
        # linear scan this used to do per plot, and gives a real error rather
        # than a NameError when a name is missing.
        by_shorthand = {
            str(v): str(k) for k, v in self.metadata.shorhand_dict.items()
        }
        if str(short_name) not in by_shorthand:
            raise KeyError(
                "No template file for {!r} in the bank metadata; the results "
                "were produced against a different bank.".format(short_name)
            )
        sn_best_fullname = by_shorthand[str(short_name)]
        subtype = str(short_name)[: str(short_name).rfind("/")]

        sn_path = os.path.join(sne_dir(10), subtype, sn_best_fullname)
        hg_path = os.path.join(gal_dir(10), row["GALAXY"])

        nova = kill_header(sn_path)
        nova[:, 1] = nova[:, 1] / np.nanmedian(nova[:, 1])

        host = np.loadtxt(hg_path)
        host[:, 1] = host[:, 1] / np.nanmedian(host[:, 1])

        # Reconstruct the model the fitter scored, so the reddening has to
        # match it: the law is evaluated at nova[:, 0], the template's REST
        # wavelength, exactly as the fit evaluates it at lam / (1 + z).
        redshifted_nova = nova[:, 0] * (z + 1)
        extinct_nova = (
            nova[:, 1]
            * 10 ** (-0.4 * extmag * Alam(nova[:, 0], R_v=parameters.R_v))
            / (z + 1)
        )

        nova_int = interpolate.interp1d(
            redshifted_nova, extinct_nova, bounds_error=False, fill_value="nan"
        )
        host_int = interpolate.interp1d(
            host[:, 0] * (z + 1), host[:, 1] / (z + 1),
            bounds_error=False, fill_value="nan",
        )
        host_nova = bb * nova_int(parameters.lam) + dd * host_int(parameters.lam)

        sn_type = short_name[: short_name.find("/")]
        subclass = short_name[short_name.find("/") + 1 : short_name.rfind("/")]
        phase = str(short_name[short_name.rfind(":") + 1 : -1])

        path = self.output.plot(j + 1, png=parameters.show_plot_png)

        # Held explicitly rather than left on pyplot's global stack: a batch
        # of a few hundred fits used to accumulate every figure it had ever
        # drawn, because nothing ever closed them.
        figure = plt.figure(figsize=(8 * np.sqrt(2), 8))
        try:
            plt.plot(
                parameters.lam, self.int_obj, "r",
                label="Input object: " + self.name,
            )
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
                + str(os.path.basename(hg_path))
                + "\nSN contrib: {0: .1f}%".format(100 * sn_cont),
            )
            plt.legend(framealpha=1, frameon=True, fontsize=12)
            plt.ylabel("Flux arbitrary", fontsize=14)
            plt.xlabel("Lamda", fontsize=14)
            plt.title("Best fit for z = " + str(z), fontsize=15, fontweight="bold")

            plt.savefig(path)

            if parameters.show:
                plt.show()
        finally:
            plt.close(figure)

        return path

    def any_result(self, j):
        """Backwards-compatible alias for :meth:`plot_rank`."""

        return self.plot_rank(j)

    def convolution(self):

        obj_res = 100

        obj_med = np.median(self.lamda)
        width = obj_med / obj_res
        sig = width / (2 * np.sqrt(2 * np.log(2)))

        filtered = gaussian_filter1d(self.flux, sig)
        plt.plot(self.lamda, filtered)
