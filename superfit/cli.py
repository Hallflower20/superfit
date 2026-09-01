#!/usr/bin/env python
"""Command line for superfit.

Installed as the ``superfit`` command; also reachable as ``python run.py``
from a source checkout.

The common case is one spectrum and one redshift, and it should take one
line to say so::

    superfit fit spectrum.flm --z 0.127

Everything else has a default. A JSON parameter file is still accepted, for
runs that need to be reproducible or that set more than a handful of things::

    superfit fit spectrum.flm --config parameters.json

The setup commands exist because getting to that first line used to involve
a curl, an unzip and an environment variable::

    superfit bank install
    superfit bank status
    superfit doctor
"""

import argparse
import json
import os
import sys

SUBCOMMANDS = ("fit", "bank", "doctor", "config")

EPILOG = """\
examples:
  superfit fit spectrum.flm --z 0.127
  superfit fit spectrum.flm --scan-z 0.0 0.3 --z-step 0.005
  superfit fit spectrum.flm --z 0.127 --plots 3 --output results/
  superfit fit spectrum.flm --config parameters.json --overwrite

  superfit bank install          download and unpack the template bank
  superfit bank pack             pack it so a fit opens two files, not a thousand
  superfit bank status           say which bank a fit would use
  superfit doctor                check this installation end to end
  superfit config create         write a documented parameter file
"""


# -- argument parsing ------------------------------------------------------


def build_parser(prog="superfit"):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Classify a supernova spectrum against the superfit template bank.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    from superfit import __version__

    parser.add_argument("--version", action="version", version="superfit " + __version__)

    subparsers = parser.add_subparsers(dest="command", metavar="<command>")

    _add_fit_parser(subparsers)
    _add_bank_parser(subparsers)
    _add_doctor_parser(subparsers)
    _add_config_parser(subparsers)

    return parser


def _add_fit_parser(subparsers):
    from superfit.config import PROFILES

    fit = subparsers.add_parser(
        "fit",
        help="fit one spectrum against the template bank",
        description="Fit one spectrum against the template bank.",
    )
    fit.set_defaults(handler=run_fit)

    fit.add_argument(
        "spectrum",
        nargs="?",
        help="Spectrum to fit: two or three columns of wavelength, flux and "
        "optionally error. Omit only when --config names one.",
    )
    fit.add_argument(
        "--config",
        help="JSON parameter file (or a JSON string) to take settings from. "
        "Command-line options override it.",
    )
    fit.add_argument("--name", help="Name for the spectrum; defaults to the filename.")

    reading = fit.add_argument_group("reading the spectrum")
    reading.add_argument(
        "--columns",
        metavar="W,F[,E]",
        help="Which columns hold wavelength, flux and error. Names for csv "
        "and FITS ('wave,flux,ivar'), positions for ascii ('0,1,3'). "
        "Left out, they are identified from the file.",
    )
    reading.add_argument(
        "--wavelength-unit",
        metavar="UNIT",
        help="Unit the wavelengths are in: AA, nm, micron, log10(AA). Left "
        "out, taken from the file, and Angstroms assumed.",
    )
    reading.add_argument(
        "--hdu", help="Which FITS HDU to read, by index or by name."
    )

    redshift = fit.add_argument_group("redshift")
    redshift.add_argument(
        "--z", type=float, help="Fit at this exact redshift. The usual case."
    )
    redshift.add_argument(
        "--scan-z",
        nargs=2,
        type=float,
        metavar=("BEGIN", "END"),
        help="Search this redshift range instead of fixing z.",
    )
    redshift.add_argument(
        "--z-step",
        type=float,
        metavar="STEP",
        help="Redshift step for --scan-z (default 0.01).",
    )

    fitting = fit.add_argument_group("the fit")
    fitting.add_argument(
        "--bank",
        metavar="NAME",
        help="Which template bank to fit against: legacy, modern-curated, "
        "modern, or any name registered with `superfit bank install --from`. "
        "See `superfit bank list`.",
    )
    fitting.add_argument(
        "--resolution", type=int, metavar="A", help="Binning resolution in Angstroms."
    )
    fitting.add_argument(
        "--error-model",
        choices=("sg", "linear", "included"),
        help="How the error spectrum is obtained: a Savitzky-Golay estimate "
        "(sg, the default), a linear estimate, or the observation's own "
        "error column (included).",
    )
    fitting.add_argument(
        "--profile",
        choices=sorted(PROFILES),
        help="Named set of scientific defaults. 'legacy' (the default) "
        "reproduces published superfit results; 'modern' weights the "
        "amplitude solve consistently with the chi2.",
    )
    fitting.add_argument(
        "--rv", type=float, metavar="R_V", help="Total-to-selective extinction ratio."
    )
    fitting.add_argument(
        "--av-range",
        nargs=2,
        type=float,
        metavar=("LOW", "HIGH"),
        help="Extinction range to search, in magnitudes at V.",
    )
    fitting.add_argument(
        "--epochs",
        nargs=2,
        type=float,
        metavar=("LOW", "HIGH"),
        help="Restrict templates to this phase window, in days. Equal values "
        "mean no restriction.",
    )
    fitting.add_argument(
        "--no-stars",
        action="store_true",
        help="Do not fit the bank's stellar templates. By default a bank "
        "that carries them (the modern banks) has each star fit on its own, "
        "at z = 0, with no host galaxy.",
    )
    fitting.add_argument(
        "--no-qsos",
        action="store_true",
        help="Do not fit the bank's QSO templates. By default a bank that "
        "carries them has each QSO fit on its own over the redshift grid; "
        "QSOs are never offered as hosts for transients either way.",
    )
    fitting.add_argument(
        "--no-mask-galaxy-lines",
        action="store_true",
        help="Do not mask host galaxy emission lines.",
    )
    fitting.add_argument(
        "--no-mask-telluric",
        action="store_true",
        help="Do not mask the telluric A band.",
    )
    fitting.add_argument(
        "--n-cores",
        type=int,
        dest="n_cores",
        metavar="N",
        help="Worker processes for the (z, A_v) grid. 0 means auto.",
    )
    fitting.add_argument(
        "--weighted-solve",
        action="store_true",
        help="Shorthand for --profile modern.",
    )

    output = fit.add_argument_group("output")
    output.add_argument(
        "--output",
        "-o",
        metavar="DIR",
        help="Directory to create the run directory in (default: here).",
    )
    output.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing run of this spectrum.",
    )
    output.add_argument(
        "--plots",
        type=int,
        metavar="N",
        help="Plot the N best fits (default 0).",
    )
    output.add_argument(
        "--show", action="store_true", help="Display the plots as well as saving them."
    )
    output.add_argument("--png", action="store_true", help="Save plots as PNG, not PDF.")
    output.add_argument(
        "--quiet", action="store_true", help="Only print the best match."
    )

    # Accepted, not advertised: the names the previous CLI used.
    fit.add_argument("--object", dest="spectrum_legacy", help=argparse.SUPPRESS)
    fit.add_argument("--out", dest="output_legacy", help=argparse.SUPPRESS)
    fit.add_argument("--no-plots", action="store_true", help=argparse.SUPPRESS)

    return fit


def _add_bank_parser(subparsers):
    bank = subparsers.add_parser(
        "bank",
        help="install and inspect the template bank",
        description="Install and inspect the template bank.",
    )
    actions = bank.add_subparsers(dest="bank_command", metavar="<action>")

    install = actions.add_parser(
        "install",
        help="install a template bank, by name",
        description="Download, verify and unpack a template bank into a "
        "standard per-user data directory, which superfit searches "
        "automatically. No environment variable needed. `superfit bank list` "
        "shows the names.",
    )
    install.set_defaults(handler=run_bank_install)
    install.add_argument(
        "name",
        nargs="?",
        help="Which bank to install (default: legacy). See `superfit bank list`.",
    )
    install.add_argument(
        "--from",
        dest="from_dir",
        metavar="DIR",
        help="Install from a bank directory already on disk instead of "
        "downloading. Records where it is rather than copying it, which is "
        "what you want for a bank on shared storage.",
    )
    install.add_argument(
        "--copy",
        action="store_true",
        help="With --from, copy the bank into the per-user data directory "
        "rather than pointing at it where it is.",
    )
    install.add_argument("--dir", help="Install here instead of the default location.")
    install.add_argument(
        "--archive",
        help="Use a zip already on disk instead of downloading. For machines "
        "with no direct internet access.",
    )
    install.add_argument("--url", help="Download from here instead of the default.")
    install.add_argument(
        "--overwrite", action="store_true", help="Replace an existing installation."
    )
    install.add_argument(
        "--no-verify",
        action="store_true",
        help="Skip the checksum check. Only for installing a bank this "
        "version of superfit does not know the checksum of.",
    )
    install.add_argument("--quiet", action="store_true", help="No progress bars.")
    install.add_argument(
        "--no-pack",
        action="store_true",
        help="Skip building the packed bank. Fits then read the template text "
        "files, which is several seconds slower per run.",
    )

    pack = actions.add_parser(
        "pack",
        help="pack the template bank into one array per directory",
        description="Concatenate each template directory into a single .npy "
        "so a fit opens two files instead of a thousand. `bank install` does "
        "this already; run it by hand for a bank installed some other way, or "
        "after editing one. Packing is optional -- without it fits read the "
        "text files, just more slowly.",
    )
    pack.set_defaults(handler=run_bank_pack)
    pack.add_argument("--dir", help="Pack this bank instead of the one a fit would use.")
    pack.add_argument("--bank", help="Pack the bank with this name.")
    pack.add_argument(
        "--jobs",
        type=int,
        metavar="N",
        help="Processes to parse templates with (default 8). The largest "
        "bank is ten minutes of parsing on one core.",
    )
    pack.add_argument(
        "--quiet", action="store_true", help="Do not list what was packed."
    )

    verify = actions.add_parser(
        "verify",
        help="prove a pack still describes the bank's text",
        description="Re-read every template and compare against the digest "
        "recorded when the bank was packed. A fit checks names, sizes and "
        "modification times, which is cheap; this checks the bytes, which "
        "costs a full read of the bank. Use it when a result has to be "
        "provably reproducible.",
    )
    verify.set_defaults(handler=run_bank_verify)
    verify.add_argument("--dir", help="Verify this bank instead of the default.")
    verify.add_argument("--bank", help="Verify the bank with this name.")
    verify.add_argument("--jobs", type=int, metavar="N", help="Processes to read with.")
    verify.add_argument("--quiet", action="store_true", help="No progress bars.")

    listing = actions.add_parser(
        "list",
        help="list the template banks superfit knows about",
        description="List the banks superfit knows by name, which of them "
        "are installed on this machine, and where.",
    )
    listing.set_defaults(handler=run_bank_list)

    status = actions.add_parser(
        "status",
        help="say which template bank a fit would use",
        description="Report which template bank a fit would use, where it "
        "came from, and whether it is complete.",
    )
    status.set_defaults(handler=run_bank_status)
    status.add_argument("--dir", help="Inspect this directory instead of searching.")
    status.add_argument("--bank", help="Inspect the bank with this name.")

    bank.set_defaults(handler=lambda args: _needs_action(bank))
    return bank


def _add_doctor_parser(subparsers):
    doctor = subparsers.add_parser(
        "doctor",
        help="check this installation end to end",
        description="Check the Python version, the dependencies, the template "
        "bank and the output locations, and say what is wrong with any of "
        "them.",
    )
    doctor.set_defaults(handler=run_doctor)
    return doctor


def _add_config_parser(subparsers):
    from superfit.config import PROFILES

    config = subparsers.add_parser(
        "config",
        help="create and inspect parameter files",
        description="Create and inspect parameter files.",
    )
    actions = config.add_subparsers(dest="config_command", metavar="<action>")

    create = actions.add_parser(
        "create",
        help="write a documented parameter file",
        description="Write a parameter file with the settings most runs "
        "touch, each one explained.",
    )
    create.set_defaults(handler=run_config_create)
    create.add_argument(
        "path", nargs="?", default="parameters.json", help="Where to write it."
    )
    create.add_argument(
        "--full",
        action="store_true",
        help="Include every setting, not just the common ones.",
    )
    create.add_argument("--profile", choices=sorted(PROFILES))
    create.add_argument(
        "--force", action="store_true", help="Overwrite an existing file."
    )

    show = actions.add_parser(
        "show",
        help="print the effective configuration",
        description="Print the configuration a fit would run with, after "
        "defaults, profile and overrides are merged.",
    )
    show.set_defaults(handler=run_config_show)
    show.add_argument("source", nargs="?", help="A parameter file to merge in.")
    show.add_argument("--profile", choices=sorted(PROFILES))

    config.set_defaults(handler=lambda args: _needs_action(config))
    return config


def _needs_action(parser):
    parser.print_help()
    return 2


# -- fit -------------------------------------------------------------------


def _int_or_name(text):
    """A FITS HDU is addressed by index or by name; the shell gives us text."""

    try:
        return int(text)
    except ValueError:
        return text


def parse_columns(text):
    """``--columns wave,flux,ivar`` or ``--columns 0,1,3``.

    Positions stay integers so the ascii reader takes them as positions;
    anything else is a column name.
    """

    names = [part.strip() for part in text.split(",") if part.strip()]
    if not 2 <= len(names) <= 3:
        raise ValueError(
            "--columns takes two or three comma-separated names or positions "
            "(wavelength, flux, and optionally error), got {!r}.".format(text)
        )

    resolved = [_int_or_name(name) for name in names]

    # `ivar` is inverse variance, not an uncertainty, and has to be converted.
    if len(resolved) == 3 and str(resolved[2]).lower() in ("ivar", "invvar"):
        return {"wavelength": resolved[0], "flux": resolved[1], "ivar": resolved[2]}
    return resolved


def fit_overrides(args):
    """Turn parsed ``fit`` arguments into configuration overrides.

    Separated from :func:`run_fit` so it can be tested without a template
    bank: this is where a misread flag turns into a different fit.
    """

    notes = []
    overrides = {}

    if args.z is not None and args.scan_z:
        raise ValueError("Pass either --z or --scan-z, not both.")
    if args.z_step is not None and not args.scan_z:
        raise ValueError("--z-step only means something with --scan-z.")

    if args.z is not None:
        overrides["z"] = args.z
    elif args.scan_z:
        begin, end = args.scan_z
        overrides["use_exact_z"] = False
        overrides["z_range_begin"] = begin
        overrides["z_range_end"] = end
        if args.z_step is not None:
            overrides["z_int"] = args.z_step
        if not args.no_mask_galaxy_lines:
            # Host lines are placed at the object's redshift, so there has to
            # be one. Refusing the whole command over this would just teach
            # people to type --no-mask-galaxy-lines without reading it.
            overrides["mask_galaxy_lines"] = False
            notes.append(
                "--scan-z turns off host-line masking: the lines have no "
                "single redshift to sit at."
            )

    if args.columns:
        overrides["spectrum_columns"] = parse_columns(args.columns)
    if args.wavelength_unit:
        overrides["wavelength_unit"] = args.wavelength_unit
    if args.hdu is not None:
        overrides["spectrum_hdu"] = _int_or_name(args.hdu)

    if args.bank:
        overrides["bank"] = args.bank
    if args.resolution is not None:
        overrides["resolution"] = args.resolution
    if args.error_model is not None:
        overrides["error_model"] = args.error_model
    if args.rv is not None:
        overrides["R_v"] = args.rv
    if args.av_range:
        overrides["Alam_low"], overrides["Alam_high"] = args.av_range
    if args.epochs:
        overrides["epoch_low"], overrides["epoch_high"] = args.epochs
    if args.n_cores is not None:
        overrides["n_cores"] = args.n_cores

    if args.no_stars:
        overrides["fit_stars"] = False
    if args.no_qsos:
        overrides["fit_qsos"] = False

    if args.no_mask_galaxy_lines:
        overrides["mask_galaxy_lines"] = False
    if args.no_mask_telluric:
        overrides["mask_telluric"] = False

    output = args.output or args.output_legacy
    if output is not None:
        overrides["output_dir"] = output
    if args.overwrite:
        overrides["overwrite"] = True

    if args.plots is not None:
        overrides["n_plots"] = args.plots
    if args.no_plots:
        overrides["n_plots"] = 0
    if args.show:
        overrides["show_plot"] = True
        # Asking to see the plots and being shown none is not what anyone
        # means; plotting is off by default.
        overrides.setdefault("n_plots", 1)
    if args.png:
        overrides["show_plot_png"] = True

    return overrides, notes


def run_fit(args):
    from superfit.sf_class import Superfit

    spectrum = args.spectrum or args.spectrum_legacy
    if spectrum is None and not args.config:
        raise ValueError(
            "No spectrum to fit. Give one as `superfit fit spectrum.flm`, or "
            "name it in a --config file."
        )

    overrides, notes = fit_overrides(args)

    profile = args.profile
    if args.weighted_solve:
        if profile and profile != "modern":
            raise ValueError(
                "--weighted-solve is --profile modern; it contradicts "
                "--profile {}.".format(profile)
            )
        profile = "modern"

    for note in notes:
        print("note: {}".format(note))

    result = Superfit(
        spectrum, config=args.config, name=args.name, profile=profile, **overrides
    ).run()

    print_summary(result, quiet=args.quiet)
    return 0


def print_summary(result, quiet=False, top=3):
    """Say what won, and where the run was written."""

    if len(result) == 0:
        print("\nNo template survived the overlap cut. Nothing to report.")
        return

    best = result.best
    print("\nBest match: {}".format(best["SN"]))
    print(
        "  z = {:g}   A_v = {:g}   chi2/dof = {:.4g}   SN contributes "
        "{:.0f}%".format(
            best["Z"], best["A_v"], best["CHI2/dof"], 100 * best["Frac(SN)"]
        )
    )

    if not quiet and len(result) > 1:
        print("\nNext best:")
        for rank in range(1, min(top, len(result))):
            row = result.results.iloc[rank]
            print(
                "  {}. {:<44} chi2/dof = {:.4g}".format(
                    rank + 1, str(row["SN"])[:44], row["CHI2/dof"]
                )
            )

    print("\nWritten to {}".format(result.directory))
    if not quiet:
        for path in result.artifacts:
            print("  {}".format(path.name))


# -- bank ------------------------------------------------------------------


def run_bank_list(args):
    from superfit import bank, packed

    installed = bank.installed_banks()

    print("Template banks superfit knows by name:\n")
    for name, entry in bank.KNOWN_BANKS.items():
        directory = installed.get(name)
        mark = "installed" if directory else "not installed"
        print("  {:<16} {}  [{}]".format(name, entry["summary"], mark))
        print("      {}".format(entry["description"]))
        if directory:
            pack = packed.open_pack(directory, "binnings/10A/sne")
            print("      at {}{}".format(
                directory, "" if pack else "   (not packed -- `superfit bank pack`)"))
        elif entry["url"] is None:
            print("      Not published; install with --from <directory>.")
        print()

    extra = {n: d for n, d in installed.items() if n not in bank.KNOWN_BANKS}
    if extra:
        print("Also registered on this machine:\n")
        for name, directory in sorted(extra.items()):
            print("  {:<16} {}".format(name, directory))
        print()

    print("Fit against one with:  superfit fit spectrum.flm --z 0.1 --bank <name>")
    return 0


def run_bank_install(args):
    from superfit import bank, packed, paths

    name = args.name or bank.DEFAULT_BANK

    if args.from_dir:
        directory = bank.install_from_directory(
            name, args.from_dir, copy=args.copy, quiet=args.quiet
        )
        if not args.no_pack:
            _pack_installed(str(directory), args)
        print(
            "\nBank {!r} ready. Fit against it with:\n\n"
            "    superfit fit spectrum.flm --z 0.1 --bank {}\n".format(name, name)
        )
        return 0

    entry = bank.KNOWN_BANKS.get(name)
    if entry is None:
        print("No bank named {!r}. Try `superfit bank list`.".format(name))
        return 1
    if entry["url"] is None and not args.url:
        print(
            "The bank {!r} is not published, so there is nothing to download.\n"
            "Install it from a directory you already have:\n\n"
            "    superfit bank install {} --from <directory>\n".format(name, name)
        )
        return 1

    destination = args.dir or bank.bank_install_dir(name)

    manifest = bank.install(
        destination=destination,
        url=args.url or entry["url"],
        archive=args.archive,
        overwrite=args.overwrite,
        expected_sha256=None if args.no_verify else entry["sha256"],
        quiet=args.quiet,
    )

    bank.register(name, destination)

    print("\nTemplate bank {!r} installed in {}".format(name, destination))
    print("  sha256 {}".format(manifest["sha256"]))

    if not args.no_pack:
        _pack_installed(str(destination), args)

    print(
        "\nNothing else to set up. Try:\n\n"
        "    superfit fit spectrum.flm --z 0.1 --bank {}\n".format(name)
    )
    return 0


def _pack_installed(directory, args):
    """Pack a bank just installed, reporting rather than raising on failure.

    A fit that reads the text files opens a thousand of them, which on a
    parallel filesystem is most of a cold run. Doing it at install time means
    the first fit is already fast, and a bank that cannot be packed says so
    now rather than costing every run a few seconds in silence.
    """

    from superfit import packed

    if not getattr(args, "quiet", False):
        print("\nPacking the bank so fits open two files instead of a thousand")
    try:
        packed.pack(
            directory, quiet=getattr(args, "quiet", False), jobs=getattr(args, "jobs", None)
        )
    except (packed.PackError, OSError) as exc:
        print("  WARNING: could not pack the bank: {}".format(exc))
        print("  Fits will read the template files directly, just slower.")


def _bank_directory(args):
    """The bank a `superfit bank` subcommand should act on."""

    from superfit import paths

    if getattr(args, "dir", None):
        return args.dir
    return paths.find_bank_dir(name=getattr(args, "bank", None) or None)


def run_bank_pack(args):
    from superfit import packed

    directory = _bank_directory(args)

    try:
        indexes = packed.pack(directory, quiet=args.quiet, jobs=args.jobs)
    except (packed.PackError, OSError) as exc:
        print("Could not pack {}: {}".format(directory, exc))
        return 1

    print(
        "Packed {} templates from {} directories.".format(
            sum(index["n_templates"] for index in indexes), len(indexes)
        )
    )
    return 0


def run_bank_verify(args):
    from superfit import packed

    directory = _bank_directory(args)
    print("Verifying {} against its packs\n".format(directory))

    results = packed.verify(directory, quiet=args.quiet, jobs=args.jobs)

    failed = False
    for relative, ok, detail in results:
        mark = "ok  " if ok else ("--  " if ok is None else "FAIL")
        if ok is False:
            failed = True
        print("  [{}] {:<28} {}".format(mark, relative, detail))

    if failed:
        print("\nAt least one pack no longer matches the text it was built "
              "from. Rebuild with `superfit bank pack`.")
        return 1

    print("\nEvery pack matches the text it was built from.")
    return 0


def run_bank_status(args):
    from superfit import bank

    report = bank.status(_bank_directory(args) if getattr(args, "bank", None) else args.dir)

    if not report["found"]:
        print("No template bank found.\n")
        print(report["error"])
        print("\nInstall one with:\n\n    superfit bank install\n")
        print("It would go in {}".format(report["install_dir"]))
        return 1

    print("Template bank: {}".format(report["path"]))
    print("  located via: {}".format(report["source"]))

    if not report["complete"]:
        print("  INCOMPLETE: missing {}".format(", ".join(report["missing"])))
        print("\nReinstall with:\n\n    superfit bank install --overwrite\n")
        return 1

    print(
        "  contents:    {} supernova types, {} galaxy templates".format(
            report["n_sn_types"], report["n_galaxy_templates"]
        )
    )

    pack = report["packed"] or {}
    if pack.get("packed") and not pack.get("unpacked"):
        print("  packed:      all {} directories, in {}".format(
            len(pack["packed"]), pack["root"]))
    elif pack.get("packed"):
        print("  packed:      {} of {} directories, in {}".format(
            len(pack["packed"]),
            len(pack["packed"]) + len(pack["unpacked"]),
            pack["root"],
        ))
        print("               not packed: {}".format(", ".join(pack["unpacked"])))
        print("               run `superfit bank pack` to finish")
    else:
        print("  packed:      no -- fits read the template text files")
        print("               `superfit bank pack` takes several seconds off "
              "every run")

    manifest = report["manifest"]
    if manifest:
        print("  installed:   {}".format(manifest.get("installed_at", "unknown")))
        print("  sha256:      {}".format(manifest.get("sha256", "unknown")))
        if manifest.get("sha256") != _known_checksum():
            print(
                "  NOTE: this is not the archive this version of superfit "
                "knows about."
            )
    else:
        print("  installed:   by hand (no install manifest)")

    return 0


def _known_checksum():
    from superfit import bank

    return bank.BANK_SHA256


# -- doctor ----------------------------------------------------------------

# (name, import name, why it is needed)
REQUIREMENTS = [
    ("numpy", "numpy", "arrays"),
    ("scipy", "scipy", "interpolation and binning"),
    ("astropy", "astropy", "result tables"),
    ("pandas", "pandas", "reading results and bank metadata"),
    ("matplotlib", "matplotlib", "plots"),
    ("extinction", "extinction", "the CCM89 extinction law"),
    ("PyAstronomy", "PyAstronomy", "spectral utilities"),
    ("tqdm", "tqdm", "progress bars"),
]

OPTIONAL = [
    ("threadpoolctl", "threadpoolctl", "stops worker processes fighting over BLAS threads"),
]


def run_doctor(args):
    """Check everything a fit needs, and report each item separately."""

    failures = 0
    warnings = 0

    def report(ok, label, detail=""):
        mark = "ok  " if ok else "FAIL"
        print("[{}] {}{}".format(mark, label, "  -- " + detail if detail else ""))

    print("superfit doctor\n")

    from superfit import __version__

    print("superfit {} from {}".format(__version__, _package_location()))
    print("python {}\n".format(sys.version.split()[0]))

    if sys.version_info < (3, 9):
        report(False, "python 3.9 or newer", "found {}".format(sys.version.split()[0]))
        failures += 1
    else:
        report(True, "python version")

    for label, module, why in REQUIREMENTS:
        version = _import_version(module)
        if version is None:
            report(False, label, "not importable; needed for {}".format(why))
            failures += 1
        else:
            report(True, label, version)

    for label, module, why in OPTIONAL:
        version = _import_version(module)
        if version is None:
            print("[note] {} is not installed  -- {}".format(label, why))
            warnings += 1
        else:
            report(True, label, version)

    print()

    from superfit import bank

    status = bank.status()
    if not status["found"]:
        report(False, "template bank", "not found; run `superfit bank install`")
        failures += 1
    elif not status["complete"]:
        report(
            False,
            "template bank",
            "{} is missing {}".format(status["path"], ", ".join(status["missing"])),
        )
        failures += 1
    else:
        report(
            True,
            "template bank",
            "{} ({} SN types)".format(status["path"], status["n_sn_types"]),
        )

    writable, detail = _writable(bank.default_install_dir())
    if writable:
        report(True, "install directory", detail)
    else:
        print("[note] install directory not writable  -- {}".format(detail))
        warnings += 1

    backend, interactive = _matplotlib_backend()
    report(True, "matplotlib backend", backend)
    if not interactive:
        print(
            "[note] {} cannot display plots; --show will do nothing, but "
            "saved plots are unaffected".format(backend)
        )
        warnings += 1

    print()
    if failures:
        print("{} problem(s) found.".format(failures))
        return 1
    print("Everything checks out{}.".format(
        "" if not warnings else ", with {} note(s)".format(warnings)
    ))
    return 0


def _package_location():
    import superfit

    return os.path.dirname(superfit.__file__)


def _import_version(name):
    """The installed version of ``name``, or None if it cannot be imported."""

    try:
        module = __import__(name)
    except Exception:
        return None
    return str(getattr(module, "__version__", "installed"))


def _writable(directory):
    import tempfile

    probe = directory if os.path.isdir(directory) else os.path.dirname(str(directory))
    while probe and not os.path.isdir(probe):
        parent = os.path.dirname(probe)
        if parent == probe:
            break
        probe = parent

    try:
        with tempfile.NamedTemporaryFile(dir=probe or "."):
            pass
    except OSError as exc:
        return False, "{}: {}".format(probe, exc)
    return True, str(directory)


def _matplotlib_backend():
    try:
        import matplotlib

        backend = matplotlib.get_backend()
    except Exception as exc:  # pragma: no cover - matplotlib is a hard dep
        return "unavailable ({})".format(exc), False
    return backend, backend.lower() not in ("agg", "pdf", "ps", "svg", "template")


# -- config ----------------------------------------------------------------

# The settings most runs actually touch, and what each one is for. Written
# into a generated parameter file as underscore-prefixed sibling keys, which
# load_config ignores -- JSON has no comments, and a configuration nobody can
# read is how you end up with a CLI flag for everything.
COMMON_KEYS = [
    ("object_to_fit", "Path to the spectrum. Two or three columns: wavelength, flux, error."),
    ("z_exact", "The object's redshift."),
    ("use_exact_z", "true to fit at z_exact; false to scan z_range_begin..z_range_end."),
    ("z_range_begin", "First redshift of the scan; ignored when use_exact_z is true."),
    ("z_range_end", "Last redshift of the scan."),
    ("z_int", "Redshift step of the scan."),
    ("resolution", "Binning resolution in Angstroms. 10 and 30 use the pre-binned bank."),
    ("fit_stars", "Fit the bank's stars, alone at z=0. 'auto' means when the bank has them."),
    ("fit_qsos", "Fit the bank's QSOs, alone over the z grid. 'auto' means when the bank has them."),
    ("error_spectrum", "'sg', 'linear', or 'included' to use the file's own error column."),
    ("mask_galaxy_lines", "Mask host emission lines. Needs a single redshift."),
    ("mask_telluric", "Mask the telluric A band, 7594-7680 A observed."),
    ("minimum_overlap", "Least fraction of the spectrum a template must cover, 0 to 1."),
    ("saving_results_path", "Where run directories are created. Empty means here."),
    ("how_many_plots", "Plot this many of the best fits."),
    ("profile", "'legacy' reproduces published results; 'modern' weights the solve."),
]


def run_config_create(args):
    from superfit.config import DEFAULT_CONFIG, load_config

    if os.path.exists(args.path) and not args.force:
        raise ValueError(
            "{} already exists. Pass --force to overwrite it.".format(args.path)
        )

    effective = load_config(profile=args.profile)

    if args.full:
        described = dict(COMMON_KEYS)
        keys = [(k, described.get(k, "")) for k in sorted(DEFAULT_CONFIG)]
    else:
        keys = COMMON_KEYS

    document = {
        "_about": (
            "superfit parameters. Keys starting with an underscore are "
            "comments and are ignored. Every setting not named here takes "
            "its default; run `superfit config show` to see them all."
        )
    }
    for key, description in keys:
        if description:
            document["_" + key] = description
        document[key] = effective[key]

    with open(args.path, "w") as handle:
        json.dump(document, handle, indent=4)
        handle.write("\n")

    print("Wrote {}".format(args.path))
    print("\nRun it with:\n\n    superfit fit --config {}\n".format(args.path))
    return 0


def run_config_show(args):
    from superfit.config import load_config

    print(json.dumps(load_config(args.source, profile=args.profile), indent=2))
    return 0


# -- entry point -----------------------------------------------------------


def rewrite_legacy(argv):
    """Accept ``superfit parameters.json``, which predates the subcommands."""

    if not argv or argv[0].startswith("-") or argv[0] in SUBCOMMANDS:
        return argv

    first = argv[0]
    if first.strip().startswith("{") or first.lower().endswith(".json"):
        return ["fit", "--config", first] + argv[1:]

    return argv


def main(argv=None, prog="superfit"):
    argv = list(sys.argv[1:] if argv is None else argv)

    parser = build_parser(prog)
    args = parser.parse_args(rewrite_legacy(argv))

    if getattr(args, "handler", None) is None:
        parser.print_help()
        return 2

    from superfit.bank import BankError
    from superfit.config import ConfigError
    from superfit.output import OutputExistsError

    try:
        return args.handler(args)
    except (BankError, ConfigError, OutputExistsError, OSError, ValueError) as exc:
        # OSError covers the output directory being unwritable or out of
        # quota. It used to reach the user as a traceback, or -- worse --
        # disguised as "already exists", which sent them after an --overwrite
        # flag that could not have helped.
        print("error: {}".format(exc), file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
