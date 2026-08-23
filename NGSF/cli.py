#!/usr/bin/env python
"""Command-line entry point for NGSF.

Installed as the ``superfit`` command; also reachable as ``python run.py``
from a source checkout.

    superfit parameters.json
    superfit parameters.json --object spectrum.flm --out results/
"""

import argparse
import sys


def build_parser(prog="superfit"):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Fit a supernova spectrum against the NGSF template bank.",
    )
    parser.add_argument(
        "config",
        help="Path to a JSON parameter file (or a JSON string).",
    )
    parser.add_argument(
        "--object",
        dest="object_to_fit",
        help="Override object_to_fit from the config.",
    )
    parser.add_argument(
        "--out",
        dest="saving_results_path",
        help="Override saving_results_path from the config.",
    )
    parser.add_argument(
        "--resolution",
        type=int,
        help="Override the binning resolution in Angstroms.",
    )
    parser.add_argument(
        "--n-cores",
        type=int,
        dest="n_cores",
        help="Worker processes for the (z, A_v) grid. 0 means auto.",
    )
    parser.add_argument(
        "--weighted-solve",
        action="store_true",
        help=(
            "Solve for the template amplitudes with the same 1/sigma^2 "
            "weights the chi2 uses. Changes results; see the README."
        ),
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plotting entirely (sets how_many_plots to 0).",
    )
    return parser


def apply_overrides(config, args):
    """Fold command-line overrides into a loaded config dict."""

    config = dict(config)

    if args.object_to_fit is not None:
        config["object_to_fit"] = args.object_to_fit
    if args.saving_results_path is not None:
        path = args.saving_results_path
        # Downstream code concatenates this prefix onto file names.
        config["saving_results_path"] = path if path.endswith("/") else path + "/"
    if args.resolution is not None:
        config["resolution"] = args.resolution
    if args.n_cores is not None:
        config["n_cores"] = args.n_cores
    if args.weighted_solve:
        config["weighted_solve"] = 1
    if args.no_plots:
        config["how_many_plots"] = 0
        config["show_plot"] = 0

    return config


def main(argv=None, prog="superfit"):
    args = build_parser(prog).parse_args(argv)

    from NGSF.params import load_config

    config = apply_overrides(load_config(args.config), args)

    from NGSF.sf_class import Superfit

    supernova = Superfit(config)
    supernova.superfit()
    return 0


if __name__ == "__main__":
    sys.exit(main())
