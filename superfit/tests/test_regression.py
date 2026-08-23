"""End-to-end regression against a committed golden result.

This is the test that answers "did anything major change?". It runs the real
pipeline over the real template bank and compares the ranked result table
against ``golden/SN2021urb_exact_z_10A.csv``, which was captured from the
code as it stood before any optimisation work.

Three levels of strictness, so a failure says *what* changed:

* the identity and order of the best-fit templates -- exact match required,
  because a reordering means the science answer moved;
* the chi2 values and fitted scalings -- compared with a tolerance, because
  a legitimate refactor can shift the last bits of a float;
* the top-1 match -- checked on its own, since that is the number anyone
  actually quotes.

Regenerate the golden file only when a change to the science is intended,
and say so in the commit message:

    python superfit/tests/test_regression.py --regenerate
"""

import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

if __name__ == "__main__" and __package__ is None:
    # Allow `python superfit/tests/test_regression.py --regenerate` from anywhere.
    sys.path.insert(
        0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    )

from superfit.tests.conftest import GOLDEN, REPO, TEST_SPECTRUM, needs_bank

GOLDEN_CSV = os.path.join(GOLDEN, "SN2021urb_exact_z_10A.csv")

# The configuration the golden file was produced with, stated here rather
# than read from parameters.json. The golden result then cannot be
# invalidated by someone editing the example config, and this test works
# against an installed wheel, where parameters.json is not present.
GOLDEN_CONFIG = {
    "use_exact_z": 1,
    "z_exact": 0.127,
    "resolution": 10,
    "lower_lam": 0,
    "upper_lam": 0,
    "error_spectrum": "sg",
    "mask_galaxy_lines": 1,
    "mask_telluric": 1,
    "minimum_overlap": 0.7,
    "epoch_low": 0,
    "epoch_high": 0,
    "Alam_low": -2,
    "Alam_high": 2,
    "Alam_interval": 0.2,
    "R_v": 3.1,
    "iterations": 10,
    "weighted_solve": 0,
    "show_plot": 0,
    "how_many_plots": 0,
}

# Columns that must match exactly -- these identify *which* template won.
IDENTITY_COLUMNS = ["GALAXY", "SN", "Z", "A_v", "Phase", "Band"]

# Columns compared numerically.
NUMERIC_COLUMNS = [
    "CONST_SN",
    "CONST_GAL",
    "Frac(SN)",
    "Frac(gal)",
    "CHI2/dof",
    "CHI2/dof2",
]


def run_pipeline(config, out_dir):
    """Run the fit in a subprocess and return the results DataFrame."""

    import json

    os.makedirs(out_dir, exist_ok=True)
    config = dict(config)
    config["object_to_fit"] = TEST_SPECTRUM
    config["saving_results_path"] = out_dir.rstrip("/") + "/"
    config["show_plot"] = 0
    config["how_many_plots"] = 0

    config_path = os.path.join(out_dir, "params.json")
    with open(config_path, "w") as fh:
        json.dump(config, fh)

    env = dict(os.environ)
    env["MPLBACKEND"] = "Agg"

    # Drive the installed entry point rather than run.py, so this works
    # against a wheel as well as a checkout.
    proc = subprocess.run(
        [sys.executable, "-m", "superfit.cli", config_path],
        cwd=out_dir,
        env=env,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise AssertionError(
            "superfit.cli failed ({}):\n{}\n{}".format(
                proc.returncode, proc.stdout[-4000:], proc.stderr[-4000:]
            )
        )

    stem = os.path.basename(TEST_SPECTRUM).rsplit(".", 1)[0]
    return pd.read_csv(os.path.join(out_dir, stem + ".csv"))


@needs_bank
@pytest.mark.endtoend
class TestGoldenResult:
    @pytest.fixture(scope="class")
    def results(self, tmp_path_factory):
        out = str(tmp_path_factory.mktemp("regression"))
        return run_pipeline(GOLDEN_CONFIG, out)

    @pytest.fixture(scope="class")
    def golden(self):
        return pd.read_csv(GOLDEN_CSV)

    def test_same_number_of_results(self, results, golden):
        assert len(results) == len(golden), (
            "result count changed: {} rows now vs {} in the golden file. "
            "A different number of surviving templates usually means the "
            "overlap cut, the dedup step, or the template selection moved."
        ).format(len(results), len(golden))

    def test_best_fit_is_unchanged(self, results, golden):
        """The headline answer: top-ranked template and its redshift."""

        got = results.iloc[0]
        want = golden.iloc[0]
        for col in IDENTITY_COLUMNS:
            assert got[col] == want[col], "top-1 {} changed: {!r} -> {!r}".format(
                col, want[col], got[col]
            )

    def test_ranking_is_unchanged(self, results, golden):
        for col in IDENTITY_COLUMNS:
            assert list(results[col]) == list(golden[col]), (
                "the ranked list of {} changed:\n  golden: {}\n  now:    {}".format(
                    col, list(golden[col]), list(results[col])
                )
            )

    @pytest.mark.parametrize("col", NUMERIC_COLUMNS)
    def test_numbers_match(self, results, golden, col):
        np.testing.assert_allclose(
            results[col].to_numpy(),
            golden[col].to_numpy(),
            rtol=1e-5,
            atol=1e-8,
            err_msg="column {} drifted beyond float noise".format(col),
        )

    def test_chi2_is_finite_and_ordered(self, results):
        chi2 = results["CHI2/dof2"].to_numpy()
        assert np.isfinite(chi2).all(), "NaN or inf chi2 in the results"
        assert (np.diff(chi2) >= 0).all(), "results are not sorted by CHI2/dof2"

    def test_redshift_is_the_requested_one(self, results, golden):
        assert np.allclose(results["Z"], golden["Z"].iloc[0])

    def test_fractions_sum_to_one(self, results):
        """Frac(SN) and Frac(gal) are a partition of the model flux."""

        total = results["Frac(SN)"] + results["Frac(gal)"]
        np.testing.assert_allclose(total, 1.0, rtol=1e-5)


def _regenerate():
    """Rewrite the golden file from a fresh run. Intentional changes only."""

    import tempfile

    tmp = tempfile.mkdtemp(prefix="superfit-golden-")
    try:
        results = run_pipeline(GOLDEN_CONFIG, tmp)
        os.makedirs(GOLDEN, exist_ok=True)
        results.to_csv(GOLDEN_CSV, index=False)
        print("wrote {} ({} rows)".format(GOLDEN_CSV, len(results)))
        print(results.to_string())
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    if "--regenerate" in sys.argv:
        _regenerate()
    else:
        print(__doc__)
