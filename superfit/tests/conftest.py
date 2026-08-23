"""Shared fixtures for the superfit test suite.

The numerical-kernel tests deliberately avoid the template bank so that they
run anywhere in a second or two. Only the end-to-end regression test needs
the bank, and it skips itself when the bank is absent.
"""

import json
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(HERE, "data")
GOLDEN = os.path.join(HERE, "golden")

TEST_SPECTRUM = os.path.join(
    DATA, "SN2021urb_2021-08-06_00-00-00_Keck1_LRIS_TNS.flm"
)


def bank_is_available():
    """True when a template bank with both layouts can be resolved."""

    try:
        from superfit import paths

        return os.path.isdir(paths.original_resolution_dir()) and os.path.isdir(
            paths.binnings_dir()
        )
    except Exception:
        return False


needs_bank = pytest.mark.skipif(
    not bank_is_available(),
    reason="template bank not found; set SUPERFIT_BANK_DIR or unzip supyfit_bank.zip",
)


@pytest.fixture
def rng():
    """Deterministic random numbers -- never seed from the clock in tests."""

    return np.random.default_rng(20240822)


@pytest.fixture
def smooth_spectrum():
    """A clean, well-behaved 2-column spectrum: (wavelength, flux)."""

    lam = np.linspace(4000.0, 8000.0, 600)
    flux = 1.0 + 0.1 * np.sin(lam / 200.0)
    return np.column_stack([lam, flux])


@pytest.fixture
def noisy_spectrum(rng):
    """A smooth continuum plus a realistic amount of white noise."""

    lam = np.linspace(4000.0, 8000.0, 600)
    flux = 1.0 + 0.1 * np.sin(lam / 200.0) + rng.normal(0.0, 0.02, lam.size)
    return np.column_stack([lam, flux])


@pytest.fixture
def toy_fit_problem(rng):
    """A small stand-in for the real fit: observed flux, SN and galaxy banks.

    Shapes match what ``core`` expects internally: ``sn`` is (1, n_sn, n_lam)
    and ``gal`` is (n_gal, 1, n_lam), so the two broadcast into the
    (n_gal, n_sn, n_lam) grid the chi2 is evaluated over.
    """

    n_lam, n_sn, n_gal = 120, 7, 4
    lam = np.linspace(4000.0, 8000.0, n_lam)

    sn = np.abs(rng.normal(1.0, 0.3, (1, n_sn, n_lam)))
    gal = np.abs(rng.normal(1.0, 0.3, (n_gal, 1, n_lam)))

    # Build the observation from a known SN/galaxy pair so there is a true answer.
    truth_sn, truth_gal, b_true, d_true = 3, 2, 0.7, 0.4
    int_obj = b_true * sn[0, truth_sn] + d_true * gal[truth_gal, 0]
    int_obj = int_obj + rng.normal(0.0, 0.01, n_lam)

    sigma = np.full(n_lam, 0.05)

    return {
        "lam": lam,
        "sn": sn,
        "gal": gal,
        "int_obj": int_obj,
        "sigma": sigma,
        "truth_sn": truth_sn,
        "truth_gal": truth_gal,
        "b_true": b_true,
        "d_true": d_true,
    }


@pytest.fixture
def base_parameters():
    """The shipped parameters.json, as a plain dict."""

    with open(os.path.join(REPO, "parameters.json")) as fh:
        return json.load(fh)
