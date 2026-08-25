"""superfit -- spectral classification of supernovae against a template bank.

Fit an observed spectrum against a library of supernova and host-galaxy
templates over a grid of redshift and extinction, and rank the matches by
chi2.

Quick start, from arrays already in memory::

    from superfit import Superfit

    fit = Superfit(wavelength=lam, flux=flux, error=err, z=0.127)
    results = fit.run()          # a pandas DataFrame, best match first

or from a file::

    from superfit import Spectrum, Superfit

    spectrum = Spectrum.from_file("SN2021urb.flm")
    results = Superfit(spectrum, z=0.127).run()

For many spectra in one process, a ``Session`` prepares the bank once instead
of once per fit::

    from superfit import Session

    with Session(bank="modern-curated", lower_lam=3500, upper_lam=9000) as s:
        for path in paths:
            print(s.fit(path, z=0.1).results.iloc[0]["SN"])

Everything else -- the template lists, the A_v grid, the error model -- has a
default, and any of it can be overridden by passing ``config=`` a dict or by
keyword. See ``superfit.config.DEFAULT_CONFIG``.
"""

__version__ = "0.1.0"

__all__ = [
    "DEFAULT_CONFIG",
    "FitResult",
    "OutputExistsError",
    "Session",
    "Spectrum",
    "Superfit",
    "__version__",
    "bank_is_available",
    "load_config",
    "set_bank_dir",
]

# Imported lazily (PEP 562): pulling in Superfit costs matplotlib, astropy and
# scipy, which is a lot to pay for `import superfit` when all the caller wants
# is the version number.
_LAZY = {
    "Superfit": ("superfit.sf_class", "Superfit"),
    "Session": ("superfit.session", "Session"),
    "Spectrum": ("superfit.spectrum", "Spectrum"),
    "FitResult": ("superfit.output", "FitResult"),
    "OutputExistsError": ("superfit.output", "OutputExistsError"),
    "DEFAULT_CONFIG": ("superfit.config", "DEFAULT_CONFIG"),
    "load_config": ("superfit.config", "load_config"),
    "set_bank_dir": ("superfit.paths", "set_bank_dir"),
    "bank_is_available": ("superfit.paths", "bank_is_available"),
}


def __getattr__(name):
    if name in _LAZY:
        import importlib

        module_name, attr = _LAZY[name]
        return getattr(importlib.import_module(module_name), attr)
    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name))


def __dir__():
    return sorted(list(globals()) + list(_LAZY))
