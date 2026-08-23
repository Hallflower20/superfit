"""Backwards-compatible alias for the package now called ``superfit``.

``import NGSF`` and ``from NGSF.sf_class import Superfit`` keep working, so
existing scripts and notebooks do not have to change at once. Both names
refer to the *same* module objects -- this registers the submodules under
their old names in ``sys.modules`` rather than re-importing them, so state
such as the loaded configuration is shared and ``isinstance`` checks hold
across the two spellings.

New code should import from ``superfit``.
"""

import importlib
import sys
import warnings

_SUBMODULES = [
    "auxiliary",
    "cli",
    "error_routines",
    "get_metadata",
    "Header_Binnings",
    "params",
    "paths",
    "SF_functions",
    "sf_class",
    "spectrum",
]

warnings.warn(
    "The 'NGSF' package has been renamed to 'superfit'. 'import NGSF' still "
    "works but will be removed in a future release; import from 'superfit' "
    "instead.",
    DeprecationWarning,
    stacklevel=2,
)

_superfit = importlib.import_module("superfit")


def __getattr__(name):
    """Resolve NGSF.<name> to superfit.<name>, submodule or attribute alike.

    Lazy for the same reason superfit itself is: touching the alias should
    not drag in matplotlib and scipy, and it must not fail on a submodule
    that happens to be unimportable for unrelated reasons.
    """

    if name in _SUBMODULES:
        module = importlib.import_module("superfit.{}".format(name))
        # Register under the old dotted name too, so `from NGSF.x import y`
        # and `import superfit.x` yield the identical module object rather
        # than two copies with independent state.
        sys.modules["{}.{}".format(__name__, name)] = module
        globals()[name] = module
        return module

    return getattr(_superfit, name)


def __dir__():
    return sorted(set(list(globals()) + _SUBMODULES + list(_superfit.__all__)))


__all__ = list(_superfit.__all__)
__path__ = list(getattr(_superfit, "__path__", []))
__version__ = _superfit.__version__
