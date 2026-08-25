"""Fitting many spectra in one process, without re-preparing the bank each time.

A single fit spends most of its time on work that has nothing to do with the
spectrum. On the largest bank: a directory listing to validate the pack,
twelve thousand template reads, host-line masking, and the resample onto the
log grid -- about 2.6 s warm, against 2.0 s of actual fitting. Run a hundred
spectra as a hundred fits and that preparation is paid a hundred times.

A Session pays it once:

    from superfit import Session

    with Session(bank="modern-curated", z=0.1, resolution=10,
                 lower_lam=3500, upper_lam=9000) as session:
        for path in spectra:
            result = session.fit(path)
            print(result.results.iloc[0]["SN"])

Two things are worth knowing about what a Session pins.

The prepared bank is reused only when the *observed grid* matches, because a
redshift is a shift along that grid and templates resampled onto one grid mean
nothing on another. The grid comes from the observation's own wavelength range
unless ``lower_lam`` and ``upper_lam`` are set, so a batch of spectra with
different ranges gets a different grid each and reuses nothing. Setting an
explicit range is what makes a batch share the preparation -- and it is also
the honest thing to do when comparing classifications across spectra, since
otherwise each was fitted over a different wavelength span.

The bank itself is pinned for the session's lifetime: it is validated once, at
the first fit, and not re-listed for each one after. That is a deliberate
trade. Re-listing the largest bank is 15561 stats, and a batch is a thing you
want fitted against *one* bank anyway -- a bank changing halfway through is a
problem to notice, not a change to follow silently. Call
:meth:`revalidate` to check again.
"""

from superfit import packed
from superfit.sf_class import Superfit


class Session:
    """Settings shared by many fits, and the prepared bank they reuse.

    Every keyword is a superfit setting, applied to each fit as a default.
    :meth:`fit` overrides them per spectrum.
    """

    def __init__(self, **settings):
        self._settings = dict(settings)
        self._fits = 0
        self._closed = False

        # The bank is validated on the first fit and pinned after it.
        self._validated = False

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def __repr__(self):
        return "<superfit.Session: {} fit(s), bank={!r}>".format(
            self._fits, self._settings.get("bank") or "default"
        )

    @property
    def settings(self):
        """The shared settings, as a dict. A copy; edit through the session."""

        return dict(self._settings)

    def set(self, **settings):
        """Change the shared settings for subsequent fits."""

        self._settings.update(settings)
        return self

    def revalidate(self):
        """Check the bank against the filesystem again on the next fit.

        A session pins the bank, which is what makes it fast. This is how to
        pick up a bank that has been repacked or edited without starting a new
        process.
        """

        self._validated = False
        packed.begin_fit(revalidate=True)
        return self

    def fit(self, spectrum=None, **overrides):
        """Fit one spectrum with the session's settings.

        Returns the same FitResult a standalone ``Superfit(...).run()`` does.
        """

        if self._closed:
            raise RuntimeError(
                "This Session has been closed; build a new one to fit again."
            )

        config = dict(self._settings)
        config.update(overrides)
        if spectrum is not None:
            config["object_to_fit"] = spectrum

        # The first fit validates the bank; the rest trust it. Anything that
        # would change which bank is in play -- a different name, a repack --
        # goes through revalidate(). A caller who asks explicitly wins.
        config.setdefault("revalidate_bank", not self._validated)

        result = Superfit(config).run()

        self._validated = True
        self._fits += 1
        return result

    def fit_all(self, spectra, **overrides):
        """Fit each spectrum in turn, yielding results as they finish.

        A generator, so a long batch reports progress rather than going quiet
        until the end. Exceptions are not caught: one unreadable spectrum in a
        batch is a thing to see, and the results so far have already been
        yielded and written to disk.
        """

        for spectrum in spectra:
            yield self.fit(spectrum, **overrides)

    def close(self):
        """Release the prepared bank this session was holding.

        The prepared bank for the largest catalogue is around 90 MB of
        resampled templates, which is worth handing back when a long-lived
        process is done fitting.
        """

        from superfit import SF_functions

        SF_functions.forget_prepared_bank()
        self._closed = True
        return self
