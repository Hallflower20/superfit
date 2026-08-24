"""Reading a spectrum off disk: ascii, csv, FITS, and the units question.

A spectrum read a factor of ten off in wavelength is not an error message,
it is a misclassification, so these tests care most about the cases where a
file says something the reader could plausibly get wrong: nm instead of
Angstroms, inverse variance instead of an uncertainty, a wavelength axis
that exists only as WCS keywords.
"""

import numpy as np
import pytest

from superfit import Spectrum
from superfit.io import SpectrumReadError, read_spectrum, to_angstrom

pytest.importorskip("astropy.io.fits")

LAM = np.linspace(4000.0, 8000.0, 200)
FLUX = 1.0 + 0.1 * np.sin(LAM / 200.0)
ERR = np.full_like(FLUX, 0.02)


@pytest.fixture
def ascii_file(tmp_path):
    path = tmp_path / "spectrum.flm"
    with open(path, "w") as handle:
        handle.write("# wavelength flux fluxerr\n")
        handle.write("# instrument: a telescope\n")
        for row in zip(LAM, FLUX, ERR):
            handle.write("{} {} {}\n".format(*row))
    return path


class TestAscii:
    def test_columns_come_out_in_order(self, ascii_file):
        wavelength, flux, error = read_spectrum(ascii_file)

        np.testing.assert_allclose(wavelength, LAM)
        np.testing.assert_allclose(flux, FLUX)
        np.testing.assert_allclose(error, ERR)

    def test_a_third_column_is_not_thrown_away(self, ascii_file):
        """It used to be, which broke error_spectrum='included' silently."""

        assert read_spectrum(ascii_file)[2] is not None

    def test_two_columns_give_no_error(self, tmp_path):
        path = tmp_path / "two.flm"
        np.savetxt(path, np.column_stack([LAM, FLUX]))

        assert read_spectrum(path)[2] is None

    def test_columns_can_be_given_by_position(self, tmp_path):
        path = tmp_path / "wide.txt"
        np.savetxt(path, np.column_stack([np.arange(LAM.size), LAM, FLUX, ERR]))

        wavelength, flux, error = read_spectrum(path, columns=[1, 2, 3])

        np.testing.assert_allclose(wavelength, LAM)
        np.testing.assert_allclose(error, ERR)

    def test_columns_can_be_given_by_header_name(self, ascii_file):
        wavelength, flux, _ = read_spectrum(ascii_file, columns=["wavelength", "flux"])

        np.testing.assert_allclose(wavelength, LAM)

    def test_an_unknown_name_lists_what_the_file_has(self, ascii_file):
        with pytest.raises(SpectrumReadError, match="wavelength, flux, fluxerr"):
            read_spectrum(ascii_file, columns=["nope", "flux"])

    def test_a_column_past_the_end_is_reported(self, tmp_path):
        path = tmp_path / "two.flm"
        np.savetxt(path, np.column_stack([LAM, FLUX]))

        with pytest.raises(SpectrumReadError, match="the file has 2 columns"):
            read_spectrum(path, columns=[0, 7])

    def test_a_file_with_no_data_is_reported(self, tmp_path):
        path = tmp_path / "empty.flm"
        path.write_text("# only a header\n")

        with pytest.raises(SpectrumReadError, match="no data rows"):
            read_spectrum(path)


class TestAsciiIsNotSilentlyLenient:
    def test_a_ragged_file_is_refused_rather_than_truncated(self, tmp_path):
        """One short line used to drop a column from the whole file."""

        path = tmp_path / "ragged.flm"
        with open(path, "w") as handle:
            for row in zip(LAM, FLUX, ERR):
                handle.write("{} {} {}\n".format(*row))
            handle.write("9000.0 1.0\n")

        with pytest.raises(SpectrumReadError, match="line 201 has 2 columns"):
            read_spectrum(path)

    def test_an_out_of_range_error_column_is_refused(self, tmp_path):
        """`--columns 0,1,3` on a 3-column file used to yield no error."""

        path = tmp_path / "three.flm"
        np.savetxt(path, np.column_stack([LAM, FLUX, ERR]))

        with pytest.raises(SpectrumReadError, match="asks for error in column 3"):
            read_spectrum(path, columns=[0, 1, 3])


class TestUnits:
    @pytest.mark.parametrize(
        "unit,factor",
        [("nm", 10.0), ("micron", 1e4), ("um", 1e4), ("AA", 1.0), ("angstrom", 1.0)],
    )
    def test_wavelengths_are_converted_to_angstroms(self, tmp_path, unit, factor):
        path = tmp_path / "s.txt"
        np.savetxt(path, np.column_stack([LAM / factor, FLUX]))

        wavelength, _, _ = read_spectrum(path, wavelength_unit=unit)

        np.testing.assert_allclose(wavelength, LAM)

    def test_log_wavelengths_are_understood(self, tmp_path):
        """SDSS and friends store log10(lambda / Angstrom)."""

        path = tmp_path / "s.txt"
        np.savetxt(path, np.column_stack([np.log10(LAM), FLUX]))

        wavelength, _, _ = read_spectrum(path, wavelength_unit="log10(AA)")

        np.testing.assert_allclose(wavelength, LAM)

    def test_an_unknown_unit_lists_the_known_ones(self):
        with pytest.raises(SpectrumReadError, match="nm"):
            to_angstrom([1.0], "furlongs")

    def test_no_unit_means_angstroms(self):
        np.testing.assert_allclose(to_angstrom([5000.0], None), [5000.0])


class TestCsv:
    def test_named_columns(self, tmp_path):
        path = tmp_path / "s.csv"
        path.write_text(
            "wavelength,flux,fluxerr\n"
            + "\n".join("{},{},{}".format(*row) for row in zip(LAM, FLUX, ERR))
        )

        wavelength, flux, error = read_spectrum(path)

        np.testing.assert_allclose(wavelength, LAM)
        np.testing.assert_allclose(error, ERR)

    def test_unusual_names_can_be_pointed_at(self, tmp_path):
        path = tmp_path / "s.csv"
        path.write_text(
            "pixel,my_lambda,my_flux\n"
            + "\n".join(
                "{},{},{}".format(i, w, f) for i, (w, f) in enumerate(zip(LAM, FLUX))
            )
        )

        wavelength, flux, _ = read_spectrum(path, columns=["my_lambda", "my_flux"])

        np.testing.assert_allclose(wavelength, LAM)

    def test_columns_it_cannot_identify_are_reported(self, tmp_path):
        path = tmp_path / "s.csv"
        path.write_text("a,b\n1,2\n3,4\n")

        with pytest.raises(SpectrumReadError, match="Could not find the wavelength"):
            read_spectrum(path)


def write_fits_table(path, columns, units=None):
    from astropy.io import fits

    units = units or {}
    cols = [
        fits.Column(name=name, format="D", unit=units.get(name), array=values)
        for name, values in columns.items()
    ]
    hdu = fits.BinTableHDU.from_columns(cols)
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)
    return path


class TestFitsTables:
    def test_a_binary_table_is_read_by_column_name(self, tmp_path):
        path = write_fits_table(
            tmp_path / "s.fits", {"wave": LAM, "flux": FLUX, "err": ERR}
        )

        wavelength, flux, error = read_spectrum(path)

        np.testing.assert_allclose(wavelength, LAM)
        np.testing.assert_allclose(flux, FLUX)
        np.testing.assert_allclose(error, ERR)

    def test_units_are_taken_from_the_header(self, tmp_path):
        path = write_fits_table(
            tmp_path / "s.fits",
            {"wave": LAM / 10.0, "flux": FLUX},
            units={"wave": "nm"},
        )

        wavelength, _, _ = read_spectrum(path)

        np.testing.assert_allclose(wavelength, LAM)

    def test_an_explicit_unit_overrides_the_header(self, tmp_path):
        """The caller correcting a mislabelled file has to win."""

        path = write_fits_table(
            tmp_path / "s.fits", {"wave": LAM, "flux": FLUX}, units={"wave": "nm"}
        )

        wavelength, _, _ = read_spectrum(path, wavelength_unit="AA")

        np.testing.assert_allclose(wavelength, LAM)

    def test_loglam_is_recognised_without_a_unit(self, tmp_path):
        path = write_fits_table(
            tmp_path / "s.fits", {"loglam": np.log10(LAM), "flux": FLUX}
        )

        wavelength, _, _ = read_spectrum(path)

        np.testing.assert_allclose(wavelength, LAM)

    def test_inverse_variance_becomes_an_uncertainty(self, tmp_path):
        ivar = 1.0 / ERR**2
        path = write_fits_table(
            tmp_path / "s.fits", {"wave": LAM, "flux": FLUX, "ivar": ivar}
        )

        _, _, error = read_spectrum(path)

        np.testing.assert_allclose(error, ERR)

    def test_zero_inverse_variance_is_not_an_infinite_error(self, tmp_path):
        ivar = 1.0 / ERR**2
        ivar[:5] = 0.0
        path = write_fits_table(
            tmp_path / "s.fits", {"wave": LAM, "flux": FLUX, "ivar": ivar}
        )

        _, _, error = read_spectrum(path)

        assert np.isnan(error[:5]).all()
        assert np.isfinite(error[5:]).all()

    def test_a_named_error_column_beats_a_guessed_ivar(self, tmp_path):
        """The caller's explicit choice must not lose to name-guessing.

        A file carrying both, where the ivar is mostly masked, gave an
        all-NaN error -- and an all-NaN sigma makes every chi2 infinite,
        which sorts into memory order and reports the first template in the
        bank as a confident match.
        """

        masked_ivar = np.zeros_like(ERR)
        path = write_fits_table(
            tmp_path / "s.fits",
            {"wave": LAM, "flux": FLUX, "fluxerr": ERR, "ivar": masked_ivar},
        )

        _, _, error = read_spectrum(path, columns=["wave", "flux", "fluxerr"])

        np.testing.assert_allclose(error, ERR)

    def test_ivar_can_still_be_asked_for_explicitly(self, tmp_path):
        ivar = 1.0 / ERR**2
        path = write_fits_table(
            tmp_path / "s.fits",
            {"wave": LAM, "flux": FLUX, "fluxerr": ERR * 10, "ivar": ivar},
        )

        _, _, error = read_spectrum(
            path, columns={"wavelength": "wave", "flux": "flux", "ivar": "ivar"}
        )

        np.testing.assert_allclose(error, ERR)

    def test_a_named_hdu_can_be_asked_for(self, tmp_path):
        from astropy.io import fits

        good = fits.BinTableHDU.from_columns(
            [
                fits.Column(name="wave", format="D", array=LAM),
                fits.Column(name="flux", format="D", array=FLUX),
            ],
            name="SPECTRUM",
        )
        decoy = fits.BinTableHDU.from_columns(
            [fits.Column(name="wave", format="D", array=LAM * 2)], name="JUNK"
        )
        path = tmp_path / "s.fits"
        fits.HDUList([fits.PrimaryHDU(), decoy, good]).writeto(path)

        wavelength, _, _ = read_spectrum(path, hdu="SPECTRUM")

        np.testing.assert_allclose(wavelength, LAM)

    def test_later_hdus_are_tried_when_the_first_is_unusable(self, tmp_path):
        from astropy.io import fits

        decoy = fits.BinTableHDU.from_columns(
            [fits.Column(name="unrelated", format="D", array=LAM)]
        )
        good = fits.BinTableHDU.from_columns(
            [
                fits.Column(name="wave", format="D", array=LAM),
                fits.Column(name="flux", format="D", array=FLUX),
            ]
        )
        path = tmp_path / "s.fits"
        fits.HDUList([fits.PrimaryHDU(), decoy, good]).writeto(path)

        wavelength, _, _ = read_spectrum(path)

        np.testing.assert_allclose(wavelength, LAM)

    def test_a_file_with_nothing_readable_says_what_it_tried(self, tmp_path):
        path = write_fits_table(tmp_path / "s.fits", {"unrelated": LAM})

        with pytest.raises(SpectrumReadError, match="No HDU"):
            read_spectrum(path)


class TestFitsImages:
    def _image(self, tmp_path, header):
        from astropy.io import fits

        hdu = fits.PrimaryHDU(FLUX)
        hdu.header.update(header)
        path = tmp_path / "s.fits"
        hdu.writeto(path, overwrite=True)
        return path

    def test_a_linear_wavelength_wcs(self, tmp_path):
        step = LAM[1] - LAM[0]
        path = self._image(
            tmp_path, {"CRVAL1": LAM[0], "CDELT1": step, "CRPIX1": 1}
        )

        wavelength, flux, _ = read_spectrum(path)

        np.testing.assert_allclose(wavelength, LAM)
        np.testing.assert_allclose(flux, FLUX)

    def test_crpix_is_one_based(self, tmp_path):
        """CRVAL holds at pixel CRPIX, counting from 1 as FITS does."""

        step = LAM[1] - LAM[0]
        path = self._image(
            tmp_path, {"CRVAL1": LAM[10], "CDELT1": step, "CRPIX1": 11}
        )

        wavelength, _, _ = read_spectrum(path)

        np.testing.assert_allclose(wavelength, LAM)

    def test_cd1_1_is_accepted_as_the_step(self, tmp_path):
        step = LAM[1] - LAM[0]
        path = self._image(tmp_path, {"CRVAL1": LAM[0], "CD1_1": step, "CRPIX1": 1})

        wavelength, _, _ = read_spectrum(path)

        np.testing.assert_allclose(wavelength, LAM)

    def test_a_log_linear_axis(self, tmp_path):
        # Evenly spaced in log10(lambda), which is what DC-FLAG announces.
        log = np.linspace(np.log10(4000.0), np.log10(8000.0), FLUX.size)
        path = self._image(
            tmp_path,
            {"CRVAL1": log[0], "CDELT1": log[1] - log[0], "CRPIX1": 1, "DC-FLAG": 1},
        )

        wavelength, _, _ = read_spectrum(path)

        np.testing.assert_allclose(wavelength, 10.0**log, rtol=1e-10)

    def test_ctype_log_is_understood_as_well_as_dc_flag(self, tmp_path):
        log = np.linspace(np.log10(4000.0), np.log10(8000.0), FLUX.size)
        path = self._image(
            tmp_path,
            {
                "CRVAL1": log[0],
                "CDELT1": log[1] - log[0],
                "CRPIX1": 1,
                "CTYPE1": "AWAV-LOG",
            },
        )

        wavelength, _, _ = read_spectrum(path)

        np.testing.assert_allclose(wavelength, 10.0**log, rtol=1e-10)

    def test_a_named_multispec_stack_uses_the_first_band(self, tmp_path):
        """Naming the HDU is how you say "yes, row 0 is the flux"."""

        from astropy.io import fits

        step = LAM[1] - LAM[0]
        stack = np.vstack([FLUX, FLUX * 99, ERR])
        hdu = fits.PrimaryHDU(stack)
        hdu.header.update({"CRVAL1": LAM[0], "CDELT1": step, "CRPIX1": 1})
        path = tmp_path / "s.fits"
        hdu.writeto(path)

        wavelength, flux, _ = read_spectrum(path, hdu=0)

        np.testing.assert_allclose(flux, FLUX)
        np.testing.assert_allclose(wavelength, LAM)

    def test_an_image_with_no_wcs_is_reported(self, tmp_path):
        path = self._image(tmp_path, {})

        with pytest.raises(SpectrumReadError, match="No HDU"):
            read_spectrum(path)

    def test_a_missing_dispersion_keyword_is_refused(self, tmp_path):
        """It used to default to 1 A/pixel, silently mis-registering the file."""

        path = self._image(tmp_path, {"CRVAL1": 3800.0})

        with pytest.raises(SpectrumReadError, match="wavelength step is unknown"):
            read_spectrum(path)

    def test_a_zero_dispersion_falls_through_to_the_next_keyword(self, tmp_path):
        """Some IRAF writers set CDELT1 = 0 alongside a real CD1_1."""

        step = LAM[1] - LAM[0]
        path = self._image(
            tmp_path,
            {"CRVAL1": LAM[0], "CDELT1": 0.0, "CD1_1": step, "CRPIX1": 1},
        )

        wavelength, _, _ = read_spectrum(path)

        np.testing.assert_allclose(wavelength, LAM)

    def test_a_two_dimensional_image_is_not_guessed_at(self, tmp_path):
        """Row 0 of a long-slit frame is sky, not the spectrum."""

        from astropy.io import fits

        step = LAM[1] - LAM[0]
        hdu = fits.PrimaryHDU(np.vstack([FLUX * 99, FLUX]))
        hdu.header.update({"CRVAL1": LAM[0], "CDELT1": step, "CRPIX1": 1})
        path = tmp_path / "s.fits"
        hdu.writeto(path)

        with pytest.raises(SpectrumReadError, match="2-dimensional"):
            read_spectrum(path)


class TestThroughSpectrum:
    def test_from_file_reads_fits(self, tmp_path):
        path = write_fits_table(
            tmp_path / "SN2011fe.fits", {"wave": LAM, "flux": FLUX, "err": ERR}
        )

        spectrum = Spectrum.from_file(path)

        assert len(spectrum) == LAM.size
        assert spectrum.error is not None

    def test_the_fits_suffix_is_not_kept_in_the_name(self, tmp_path):
        path = write_fits_table(tmp_path / "SN2011fe.fits", {"wave": LAM, "flux": FLUX})

        assert Spectrum.from_file(path).name == "SN2011fe"

    def test_reading_options_are_refused_for_in_memory_input(self):
        with pytest.raises(TypeError, match="read from a file"):
            Spectrum.coerce(
                np.column_stack([LAM, FLUX]), wavelength_unit="nm"
            )

    def test_an_array_still_coerces(self):
        assert len(Spectrum.coerce(np.column_stack([LAM, FLUX]))) == LAM.size
