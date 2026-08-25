import glob
import numpy as np
import os
import pandas as pd
import csv
import threading
from superfit.paths import MJD_MAX_BRIGHTNESS_CSV, mjd_max_brightness_csv, sne_dir


def JD(mjd):
    return float(mjd) + 2400000.5


def _as_float(value):
    """``value`` as a float, or None if it is not one.

    Phase arithmetic needs a real number at both ends. A blank cell, a NaN,
    a stray word in a hand-edited table: all of them mean the same thing here
    -- there is no epoch to measure from -- and all of them used to raise or,
    worse, propagate NaN into a comparison that silently comes out False.
    """

    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if np.isnan(number) else number

def list_folders(path):
    """Immediate subdirectories of ``path``.

    Built with os.path.join rather than by appending "/" to the string: on
    Windows the separator is a backslash, and a path assembled by hand is one
    the rest of this module then fails to split apart again.
    """

    return sorted(
        entry
        for entry in glob.glob(os.path.join(path, "*"))
        if os.path.isdir(entry)
    )

class Metadata(object):

    def __init__(self, parameters):

        # This fit's bank, and this fit's phase table -- both settled by
        # Parameters, not looked up from process-wide state that another
        # Superfit's construction may since have moved. A bank built
        # elsewhere knows its own objects; the copy inside the package only
        # knows the 189 the legacy bank was built from.
        bank_dir = getattr(parameters, "bank_dir", None)
        mjd_max_brightness = getattr(parameters, "phase_table", None)
        if mjd_max_brightness is None:
            mjd_max_brightness = mjd_max_brightness_csv(bank_dir)

        sne_root = sne_dir(bank_dir=bank_dir)



        with open(mjd_max_brightness, mode='r') as inp:
            reader = csv.reader(inp)
            band_dictionary = {rows[0]:rows[2] for rows in reader}


        with open(mjd_max_brightness, mode='r') as inp:
            reader = csv.reader(inp)
            MJD_dictionary = {rows[0]:rows[1] for rows in reader}




        folders = [os.path.join(sne_root, x) for x in parameters.temp_sn_tr]
        have_wiserep=[]
        no_wiserep=[]
        z_dic={}
        path_dic={}
        dictionary_all_trunc_objects ={}
        JD_dic={}
        coord_dic={}
        spec_file_dic={}
        inst_dic={}
        obs_date_dict={}
        shorhand_dict={}
        Type_dic={}
        subfolders=[]
        short_path_dict={}
        for folder in folders:
            subs=list_folders(folder)
            for sub in subs:
                subpath=sub
                # basename/dirname rather than rfind('/'): the bank lives
                # under a Windows path as readily as a POSIX one, and there
                # the separator these used to search for never appears.
                sub=os.path.basename(subpath)
                subfolders.append(subpath)
                sn_type=os.path.basename(os.path.dirname(subpath))
                Type_dic[sub]=sn_type
                wiserep_csv=os.path.join(subpath, 'wiserep_spectra.csv')
                if os.path.exists(wiserep_csv):
                    have_wiserep.append(subpath)
                    # pandas parses these ~190 files about 5x faster than
                    # astropy.io.ascii, which was the single largest cost in
                    # building the metadata.
                    wise=pd.read_csv(wiserep_csv)

                    # A header with no rows is a bank saying "no spectra
                    # survived curation for this object", which is a
                    # statement, not a fault. Reading .iloc[0] off it raised
                    # IndexError and took the whole fit down; one bank has 51
                    # such objects. Treated as an object with no metadata,
                    # which is what the else-branch below already handles.
                    if len(wise) == 0:
                        no_wiserep.append(subpath)
                        continue

                    path_dic[sub]=subpath
                    z_dic[sub]=wise['Redshift'].iloc[0]
                    coord_dic[sub]=np.array(list(wise[['Obj. RA','Obj. DEC']].iloc[0]))



                    JD_dic[sub]=np.array(wise['JD'])
                    obs_date_dict[sub]=np.array(wise['Obs-date'])
                    spec_file_dic[sub]=np.array(wise['Ascii file'])
                    inst_dic[sub]=np.array(wise['Instrument'])
                    lis=[]
                    for i,spec_file in enumerate(spec_file_dic[sub]):



                        # An object the phase table does not list, or lists
                        # with the -1 sentinel, has an unknown phase -- the
                        # same answer, and one this already understood. It
                        # used to be a bare subscript, so a bank whose table
                        # does not cover every object it ships raised
                        # KeyError instead; the larger banks leave thousands
                        # of objects uncovered. Same for a maximum that will
                        # not parse, and for a spectrum with no epoch of its
                        # own to subtract it from.
                        mjd_peak = _as_float(MJD_dictionary.get(sub))
                        observed_jd = _as_float(wise['JD'].iloc[i])

                        if mjd_peak is None or mjd_peak == -1 or observed_jd is None:

                            phase = 'u'

                        else:

                            phase = round(observed_jd - JD(mjd_peak), 2)


                        if parameters.epoch_high == parameters.epoch_low:

                            band = band_dictionary.get(sub, '')

                            shorhand_dict[spec_file]=sn_type + '/' + sub + '/' + wise['Instrument'].iloc[i]+' phase-band : '+ str(phase) + str(band)

                            short_path_dict[shorhand_dict[spec_file]]=spec_file

                            dictionary_all_trunc_objects[spec_file] = os.path.join(sne_root, sn_type, sub, spec_file)



                        else:

                            if phase!='u' and phase >= parameters.epoch_low and phase <= parameters.epoch_high:

                                band = band_dictionary.get(sub, '')

                                shorhand_dict[spec_file]=sn_type + '/' + sub + '/' + wise['Instrument'].iloc[i]+' phase-band : '+ str(phase) + str(band)

                                short_path_dict[shorhand_dict[spec_file]]=spec_file

                                dictionary_all_trunc_objects[spec_file] = os.path.join(sne_root, sn_type, sub, spec_file)



                else:
                    no_wiserep.append(subpath)

        self.shorhand_dict = shorhand_dict
        self.no_wiserep = no_wiserep
        self.dictionary_all_trunc_objects = dictionary_all_trunc_objects


_cached_metadata = None
_cached_key = None

# Held across the check, the build and the store, so that the three cannot
# be interleaved with another caller's. Returning the built object from a
# local (below) already narrows the window to a single store instruction,
# but "narrow" is not "closed", and a scan that hands back the wrong
# template list is not the kind of bug worth leaving a window for. The lock
# costs nothing: it is taken once per fit, around work that takes a second.
_cache_lock = threading.Lock()


def get_metadata(parameters):
    """Return the bank metadata for ``parameters``, scanning the bank rarely.

    Building this walks every object directory and parses ~190 wiserep CSVs.
    It used to be done twice per run -- once in Superfit.__init__ and once in
    all_parameter_space -- for identical results.

    The cache is keyed on what the scan actually reads (the SN types and the
    epoch window), not on the identity of the Parameters object, so a second
    fit that differs only in redshift reuses the first one's scan.
    """

    global _cached_metadata, _cached_key

    key = parameters.metadata_key

    with _cache_lock:
        if _cached_metadata is not None and _cached_key == key:
            return _cached_metadata

        # Built into a local and returned from the local. Returning the
        # global instead could hand back whatever another caller had stored
        # in between, and leave the cache holding one fit's metadata under
        # another fit's key for the rest of the process.
        metadata = Metadata(parameters)
        _cached_metadata = metadata
        _cached_key = key
        return metadata
