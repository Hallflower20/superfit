import glob
import numpy as np
import os
import pandas as pd
import csv
from superfit.paths import MJD_MAX_BRIGHTNESS_CSV, sne_dir


def JD(mjd):
    return float(mjd) + 2400000.5

def list_folders(path):
    if path[-1] != '/':
        path=path+'/'

    folders=[]
    dirs=glob.glob(path+'*')
    for dir in dirs:
        if os.path.isdir(dir):
            folders.append(dir)

    return folders

class Metadata(object):

    def __init__(self, parameters):

        mjd_max_brightness = MJD_MAX_BRIGHTNESS_CSV



        with open(mjd_max_brightness, mode='r') as inp:
            reader = csv.reader(inp)
            band_dictionary = {rows[0]:rows[2] for rows in reader}


        with open(mjd_max_brightness, mode='r') as inp:
            reader = csv.reader(inp)
            MJD_dictionary = {rows[0]:rows[1] for rows in reader}




        folders = [os.path.join(sne_dir(), x) for x in parameters.temp_sn_tr]
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
                idx=subpath.rfind('/')
                sub=subpath[(idx+1):]
                subfolders.append(subpath)
                idx2=subpath[0:idx].rfind('/')
                sn_type=subpath[idx2+1:idx]
                Type_dic[sub]=sn_type
                if os.path.exists(subpath+'/wiserep_spectra.csv'):
                    have_wiserep.append(subpath)
                    # pandas parses these ~190 files about 5x faster than
                    # astropy.io.ascii, which was the single largest cost in
                    # building the metadata.
                    wise=pd.read_csv(subpath+'/wiserep_spectra.csv')
                    path_dic[sub]=subpath
                    z_dic[sub]=wise['Redshift'].iloc[0]
                    coord_dic[sub]=np.array(list(wise[['Obj. RA','Obj. DEC']].iloc[0]))



                    JD_dic[sub]=np.array(wise['JD'])
                    obs_date_dict[sub]=np.array(wise['Obs-date'])
                    spec_file_dic[sub]=np.array(wise['Ascii file'])
                    inst_dic[sub]=np.array(wise['Instrument'])
                    lis=[]
                    for i,spec_file in enumerate(spec_file_dic[sub]):



                        if float(MJD_dictionary[sub]) == -1:

                            phase = 'u'

                        else:

                            phase = float(wise['JD'].iloc[i]) - JD(float(MJD_dictionary[sub]))

                            phase = round(phase,2)


                        if parameters.epoch_high == parameters.epoch_low:

                            band = band_dictionary[sub]

                            shorhand_dict[spec_file]=sn_type + '/' + sub + '/' + wise['Instrument'].iloc[i]+' phase-band : '+ str(phase) + str(band)

                            short_path_dict[shorhand_dict[spec_file]]=spec_file

                            dictionary_all_trunc_objects[spec_file] = os.path.join(sne_dir(), sn_type, sub, spec_file)



                        else:

                            if phase!='u' and phase >= parameters.epoch_low and phase <= parameters.epoch_high:

                                band = band_dictionary[sub]

                                shorhand_dict[spec_file]=sn_type + '/' + sub + '/' + wise['Instrument'].iloc[i]+' phase-band : '+ str(phase) + str(band)

                                short_path_dict[shorhand_dict[spec_file]]=spec_file

                                dictionary_all_trunc_objects[spec_file] = os.path.join(sne_dir(), sn_type, sub, spec_file)



                else:
                    no_wiserep.append(subpath)

        self.shorhand_dict = shorhand_dict
        self.no_wiserep = no_wiserep
        self.dictionary_all_trunc_objects = dictionary_all_trunc_objects


_cached_metadata = None
_cached_key = None


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
    if _cached_metadata is not None and _cached_key == key:
        return _cached_metadata

    # Built into a local and returned from the local. Returning the global
    # instead handed back whatever a concurrent caller had just cached --
    # a different set of templates, silently -- and left the cache holding
    # one fit's metadata under another fit's key for the rest of the process.
    metadata = Metadata(parameters)
    _cached_metadata = metadata
    _cached_key = key
    return metadata
