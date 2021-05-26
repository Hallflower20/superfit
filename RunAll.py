#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys 
from SF_functions import *
from Header_Binnings import *
from params import *
import numpy as np 
import glob
from astropy.table import Table
# Enter path of object of interest, can also be specified as input



#original=sys.argv[1]


# In[2]:


objects_path = '/mnt/c/Users/20xha/Documents/Caltech/Research/superfit/spectra_nonan/'
objects_list = glob.glob(objects_path+'*.ascii') 


# In[3]:


sample = Table.read("/mnt/c/Users/20xha/Documents/Caltech/Research/superfit/ML_sample.ascii", format = "ascii")
sample.rename_column('col1', 'ZTF_Name')
sample.rename_column('col2', "Class")
sample.rename_column('col3', "redshift")
sample.rename_column('col8', "Version")


# In[5]:


counter = 820
print(counter)
for row in sample[820:]:
    try:
        with open("parameters.json", "r") as read_file:
            data = json.load(read_file)


        # Choose saving paths for binned data and results 

        path = data['path']
        sys.path.insert(1,path)
        save_bin_path     = path + "bin/"
        save_results_path = path + "results/"




        # Path where original bank is located for metadata

        original_bank_path = path + 'bank/original_resolution/sne/'

        #--------------------------------------------------------------------------------------------------

        # Select a range and number of steps for z

        z_start  = data['z_start'] 
        z_end    = data['z_end']
        z_num    = data['z_num']



        redshift      =    np.linspace(z_start, z_end,z_num)



        # Number of steps for A_v (do not change)

        alam_num = 21
        extconstant   =    np.linspace(-2,2,alam_num)



        # What part of the library do you want to look at?  

        temp_gal_tr = data['temp_gal_tr']
        temp_sn_tr  = data['temp_sn_tr']


        # Select a wavelength range and resolution
        resolution = data['resolution']
        upper      = data['upper']
        lower      = data['lower']
        interval   = int((upper - lower)/resolution)
        lam        =     np.linspace(lower, upper, interval)
        # Kind of error spectrum ('SG', 'linear' or 'included')
        kind = data['kind']
        # To show plot? 
        show = data['show']   
        # How many top results to plot? 
        n = data['n']
        
        if (counter % 250 == 0):
            print("{} / {}".format(counter, len(sample)))
        counter += 1
        
        filename = row["Version"]
        original = objects_path + filename
        if(original in objects_list):
            n = 3
            resolution = data['resolution']
            upper      = data['upper']
            lower      = data['lower']
            
            maxmin_lambda = Table.read(original, format = "ascii")["col1"]
            min_lambda = min(maxmin_lambda)
            max_lambda = max(maxmin_lambda)
            
            if(upper < max_lambda):
                upper = max_lambda + 100
            if(lower > min_lambda):
                lower = min_lambda - 100

            interval   = int((upper - lower)/resolution)
            lam        =     np.linspace(lower, upper, interval)
            
            redshift_0 = redshifts_all.iloc[np.where(redshifts_all["ZTFID"] == row["ZTF_Name"])[0]]["redshift"].values[0]
            redshift_1 = row["redshift"]
            if(redshift_0 != '-'):
                redshift_0 = float(redshift_0)
                redshift = np.linspace(redshift_0 - 0.05 * redshift_0, redshift_0 + 0.05 * redshift_0, 3)
            elif(not(np.isnan(redshift_1))):
                redshift = np.linspace(redshift_1 - 0.05 * redshift_1, redshift_1 + 0.05 * redshift_1, 3)
            else:
                redshift = np.linspace(0, 0.2, 21)

            try:
                resolution=10
                binned_name= obj_name_int(original, lam, resolution)[3]
                print('Running optimization for spectrum file: {0} with resolution = {1} Ang'.format(binned_name,resolution))
                #Obtaining the binned file name (obj to be analyzed)
                save_bin = save_bin_path + binned_name
                #Calling the original file, getting rid of the header and binning it (default set to 20A)
                kill_header_and_bin(original,resolution, save_bin = save_bin)
                #Core superfit function on the binned file, default to plot and save n fits
                all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
                lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show)
            except:
                try:
                    resolution=20
                    print('Failed. Retrying with resolution = {0} Ang'.format(resolution))

                    #Obtaining the binned file name (obj to be analyzed)
                    save_bin = save_bin_path + binned_name
                    #Calling the original file, getting rid of the header and binning it (default set to 20A)
                    kill_header_and_bin(original,resolution, save_bin = save_bin)
                    #Core superfit function on the binned file, default to plot and save n fits
                    all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
                    lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show)
                except:
                    resolution=30
                    print('Failed. Retrying with resolution = {0} Ang'.format(resolution))

                    #Obtaining the binned file name (obj to be analyzed)
                    try:
                        save_bin = save_bin_path + binned_name
                    except: 
                        import ipdb; ipdb.set_trace()
                    #Calling the original file, getting rid of the header and binning it (default set to 20A)
                    kill_header_and_bin(original,resolution, save_bin = save_bin)
                    #Core superfit function on the binned file, default to plot and save n fits
                    all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
                    lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show)
    except:
        print('An error has occured when trying to optimize for spectrum file {}. Inspect input spectrum and parameters. Stopping at {}'.format(binned_name, counter - 1))


# In[ ]:




