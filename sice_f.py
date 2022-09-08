# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:35:10 2019

tip list:
    %matplotlib inline
    %matplotlib qt

@author: bav@geus.dk
"""

import numpy as np
import rasterio as rio
import time
import os
import sys
#import bav_lib as bl
import pandas as pd
import rioxarray
import shutil
import glob
import subprocess
InputFolder = './data/2019-06-14/'
OutputFolder = './data/2019-06-14/fortran/'
ProcessingFolder = './fortran/'
# %% 
start_time = time.time()
InputFolder = InputFolder + "/"
Oa01 = rio.open(InputFolder + "r_TOA_01.tif")
meta = Oa01.meta

def WriteOutput(var, var_name, in_folder):
    # this functions write tif files based on a model file, here "Oa01"
    # opens a file for writing
    with rio.open(in_folder + var_name + ".tif", "w+", **meta) as dst:
        dst.write(var.astype("float32"), 1)

toa = np.tile(Oa01.read(1) * np.nan, (21, 1, 1))

for i in range(21):
    dat = rio.open((InputFolder + "r_TOA_" + str(i + 1).zfill(2) + ".tif"))
    toa[i, :, :] = dat.read(1)

ozone = rio.open(InputFolder + "O3.tif").read(1)
water = rio.open(InputFolder + "WV.tif").read(1)
sza = rio.open(InputFolder + "SZA.tif").read(1)
saa = rio.open(InputFolder + "SAA.tif").read(1)
vza = rio.open(InputFolder + "OZA.tif").read(1)
vaa = rio.open(InputFolder + "OAA.tif").read(1)
height = rio.open(InputFolder + "height.tif").read(1)

sza[np.isnan(toa[0, :, :])] = np.nan
saa[np.isnan(toa[0, :, :])] = np.nan
vza[np.isnan(toa[0, :, :])] = np.nan
vaa[np.isnan(toa[0, :, :])] = np.nan

olci_toa = np.vstack(
    (
        np.arange(1, len(sza.flatten()) + 1),  # pixel number_x
        np.arange(1, len(sza.flatten()) + 1),  # pixel number_y
        np.arange(1, len(sza.flatten()) + 1),  # latitude
        np.arange(1, len(sza.flatten()) + 1),  # longitude
        sza.flatten(),  # solar zenith angle
        saa.flatten(),  # soalr azimuthal angle
        vza.flatten(),  # viewing zenith angle
        vaa.flatten(),  # viewing azimuth angle
    )
)
for i in range(21):
    olci_toa = np.vstack((olci_toa, toa[i, :, :].flatten()))
olci_toa = np.vstack((olci_toa, height.flatten())) 
olci_toa = np.vstack((olci_toa, ozone.flatten()))

olci_toa = olci_toa.T
olci_toa_save = olci_toa
ind_good_pixels = np.logical_not(np.isnan(toa[0, :, :]).flatten())
olci_toa = olci_toa[ind_good_pixels, :]
olci_toa[np.isnan(olci_toa)] = 999

# os.mkdir(OutputFolder)
print("\nInput file saved: " + OutputFolder + "input.dat")
np.savetxt(
    OutputFolder + "input.dat",
    X=olci_toa,
    delimiter="\t",
    fmt="%i "*2 + "%10.5f " * 29,
)

#%%  You can now run sice.f
shutil.copy(OutputFolder+'input.dat', ProcessingFolder[:-1])
os.chdir('./fortran')
print('Running sice.f')
#subprocess.run('gfortran ./version\ 6.1\ 2022/s.f -o ./sice.exe', shell=True)
subprocess.check_call('./sice.exe')
for file in glob.glob(r'*.dat'):
    if file == 'thv.dat':
        continue
    if file == 'input.dat':
        continue
    
    print('Moving '+file)
    shutil.move(file, '../'+OutputFolder+'/'+file)
os.chdir('..')


#%% Converting text output to geotiff
Oa01 = rio.open(InputFolder + "r_TOA_01.tif")
meta = Oa01.meta
with rio.Env():
    meta.update(compress="DEFLATE")

def output_sice_f(file_name, var_name, var_id):
    sice_out = pd.read_csv(file_name, delim_whitespace=True, header=None).values
    ind_orig = np.arange(1, len(Oa01.read(1).flatten()) + 1)
    var = ind_orig * np.nan
    ind_pix = sice_out[:, 1].astype(int)
    var[ind_pix - 1] = sice_out[:, var_id]
    var_mat = np.reshape(var, np.shape(Oa01.read(1)))
    var_mat[var_mat == 999] = np.nan
    with rio.open(OutputFolder + var_name + ".tif", "w+", **meta) as dst:
        dst.write(var_mat.astype("float32"), 1)

output_sice_f(OutputFolder + "snow_parameters.dat", "isnow", 3)
output_sice_f(OutputFolder + "snow_parameters.dat", "factor", 4)
output_sice_f(OutputFolder + "snow_parameters.dat", "grain_diameter", 5)
output_sice_f(OutputFolder + "snow_parameters.dat", "snow_specific_area", 6)
output_sice_f(OutputFolder + "snow_parameters.dat", "al", 7)
output_sice_f(OutputFolder + "snow_parameters.dat", "r0", 8)
output_sice_f(OutputFolder + "snow_parameters.dat", "bm", 10)
output_sice_f(OutputFolder + "snow_parameters.dat", "polut", 11)
output_sice_f(OutputFolder + "snow_parameters.dat", "albedo_bb_planar_sw", 16)
output_sice_f(OutputFolder + "snow_parameters.dat", "albedo_bb_spherical_sw", 19)

i = 0
output_sice_f(
    OutputFolder + "spectral_spherical_albedo.dat",
    "alb_sph_" + str(i + 1).zfill(2),
    5 + i,
)
output_sice_f(
    OutputFolder + "spectral_spherical_albedo_solved.dat",
    "alb_sph_" + str(i + 1).zfill(2) + '_solved',
    5 + i,
)
output_sice_f(
    OutputFolder + "spectral_plane_albedo.dat",
    "alb_pl_" + str(i + 1).zfill(2),
    5 + i,
)

output_sice_f(
    OutputFolder + "spectral_plane_albedo_solved.dat",
    "alb_pl_" + str(i + 1).zfill(2)+ '_solved',
    5 + i,
)
output_sice_f(
    OutputFolder + "spectral_bOAR.dat",
    "BOAR_" + str(i + 1).zfill(2),
    5 + i,
)

output_sice_f(
    OutputFolder + "spectral_BOAR_solved.dat",
    "bOAR_" + str(i + 1).zfill(2)+ '_solved',
    5 + i,
)


#    for i in range(21):
#        output_sice_f(
#            OutputFolder + "spherical_albedo.dat",
#            "albedo_spectral_spherical_" + str(i + 1).zfill(2),
#            4 + i,
#        )
#        output_sice_f(
#            OutputFolder + "planar_albedo.dat",
#            "albedo_spectral_planar_" + str(i + 1).zfill(2),
#            4 + i,
#        )
# format snow_parameters.dat:
#            "j", "alat", "alon", "NCLASS",
#            "factor", "diam", "ssa", "dlina", "rv", "aload1",
#            "powe", "polut", "eff","absor1","absef660","absor1000",
#            "rsw", "rvis","rnir", "rsws", "rviss","rnirs",
#            "andbi","andsi","ratka", "NPOLE",
#            "NBARE", "NSNOW", "sza","vza","raa","toa(1)","toa(21)",
#            "tocos","akozon","difka","cv1","cv2"
# format spectral_spherical_albedo.dat:
# WRITE(1001,*) j,alat,alon, NCLASS,factor,(albs(ir),ir=1,21)

#%%   plotting oulooks
import xarray as xr
import rioxarray
import matplotlib.pyplot as plt
#%matplotlib qt
import numpy as np
data_folder = 'data/2019-06-14/'
code_ver1 = 'pySICEv1.6'
code_ver2 = 'pySICEv2.0'
folder1 = data_folder+code_ver1+'/'
folder2 =  data_folder+code_ver2+'/'
var_list = ['isnow', "grain_diameter", "snow_specific_area",   
            "r0", "albedo_bb_planar_sw",  "albedo_bb_spherical_sw"]
# var_list = ["isnow", "alb_sph_01", "alb_sph_01_solved", "alb_pl_01", "alb_pl_01_solved", "BOAR_01", "BOAR_01_solved"]
plt.close('all')
for var in var_list:
    ds_f = rioxarray.open_rasterio(folder1+var+'.tif').squeeze()
    ds_p = rioxarray.open_rasterio(folder2+var+'.tif').squeeze().interp_like(ds_f, method ='nearest')
    
    if var == 'polut':
        ds_f = np.log10(ds_f)
        ds_p = np.log10(ds_p)
        
    if var in ['albedo_bb_planar_sw', 'albedo_bb_spherical_sw']:
        ds_f = ds_f.where(ds_f>0).where(ds_f<2)
    if var == 'isnow':
        ds_f = ds_f.where(ds_f>=0).where(ds_f<4)

        
    vmin = min(ds_f.min(), ds_p.min())
    vmax = max(ds_f.max(), ds_p.max())
    fig, ax = plt.subplots(2,2, figsize=(15,15))
    ax=ax.flatten()

        
    ds_f.dropna('x', 'all').dropna('y', 'all').plot(ax=ax[0], vmin=vmin, vmax=vmax)
    ax[0].set_title(code_ver1)
    ds_p.dropna('x', 'all').dropna('y', 'all').plot(ax=ax[1], vmin=vmin, vmax=vmax)
    ax[1].set_title(code_ver2)
    
    (ds_p-ds_f).dropna('x', 'all').dropna('y', 'all').plot(ax=ax[2])
    ax[2].set_title(code_ver2+' - '+code_ver1)

    ax[3].plot(ds_f.values.flatten(),
               ds_p.values.flatten(),
               marker ='.', linestyle='None')
    ax[3].plot([ds_f.min(), ds_f.max()], [ds_f.min(), ds_f.max()],color='k')
    ax[3].set_xlabel(code_ver1)
    ax[3].set_ylabel(code_ver2)
    

    
