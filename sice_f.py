# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:35:10 2019

This script reads a set of input tif files and outputs the csv needed to run sice.f

You will need to update 
    InputFoder 
    the path to water_vod and ozone_vod
    the output path
    
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

start_time = time.time()

#InputFolder = sys.argv[1] + "/"
InputFolder = './data/5_km_res/'
OutputFolder = './data/5_km_res/fortran/'
ProcessingFolder = './fortran/'

# input text file
if os.path.isfile(InputFolder):
    InputFolder = os.path.dirname(InputFolder) + "/"
    # data_in = pd.read_csv('validation/data/S3_PROMICE.csv')
    data_in = pd.read_csv(InputFolder)
    toa = np.expand_dims(
        data_in[[c for c in data_in.columns if c.find("reflec") >= 0]]
        .to_numpy()
        .transpose(),
        axis=2,
    )

    ozone = np.expand_dims(data_in["total_ozone"], axis=1)
    water = np.expand_dims(data_in["total_columnar_water_vapour"], axis=1)
    sza = np.expand_dims(data_in["sza"], axis=1)
    saa = np.expand_dims(data_in["saa"], axis=1)
    vza = np.expand_dims(data_in["vza"], axis=1)
    vaa = np.expand_dims(data_in["vaa"], axis=1)
    height = np.expand_dims(data_in["altitude"], axis=1)

    sza[np.isnan(toa[0, :, :])] = np.nan
    saa[np.isnan(toa[0, :, :])] = np.nan
    vza[np.isnan(toa[0, :, :])] = np.nan
    vaa[np.isnan(toa[0, :, :])] = np.nan

#  input tif 
elif os.path.isdir(InputFolder):
    print("\n tiff input")
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

#  New format

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
# 21 OLCI TOA reflectances (pi*I/cos(sza)/F0,
for i in range(21):
    olci_toa = np.vstack((olci_toa, toa[i, :, :].flatten()))
# height of the surface (m)
olci_toa = np.vstack((olci_toa, height.flatten())) 
# TOTAL OZONE load (ECMWF)
olci_toa = np.vstack((olci_toa, ozone.flatten()))

olci_toa = olci_toa.T
olci_toa_save = olci_toa
ind_good_pixels = np.logical_not(np.isnan(toa[0, :, :]).flatten())
olci_toa = olci_toa[ind_good_pixels, :]
olci_toa[np.isnan(olci_toa)] = 999
OutputFolder = InputFolder + "fortran/"
# os.mkdir(OutputFolder)
print("\nInput file saved: " + OutputFolder + "input.dat")
np.savetxt(
    OutputFolder + "input.dat",
    X=olci_toa,
    delimiter="\t",
    fmt="%i "*2 + "%10.5f " * 29,
)

#%%  You can now run sice.f
import shutil
import glob

shutil.copy(OutputFolder+'input.dat', ProcessingFolder[:-1])
os.chdir('./fortran')
subprocess.check_call('./sice.exe')
for file in glob.glob(r'*.dat'):
    if file == 'thv.dat':
        continue
    if file == 'input.dat':
        continue
    
    print('Moving '+file)
    shutil.move(file, '../'+OutputFolder)
os.chdir('..')


#%% Loading result from sice_debug.f


#WRITE(1004,*)
#     c               j,alat,alon,
#     c       NCLASS,factor,diam,ssa,dlina,rv,
#     c       aload1,powe,polut,
#     c       deff,absor1,absef660,absor1000,
#     c rsw,
#     c  rvis,rnir,
#     c  rsws,
#     c  rviss,rnirs,
#     c  andbi,andsi,ratka,
#     c      NPOLE,
#     c      NBARE,
#     c       NSNOW,
#     c sza,vza,raa,toa(1),toa(21),
#     c       tocos,akozon,difka,cv1,cv2
if os.path.isfile(InputFolder):
    snow_parameters = pd.read_csv(
        OutputFolder + "snow_parameters.dat",
        names=[
            "j", "alat", "alon", "NCLASS",
            "factor", "diam", "ssa", "dlina", "rv", "aload1",
            "powe", "polut", "eff","absor1","absef660","absor1000",
            "rsw", "rvis","rnir", "rsws", "rviss","rnirs",
            "andbi","andsi","ratka", "NPOLE",
            "NBARE", "NSNOW", "sza","vza","raa","toa(1)","toa(21)",
            "tocos","akozon","difka","cv1","cv2"
        ],
        sep=" ",
        skipinitialspace=True,
        header=None,
    )

    data_out = data_in
    data_out["grain_diameter"] = np.nan
    data_out["snow_specific_area"] = np.nan
    data_out["al"] = np.nan
    data_out["r0"] = np.nan
    data_out["diagnostic_retrieval"] = np.nan
    data_out["conc"] = np.nan
    data_out["albedo_bb_planar_sw"] = np.nan
    data_out["albedo_bb_spherical_sw"] = np.nan

    data_out["grain_diameter"][size.ns - 1] = size.D
    # data_out.loc[size.ns-1, ['grain_diameter']]=size.D
    data_out.loc[size.ns - 1, ["snow_specific_area"]] = size.area
    data_out.loc[size.ns - 1, ["al"]] = size.al
    data_out.loc[size.ns - 1, ["r0"]] = size.r0
    data_out.loc[size.ns - 1, ["diagnostic_retrieval"]] = bba.isnow
    data_out.loc[size.ns - 1, ["conc"]] = impurity.conc
    data_out["albedo_bb_planar_sw"][size.ns - 1] = bba.rp3
    data_out["albedo_bb_spherical_sw"] = np.nan
    data_out["albedo_bb_spherical_sw"][size.ns - 1] = bba.rs3

    # data_out.iloc[size.ns-1, data_out.columns.get_loc('grain_diameter')]=size.D
    # data_out.iloc[size.ns-1, data_out.columns.get_loc('snow_specific_area')]=size.area
    # data_out.iloc[size.ns-1, data_out.columns.get_loc('al')]=size.al
    # data_out.iloc[size.ns-1, data_out.columns.get_loc('r0')]=size.r0
    # data_out.iloc[size.ns-1, data_out.columns.get_loc('diagnostic_retrieval')]=bba.isnow
    # data_out.iloc[size.ns-1, data_out.columns.get_loc('conc')]=impurity.conc
    # data_out.iloc[size.ns-1, data_out.columns.get_loc('albedo_bb_planar_sw')]=bba.rp3
    # data_out.iloc[size.ns-1, data_out.columns.get_loc('albedo_bb_spherical_sw')]=bba.rs3

    for i in np.append(np.arange(11), np.arange(15, 21)):
        # for i in np.arange(21):
        data_out["albedo_spectral_spherical_" + str(i + 1).zfill(2)] = np.nan
        data_out.loc[
            size.ns - 1, ["albedo_spectral_spherical_" + str(i + 1).zfill(2)]
        ] = spherical_albedo[:, 4 + i]

    for i in np.append(np.arange(11), np.arange(15, 21)):
        data_out["rBRR_" + str(i + 1).zfill(2)] = np.nan
        data_out.loc[size.ns - 1, ["rBRR_" + str(i + 1).zfill(2)]] = planar_albedo[
            :, 4 + i
        ]

    data_out.to_csv(InputFolder[:-4] + "_fortran_out.csv")
    print("\nOutput: " + InputFolder[:-4] + "_fortran_out.csv")

# ========= input tif ===============
elif os.path.isdir(InputFolder):
    Oa01 = rio.open(InputFolder + "r_TOA_01.tif")
    meta = Oa01.meta
    with rio.Env():
        meta.update(compress="DEFLATE")

    def output_sice_f(file_name, var_name, var_id):
        sice_out = np.loadtxt(file_name)
        ind_orig = np.arange(1, len(Oa01.read(1).flatten()) + 1)
        var = ind_orig * np.nan
        ind_pix = sice_out[:, 0].astype(int)
        var[ind_pix - 1] = sice_out[:, var_id]
        var_mat = np.reshape(var, np.shape(Oa01.read(1)))
        var_mat[var_mat == 999] = np.nan
        with rio.open(OutputFolder + var_name + ".tif", "w+", **meta) as dst:
            dst.write(var_mat.astype("float32"), 1)

    output_sice_f(OutputFolder + "bba.dat", "albedo_bb_planar_sw", 4)
    output_sice_f(OutputFolder + "bba.dat", "albedo_bb_spherical_sw", 5)
    output_sice_f(OutputFolder + "size.dat", "D", 4)
    output_sice_f(OutputFolder + "size.dat", "area", 5)
    output_sice_f(OutputFolder + "size.dat", "al", 6)
    output_sice_f(OutputFolder + "size.dat", "r0", 7)
    output_sice_f(OutputFolder + "bba_alex_reduced.dat", "diagnostic_retrieval", 3)

    for i in range(21):
        output_sice_f(
            OutputFolder + "spherical_albedo.dat",
            "albedo_spectral_spherical_" + str(i + 1).zfill(2),
            4 + i,
        )
        output_sice_f(
            OutputFolder + "planar_albedo.dat",
            "albedo_spectral_planar_" + str(i + 1).zfill(2),
            4 + i,
        )
        output_sice_f(OutputFolder + "boar.dat", "rBRR_" + str(i + 1), 4 + i)

    #%%   plotting oulooks
    try:
        os.mkdir(OutputFolder + "plots")
    except:
        print("folder exist")

    import matplotlib.pyplot as plt

    # fig,ax=bl.heatmap_discrete(rio.open(OutputFolder+'diagnostic_retrieval.tif').read(1),
    #                         'diagnostic_retrieval ')
    # ax.set_title(OutputFolder)
    # fig.savefig(OutputFolder+'plots/diagnostic_retrieval.png',bbox_inches='tight')

    var_list = ("albedo_bb_planar_sw", "albedo_bb_spherical_sw")
    for i in range(len(var_list)):
        var_1 = rio.open(OutputFolder + var_list[i] + ".tif").read(1)
        plt.figure(figsize=(10, 15))
        bl.heatmap(var_1, var_list[i], col_lim=(0, 1), cmap_in="jet")
        plt.title(OutputFolder)
        plt.savefig(
            OutputFolder + "plots/" + var_list[i] + "_diff.png", bbox_inches="tight"
        )
    plt.ioff()
    for i in np.append(np.arange(11), np.arange(21)):
        var_name = "albedo_spectral_spherical_" + str(i + 1).zfill(2)
        var_1 = rio.open(OutputFolder + var_name + ".tif").read(1)
        plt.figure(figsize=(10, 15))
        bl.heatmap(var_1, var_name, col_lim=(0, 1), cmap_in="jet")
        plt.title(OutputFolder)
        plt.savefig(OutputFolder + "plots/" + var_name + ".png", bbox_inches="tight")
    plt.ion()
