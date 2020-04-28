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
from shutil import copyfile
import os

# Folder where tifs are located
InputFolder = './out/sice_2020_fortran/'
OutFolder = 'out/SICE_2020_py/'
#
#os.mkdir(OutFolder)

#%% Copying files into processing folder   
def transfer_to_out_folder(old_name,new_name):
    copyfile(InputFolder+old_name, 
             OutFolder+new_name )

for i in range(21):
    copyfile(InputFolder+'olci_toa_toa'+str(i+1)+'.tif', 
             OutFolder+'r_TOA_'+str(i+1).zfill(2)+'.tif' )
#    bv.OutlookRaster(dat,('Oa'+str(i+1).zfill(2)))

transfer_to_out_folder('olci_toa_ozone.tif','O3.tif')
transfer_to_out_folder('olci_toa_water.tif','WV.tif')
transfer_to_out_folder('olci_toa_sza.tif','SZA.tif')
transfer_to_out_folder('olci_toa_saa.tif','SAA.tif')
transfer_to_out_folder('olci_toa_vza.tif','OZA.tif')
transfer_to_out_folder('olci_toa_vaa.tif','OAA.tif')
transfer_to_out_folder('olci_toa_height.tif','height.tif')
copyfile('masks/greenland.tif',  OutFolder+'mask.tif' )

#%% turning geotiffs into sice input
print("Reading input ")

import rasterio as rio
import numpy as np
from numpy import genfromtxt

InputFolder = 'out/SICE_test_2020_v2/'
Oa01 = rio.open(InputFolder+'r_TOA_01.tif')
meta = Oa01.meta

def WriteOutput(var,var_name,in_folder):
    # this functions write tif files based on a model file, here "Oa01"
    # opens a file for writing
    with rio.open(in_folder+var_name+'.tif', 'w+', **meta) as dst:
        dst.write(var.astype('float32'),1)
    
toa = np.tile(Oa01.read(1)*np.nan, (21,1,1))

for i in range(21):
    dat = rio.open((InputFolder+'r_TOA_'+str(i+1).zfill(2)+'.tif'))
    toa[i,:,:] = dat.read(1)
    
ozone = rio.open(InputFolder+'O3.tif').read(1)
water = rio.open(InputFolder+'WV.tif').read(1)
sza = rio.open(InputFolder+'SZA.tif').read(1)
saa = rio.open(InputFolder+'SAA.tif').read(1)
vza = rio.open(InputFolder+'OZA.tif').read(1)
vaa = rio.open(InputFolder+'OAA.tif').read(1)
height = rio.open(InputFolder+'height.tif').read(1)

sza[np.isnan(toa[0,:,:])] = np.nan
saa[np.isnan(toa[0,:,:])] = np.nan
vza[np.isnan(toa[0,:,:])] = np.nan
vaa[np.isnan(toa[0,:,:])] = np.nan

altitude = genfromtxt(InputFolder+'altitude.csv',  delimiter='   ')
print(len(altitude))
len(olci_toa[:,0])

#c                                           READING OLCI DATA:           
#          read(1,*,end=87) ns,alat,alon,sza,vza,saa,vaa,height,
#     c                   (toa(iks),iks=1,21),OZON,WATER
olci_toa = np.vstack((np.arange(1,len(sza.flatten())+1), 
                     sza.flatten()*np.nan, 
                     sza.flatten()*np.nan, 
                     sza.flatten(),
                     vza.flatten(),
                     saa.flatten(),
                     vaa.flatten(),
                     height.flatten()))

for i in range(21):
    olci_toa = np.vstack((olci_toa, toa[i,:,:].flatten()))

olci_toa=olci_toa.T
olci_toa = olci_toa[np.logical_not(np.isnan(toa[0,:,:]).flatten()),:]
olci_toa[np.isnan(olci_toa)] = 999
OutputFolder='out/sice_2020_fortran/'
np.savetxt(OutputFolder+'olci_toa.dat', X=olci_toa, delimiter='\t', fmt='%i '+ '%10.5f '*5 + '%i '+ '%10.5f'*22)


#%% 
folder = 'out/SICE_2020_py/'
import time
from sice import sice

start_time = time.time()
sice(folder)
time_new = time.time() - start_time

#%%
from sice import sice_old

start_time = time.time()
sice_old(folder)
time_old = time.time() - start_time

#%%  input tif (old format)
#print("Reading input ")
#
#Oa01 = rio.open(InputFolder+'Oa01_reflectance.tif')
#meta = Oa01.meta  
#    
#olci_data = np.tile(Oa01.read(1).flatten()*np.nan, (21,1)).transpose()
#
#for i in range(21):
#    dat = rio.open((InputFolder+'Oa'+str(i+1).zfill(2)+'_reflectance.tif'))
##    bv.OutlookRaster(dat,('Oa'+str(i+1).zfill(2)))
#    olci_data[:,i] = dat.read(1).flatten()
#    
#ozone_dat = rio.open(InputFolder+'ozone.tif')
#ozone = ozone_dat.read(1).flatten()
#water_dat = rio.open(InputFolder+'water.tif')
#water = water_dat.read(1).flatten()
#ozone = ozone_dat.read(1).flatten()
#sza_dat = rio.open(InputFolder+'SZA.tif')
#sza = sza_dat.read(1).flatten()
#saa_dat = rio.open(InputFolder+'SAA.tif')
#saa = saa_dat.read(1).flatten()
#vza_dat = rio.open(InputFolder+'OZA.tif')
#vza = vza_dat.read(1).flatten()
#vaa_dat = rio.open(InputFolder+'OAA.tif')
#vaa = vaa_dat.read(1).flatten()
#height_dat = rio.open(InputFolder+'height.tif')
#height = height_dat.read(1).flatten()
#
#water_vod = genfromtxt('tg_water_vod.dat', delimiter='   ')
#voda = water_vod[range(21),1]
#
#ozone_vod = genfromtxt('tg_vod.dat', delimiter='   ',skip_header=2)
#tozon = ozone_vod[range(21),1]
#aot = 0.1

#% Generating sice intput (old format)
#tmp =np.arange(0, sza.size, 1)
#ind = np.array(np.where(~np.isnan(olci_data[:,1])))
#olci_data_newformat =  np.vstack((tmp, tmp,tmp, sza, vaa, saa,vaa,height,olci_data.transpose(), ozone, water)).transpose()
##                                ns, alat alon,  sza,vza,saa,vaa,height, (toa(iks),iks=1,21),   ozone,water
#olci_data_newformat=olci_data_newformat[ind[0,:],:]
#np.savetxt("out/SICE_try/olci_toa_newformat.dat", olci_data_newformat, delimiter="    ", fmt='%06.6f')

#%% new format


#%%  You can now run sice.f


#%% Loading result from sice_debug.f

#def output_sice_f(file_name,var_name,var_id):
#    sice_out = np.loadtxt(file_name)
#    var = sza*0+999
#    var[np.array(np.where(np.isnan(olci_data[:,1])))] = np.nan
#    var[np.array(np.where(~np.isnan(olci_data[:,1])))] = sice_out[:,var_id]
#    x,y,z = WriteOutput(var,(var_name+'_f'),InputFolder)
#    return x,y,z

#x, y, rp1_f = output_sice_f("../SnowProcessor/2.2/bba.dat",'rp1',5)
#x, y, rp2_f = output_sice_f("../SnowProcessor/2.2/bba.dat",'rp2',6)
#x, y, rs3_f = output_sice_f("../SnowProcessor/2.2/bba.dat",'rs3',7)
#x, y, rs1_f = output_sice_f("../SnowProcessor/2.2/bba.dat",'rs1',8)
#x, y, rs2_f = output_sice_f("../SnowProcessor/2.2/bba.dat",'rs2',9)
#x, y, isnow_f = output_sice_f("../SnowProcessor/2.2/bba.dat",'isnow',10)
#output_sice_f("../SnowProcessor/2.2/size.dat",'D',4)
#output_sice_f("../SnowProcessor/2.2/size.dat",'area',5)
#output_sice_f("../SnowProcessor/2.2/size.dat",'al',6)
#output_sice_f("../SnowProcessor/2.2/size.dat",'r0',7)
    