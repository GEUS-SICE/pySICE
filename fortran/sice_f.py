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
import bav_lib as bl

start_time = time.time()

InputFolder = sys.argv[1] + '/'
# InputFolder = 'out/SICE_2020_py1.4/'

print(os.getcwd())


#%% turning geotiffs into sice input
print("Reading input ")

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

# ns,alat,alon,sza,vza,saa,vaa,height,  (toa(iks),iks=1,21),OZON,WATER
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

olci_toa = np.vstack((olci_toa, ozone.flatten()))
olci_toa = np.vstack((olci_toa, water.flatten()))

olci_toa=olci_toa.T
olci_toa_save = olci_toa
ind_good_pixels = np.logical_not(np.isnan(toa[0,:,:]).flatten())
olci_toa = olci_toa[ind_good_pixels,:]
olci_toa[np.isnan(olci_toa)] = 999
OutputFolder= InputFolder+'fortran/'
#os.mkdir(OutputFolder)
np.savetxt(OutputFolder+'olci_toa_newformat.dat', X=olci_toa, delimiter='\t', \
           fmt='%i '+ '%10.5f '*6 + '%i '+ '%10.5f'*23)

#%%  You can now run sice.f
print('bash ./fortran/sice_f.sh -i '+OutputFolder)
os.system('bash ./fortran/sice_f.sh -i '+OutputFolder)


#%% Loading result from sice_debug.f
Oa01 = rio.open(InputFolder+'r_TOA_01.tif')
meta = Oa01.meta
with rio.Env():    
    meta.update(compress='DEFLATE')
def output_sice_f(file_name,var_name,var_id):
    sice_out = np.loadtxt(file_name)
    ind_orig = np.arange(1,len(Oa01.read(1).flatten())+1)
    var = ind_orig*np.nan
    ind_pix = sice_out[:,0].astype(int)
    var[ind_pix-1] = sice_out[:,var_id]
    var_mat = np.reshape(var,np.shape(Oa01.read(1)))
    var_mat[var_mat==999] = np.nan
    with rio.open(OutputFolder+var_name+'.tif', 'w+', **meta) as dst: 
        dst.write(var_mat.astype('float32'),1)
    
output_sice_f(OutputFolder+"bba.dat",'albedo_bb_planar_sw',4)
output_sice_f(OutputFolder+"bba.dat",'albedo_bb_spherical_sw',5)
output_sice_f(OutputFolder+"size.dat",'D',4)
output_sice_f(OutputFolder+"size.dat",'area',5)
output_sice_f(OutputFolder+"size.dat",'al',6)
output_sice_f(OutputFolder+"size.dat",'r0',7)
output_sice_f(OutputFolder+"bba_alex_reduced.dat",'diagnostic_retrieval',3)

for i in range(21):
    output_sice_f(OutputFolder+"spherical_albedo.dat",'albedo_spectral_spherical_'+str(i+1).zfill(2),4+i)
    output_sice_f(OutputFolder+"planar_albedo.dat",'albedo_spectral_planar_'+str(i+1).zfill(2),4+i)
    output_sice_f(OutputFolder+"boar.dat",'rBRR_'+str(i+1),4+i)

# OUTPUT files:
#  file= 'spherical_albedo.dat' )           ns,ndate(3),alat,alon,(answer(i),i=1,21),isnow
#  file= 'planar_albedo.dat'     )          ns,ndate(3),alat,alon,(rp(i),i=1,21),isnow
#  file= 'boar.dat'                    )    ns,ndate(3),alat,alon,(refl(i),i=1,21),isnow
#  file= 'size.dat'                     )   ns,ndate(3),alat,alon,D,area,al,r0, andsi,andbi,indexs,indexi,indexd,isnow
#  file=   'impurity.dat'              )    ns,ndate(3),alat,alon,ntype,conc,bf,bm,thv, toa(1),isnow
#  file=   'bba.dat'                     )  ns,ndate(3),alat,alon,rp3,rs3, isnow
#  file=   'bba_alex_reduced.dat')          ns,ndate(3),rp3,isnow   
#  file=   'notsnow.dat')                   ns,ndate(3),icloud,iice
#  file='retrieved_O3.dat')                 ns,alat,alon,BXXX,totadu,deltak,sza,vza,amf

#%%   plotting oulooks  
try:
    os.mkdir(OutputFolder+'plots')
except:
    print('folder exist')
    
import matplotlib.pyplot as plt
   
# fig,ax=bl.heatmap_discrete(rio.open(OutputFolder+'diagnostic_retrieval.tif').read(1),
#                         'diagnostic_retrieval ')
# ax.set_title(OutputFolder)
# fig.savefig(OutputFolder+'plots/diagnostic_retrieval.png',bbox_inches='tight')

var_list = ('albedo_bb_planar_sw','albedo_bb_spherical_sw')
for i in range(len(var_list)):
    var_1 = rio.open(OutputFolder+var_list[i]+'.tif').read(1)
    plt.figure(figsize=(10,15))
    bl.heatmap(var_1,var_list[i], col_lim=(0, 1) ,cmap_in='jet')
    plt.title(OutputFolder)
    plt.savefig(OutputFolder+'plots/'+var_list[i]+'_diff.png',bbox_inches='tight')
plt.ioff()  
for i in np.append(np.arange(11), np.arange(21)):
    var_name = 'albedo_spectral_spherical_'+ str(i+1).zfill(2)
    var_1 = rio.open(OutputFolder+var_name+'.tif').read(1)
    plt.figure(figsize=(10,15))
    bl.heatmap(var_1,var_name, col_lim=(0, 1) ,cmap_in='jet')
    plt.title(OutputFolder)
    plt.savefig(OutputFolder+'plots/'+var_name+'.png',bbox_inches='tight')
plt.ion()