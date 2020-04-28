# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 07:26:21 2020

@author: bav
"""
#%matplotlib qt
import matplotlib.pyplot as plt
plt.close('all')
plt.rcParams.update({'font.size': 18})

import sice_lib as sl
import bav_lib as bl
import numpy as np
from numpy import genfromtxt
import sice_lib as sl
import rasterio as rio
from constants import w, bai, sol1_clean, sol2, sol3_clean, sol1_pol, sol3_pol, asol

InputFolder = 'out/SICE_2020_py/'


#%%

#% ========= input tif ================
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

water_vod = genfromtxt('./tg_water_vod.dat', delimiter='   ')
voda = water_vod[range(21),1]

ozone_vod = genfromtxt('./tg_vod.dat', delimiter='   ',skip_header=2)
tozon = ozone_vod[range(21),1]
aot = 0.1

# declaring variables
BXXX, isnow, D, area, al, r0, isnow, conc, ntype, rp1, rp2, rp3, rs1, rs2, rs3 =  \
vaa*np.nan, vaa*np.nan, vaa*np.nan, vaa*np.nan, vaa*np.nan, vaa*np.nan, \
vaa*np.nan, vaa*np.nan, vaa*np.nan, vaa*np.nan, vaa*np.nan, vaa*np.nan, \
vaa*np.nan, vaa*np.nan, vaa*np.nan

alb_sph, rp, refl =  toa*np.nan, toa*np.nan, toa*np.nan

# =========== ozone scattering  ====================================
BXXX, toa_cor_o3 = sl.ozone_scattering(ozone,tozon, sza, vza,toa)

# Filtering pixels unsuitable for retrieval
isnow[sza>75] = 100
isnow[toa_cor_o3[20, :,:] < 0.1] = 102
for i_channel in range(21):
    toa_cor_o3[i_channel, np.logical_not(np.isnan(isnow))] = np.nan

vaa[np.logical_not(np.isnan(isnow))] = np.nan
saa[np.logical_not(np.isnan(isnow))] = np.nan
sza[np.logical_not(np.isnan(isnow))] = np.nan
vza[np.logical_not(np.isnan(isnow))] = np.nan
height = height.astype(float)
height[np.logical_not(np.isnan(isnow))] = np.nan

# =========== view geometry and atmosphere propeties  ==============
raa, am1, am2, ak1, ak2, amf, co = sl.view_geometry(vaa, saa, sza, vza, aot, height)
tau, p, g = sl.aerosol_properties(aot, height, co)
        
# =========== snow properties  ====================================
D, area, al, r0, bal = sl.snow_properties(toa_cor_o3, ak1, ak2)

# =========== snow properties  ====================================
D, area, al, r0, bal = sl.snow_properties(toa_cor_o3, ak1, ak2)
# filtering small D
D_thresh = 0.1
isnow[D<D_thresh] = 104

for i in range(21):
    toa_cor_o3[i,D<D_thresh] = np.nan
area[D<D_thresh] = np.nan
al[D<D_thresh] = np.nan
r0[D<D_thresh] = np.nan
bal[D<D_thresh] = np.nan
am1[D<D_thresh] = np.nan
am2[D<D_thresh] = np.nan
#D[D<D_thresh] = np.nan

# =========== clean snow  ====================================
# for that we calculate the theoretical reflectance at band 1 of a surface with:
# r0 = 1, a (albedo) = 1, ak1 = 1, ak2 = 1
t1, t2, ratm, r, astra, rms = sl.prepare_coef(tau, g, p, am1, am2, amf)
rs_1 = sl.alb2rtoa(1, t1[0,:,:], t2[0,:,:], np.ones_like(r0), np.ones_like(ak1), 
                   np.ones_like(ak2), ratm[0,:,:], r[0,:,:])

# we then compare it to the observed toa[0] value
ind_clean = toa_cor_o3[0,:,:] >= rs_1
isnow[ind_clean] = 0

ind_all_clean = np.logical_or(isnow == 0, isnow == 7)

# ============== CalCULATION OF BBA of clean snow
BBA_v = np.vectorize(sl.BBA_calc_clean)
p1,p2,s1,s2 = BBA_v(al[ind_all_clean], ak1[ind_all_clean])
# shortwave(0.3-2.4 micron)
rp3[ind_all_clean]=(p1+p2)/sol3_clean
rs3[ind_all_clean]=(s1+s2)/sol3_clean
rs3[~ind_all_clean]=np.nan
rs3[~ind_all_clean]=np.nan
# approximation
# planar albedo
#rp1 and rp2 not derived anymore1
rp3_2 = rp3
rp3_2[ind_all_clean]=sl.plane_albedo_sw_approx(D[ind_all_clean],am1[ind_all_clean])
#     spherical albedo
#rs1 and rs2 not derived anymore
rs3_2 = rs3
rs3_2[ind_all_clean]= sl.spher_albedo_sw_approx(D[ind_all_clean])

#%% Check of planar and spherical shortwave albedo
fig,ax = plt.subplots(2,2,figsize=(15, 8))
plt.subplots_adjust(wspace = 0.4, hspace=0.4)
ax[0,0].scatter(D[ind_all_clean],rs3[ind_all_clean],label='integrated')
ax[0,0].set_xlabel('D')
ax[0,0].set_ylabel('SW spherical albedo')
x = np.linspace(np.nanmin(D),np.nanmax(D),100)
ax[0,0].plot(x,sl.spher_albedo_sw_approx(x),label='approximated')
ax[0,0].legend()

ax[0,1].scatter(rs3[ind_all_clean],rs3_2[ind_all_clean],label='integrated')
ax[0,1].set_xlabel('integrated')
ax[0,1].set_ylabel('approximated')
ax[0,1].set_title('SW spherical albedo')

ax[1,0].scatter(D[ind_all_clean],rp3[ind_all_clean],label='integrated')
ax[1,0].set_xlabel('D')
ax[1,0].set_ylabel('SW planar albedo')
x = np.linspace(np.nanmin(D),np.nanmax(D),100)
ax[1,0].plot(x,sl.plane_albedo_sw_approx(x,np.nanmean(am1)),label='approximated')
ax[1,0].legend()
ax[1,0].text(1.5,0.8,'am1 = '+str(np.nanmean(am1)))

ax[1,1].scatter(rp3[ind_all_clean],rp3_2[ind_all_clean],label='integrated')
ax[1,1].set_xlabel('integrated')
ax[1,1].set_ylabel('approximated')
ax[1,1].set_title('SW planar albedo')
fig.savefig('clean_SW_BBA_approx.png',bbox_inches='tight')

plt.figure()
plt.hist(D.flatten(),bins= 50,label='D')
#plt.hist(al.flatten(),bins= 50,label='al')
plt.hist(am1.flatten(),bins= 50,label='am1')
# plt.xlim([0, 1])
plt.ylim([0, 3500])
plt.xlim([0, 2])
plt.legend()
plt.ylabel('count')
plt.xlabel('value')

#%% Comparison output
InputFolder_py = 'out/SICE_2020_py/'
InputFolder_fortran = 'out/SICE_2020_py - Kopi/'
#InputFolder_fortran = 'out/SICE_2020_py//'
tmp=bl.heatmap(rio.open(InputFolder_py+'diagnostic_retrieval.tif').read(1),'isnow_py')
isnow = rio.open(InputFolder_py+'diagnostic_retrieval.tif').read(1)
ind_all_clean = np.logical_or(isnow == 0, isnow == 7)

var_list = ('albedo_bb_planar_sw','albedo_bb_spherical_sw')
for i in range(len(var_list)):
    var_py = rio.open(InputFolder_py+var_list[i]+'.tif').read(1)
    var_f = rio.open(InputFolder_fortran+var_list[i]+'.tif').read(1)
    diff = var_f-var_py
    diff[~ind_all_clean] = np.nan
    plt.figure(figsize=(10,15))
    bl.heatmap(diff,'Approximation - Integration',col_lim=(-0.01, 0.01),cmap_in='seismic')
    plt.title(var_list[i])
    plt.savefig(var_list[i]+'_diff.png',bbox_inches='tight')

#%%
plt.close('all')
plt.figure()
plt.scatter(al,D)
    
