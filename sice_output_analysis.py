# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 13:35:44 2020
%matplotlib qt
@author: bav
"""
from IPython import get_ipython
ipython = get_ipython()
ipython.magic("matplotlib qt")

import numpy as np
import rasterio as rio
import bav_lib as bl
import matplotlib.pyplot as plt
import os
plt.close('all')

#%% Running sice
InputFolder = 'out/SICE_2020_py1.4.1/'
InputFolder = 'C:/Users/bav/OneDrive - Geological survey of Denmark and Greenland/Data/Cook_data/Cook_pySICEv1.4_output/20170722T141728/'
ipython.magic('run sice.py '+InputFolder)

#% Plotting output
ipython.magic("matplotlib inline")

try:
    os.mkdir(InputFolder+'plots')
except:
    print('folder exist')
    
fig,ax=bl.heatmap_discrete(rio.open(InputFolder+'diagnostic_retrieval.tif').read(1),  'diagnostic_retrieval ')
ax.set_title(InputFolder)
fig.savefig(InputFolder+ 'plots/diagnostic_retrieval.png', bbox_inches='tight')

var_list = ('albedo_bb_planar_sw','albedo_bb_spherical_sw','r0')
for i in range(len(var_list)):
    var_1 = rio.open(InputFolder+var_list[i]+'.tif').read(1)
    plt.figure(figsize=(10,15))
    bl.heatmap(var_1,var_list[i], col_lim=(0, 1) ,cmap_in='jet')
    plt.title(InputFolder)
    plt.savefig(InputFolder+'plots/'+var_list[i]+'.png',bbox_inches='tight')
    plt.close()
    
var_list = ('O3_SICE', 'grain_diameter', 'snow_specific_area',
            'al','conc')
for i in range(len(var_list)):
    var_1 = rio.open(InputFolder+var_list[i]+'.tif').read(1)
    plt.figure(figsize=(10,15))
    bl.heatmap(var_1,var_list[i],cmap_in='jet')
    plt.title(InputFolder)
    plt.savefig(InputFolder+'plots/'+var_list[i]+'.png',bbox_inches='tight')
    plt.close()

for i in np.append(np.arange(11), np.arange(15,21)):
    var_name = 'albedo_spectral_spherical_'+ str(i+1).zfill(2)
    var_1 = rio.open(InputFolder+var_name+'.tif').read(1)
    plt.figure(figsize=(10,15))
    bl.heatmap(var_1,var_name, col_lim=(0, 1) ,cmap_in='jet')
    plt.title(InputFolder)
    plt.savefig(InputFolder+'plots/'+var_name+'.png',bbox_inches='tight')
    plt.close()
    var_name = 'albedo_spectral_planar_'+ str(i+1).zfill(2)
    var_1 = rio.open(InputFolder+var_name+'.tif').read(1)
    plt.figure(figsize=(10,15))
    bl.heatmap(var_1,var_name, col_lim=(0, 1) ,cmap_in='jet')
    plt.title(InputFolder)
    plt.savefig(InputFolder+'plots/'+var_name+'.png',bbox_inches='tight')
    plt.close()
    var_name = 'rBRR_'+ str(i+1).zfill(2)
    var_1 = rio.open(InputFolder+var_name+'.tif').read(1)
    plt.figure(figsize=(10,15))
    bl.heatmap(var_1,var_name, col_lim=(0, 1) ,cmap_in='jet')
    plt.title(InputFolder)
    plt.savefig(InputFolder+'plots/'+var_name+'.png',bbox_inches='tight')
    plt.close()
ipython.magic("matplotlib qt")

bl.heatmap(isnow)
plt.figure()
bl.heatmap(toa_cor_o3[20, :,:])
plt.figure()
bl.heatmap(toa[20, :,:])
#%% ========== Compare two folders ==================

plt.close('all')
plt.rcParams.update({'font.size': 18})

InputFolder_1 = 'SICE_2020_py1.4/'
InputFolder_2 = 'SICE_2020_py1.4.1/'

tmp=bl.heatmap_discrete(rio.open('out/'+InputFolder_1+'diagnostic_retrieval.tif').read(1),
                        'isnow_py '+InputFolder_1)
try: tmp=bl.heatmap_discrete(rio.open('out/'+InputFolder_2+'diagnostic_retrieval.tif').read(1), 
                        'isnow_py '+InputFolder_2)
except:
    print('invalid diagnostic_retrieval in '+InputFolder_2)
    
plt.figure()
tmp=bl.heatmap(rio.open('out/'+InputFolder_1+'albedo_bb_planar_sw.tif').read(1),'albedo_bb_planar_sw')

var_list = ('albedo_bb_planar_sw','albedo_bb_spherical_sw')
for i in range(len(var_list)):
    var_1 = rio.open('out/'+InputFolder_1+var_list[i]+'.tif').read(1)
    var_2 = rio.open('out/'+InputFolder_2+var_list[i]+'.tif').read(1)
    diff = var_1 - var_2
    plt.figure(figsize=(10,15))
    bl.heatmap(diff,InputFolder_1 + '  -  ' + InputFolder_2)
    plt.title(var_list[i])
    plt.savefig('Plots/'+var_list[i]+'_diff.png',bbox_inches='tight')
    
for i in np.append(np.arange(11), np.arange(15,21)):
    var_name = 'albedo_spectral_planar_'+ str(i+1).zfill(2)
    var_1 = rio.open('out/'+InputFolder_1+var_name+'.tif').read(1)
    var_2 = rio.open('out/'+InputFolder_2+var_name+'.tif').read(1)
    diff = var_1 - var_2
    plt.figure(figsize=(10,15))
    bl.heatmap(diff,InputFolder_1 + '  -  ' + InputFolder_2)
    plt.title(var_name)
    plt.savefig('Plots/'+var_name+'_diff.png',bbox_inches='tight')

#%% Comparing with Fortran output

# bl.heatmap(rio.open('out/'+InputFolder+'diagnostic_retrieval.tif').read(1),'isnow_py')
# bl.heatmap(rio.open('out/'+InputFolder_2+'isnow.tif').read(1)-rio.open('out/'+InputFolder_1+'isnow.tif').read(1),'isnow_f-isnow_py')


#%% fixing diagnostic plots
path = 'C:/Users/bav/OneDrive - Geological survey of Denmark and Greenland/Data/Cook_data/Cook_pySICEv1.4_output'
ipython.magic("matplotlib inline")

import glob
path_files = glob.glob(path+'/*/diagnostic_retrieval.tif')
for i in range(len(path_files)):
    var_1 = rio.open(path_files[i]).read(1)

    fig,ax = bl.heatmap_discrete(var_1,cmap_in='tab20')
    ax.set_title('diagnostic_retrieval   '+path_files[i][105:len(path_files[i])-25])
    plt.savefig(path_files[i][0:len(path_files[i])-25]+'/plots/diagnostic_retrieval.png',bbox_inches='tight')
    plt.close()