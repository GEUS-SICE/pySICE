# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 13:35:44 2020
%matplotlib qt
@author: bav
"""
#%matplotlib qt
#%% Running dice
InputFolder = 'SICE_2020_py'

exec('%run sice.py out/'+InputFolder)

#%% Plotting output

bl.heatmap(rio.open('out/'+InputFolder_2+'isnow.tif').read(1),'isnow_f')
bl.heatmap(rio.open('out/'+InputFolder_1+'diagnostic_retrieval.tif').read(1),'isnow_py')
bl.heatmap(rio.open('out/'+InputFolder_2+'isnow.tif').read(1)-rio.open('out/'+InputFolder_1+'isnow.tif').read(1),'isnow_f-isnow_py')


#%% ========== Compare two folders ==================
import rasterio as rio
import bav_lib as bl
import matplotlib.pyplot as plt
plt.close('all')
plt.rcParams.update({'font.size': 18})


InputFolder_1 = 'SICE_2020_py1.3/'
InputFolder_2 = 'SICE_2020_py1.4/'

tmp=bl.heatmap_discrete(rio.open('out/'+InputFolder_1+'diagnostic_retrieval.tif').read(1),
                        'isnow_py '+InputFolder_1)
tmp=bl.heatmap_discrete(rio.open('out/'+InputFolder_2+'diagnostic_retrieval.tif').read(1), 
                        'isnow_py '+InputFolder_2)
plt.figure()
tmp=bl.heatmap(rio.open('out/'+InputFolder_1+'albedo_bb_planar_sw.tif').read(1),'albedo_bb_planar_sw')
#%%
plt.close('all')


isnow = rio.open('out/'+InputFolder_1+'diagnostic_retrieval.tif').read(1)
ind_all_clean = np.logical_or(isnow == 0, isnow == 7)

var_list = ('albedo_bb_planar_sw','albedo_bb_spherical_sw')
for i in range(len(var_list)):
    var_1 = rio.open('out/'+InputFolder_1+var_list[i]+'.tif').read(1)
    var_2 = rio.open('out/'+InputFolder_2+var_list[i]+'.tif').read(1)
    diff = var_1 - var_2
    plt.figure(figsize=(10,15))
    bl.heatmap(diff,'sice_f_5.0  -  sice_f_5.1')
#    bl.heatmap(diff,'Approximation - Integration',col_lim=(-0.01, 0.01),cmap_in='seismic')

    plt.title(var_list[i])
    plt.savefig(var_list[i]+'_diff.png',bbox_inches='tight')

    