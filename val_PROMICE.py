# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 08:58:21 2020
Data processing for the evaluation of SICE against PROMICE measurements
@author: bav
"""

import numpy as np
import pandas as pd
import os
import errno
from pathlib import Path
import bav_lib as bl
import matplotlib.pyplot as plt

path_to_file = 'validation/data/S3_PROMICE.csv'

# Read merged csv
data_all = pd.read_csv(path_to_file) 

#%% Running pySICE
cmd = '%run sice_f.py '+path_to_file
# %run sice.py validation/data/data_in.csv
print(cmd)
eval(cmd)
#%% Analysis
#%matplotlib qt
#%matplotlib inline  
data_out = pd.read_csv('validation/data/S3_PROMICE_fortran_out.csv')

# Removing cloudy measurements
cloud_thd = 0.3
data_cloud = data_out[data_out['cloud']>cloud_thd]
data_out=data_out[data_out['cloud']<=cloud_thd]
print('\nRemoving '+str(data_cloud.shape[0]) +'/'+
      str(data_cloud.shape[0]+data_out.shape[0]) + ' due to clouds')

#%% Plot results ##
fs = 15
ss=5
f1, ax = plt.subplots(1,3,figsize=(15, 7))
f1.subplots_adjust(hspace=0.2, wspace=0.1,
                   left = 0.08 , right = 0.95 ,
                   bottom = 0.2 , top = 0.9)
# ax[0] = bl.density_scatter(data_out['albedo_bb_planar_sw'], data_out['PROMICE_alb'],ax[0],ss)
ax[0].scatter(data_out['albedo_bb_planar_sw'], data_out['PROMICE_alb'], c= "black",s=5)
plt.tight_layout()
ax[0].set_xlabel("Albedo from pySICEv1.4 (-)",fontsize=fs)
ax[0].set_ylabel("PROMICE albedo (-)",fontsize=fs)

ax[1] = bl.density_scatter(data_out['BBA_emp'], data_out['PROMICE_alb'],ax[1],ss)
ax[1].scatter(data_out['BBA_emp'], data_out['PROMICE_alb'],s=ss,c="black")
plt.tight_layout()
ax[1].set_xlabel("Empirical (-)",fontsize=fs)
ax[1].set_ylabel("PROMICE albedo (-)",fontsize=fs)

var= 'rBRR_20'
# ax[2] = bl.density_scatter(data_out[var], data_out['PROMICE_alb'],ax[2],ss)
ax[2].scatter(data_out[var], data_out['PROMICE_alb'],s=ss,c="black")
plt.tight_layout()
ax[2].set_xlabel(var,fontsize=fs)
ax[2].set_ylabel("PROMICE albedo (-)",fontsize=fs)

f1.savefig('validation/figures/scatter.png')

#%%
data_out['date'] = pd.to_datetime(data_out[['year','month','day','hour','minute','second']])
data_out.set_index(['date', 'station'])
data_out.sort_values(by=['station','date'],inplace=True)
data_out.head(5)
data_out.station.unique()

bl.multi_plot(data_out, title = 'Albedo (-)',
              sites =['KAN_L','KAN_M','KAN_U','KPC_L','KPC_U','SCO_L','SCO_U','EGP'],
              OutputFolder = 'validation/figures/', filename_out='PROMICE_comp_1')
bl.multi_plot(data_out, title='Albedo (-)',
              sites =['QAS_L','QAS_U','TAS_L','TAS_A','THU_L','THU_U','UPE_L','UPE_U'],
              OutputFolder = 'validation/figures/', filename_out='PROMICE_comp_2')

