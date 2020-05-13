# -*- coding: utf-8 -*-
"""
Created on Fri May  1 09:50:56 2020

@author: bav
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import bav_lib as bl
import seaborn as sns
from numpy import genfromtxt
plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 1
plt.rcParams["font.size"] = 16
plt.close("all")


#%%  top and bottom of the atmosphere solar irradiance

fig, ax = plt.subplots(1,1, sharey='row',figsize=(15,7))
fig.subplots_adjust(hspace=0.4, wspace=0.15, left = 0.15, top=0.7)

#  ASTM G173-03 Reference Spectra Derived from SMARTS v. 2.9.2
df_sol = pd.read_excel (r'misc/astmg173.xls',header=0, nrows=None)
#Wvlgth nm	Etr W*m-2*nm-1	Global tilt  W*m-2*nm-1	Direct+circumsolar W*m-2*nm-1
ax.plot(df_sol['Wvlgth nm'],df_sol['Etr W*m-2*nm-1'],color='black')
ax.plot(df_sol['Wvlgth nm'],df_sol['Direct+circumsolar W*m-2*nm-1'],color='red')
ax.set_ylabel('Solar irradiance \n($\mathregular{w m^{-2} nm^{-1}}$)')
ax.legend({'Top of the atmosphere','Bottom of the atmosphere','OLCI bands'},
          loc='upper center',
          bbox_to_anchor=(0.5, 1.4))
ax.set_ylim(0, 3)

ax = bl.plot_OLCI_bands(ax)
fig.savefig('./atbd/solar_irradiance.png',bbox_inches='tight')

#%%  Ozone optical depth
ozone_vod = genfromtxt('./tg_vod.dat', delimiter='   ')

fig, ax = plt.subplots(1,1, sharey='row',figsize=(15,7))
fig.subplots_adjust(hspace=0.4, wspace=0.15, left = 0.15, top=0.7)
ax.plot(ozone_vod[:,0],ozone_vod[:,1],color='black',label=None)
ax.set_ylabel('Effective ozone optical depth (-)')
ax.set_ylim(0, 0.07)
ax = bl.plot_OLCI_bands(ax)
ax.legend(loc='upper right')
fig.savefig('./atbd/ozone_optical_depth.png',bbox_inches='tight')

