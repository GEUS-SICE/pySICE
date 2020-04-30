# -*- coding: utf-8 -*-
"""

"""
ly='p' # x for plot window, p for .png
do_gif=0


import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['font.sans-serif'] = ['Georgia']
plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['axes.edgecolor'] = 'w'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 1
plt.rcParams['grid.color'] = "grey"
plt.rcParams["font.size"] = 16

version='20190807'
SICE_version_name='2.0'

#version='20190906'
#SICE_version_name='3.0'

version='20190830'
SICE_version_name='2.1'

version='20190926'
SICE_version_name='2.2'

figpath="/Users/jason/Dropbox/S3/validation_source_data/Cook_data/Figs/"

fn='/Users/jason/Dropbox/S3/ancil/band_centers.txt'
df = pd.read_csv(fn,header=None)
central_wavelengths=df.iloc[:,0]

print(central_wavelengths)
central_wavelengths=central_wavelengths.to_numpy().flatten()

SICE_path="/Users/jason/Dropbox/S3/AK_FORTRAN/"+version+"_SICE/output/"

inv= np.where((central_wavelengths < 1020) & (central_wavelengths >= 900)) 
print(central_wavelengths[inv])

#for i in range(10,40):
# for i in range(10,11):
for i in range(11,12):
    
    if ly == 'p' :
        plt.close("all")

    idx="{:02d}".format(i+1)
    
    print(idx)
    
    figname=idx+"_comp_Cook_OLCI_"+version

    fn=SICE_path+idx+".albedo_spectral_planar.csv"
    print(fn)
    df = pd.read_csv(fn)
    SICE=df.loc[0,'r1':'r21']
    SICE=SICE.to_numpy().flatten()
#    SICE[inv]=np.nan
    
    fn=SICE_path+idx+".BOAR.csv"
    print(fn)
    df = pd.read_csv(fn)
    BOAR=df.loc[0,'r1':'r21']
    BOAR=BOAR.to_numpy().flatten()
    BOAR[inv]=np.nan

    fn=SICE_path+idx+".TOAR.csv"
    df = pd.read_csv(fn)
    TOAR=df.loc[0,'r1':'r21']
    TOAR=TOAR.to_numpy().flatten()
    
    data_fn = "/Users/jason/Dropbox/S3/validation_source_data/Cook_data/stats_lo.csv"
    df = pd.read_csv(data_fn)
    
    fn="/Users/jason/Dropbox/S3/validation_source_data/Cook_data/output_S3_extract/"+idx+".csv"
    df_OLCI=pd.read_csv(fn)
    #OLCI = df2["albedo_spectral_planar_400"]
    #albedo_spectral_planar_400	albedo_spectral_planar_412	albedo_spectral_planar_442	albedo_spectral_planar_490	albedo_spectral_planar_510	albedo_spectral_planar_560	albedo_spectral_planar_620	albedo_spectral_planar_665	albedo_spectral_planar_673	albedo_spectral_planar_681	albedo_spectral_planar_708	albedo_spectral_planar_753	albedo_spectral_planar_761	albedo_spectral_planar_764	albedo_spectral_planar_767	albedo_spectral_planar_778	albedo_spectral_planar_865	albedo_spectral_planar_885	albedo_spectral_planar_900	albedo_spectral_planar_940	albedo_spectral_planar_1020
    v= np.where(df_OLCI["day"] == 22) 
    

#    print(v)
    
    OLCI=df_OLCI.loc[v[0],'albedo_spectral_planar_400':'albedo_spectral_planar_1020']
    OLCI=OLCI.to_numpy().flatten()
    
        # Assign variables
    X = df["lamda"].values
    Y = df["mean"].values
    Y_std = df["std"]  .values  
    
    # Convert to numpy arrays and reshape X vector (to comply with the convetion of the 'linear_model' function)
    # .values takes input from a "series" to an array

    
    #print("X min: {:.2f}".format(X.min()))
    #print("X max: {:.2f}".format(X.max()))
    #print("Y min: {:.2f}".format(Y.min()))
    #print("Y max: {:.2f}".format(Y.max()))
    
    v= np.where( (df["lamda"] < 1500)) 
    
    X=X[v]
    Y=Y[v]
    Y_std=Y_std[v]
    
    th=1 # line thickness
    plt.plot(X, Y, color='navy', linewidth=th, label='obs. J. Cook',linestyle='-')
    plt.plot(X, Y+Y_std*1.96, color='navy', linewidth=th, label='obs. 95% interval',linestyle='--')
    plt.plot(X, Y-Y_std*1.96, color='navy', linewidth=th,linestyle='--')
    
    ss=10
    th=0.5
    # plt.plot(central_wavelengths,OLCI,marker="x",markersize=ss/2.,color='m',label='planar alb, S3Snow v2.3',linewidth=th)
    #plt.plot(central_wavelengths,OLCI,':',color='m')
    plt.plot(central_wavelengths,SICE,marker="+",markersize=ss,color='r', label='planar alb, SICE '+SICE_version_name,linewidth=th)
    #plt.plot(central_wavelengths,SICE,'-',color='r')
    plt.plot(central_wavelengths,TOAR,marker="|",markersize=ss,color='g', label='rTOA',linewidth=th)
    plt.plot(central_wavelengths,TOAR,'-',color='g')
    plt.plot(central_wavelengths,BOAR,marker="|",markersize=ss,color='b', label='rBRR',linewidth=th)
    plt.plot(central_wavelengths,BOAR,'-',color='b')
    
    plt.xlabel('wavelength, nm') ; plt.ylabel('spectral albedo')
    plt.title('SW Greenland bare ice, SICE v'+SICE_version_name+', case '+idx)
    
    plt.legend(loc='upper right')
    
    plt.ylim(-0.035,0.7) 

    if ly == 'x':
        plt.show()
    if ly == 'p':
        DPI=300
        DPI=150
        fignam=figpath+figname
        plt.savefig(fignam+'.png', bbox_inches = 'tight',figsize = (13,8), dpi = DPI)
        
    #  do_gif  plt.savefig(fignam+'.eps')

if do_gif == 1:
    os.system(opt+'/local/bin/convert -delay 10 -loop 0 '+figpath+'*'+version+'*.png '+figpath+'anim'+version+'.gif')
