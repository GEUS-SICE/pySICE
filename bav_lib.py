# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:03:28 2019

@author: bav
"""

#%matplotlib qt  

import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable  
import pandas as pd

def OutlookRaster(var,title):
    l,b,r,t = var.bounds
    res = var.res
    x = np.arange(l,r, res[0])
    y = np.arange(t,b, -res[0])
    z=var.read(1)
    nan_col = ~np.all(np.isnan(z), axis=0)
    z=z[:,nan_col]
    x=x[nan_col]
    
    nan_row = ~np.all(np.isnan(z), axis=1)
    z=z[nan_row,:]
    y=y[nan_row]
     
    fig = go.Figure(
            data=go.Heatmap(x=x,
                            y=y,
                            z=z,
                            type = 'heatmap',
                            colorscale='Jet',
                            colorbar=dict(title=title),
                            showscale=True))
    fig.update_layout(
        autosize=False,
        width=500,
        height=500)    
    fig.show()
#    fig.write_image(title+".jpeg")
    return x,y,z

#%%
def heatmap(var, title='', col_lim=(np.nan, np.nan) ,cmap_in='gnuplot'):
    if np.isnan(col_lim[0]):
        col_lim=(np.nanmin(var), np.nanmax(var))
    z=var
    nan_col = ~np.all(np.isnan(z), axis=0)
    z=z[:,nan_col]
    
    nan_row = ~np.all(np.isnan(z), axis=1)
    z=z[nan_row,:]

    cmap = plt.get_cmap(cmap_in)
    cmap.set_bad(color='gray')

#    fig,ax = plt.subplots()
    im = plt.imshow(z, cmap=cmap, interpolation='nearest',vmin=col_lim[0], vmax=col_lim[1])
    cb= plt.colorbar(im)
    cb.ax.set_ylabel(title)
#    fig.show()
#    fig.write_image(title+".jpeg")
    return z
#%%
def heatmap_discrete(var, title='', col_lim=(np.nan, np.nan) ,cmap_in='tab20_r'):
    if np.isnan(col_lim[0]):
        col_lim=(np.nanmin(var), np.nanmax(var))
    z=var
    nan_col = ~np.all(np.isnan(z), axis=0)
    z=z[:,nan_col]
    
    nan_row = ~np.all(np.isnan(z), axis=1)
    z=z[nan_row,:]

    cmap = plt.get_cmap(cmap_in)
    cmap.set_bad(color='white')
    
    bounds = np.unique(var)[np.logical_not(np.isnan(np.unique(var)))]
    bounds = np.append(bounds, bounds[len(bounds)-1]+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N+1)

    fig,ax = plt.subplots(figsize=(10,15))
    im = ax.imshow(z+1e-6, cmap=cmap, norm = norm, interpolation='Nearest')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%" , pad= 0.1)
    cb = mpl.colorbar.ColorbarBase(ax=cax, cmap=cmap, norm=norm)
    
    tic_locs =bounds[0:len(bounds)-1] - (bounds[0:len(bounds)-1]-bounds[1:len(bounds)])/2
    cb.set_ticks(tic_locs)
    cb.ax.set_yticklabels(bounds[0:len(bounds)-1])
    cb.ax.set_ylabel(title)

    fig.show()
#    fig.write_image(title+".jpeg")
    return fig,ax
    
#%% Plot OLCI bands
def plot_OLCI_bands(ax):
# pots olci bands in the background

    olci_bands = pd.read_excel (r'misc/olci_bands.xlsx',header=0,thousands=' ')
    #Out[409]: Index(['Band', 'λ centre (nm)', 'Width (nm)', 'Function'], dtype='object')
    
    ax.set_xlabel('Wavelength (nm)')
    #ax.autoscale(enable=True, tight=True)
    ax.bar(olci_bands['λ centre (nm)'], ax.get_ylim()[1], olci_bands['Width (nm)'], alpha=0.5, edgecolor ='darkblue', label='OLCI bands')
    
    ax.set_xlim(350, 1050)
    
    height_text = np.array((2, 2.3, 2.3, 2.3, 2.3, 2.3, 
                                  2.3, 2, 2.4 , 
                                  1.7,2, 1.6, 2, 2.5,
                                  1.25, 1.7, 2, 1.7,2,2, 2))/3 *ax.get_ylim()[1]
    for i in range(1,22):
        ax.annotate(str(i),
                    (olci_bands['λ centre (nm)'][i-1], height_text[i-1]), 
                    ha='center',
                    bbox=dict(boxstyle="circle", fc="w"),
                    fontsize =13)
    return ax


