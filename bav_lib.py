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
def heatmap_discrete(var, title='', col_lim=(np.nan, np.nan) ,cmap_in='gnuplot'):
    if np.isnan(col_lim[0]):
        col_lim=(np.nanmin(var), np.nanmax(var))
    z=var
    nan_col = ~np.all(np.isnan(z), axis=0)
    z=z[:,nan_col]
    
    nan_row = ~np.all(np.isnan(z), axis=1)
    z=z[nan_row,:]

    cmap = plt.get_cmap(cmap_in)
    cmap.set_bad(color='gray')
    
    bounds = np.unique(var)[np.logical_not(np.isnan(np.unique(var)))]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N+1)
         
    fig,ax = plt.subplots()
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
    return z
    
    #%% =================================================================
def tmp(**kwargs):
    for arg_name in kwargs:
        return arg_name


