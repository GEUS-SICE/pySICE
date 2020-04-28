# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:03:28 2019

@author: bav
"""
import numpy as np
import plotly.graph_objects as go

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
    
    #%% =================================================================
def tmp(**kwargs):
    for arg_name in kwargs:
        return arg_name


