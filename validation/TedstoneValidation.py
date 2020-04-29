# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 18:43:46 2019

@author: bav
"""

#%% preamble
# this script will read-in and plot 3-dimensional NetCDF data in python
import os
os.environ['PROJ_LIB']=r'C:\Users\bav\AppData\Local\Continuum\anaconda3\Library\share'

from netCDF4 import Dataset
import math
import numpy as np
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import pandas as pd
import glob # for directory content
from osgeo import gdal
from osgeo import osr
from pylab import * # for diff function
from pyproj import Proj # for UTM-latlon conversion
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from osgeo import gdal
import sys

#%% Creating readable geotiffs
 #    File structure:
#    line number ( e.g.,time of measurements)   latitude,   longitude, 
#   solar zenith angle,     viewing zenith angle,       solar azimuthal angle,
#   viewing azimuthal angle,        height of underlying surface in m,
#   21 OLCI TOA reflectances   
 
# Read in NetCDF4 file. Assign directory path if necessary.
folder = ('.\\Tedstone\\UAS\\albedo_and_class')
files = [f for f in glob.glob(folder + "**\\*.nc", recursive=True)]
for f in files:
#f=files[1]    
    filename = f[len(folder):len(f)]

    print(filename)

    data = Dataset(f, mode='r')
    albedo = data.variables['albedo'][:]
    classified = data.variables['classified'][:]
    #lats=data.variables['lat'][:]
    #lons=data.variables['lon'][:]
    x=data.variables['x'][:]
    y=data.variables['y'][:]

    # find the extent and difference 
    xmin=min(x);xmax=max(x)
    ymin=min(y);ymax=max(y)
    mdx=abs(diff(x))
    mdy=abs(diff(y))

    # determine dx and dy from the median of all the non-zero difference values
    dx=median(mdx[where(mdx>0.0)[0]])
    dy=median(mdy[where(mdy>0.0)[0]])
    
    # write as 32-bit GeoTIFF using GDAL
    ny,nx = albedo.shape
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(filename[1:len(filename)-3]+'_alb.tif', nx, ny, 1, gdal.GDT_Float32)
    
    # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    ds.SetGeoTransform( [xmin, dx, 0, ymax, 0, dy ] )
    
    # set the reference info 
    srs = osr.SpatialReference()
    
    # UTM zone 22, North=1
    srs.SetUTM(22,1)
    srs.SetWellKnownGeogCS('WGS84')    
    ds.SetProjection( srs.ExportToWkt() )
    
    # write the data to a single band
    ds.GetRasterBand(1).WriteArray(albedo)
    # close
    ds = None
    
    ds = driver.Create(filename[1:len(filename)-3]+'_clas.tif', nx, ny, 1, gdal.GDT_Float32)
    
    # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    ds.SetGeoTransform( [xmin, dx, 0, ymax, 0, dy ] )
    
    # set the reference info 
    srs = osr.SpatialReference()
    
    # UTM zone 22, North=1
    srs.SetUTM(22,1)
    srs.SetWellKnownGeogCS('WGS84')    
    ds.SetProjection( srs.ExportToWkt() )
    
    # write the data to a single band
    ds.GetRasterBand(1).WriteArray(classified)
    # close
    ds = None
    
# %% Creating readable geotiffs
 #    File structure:
#    line number ( e.g.,time of measurements)   latitude,   longitude, 
#   solar zenith angle,     viewing zenith angle,       solar azimuthal angle,
#   viewing azimuthal angle,        height of underlying surface in m,
#   21 OLCI TOA reflectances   
 
# Read in NetCDF4 file. Assign directory path if necessary.
alb_2_3 = [0.62, 0.53, 0.47, 0.45 0.48, 0.39, 0.47]
alb_sicef = [nan, 0.64, 0.54, 0.56, 0.56, 0.47, 0.49]
    
folder = ('.\\Tedstone\\UAS\\albedo_and_class')
files = [f for f in glob.glob(folder + "**\\*.nc", recursive=True)]
count = 0
for f in files:
#f=files[1]    
    count = count +1
    filename = f[len(folder):len(f)]

    print(filename)
           
    data = Dataset(f, mode='r')
    albedo = data.variables['albedo'][:]
    albedo[albedo<=0]='nan'

#    albedo[albedo == 0] = 'nan'
    classified = data.variables['classified'][:]
    #lats=data.variables['lat'][:]
    #lons=data.variables['lon'][:]
    x=data.variables['x'][:]
    y=data.variables['y'][:]

    # find the extent and difference 
    xmin=min(x);xmax=max(x)
    ymin=min(y);ymax=max(y)
    mdx= np.abs(diff(x))
    mdy= np.abs(diff(y))

    # determine dx and dy from the median of all the non-zero difference values
    dx = np.median(mdx[np.where(mdx>0.0)[0]])
    dy = np.median(mdy[np.where(mdy>0.0)[0]])
    
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    
    xx, yy = np.meshgrid(x,y)
    myProj_UTM22 = Proj("+proj=utm +zone=22, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    myProj_3413 = Proj(" 	+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    lon, lat = myProj_UTM22(x_mean,y_mean, inverse=True)
    x_3413_center, y_3413_center = myProj_3413(lon, lat)
   
#        file_sice = glob.glob('.\\Tedstone v2.3_jeb\\201707' + filename[11:13] + '*', recursive=True)
#    file_sice = (file_sice[0] + '\\albedo_bb_planar_sw.tif')
#    src_ds=gdal.Open(file_sice) 
#    gt=src_ds.GetGeoTransform()
#    rb=src_ds.GetRasterBand(1)
#    gdal.UseExceptions() #so it doesn't print to screen everytime point is outside grid

#import rasterio as rio
#from rasterio.plot import plotting_extent
#import geopandas as gpd
## Rasterstats contains the zonalstatistics function that you will use to extract raster values
#import rasterstats as rs
#import pandas as pd
#import earthpy as et


#with rio.open(file_sice) as SICE_raster:
#    SICE_data = SICE_raster.read(1, masked=True)
#    SICE_meta = SICE_raster.profile
#
#df = pd.DataFrame({'Location': ['CenterPoint'],
#                   'Latitude': [lat],
#                   'Longitude': [lon]})
#gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude))
#gdf.buffer(1000).to_file( filename[1:13] +'_buffer.shp')


#    from scipy.spatial import ConvexHull
#    x_all = xx.flatten()
#    y_all = yy.flatten()
#    a = albedo.flatten()
#    
#    points = np.vstack((x_all[~np.isnan(a)], y_all[~np.isnan(a)]))
#    points = np.transpose(points)
#    hull = ConvexHull(points.__array__())
    
#    Plutting hull and barycenter
#    plt.plot(points[np.arange(0,len(points),math.floor(len(points)/1000)),0],
#         points[np.arange(0,len(points),math.floor(len(points)/1000)),1],
#         'o')
#    for simplex in hull.simplices:
#        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
#    
#    
#    plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=4)
#    plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
#    plt.scatter(x_mean,y_mean,s=300,c='r',zorder=3)
#    
#    plt.show()
    

# Plottin
    cmap = matplotlib.cm.jet
    cmap.set_bad('white',1.)
    
#fig, ax = plt.subplots(figsize=(10, 10))
#ax.imshow(SICE_data,
#          # Here you must set the spatial extent or else the data will not line up with your geopandas layer
#          extent=plotting_extent(SICE_raster),
#          cmap=cmap, 
#           interpolation='nearest',
#           vmin=0, vmax=1)
#gdf.plot(ax=ax,
#       marker='s',
#       markersize=45,
#       color='purple')
#plt.show()

    
    fig, (ax1, ax2) = plt.subplots(1, 2)
#    ax1.imshow(SICE_data,
#              # Here you must set the spatial extent or else the data will not line up with your geopandas layer
#              extent=plotting_extent(SICE_raster),
#              cmap=cmap, 
#               interpolation='nearest',
#               vmin=0, vmax=1)
    im1 = ax1.imshow(albedo, 
               cmap=cmap, 
               interpolation='nearest',
               vmin=0, vmax=1,
               extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)))
    ax1_divider = make_axes_locatable(ax1)
    # add an axes above the main axes.
    cax1 = ax1_divider.append_axes("top", size="7%", pad="10%")
    cb1 = colorbar(im1, cax=cax1, orientation="horizontal")
    # change tick position to top. Tick position defaults to bottom and overlaps
    # the image.
    cax1.xaxis.set_ticks_position("top")
    cax1.set_xlabel("Albedo (-)")
    cax1.xaxis.set_label_position('top') 
    ax1.ticklabel_format(useOffset=False)
    ax1.set_xlabel('x (m East, UTM22)')
    ax1.set_ylabel('y (m North, UTM22)')
    
    ax2 = plt.subplot(1,2,2)
    a=albedo.flatten()
    hist, bin_edges = np.histogram(a[~np.isnan(a)], bins=50) 
    plt.step(bin_edges[0:len(hist)],hist*0.05*0.05,color='black')
    plt.axvline(x=np.mean(a[~np.isnan(a)]),color='red',label='mean')
    plt.axvline(x=np.median(a[~np.isnan(a)]),color='blue',label='median')
    plt.axvline(x=alb_2_3[count],color='magenta',label='SICEv2.3')
    plt.axvline(x=alb_sicef[count],color='green',label='SICE.f')
    
    ax2.legend(loc='upper left')
    plt.title(filename[1:13])
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right') 
    ax2.set_xlabel("Albedo (-)")
    ax2.set_ylabel("Area (m^2)")
    plt.autoscale(enable=True, axis='both', tight=True)
    plt.gcf().subplots_adjust(right=0.0015)
    plt.savefig(filename[1:13] + "_plot.png", dpi=150)
    plt.show()