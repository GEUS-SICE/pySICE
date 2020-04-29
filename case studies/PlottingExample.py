# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:08:14 2019

@author: bav@geus.dk
"""
#%% preamble
from netCDF4 import Dataset
import math
import numpy as np
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap

#%% Read netcdf file
nc_f = ('./out/OLCI_scenes/2017-07-07/'
        'S3A_OL_1_EFR____20170701T113723_20170701T114023_20180505T075514_'
        '0179_019_237_1620_LR2_R_NT_002.SEN3'
        '/geo_coordinates.nc')

data = Dataset(nc_f, mode='r')

# Uncomment 'print data' line to print MERRA2 metadata. This line will print attribute and variable information.
# from the 'variables(dimensions)' list, choose which variable(s) to read in below.
# print data


# Read in 'T2M' 2-meter air temperature variable. Varible names can be printed by uncommenting 'print data' above.
lons = data.variables['longitude'][:]
lats = data.variables['latitude'][:]
alt = data.variables['altitude'][:]


#%% Plotting Data

def polar_stere(lon_w, lon_e, lat_s, lat_n, **kwargs):
    '''Returns a Basemap object (NPS/SPS) focused in a region.

    lon_w, lon_e, lat_s, lat_n -- Graphic limits in geographical coordinates.
                                  W and S directions are negative.
    **kwargs -- Aditional arguments for Basemap object.

    '''
    lon_0 = lon_w + (lon_e - lon_w) / 2.
    ref = lat_s if abs(lat_s) > abs(lat_n) else lat_n
    lat_0 = math.copysign(90., ref)
    proj = 'npstere' if lat_0 > 0 else 'spstere'
    prj = Basemap(projection=proj, lon_0=lon_0, lat_0=lat_0,
                          boundinglat=0, resolution='c')
    #prj = pyproj.Proj(proj='stere', lon_0=lon_0, lat_0=lat_0)
    lons = [lon_w, lon_e, lon_w, lon_e, lon_0, lon_0]
    lats = [lat_s, lat_s, lat_n, lat_n, lat_s, lat_n]
    x, y = prj(lons, lats)
    ll_lon, ll_lat = prj(min(x), min(y), inverse=True)
    ur_lon, ur_lat = prj(max(x), max(y), inverse=True)
    return Basemap(projection='stere', lat_0=lat_0, lon_0=lon_0,
                           llcrnrlon=ll_lon, llcrnrlat=ll_lat,
                           urcrnrlon=ur_lon, urcrnrlat=ur_lat, **kwargs)
    
%matplotlib qt

#%matplotlib inline
fig, ax = plt.subplots(figsize=(12,7))

map = polar_stere(-57, -29, 59, 84, resolution='l')
#m = Basemap(projection='npstere',boundinglat=60, lon_0=300)
map.drawcoastlines()

# format colors for elevation range
alt_min = np.min(alt)
alt_max = np.max(alt)
cmap = plt.get_cmap('gist_earth')
normalize = matplotlib.colors.Normalize(vmin=alt_min, vmax=alt_max)

# plot elevations with different colors using the numpy interpolation mapping tool
# the range [50,200] can be changed to create different colors and ranges
for ii in range(0,lats.shape[0],50):
#    print(ii/lats.shape[0])
    for jj in range(0,lats.shape[1],50):
        if lats[ii,jj] > 59 \
            and lats[ii,jj]<84 \
            and alt[ii, jj] > 1:
#            and lons[ii,jj]< -29 \
#            and lons[ii,jj]>-57 :
            x,y = map(lons[ii,jj],lats[ii,jj])
            color_interp = np.interp(alt[ii,jj],[alt_min,alt_max],[50,200])
            plt.plot(x,y,3,marker='o',color=cmap(int(color_interp)),markersize=1)

# format the colorbar 
cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,norm=normalize,label='Elevation')
plt.title("Scene coverage")
# save the figure and show it
fig.savefig('scene_coverage2.png', format='png', dpi=500,transparent=True)

#%% Transforming xyz into regular grid

# Interpolate using delaunay triangularization 
from scipy.interpolate import griddata
import time

xi = np.reshape(lons,-1)
yi = np.reshape(lats,-1)
zi = np.reshape(alt,-1)
xi,yi = map(xi,yi)

x = np.linspace(np.min(xi),np.max(xi),100)
y =  np.linspace(np.min(yi),np.max(yi),100)
X, Y = np.meshgrid(x,y)

t = time.time()
Ti = griddata((xi, yi), zi, (X, Y), method='nearest')
Ti[Ti==0]= np.nan
elapsed = time.time() - t
print(elapsed)
# Ti = np.ma.getdata(Ti)
# xi = np.ma.getdata(xi)
# yi = np.ma.getdata(yi)

#%% Plotting grid
%matplotlib inline
fig, ax = plt.subplots(figsize=(12,7))

map = polar_stere(-57, -29, 59, 84, resolution='l')
map.drawcoastlines()
cs = map.pcolor(X,Y, Ti , vmin=0.1, vmax=np.max(Ti), cmap=cm.jet)
cs.set_edgecolor('face')

# Add Grid Lines
map.drawparallels(np.arange(-90., 90., 15.), labels=[1,0,0,0], fontsize=5)
map.drawmeridians(np.arange(-180., 180., 30.), labels=[0,0,0,1], fontsize=4)

# Add Colorbar
cbar = map.colorbar(cs, location='bottom', pad="10%")
cbar.set_label('Elevation (m)')
cbar.ax.tick_params(labelsize=10)

# Add Title
plt.title('Elevation from OLCI file')

# Save figure as PDF
plt.savefig('grid test.pdf', format='pdf', dpi=360)


