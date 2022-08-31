# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:08:14 2019

@author: bav@geus.dk
"""
#%% preamble
# this script will read-in and plot 3-dimensional NetCDF data in python
import os

os.environ["PROJ_LIB"] = r"C:\Users\bav\AppData\Local\Continuum\anaconda3\Library\share"

from netCDF4 import Dataset
import math
import numpy as np
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import pandas as pd

#%% Creating the data matrix
#    File structure:
#    line number ( e.g.,time of measurements)   latitude,   longitude,
#   solar zenith angle,     viewing zenith angle,       solar azimuthal angle,
#   viewing azimuthal angle,        height of underlying surface in m,
#   21 OLCI TOA reflectances

# Read in NetCDF4 file. Assign directory path if necessary.
path = "OLCI scenes/S3A_OL_1_EFR____20170715T135546_20170715T135846_20180508T010441_0179_020_053_1620_LR2_R_NT_002.SEN3/"
filename = "geo_coordinates.nc"
filepath = path + filename
data = Dataset(filepath, mode="r")

lats = pd.DataFrame({"lat": np.reshape(data.variables["latitude"][:], -1)})
lons = pd.DataFrame({"lon": np.reshape(data.variables["longitude"][:], -1)})

rad = pd.DataFrame({"Id": np.arange(1, len(lons) + 1)})
rad = rad.join(lats)
rad = rad.join(lons)

# Loading tie information
filename = "tie_geometries.nc"
filepath = path + filename
data = Dataset(filepath, mode="r")
SZA_tie = data.variables["SZA"][:]
OZA_tie = data.variables["OZA"][:]
SAA_tie = data.variables["SAA"][:]
OAA_tie = data.variables["OAA"][:]

filename = "tie_geo_coordinates.nc"
filepath = path + filename
data = Dataset(filepath, mode="r")
lat_tie = data.variables["latitude"][:]
lon_tie = data.variables["longitude"][:]


for ii in range(1, 21):
    print(ii)
    filename = "Oa" + str(ii).zfill(2) + "_radiance.nc"
    filepath = path + filename
    data = Dataset(filepath, mode="r")
    tmp = pd.DataFrame(
        {
            "Oa"
            + str(ii).zfill(2)
            + "_radiance": np.reshape(
                data.variables["Oa" + str(ii).zfill(2) + "_radiance"][:], -1
            )
        }
    )
    rad = rad.join(tmp)

#%% removing lines with nan
rad_nan = rad[pd.isnull(rad["Oa06_radiance"])]
rad_nonan = rad[pd.notnull(rad["Oa06_radiance"])]

#%% writing to file
import csv

csv.register_dialect(
    "myDialect",
    delimiter="\t",
    quoting=csv.QUOTE_NONE,
    lineterminator="\n",
    skipinitialspace=False,
)

with open("dob.csv", "w") as f:
    writer = csv.writer(f, dialect="myDialect")
    for row in rad:
        writer.writerow(row)

f.close()


# format the colorbar
cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(
    cax, cmap=cmap, norm=normalize, label="Elevation"
)

# save the figure and show it
plt.savefig("asos_station_elevation.png", format="png", dpi=500, transparent=True)
plt.show()


# Interpolate using delaunay triangularization
import scipy

zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method="nearest")
# Plot Data
cs = map.pcolor(xi, yi, zi, vmin=np.min(alt), vmax=np.max(alt), cmap=cm.jet)
cs.set_edgecolor("face")

# Add Grid Lines
map.drawparallels(np.arange(-90.0, 90.0, 15.0), labels=[1, 0, 0, 0], fontsize=5)
map.drawmeridians(np.arange(-180.0, 180.0, 30.0), labels=[0, 0, 0, 1], fontsize=4)

# Add Coastlines, States, and Country Boundaries
map.drawcoastlines()
map.drawstates()
map.drawcountries()

# Add Colorbar
cbar = map.colorbar(cs, location="bottom", pad="10%")
cbar.set_label("K")
cbar.ax.tick_params(labelsize=10)

# Add Title
plt.title("MERRA-2 2-meter air temperature (2010-01)")

# Save figure as PDF
plt.savefig("MERRA2_2m_airTemp_TEST.pdf", format="pdf", dpi=360)
