# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 22:29:54 2020

@author: bav
"""

import numpy as np
import matplotlib.pyplot as plt
import rasterio as rio
import time
import sys
import os
import bav_lib as bl

start_time = time.process_time()

InputFolder = sys.argv[1] + "/"

#% Plotting output
try:
    os.mkdir(InputFolder + "plots")
except:
    print("folder exist")

fig, ax = bl.heatmap_discrete(
    rio.open(InputFolder + "diagnostic_retrieval.tif").read(1), "diagnostic_retrieval "
)
ax.set_title(InputFolder)
fig.savefig(InputFolder + "plots/diagnostic_retrieval.png", bbox_inches="tight")

var_list = ("albedo_bb_planar_sw", "albedo_bb_spherical_sw")
for i in range(len(var_list)):
    var_1 = rio.open(InputFolder + var_list[i] + ".tif").read(1)
    plt.figure(figsize=(10, 15))
    bl.heatmap(var_1, var_list[i], col_lim=(0, 1), cmap_in="jet")
    plt.title(InputFolder)
    plt.savefig(InputFolder + "plots/" + var_list[i] + ".png", bbox_inches="tight")
    plt.close()

var_list = ("O3_SICE", "grain_diameter", "snow_specific_area", "al", "conc", "r0")
for i in range(len(var_list)):
    var_1 = rio.open(InputFolder + var_list[i] + ".tif").read(1)
    plt.figure(figsize=(10, 15))
    bl.heatmap(var_1, var_list[i], cmap_in="jet")
    plt.title(InputFolder)
    plt.savefig(InputFolder + "plots/" + var_list[i] + ".png", bbox_inches="tight")
    plt.close()

for i in np.append(np.arange(11), np.arange(15, 21)):
    var_name = "albedo_spectral_spherical_" + str(i + 1).zfill(2)
    var_1 = rio.open(InputFolder + var_name + ".tif").read(1)
    plt.figure(figsize=(10, 15))
    bl.heatmap(var_1, var_name, col_lim=(0, 1), cmap_in="jet")
    plt.title(InputFolder)
    plt.savefig(InputFolder + "plots/" + var_name + ".png", bbox_inches="tight")
    plt.close()
    var_name = "albedo_spectral_planar_" + str(i + 1).zfill(2)
    var_1 = rio.open(InputFolder + var_name + ".tif").read(1)
    plt.figure(figsize=(10, 15))
    bl.heatmap(var_1, var_name, col_lim=(0, 1), cmap_in="jet")
    plt.title(InputFolder)
    plt.savefig(InputFolder + "plots/" + var_name + ".png", bbox_inches="tight")
    plt.close()
    var_name = "rBRR_" + str(i + 1).zfill(2)
    var_1 = rio.open(InputFolder + var_name + ".tif").read(1)
    plt.figure(figsize=(10, 15))
    bl.heatmap(var_1, var_name, col_lim=(0, 1), cmap_in="jet")
    plt.title(InputFolder)
    plt.savefig(InputFolder + "plots/" + var_name + ".png", bbox_inches="tight")
    plt.close()
    var_name = "r_TOA_" + str(i + 1).zfill(2)
    var_1 = rio.open(InputFolder + var_name + ".tif").read(1)
    plt.figure(figsize=(10, 15))
    bl.heatmap(var_1, var_name, col_lim=(0, 1), cmap_in="jet")
    plt.title(InputFolder)
    plt.savefig(InputFolder + "plots/" + var_name + ".png", bbox_inches="tight")
    plt.close()
