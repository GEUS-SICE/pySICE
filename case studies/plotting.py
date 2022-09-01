# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import matplotlib
import matplotlib.pyplot as plt

cmap = matplotlib.cm.jet
cmap.set_bad("magenta")
plt.figure()
OLCI_scene.toa.sel(band=0).unstack(dim="xy").transpose("y", "x").plot()
plt.figure()
(OLCI_scene.toa.sel(band=0) - OLCI_scene.toa_cor.sel(band=0)).unstack(
    dim="xy"
).transpose("y", "x").plot()


plt.figure()
snow.isnow.unstack(dim="xy").transpose("y", "x").plot()


plt.figure()
msk.unstack(dim="xy").transpose("y", "x").plot()


# %%
import xarray as xr

plt.close("all")
for var in [
    "grain_diameter",
    "albedo_bb_planar_sw",
    "albedo_bb_spherical_sw",
    "snow_specific_area",
]:
    gd_v2 = xr.open_rasterio("data/2019-06-14/pySICEv2.0/" + var + ".tif").squeeze()
    gd_v16 = xr.open_rasterio("data/2019-06-14/pySICEv1.5/" + var + ".tif").squeeze()

    fig, ax = plt.subplots(1, 3, figsize=(20, 12))
    gd_v2.plot(ax=ax[0], cbar_kwargs={"label": var})
    gd_v16.plot(ax=ax[1], cbar_kwargs={"label": var})
    (gd_v2 - gd_v16).plot(ax=ax[2], cbar_kwargs={"label": "Difference in " + var})
    ax[0].set_title("pySICEv2.0")
    ax[1].set_title("pySICEv1.6")
    ax[2].set_title("v2.0 - v1.6")

    plt.figure()
    plt.plot((gd_v2 - gd_v16).values.flatten(), marker=".", linestyle="None")
    plt.ylabel("Difference in " + var + "pySICEv2.0 - ")
