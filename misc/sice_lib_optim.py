# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:58:31 2019
Update 07032019

 pySICE library
 contains:
     alb2rtoa                  calculates TOA reflectance from surface albedo
     salbed                    calculates ratm for albedo correction (?)
     zbrent                    equation solver
     quad_func                 calculation of quadratic parameters
     funp                      snow spectral planar and spherical albedo function

 requires:
     constants.py               contains constants needed to run the functions below

 This code retrieves snow / ice  albedo and related snow products for clean Arctic
 atmosphere. The errors increase with the load of pollutants in air.
 Alexander  KOKHANOVSKY
 a.kokhanovsky@vitrocisetbelgium.com
 Translated to python by Baptiste Vandecrux (bav@geus.dk) 

@author: bav@geus.dk
"""

# pySICEv1.3
#
# from FORTRAN VERSION 5
# March 31, 2020
#
# Latest update of python scripts: 29-04-2020 (bav@geus.dk)
# - Fixed a bug in the indexing of the polluted pixels for which the spherical albedo equation could not be solved.  Solved the oultiers visible in bands 12-15 and 19-20 and  expended the BBA calculation to few pixels that fell out of the index.
# -compression of output
# - new backscatter fraction from Alex
# - new format for tg_vod.dat file

# **************************************************
# Inputs:
# toa_cor_o3[i_channel]            spectral OLCI TOA reflectance at 21 channels (R=pi*I_reflec/cos(SZA)/E_0)
#
# Outputs:
# snow characteristics:
# isnow                     0 = clean snow, 1 = polluted snow
# ntype                     pollutant type: 1(soot), 2( dust), 3 and 4 (other or mixture)
# conc                      pollutant concentration is defined as the volumetric concentration
#                           of pollutants devided by the volumetric concentration of ice grains
# bf                        normalized absorption coefficient of pollutants ay 1000nm ( in inverse mm)
# bm                        Angstroem absorption coefficient of pollutants ( around 1 - for soot, 3-7 for dust)
#
# alb_sph(i),i=1,21)        spherical albedo
# (rp(i),i=1,21)            planar albedo
# (refl(i),i=1,21)          relfectance (boar)
#
# D                         diamater of grains(mm)
# area                      specific surface area (kg/m/m)
# al                        effective absorption length(mm)
# r0                        reflectance of a semi-infinite non-absorbing snow layer
#
# plane BroadBand Albedo (BBA)
# rp1                       visible(0.3-0.7micron)
# rp2                       near-infrared (0.7-2.4micron)
# rp3                       shortwave(0.3-2.4 micron)shortwave(0.3-2.4 micron)
#
# spherical BBA
# rs1                       visible(0.3-0.7micron)
# rs2                       near-infrared (0.7-2.4micron)
# rs3                       shortwave(0.3-2.4 micron)shortwave(0.3-2.4 micron)
#
# Constants required:
# xa, ya                    ice refractive index ya at wavelength xa
# w                         OLCI channels
# bai                       Imaginary part of ice refrative index at OLCI channels
#
# Functions required:
# alb2rtoa                  calculates TOA reflectance from surface albedo
# salbed                    calculates ratm for albedo correction (?)
# zbrent                    equation solver
# sol                       solar spectrum
# analyt_func               calculation of surface radiance
# quad_func                 calculation of quadratic parameters
# funp                      snow spectral planar and spherical albedo function

import os
import time
from glob import glob

import xarray as xr
import netCDF4
import numpy as np

try:
    import rasterio as rio
except ImportError:
    rio = None  # make rasterio optional at this stage


os.environ['PYTROLL_CHUNK_SIZE'] = '256'

import numba

from constants_optim import wls, bai, xa, ya, f0, f1, f2, bet, gam, coef1, coef2, coef3, coef4
from constants_optim import sol1_clean, sol2, sol3_clean, sol1_pol, sol3_pol, asol, bandcoord


# %% ================================================
# tozon [i_channel]         spectral ozone vertical optical depth at the fixed ozone concentration 404.59DU ( wavelength, VOD)
# voda[i_channel]           spectral water vapour vertical optical depth at the fixed concentration 3.847e+22 molecules per square sm

# Outputs:
# Ozone retrieval:
# BXXX                      retrieved total ozone from OLCI measurements
# totadu                    ECMWF total column ozone in Dobson Unit
# toa_cor_03                       ozone-corrected OLCI toa relfectances


class SICEProcessor(object):

    def __init__(self, dirname):

        self.dirname = dirname

        if dirname.endswith(".zarr"):
            self._get_size_zarr()
            self.open = self.open_zarr
        elif os.path.exists(os.path.join(dirname, 'Oa01_radiance.nc')):
            self._get_size_satpy()
            self.open = self.open_satpy

        elif os.path.exists(os.path.join(dirname, 'r_TOA_01.tif')):
            self._get_size_tif()

            self.open = self.open_tif

    def _open_tif(self, filename):
        return rio.open(os.path.join(self.dirname, filename))

    def _get_size_tif(self):
        self.meta = self._open_tif('r_TOA_01.tif').meta
        self.original_width = self.meta['width']
        self.original_height = self.meta['height']

    def open_tif(self, x0=0, y0=0, width=None, height=None):

        # if os.path.isfile(self.dirname):
        #     raise NotImplementedError("this part needs to be cleaned")
        #     InputFolder = os.path.dirname(os.path.dirname(dirname_or_filename)) + '/'
        #     print('Text file input')
        #     # data_in = pd.read_csv(sys.argv[1])
        #     data_in = pd.read_csv(sys.argv[1])
        #     self.toa = np.expand_dims(data_in[[c for c in data_in.columns if c.find('reflec') >= 0]].to_numpy().transpose(), axis=2)

        #     self.ozone = np.expand_dims(data_in['total_ozone'], axis=1)
        #     self.water = np.expand_dims(data_in['total_columnar_water_vapour'], axis=1)
        #     self.sza = np.expand_dims(data_in['sza'], axis=1)
        #     self.saa = np.expand_dims(data_in['saa'], axis=1)
        #     self.vza = np.expand_dims(data_in['vza'], axis=1)
        #     self.vaa = np.expand_dims(data_in['vaa'], axis=1)
        #     self.height = np.expand_dims(data_in['altitude'], axis=1)

        #     self.sza[np.isnan(toa[0, :, :])] = np.nan
        #     self.saa[np.isnan(toa[0, :, :])] = np.nan
        #     self.vza[np.isnan(toa[0, :, :])] = np.nan
        #     self.vaa[np.isnan(toa[0, :, :])] = np.nan

        # # %% ========= input tif ===============
        if not os.path.isdir(self.dirname):
            raise Exception("dirname must be a directory")

        def read_tif(filename):
            chunks = None
            chunks = 'auto'
            data = xr.open_rasterio(os.path.join(self.dirname, filename), chunks=chunks).squeeze(dim='band', drop=True)

            if width is not None:
                data = data.isel(x=slice(x0, x0 + width))
            if height is not None:
                data = data.isel(y=slice(y0, y0 + height))

            return data.stack(xy=("x", "y")).compute()

        self.meta['transform'] = rio.transform.Affine(1.0, 0.0, 0.0, 0.0, -1.0, 0.0)  # to improve. This is invalid in some cases
        self.meta.update(compress='DEFLATE')

        self.toa = []
        for i in range(21):
            dat = read_tif(f'r_TOA_{i + 1:02}.tif')
            self.toa.append(dat)
        self.toa = xr.concat(self.toa, dim='band')
        print("toa=", self.toa.coords)

        self.ozone = read_tif('O3.tif')
        # self.water = read_tif('WV.tif')  # save memory, it is not used
        self.sza = read_tif('SZA.tif')
        self.saa = read_tif('SAA.tif')
        self.vza = read_tif('OZA.tif')
        self.vaa = read_tif('OAA.tif')
        self.elevation = read_tif('height.tif').astype(np.float64)

        self.latitude = read_tif('lat.tif')
        self.longitude = read_tif('lon.tif')

        mask = ~np.isnan(self.toa.sel(band=0))
        self.sza = self.sza.where(mask)
        self.saa = self.saa.where(mask)
        self.vza = self.vza.where(mask)
        self.vaa = self.vaa.where(mask)

        t = self.elevation.unstack('xy')
        self.meta['width'] = len(t.x)
        self.meta['height'] = len(t.y)

    def _get_size_satpy(self):
        filename = os.path.join(self.dirname, 'Oa01_radiance.nc')
        rootgrp = netCDF4.Dataset(filename, "r")
        self.original_width = rootgrp.dimensions['columns'].size
        self.original_height = rootgrp.dimensions['rows'].size

    def open_satpy(self, x0=0, y0=0, width=None, height=None):
        import satpy  # this is not good practice but avoid satpy to be a compulsary dependence

        filenames = glob(os.path.join(self.dirname, "*.nc"))

        scene = satpy.Scene(reader="olci_l1b", filenames=filenames)

        variables = {
            'solar_azimuth_angle': 'saa',
            'solar_zenith_angle': 'sza',
            'satellite_azimuth_angle': 'vaa',
            'satellite_zenith_angle': 'vza',
            'total_ozone': 'ozone',
            'altitude': 'elevation',
            'longitude': 'longitude',
            'latitude': 'latitude'
        }
        scene.load(list(variables.keys()))

        islice = {}
        if width is not None:
            islice['x'] = slice(x0, x0 + width)
        if height is not None:
            islice['y'] = slice(y0, y0 + height)

        def get_var(variable):
            # return the variable and remove what needs to be remove
            data = scene[variable].isel(islice).compute().stack(xy=("x", "y"))
            data.attrs = {}  # remove attributes, due to some conflict with tà_zarr being unable to serialize datatime
            if 'crs' in data.coords:
                del data.coords['crs']  # idem. zarr complains
            return data

        for variable in variables:
            setattr(self, variables[variable], get_var(variable))
        scene.unload()  # maybe useless
        coef = 1 / np.cos(np.deg2rad(self.sza)) / 100.

        bands = [f'Oa{i:02}' for i in range(1, 22)]
        scene.load(bands)

        scene.load([satpy.DataQuery(name=band, calibration='reflectance') for band in bands])
        self.toa = []
        for band in bands:
            self.toa.append(np.clip(get_var(band) * coef, 0, 1))
        self.toa = xr.concat(self.toa, dim='band')
        if 'crs' in self.toa.coords:
            del self.toa.coords['crs']  # idem. zarr complains

        scene.unload()  # probably useless

    def _get_size_zarr(self):
        ds = xr.open_zarr(self.dirname)
        self.original_width = len(ds.x)
        self.original_height = len(ds.y)

    def open_zarr(self, x0=0, y0=0, width=None, height=None):

        variables = {
            'solar_azimuth_angle': 'saa',
            'solar_zenith_angle': 'sza',
            'satellite_azimuth_angle': 'vaa',
            'satellite_zenith_angle': 'vza',
            'total_ozone': 'ozone',
            'altitude': 'elevation',
            'longitude': 'longitude',
            'latitude': 'latitude'
        }

        ds = xr.open_zarr(self.dirname)

        islice = {}
        if width is not None:
            islice['x'] = slice(x0, x0 + width)
        if height is not None:
            islice['y'] = slice(y0, y0 + height)

        def get_var(variable):
            # return the variable and remove what needs to be remove
            return ds[variable].isel(islice).stack(xy=("x", "y")).compute()

        for variable in variables:
            setattr(self, variables[variable], get_var(variable))

        bands = [f'Oa{i:02}' for i in range(1, 22)]
        self.toa = []
        for band in bands:
            self.toa.append(get_var(band))
        self.toa = xr.concat(self.toa, dim='band')

    def view_geometry(self):

        # =========== view geometry propeties  ==============
        aot = 1

        self.cos_sza, self.cos_vza, self.ak1, self.ak2, self.inv_cos_za, self.cos_sa = \
            view_geometry(self.vaa, self.saa, self.sza, self.vza, aot, self.elevation)

    def ozone_correction(self, write_ozone=False):

        # %% water and ozone spectral optical density
        # water_vod = genfromtxt('./tg_water_vod.dat', delimiter='   ')
        # self.voda = xr.DataArray(water_vod[0:21, 1], coords=[bandcoord])

        ozone_vod = np.genfromtxt('./tg_vod.dat', delimiter='   ')
        tozon = xr.DataArray(ozone_vod[0:21, 1], coords=[bandcoord])

        # %% =========== ozone scattering  ====================================
        BXXX, self.toa = molecular_absorption(self.ozone, tozon, self.sza, self.vza, self.toa)
        del self.ozone  # don't use anymore

        if write_ozone:
            write_output(BXXX, 'O3_SICE', self.dirname, self.meta)

    def prepare_processing(self):
        # Filtering pixels unsuitable for retrieval

        self.isnow = np.where(self.toa[20] < 0.1, 102, np.nan)
        self.isnow[self.sza > 75] = 100

        mask = np.isnan(self.isnow)
        self.toa = self.toa.where(mask)

        self.vaa = self.vaa.where(mask)
        self.saa = self.saa.where(mask)
        self.sza = self.sza.where(mask)
        self.vza = self.vza.where(mask)
        self.elevation = self.elevation.where(mask)

    def aerosol_properties(self):
        # =========== atmosphere propeties  ==============
        #print('compute aerosol_properties')
        aot = 1
        self.tau, self.p, self.g = aerosol_properties(aot, self.elevation, self.cos_sa)   # use to return , gaer, taumol, tauaer

    def snow_properties(self):
        # =========== snow properties  ====================================
        #print('compute snow_properties')
        self.diameter, self.al, self.r0, self.bal = snow_properties(self.toa, self.ak1, self.ak2)

        # filtering small D
        diameter_thresh = 0.1

        valid = self.diameter >= diameter_thresh

        self.isnow[~valid.values & np.isnan(self.isnow)] = 104

        self.toa = self.toa.where(valid)
        self.al = self.al.where(valid)  # absorption length
        self.r0 = self.r0.where(valid)
        self.bal = self.bal.where(valid)
        self.cos_sza = self.cos_sza.where(valid)
        self.cos_vza = self.cos_vza.where(valid)
        # D[D<D_thresh] = np.nan

    def clean_snow_albedo(self):
        # =========== clean snow  ====================================
        # for that we calculate the theoretical reflectance at band 1 of a surface with:
        # r0 = 1, a (albedo) = 1, ak1 = 1, ak2 = 1
        # t1 and t2 are the backscattering fraction

        print('prepare_coef')
        self.t1, self.t2, self.ratm, self.r = prepare_coef(self.tau, self.g, self.p, self.cos_sza, self.cos_vza, self.inv_cos_za)
        print('alb2rtoa')
        rs_1 = alb2rtoa(1, self.t1.sel(band=0), self.t2.sel(band=0), 1, 1, 1, self.ratm.sel(band=0), self.r.sel(band=0))

        # we then compare it to the observed toa[0] value
        self.ind_clean = self.toa.sel(band=0) >= rs_1
        self.isnow[self.ind_clean] = 0

        self.ind_pol = self.toa.sel(band=0) < rs_1

        # STEP 4a: clean snow retrieval
        # the spherical albedo derivation: alb_sph
        self.alb_sph = np.exp(-np.sqrt(1000. * 4. * np.pi * (bai / wls * self.al)))
        self.alb_sph = np.minimum(self.alb_sph, 1)

    def polluted_snow_albedo(self):

        # =========== polluted snow  ====================================

        if not np.any(self.ind_pol):
            return

        self.isnow[self.ind_pol] = 1
        print("#pol=", int(np.sum(self.ind_pol)))

        #  very dirty snow
        ind_very_dark = (self.toa.sel(band=20) < 0.4) & self.ind_pol
        self.isnow[ind_very_dark] = 6

        def compute_rclean(cos_sza, cos_vza, cos_sa):

            theta = np.rad2deg(np.arccos(cos_sa))

            pz = 11.1 * np.exp(-0.087 * theta) + 1.1 * np.exp(-0.014 * theta)

            sumcos = cos_sza + cos_vza
            rclean = 1.247 + 1.186 * sumcos + 5.157 * cos_sza * cos_vza + pz

            return rclean / 4. / sumcos

        self.r0 = self.r0.where(~ind_very_dark, compute_rclean(self.cos_sza, self.cos_vza, self.cos_sa))

        # approximation of the transcendental equation allowing closed-from solution
        #alb_sph[:,ind_pol] =   (toa_cor_o3[:,ind_pol] - r[:,ind_pol])/(t1[:,ind_pol]*t2[:,ind_pol]*r0[ind_pol] + ratm[:,ind_pol]*(toa_cor_o3[:,ind_pol] - r[:,ind_pol]))

        # solving iteratively the transcendental equation
        iind_pol = dict(xy=np.arange(len(self.ind_pol))[self.ind_pol])
        self.alb_sph[iind_pol] = 1

        def solver_wrapper(toa_cor_o3, tau, t1, t2, r0, ak1, ak2, ratm, r):
            # it is assumed that albedo is in the range 0.1-1.0
            # return sl.zbrent(func_solv, 0.1, 1, args=(t1, t2, r0, ak1, ak2, ratm, r, toa_cor_o3), max_iter=100, tolerance=1.e-6)
            return zbrent(0.1, 1, args=(t1, t2, r0, ak1, ak2, ratm, r, toa_cor_o3), max_iter=30, tolerance=2e-4)

        solver_wrapper_v = np.vectorize(solver_wrapper)

        # loop over all bands except band 19, 20
        for i_channel in np.append(np.arange(18), [20]):
            #print('band=', i_channel)
            self.alb_sph.sel(band=i_channel)[iind_pol] = solver_wrapper_v(
                self.toa.sel(band=i_channel)[iind_pol],
                self.tau.sel(band=i_channel)[iind_pol],
                self.t1.sel(band=i_channel)[iind_pol],
                self.t2.sel(band=i_channel)[iind_pol],
                self.r0[iind_pol],
                self.ak1[iind_pol], self.ak2[iind_pol],
                self.ratm.sel(band=i_channel)[iind_pol],
                self.r.sel(band=i_channel)[iind_pol]
            )
            ind_bad = self.alb_sph.sel(band=i_channel) == -999
            self.isnow[ind_bad] = -i_channel
        self.alb_sph = self.alb_sph.where(self.isnow >= 0)

        # INTERNal CHECK FOR CLEAN PIXELS
        # Are reprocessed as clean
        ind_clear_pol = ((self.alb_sph.sel(band=0) > 0.98) | (self.alb_sph.sel(band=2) > 0.98)) & self.ind_pol
        self.isnow[ind_clear_pol] = 7

        self.alb_sph = self.alb_sph.where(~ind_clear_pol, np.exp(-np.sqrt(4. * 1000. * np.pi * self.al * (bai / wls))))

       # re-defining polluted pixels
        self.ind_pol &= self.isnow != 7
        iind_pol = dict(xy=np.arange(len(self.ind_pol))[self.ind_pol])

        # retrieving snow impurities
        ntype, bf, conc = snow_impurities(self.alb_sph, self.bal)

        # alex   09.06.2019
        # reprocessing of albedo to remove gaseous absorption using linear polynomial approximation in the range 753-778nm.
        # Meaning: alb_sph[12],alb_sph[13] and alb_sph[14] are replaced by a linear  interpolation between alb_sph[11] and alb_sph[15]
        afirn = (self.alb_sph.sel(band=15)[iind_pol] - self.alb_sph.sel(band=11)[iind_pol]) / (wls.sel(band=15) - wls.sel(band=11))
        bfirn = self.alb_sph.sel(band=15)[iind_pol] - afirn * wls.sel(band=15)

        for b in [12, 13, 14]:
            self.alb_sph.sel(band=b)[iind_pol] = bfirn + afirn * wls.sel(band=b)

        # BAV 09-02-2020: 0.5 to 0.35
        # pixels that are clean enough in channels 18 19 20 and 21 are not affected by pollution, the analytical equation can then be used
        ind_ok = (self.toa.sel(band=20) > 0.35) & self.ind_pol
        iind_ok = dict(xy=np.arange(len(ind_ok))[ind_ok])

        for b in range(17, 21):
            self.alb_sph.sel(band=b)[iind_ok] = np.exp(-np.sqrt(4.*1000. * np.pi * self.al[iind_ok] * (bai.sel(band=b) / wls.sel(band=b))))

        # Alex, SEPTEMBER 26, 2019
        # to avoid the influence of gaseous absorption (water vapor) we linearly interpolate in the range 885-1020nm for bare ice cases only (low toa[20])
        # Meaning: alb_sph[18] and alb_sph[19] are replaced by a linear interpolation between alb_sph[17] and alb_sph[20]
        bcoef = (self.alb_sph.sel(band=20)[iind_pol] - self.alb_sph.sel(band=17)[iind_pol]) / (wls.sel(band=20) - wls.sel(band=17))
        acoef = self.alb_sph.sel(band=20)[iind_pol] - bcoef * wls.sel(band=20)

        for b in [18, 19]:
            self.alb_sph.sel(band=b)[iind_pol] = acoef + bcoef * wls.sel(band=b)

    def compute_plane_albedo(self, compute_polluted=True):
        # ========= derivation of plane albedo and reflectance ===========
        self.rp = self.alb_sph ** self.ak1
        self.refl = self.r0 * self.alb_sph ** (self.ak1 * self.ak2 / self.r0)

        ind_all_clean = self.ind_clean | (self.isnow == 7)

        # CalCULATION OF BBA of clean snow

        # old method: integrating equation
        # BBA_v = np.vectorize(sl.BBA_calc_clean)
        # p1,p2,s1,s2 = BBA_v(al[ind_all_clean], ak1[ind_all_clean])
        #
        # visible(0.3-0.7micron)
        # rp1[ind_all_clean]=p1/sol1_clean
        # rs1[ind_all_clean]=s1/sol1_clean
        # near-infrared (0.7-2.4micron)
        # rp2[ind_all_clean]=p2/sol2
        # rs2[ind_all_clean]=s2/sol2
        # shortwave(0.3-2.4 micron)
        # rp3[ind_all_clean]=(p1+p2)/sol3_clean
        # rs3[ind_all_clean]=(s1+s2)/sol3_clean

        # approximation
        # planar albedo
        # rp1 and rp2 not derived anymore

        # rp3[ind_all_clean] = sl.plane_albedo_sw_approx(D[ind_all_clean], cos_sza[ind_all_clean])
        self.rp3 = plane_albedo_sw_approx(self.diameter, self.cos_sza).where(ind_all_clean)
        #     spherical albedo
        # rs1 and rs2 not derived anymore
        self.rs3 = spher_albedo_sw_approx(self.diameter).where(ind_all_clean)


        if compute_polluted:

            # calculation of the BBA for the polluted snow
            iind_pol = dict(xy=np.arange(len(self.ind_pol))[self.ind_pol])

            # rp1[iind_pol], rp2[iind_pol], rp3[iind_pol] = sl.BBA_calc_pol(rp[iind_pol], asol, sol1_pol, sol2, sol3_pol)
            # rs1[iind_pol], rs2[iind_pol], rs3[iind_pol] = sl.BBA_calc_pol(alb_sph[iind_pol], asol, sol1_pol, sol2, sol3_pol)

            _, _, self.rp3[iind_pol] = BBA_calc_pol(self.rp[iind_pol], asol, sol1_pol, sol2, sol3_pol)
            _, _, self.rs3[iind_pol] = BBA_calc_pol(self.alb_sph[iind_pol], asol, sol1_pol, sol2, sol3_pol)

    # def to_csv(self):
    #     # %% Output

    #     print('\nText file output')
    #     # data_in = pd.read_csv(sys.argv[1])
    #     data_out = data_in
    #     data_out['grain_diameter'] = self.diameter
    #     # data_out['snow_specific_area']=area
    #     data_out['al'] = self.al
    #     data_out['r0'] = self.r0
    #     data_out['diagnostic_retrieval'] = self.isnow
    #     data_out['conc'] = self.conc
    #     data_out['albedo_bb_planar_sw'] = self.rp3
    #     data_out['albedo_bb_spherical_sw'] = self.rs3

    #     for i in np.append(np.arange(11), np.arange(15,21)):
    #     # for i in np.arange(21):
    #         data_out['albedo_spectral_spherical_' + str(i + 1).zfill(2)] = self.alb_sph[i,:,:]
    #     for i in np.append(np.arange(11), np.arange(15,21)):
    #         data_out['rBRR_'+str(i+1).zfill(2)] = self.rp[i,:,:]

    #     basename, ext = os.path.splitext(self.filename)
#        data_out.to_csv(basename + '_out.csv')

    def to_geotif(self, extended_output=False, save_spectral=False):
        # ========= input tif ===============
        # write_output(self.D, 'grain_diameter',self.dirname)
        write_output(6 / 0.917 / self.diameter, 'snow_specific_area', self.dirname, self.meta)
        write_output(self.rp3, 'albedo_bb_planar_sw', self.dirname, self.meta)
        write_output(self.rs3, 'albedo_bb_spherical_sw', self.dirname, self.meta)
        write_output(self.longitude, 'longitude', self.dirname, self.meta)
        write_output(self.latitude, 'latitude', self.dirname, self.meta)

        if isinstance(self.isnow, np.ndarray):
            self.isnow = xr.DataArray(self.isnow, coords=self.longitude.coords)
        write_output(self.isnow, 'diagnostic_retrieval', self.dirname, self.meta)

        if extended_output:
            write_output(self.al, 'al', self.dirname, self.meta)
            write_output(self.r0, 'r0', self.dirname, self.meta)
            if hasattr(conc):
                write_output(self.conc, 'conc', self.dirname, self.meta)
            else:
                print("no conc")

        if save_spectral:
            # for i in np.arange(21):
            for b in np.append(np.arange(11), np.arange(15, 21)):
                write_output(self.alb_sph.sel(band=b), f'albedo_spectral_spherical_{b+1:02}', self.dirname, self.meta)
                write_output(self.rp.sel(band=b), f'albedo_spectral_planar_{b+1:02}', self.dirname, self.meta)
                write_output(self.refl.sel(band=b), f'rBRR_{b+1:02}', self.dirname, self.meta)

    def to_zarr(self, append_dim=None):

        ds = xr.Dataset({'snow_specific_area': 6 / 0.917 / self.diameter.unstack(dim='xy'),
                         'albedo_bb_planar_sw': self.rp3.unstack(dim='xy'),
                         'albedo_bb_spherical_sw': self.rs3.unstack(dim='xy'),
                         'longitude': self.longitude.unstack(dim='xy'),
                         'latitude': self.latitude.unstack(dim='xy')})

        if append_dim:
            mode = 'a'
            encodings = None
        else:
            mode = 'w'
            encodings = {v: {"dtype": "float32"} for v in ds.variables}

        output_path = self.dirname
        if output_path.endswith('.zarr'):
            output_path, _ = os.path.splitext(output_path)
        if output_path.endswith('.SEN3'):
            output_path, _ = os.path.splitext(output_path)

        ds.to_zarr(output_path + '.OUT.zarr', mode=mode, append_dim=append_dim, encoding=encodings, consolidated=True)

    def process(self, compute_polluted=True, **kwargs):

        self.open(**kwargs)

        start_time = time.process_time()
        self.view_geometry()
        self.ozone_correction()
        self.prepare_processing()
        self.aerosol_properties()
        self.snow_properties()
        self.clean_snow_albedo()
        if compute_polluted:
            self.polluted_snow_albedo()
        self.compute_plane_albedo(compute_polluted=compute_polluted)
        self.duration = time.process_time() - start_time

    def process_by_chunk(self, chunk_size, compute_polluted=True):

        size = self.original_width * self.original_height
        nchunks = int(max(np.floor(size / chunk_size), 1))

        height = int(np.ceil(self.original_height / nchunks))

        chunks = [dict(x0=0, y0=i * height, width=None, height=height) for i in range(nchunks)]

        # self.area = []
        # self.rp3 = []
        # self.rs3 = []
        # self.longitude = []
        # self.latitude = []
        # self.isnow = []

        for i, chunk in enumerate(chunks):
            print(f"{i} / {nchunks}")
            op = SICEProcessor(self.dirname)
            op.process(compute_polluted=compute_polluted, **chunk)
            op.to_zarr(append_dim='y' if i > 0 else None)
            # self.area.append(6 / 0.917 / op.diameter)
            # self.rp3.append(op.rp3)
            # self.rs3.append(op.rs3)
            # self.longitude.append(op.longitude)
            # self.latitude.append(op.latitude)
            # self.isnow.append(op.isnow)
            del op

        # self.area = xr.concat(self.area)
        # self.rp3 = xr.concat(self.rp3)
        # self.rs3 = xr.concat(self.rs3)
        # self.longitude = xr.concat(self.longitude)
        # self.latitude = xr.concat(self.latitude)
        # self.isnow = xr.concat(self.isnow)


def write_output(var, var_name, in_folder, meta):
    # this functions write tif files based on a model file, here "Oa01"
    # opens a file for writing

    var = var.unstack(dim='xy')
    with rio.open(os.path.join(in_folder, var_name + '.tif'), 'w+', **meta) as dst:
        dst.write(var.astype('float32'), 1)


def molecular_absorption(ozone, tozon, sza, vza, toa):

    eps = 1.55
    # ecmwf ozone from OLCI file (in Kg.m-2) to DOBSON UNITS
    # 1 kg O3 / m2 = 46696.24  DOBSON Unit (DU)
    totadu = 46696.24 * ozone

    inv_cos_za = 1. / np.cos(np.deg2rad(sza)) + 1. / np.cos(np.deg2rad(vza))

    BX = (toa.sel(band=20)**(1 - eps)) * (toa.sel(band=16)**eps) / toa.sel(band=6)
    BXXX = np.log(BX) / 1.11e-4 / inv_cos_za
    BXXX = BXXX.where((BXXX >= 0) & (BXXX <= 500), 999)

    # Correcting TOA reflectance for ozone absorption
    # bav 09-02-2020: now water scattering not accounted for
    # kg/m**2. transfer to mol/cm**2
    #    roznov = 2.99236e-22  # 1 moles Ozone = 47.9982 grams
    # water vapor optical depth
    #    vap = water/roznov
    #    AKOWAT = vap/3.847e+22#    tvoda = np.exp(inv_cos_za*voda*AKOWAT)

    # tvoda = tozon * 0 + 1
    tvoda = 1

    toa_cor_o3 = toa * tvoda * np.exp(inv_cos_za * tozon * totadu / 404.59)

    return BXXX, toa_cor_o3
# %% viewing characteristics and aerosol properties
# sza                       solar zenith angle
# vza                       viewing zenith angle
# saa                       solar azimuthal angle
# vaa                       viewing azimuthal angle
# raa                   Relative azimuth angle
# aot                       threshold value on aerosol optical thickness (aot) at 500nm
# height                    height of underlying surface(meters)


def view_geometry(vaa, saa, sza, vza, aot, height):
    # transfer of OLCI relative azimuthal angle to the definition used in
    # radiative transfer code
    # raa       relative azimuth angle
    # sza       solar zenith angle
    # vza       viewing zenith angle
    # cos_sa       cosine of the scattering angle
    # ak1
    # ak2

    sin_sza = np.sin(np.deg2rad(sza))
    sin_vza = np.sin(np.deg2rad(vza))

    cos_sza = np.cos(np.deg2rad(sza))
    cos_vza = np.cos(np.deg2rad(vza))

    ak1 = 3. * (1. + 2. * cos_sza) / 7.
    ak2 = 3. * (1. + 2. * cos_vza) / 7.

    raa = 180. - (vaa - saa)
    cos_raa = np.cos(np.deg2rad(raa))
    cos_sa = -cos_sza * cos_vza + sin_sza * sin_vza * cos_raa

    inv_cos_za = 1. / cos_sza + 1. / cos_vza

    return cos_sza, cos_vza, ak1, ak2, inv_cos_za, cos_sa
# %%


def aerosol_properties(aot, height, cos_sa):
    # Atmospheric optical thickness
    tauaer = aot * (wls / 0.5)**(-1.3)

    ad = height / 7400.

    ak = np.minimum(1, np.exp(-ad))

    g0 = 0.5263
    g1 = 0.4627
    wave0 = 0.4685
    gaer = g0 + g1 * np.exp(-wls / wave0)
    pr = 0.75 * (1. + cos_sa**2)

    taumol = wls**(-4.05) * ak * 0.00877
    tau = tauaer + taumol

    # aerosol asymmetry parameter
    g = tauaer * gaer / tau

    # HG phase function for aerosol
    pa = (1 - g**2) / (1. - 2. * g * cos_sa + g**2)**1.5

    p = (taumol * pr + tauaer * pa) / tau   # the order is critical to have the right order of dims (band, xy)

    return tau, p, g  # , gaer, taumol, tauaer

# %% snow properties


def snow_properties(toa, ak1, ak2):
    # retrieval of snow properties ( R_0, size of grains from OLCI channels 865[17] and 1020nm[21]
    # assumed not influenced by atmospheric scattering and absorption processes)

    akap2 = 2.25e-6
    alpha2 = 4. * np.pi * akap2 / 1.020
    eps = 1.549559365010611

    # reflectivity of nonabsorbing snow layer
    rr1 = toa.sel(band=16)
    rr2 = toa.sel(band=20)
    r0 = (rr1**eps) * (rr2**(1. - eps))

    # effective absorption length(mm)
    bal = np.log(rr2 / r0) * np.log(rr2 / r0) / alpha2 / (ak1 * ak2 / r0)**2
    al = bal / 1000.

    # effective grain size(mm):diameter
    D = al / 16.36
    # snow specific area ( dimension: m*m/kg)
    #area = 6. / D / 0.917

    return D, al, r0, bal

# %% =================================================


def prepare_coef_orig(tau, g, p, cos_sza, cos_vza, inv_cos_za):  # , gaer, taumol, tauaer):
    #astra = tau*np.nan
    #rms = tau*np.nan
    #t1 = tau*np.nan
    #t2 = tau*np.nan

    # SOBOLEV
    oskar = 4. + 3. * (1. - g) * tau
    b1 = 1. + 1.5 * cos_sza + (1. - 1.5 * cos_sza) * np.exp(-tau / cos_sza)
    b2 = 1. + 1.5 * cos_vza + (1. - 1.5 * cos_vza) * np.exp(-tau / cos_vza)

    wa1 = 1.10363
    wa2 = -6.70122
    wx0 = 2.19777
    wdx = 0.51656
    bex = np.exp((g - wx0) / wdx)
    sssss = (wa1 - wa2) / (1. + bex) + wa2

    sumcos = cos_sza + cos_vza

    astra = (1. - np.exp(- tau * inv_cos_za)) / sumcos / 4.
    rms = 1. - b1 * b2 / oskar + (3. * (1. + g) * cos_sza * cos_vza - 2. * sumcos) * astra
    # backscattering fraction
    # t1[i,:,:] = np.exp(-(1.-g[i,:,:])*tau[i,:,:]/cos_sza/2.)
    # t2[i,:,:] = np.exp(-(1.-g[i,:,:])*tau[i,:,:]/cos_vza/2.)
    arg = - 0.5 * (1. - g) * tau / sssss
    t1 = np.exp(arg / cos_sza)
    t2 = np.exp(arg / cos_vza)

    rss = p * astra
    r = rss + rms

    # SALBED
#    ratm = salbed(tau, g)
    a_s = (.18016, -0.18229, 0.15535, -0.14223)
    bs = (.58331, -0.50662,  -0.09012, 0.0207)
    cs = (0.21475, -0.1, 0.13639, -0.21948)
    als = (0.16775, -0.06969, 0.08093, -0.08903)
    bets = (1.09188, 0.08994, 0.49647, -0.75218)

    g2 = g * g
    g3 = g2 * g

    a_cst = a_s[0] + a_s[1] * g + a_s[2] * g2 + a_s[3] * g3
    b_cst = bs[0] + bs[1] * g + bs[2] * g2 + bs[3] * g3
    c_cst = cs[0] + cs[1] * g + cs[2] * g2 + cs[3] * g3
    al_cst = als[0] + als[1] * g + als[2] * g2 + als[3] * g3
    bet_cst = bets[0] + bets[1] * g + bets[2] * g2 + bets[3] * g3

    ratm = tau * (a_cst * np.exp(-tau / al_cst) + b_cst * np.exp(-tau / bet_cst) + c_cst)
    return t1, t2, ratm, r   # used to also return: rms, astra


def prepare_coef(tau, g, p, cos_sza, cos_vza, inv_cos_za):  # , gaer, taumol, tauaer):
    args = tau, g, p, cos_sza, cos_vza, inv_cos_za  # , gaer, taumol, tauaer

    dims = [a.dims for a in args]
    dims = [('band', 'xy', )] * 3 + [('xy', )] * 3
    return xr.apply_ufunc(prepare_coef_numpy, *args, input_core_dims=dims, output_core_dims=[('band', 'xy')]*4)


#def prepare_coef_numpy(tau, g, p, cos_sza, cos_vza, inv_cos_za):  # , gaer, taumol, tauaer):
def prepare_coef(tau, g, p, cos_sza, cos_vza, inv_cos_za):  # , gaer, taumol, tauaer):

    one_g_tau = (1 - g) * tau

    # SOBOLEV
    b1 = 1. + 1.5 * cos_sza + (1. - 1.5 * cos_sza) * np.exp(-tau / cos_sza)
    b2 = 1. + 1.5 * cos_vza + (1. - 1.5 * cos_vza) * np.exp(-tau / cos_vza)

    sumcos = cos_sza + cos_vza

    astra = (1. - np.exp(- tau * inv_cos_za)) / sumcos / 4.
    oskar = 4. + 3. * one_g_tau
    rms = 1. - b1 * b2 / oskar + (3. * (1. + g) * (cos_sza * cos_vza) - 2. * sumcos) * astra

    del sumcos

    #rss = p * astra
    r = p * astra + rms

    del rms
    del p

    wa1 = 1.10363
    wa2 = -6.70122
    wx0 = 2.19777
    wdx = 0.51656
    bex = np.exp((g - wx0) / wdx)
    #sssss = (wa1 - wa2) / (1. + bex) + wa2

    # backscattering fraction
    arg = - 0.5 * one_g_tau / ((wa1 - wa2) / (1. + bex) + wa2)

    del bex
    del one_g_tau

    t1 = np.exp(arg / cos_sza)
    t2 = np.exp(arg / cos_vza)

    # SALBED
    #    ratm = salbed(tau, g)
    a_s = (.18016, -0.18229, 0.15535, -0.14223)
    bs = (.58331, -0.50662,  -0.09012, 0.0207)
    cs = (0.21475, -0.1, 0.13639, -0.21948)
    als = (0.16775, -0.06969, 0.08093, -0.08903)
    bets = (1.09188, 0.08994, 0.49647, -0.75218)

    a_cst = a_s[0] + a_s[1] * g
    b_cst = bs[0] + bs[1] * g
    c_cst = cs[0] + cs[1] * g
    al_cst = als[0] + als[1] * g
    bet_cst = bets[0] + bets[1] * g

    for n in range(2, 4):
        if n == 2:
            gg = g ** 2
        else:
            gg *= g
        a_cst += a_s[n] * gg
        b_cst += bs[n] * gg
        c_cst += cs[n] * gg
        al_cst += als[n] * gg
        bet_cst += bets[n] + gg
    del gg

    ratm = tau * (a_cst * np.exp(-tau / al_cst) + b_cst * np.exp(-tau / bet_cst) + c_cst)

    return t1, t2, ratm, r   # used to also return: , astra, rms


def prepare_coef_numpy(tau, g, p, cos_sza, cos_vza, inv_cos_za):  # , gaer, taumol, tauaer):

    one_g_tau = (1 - g) * tau

    # SOBOLEV
    b1 = 1. + 1.5 * cos_sza + (1. - 1.5 * cos_sza) * np.exp(-tau / cos_sza)
    b2 = 1. + 1.5 * cos_vza + (1. - 1.5 * cos_vza) * np.exp(-tau / cos_vza)

    sumcos = cos_sza + cos_vza

    astra = (1. - np.exp(- tau * inv_cos_za)) / sumcos / 4.
    oskar = 4. + 3. * one_g_tau
    rms = 1. - b1 * b2 / oskar + (3. * (1. + g) * (cos_sza * cos_vza) - 2. * sumcos) * astra

    del sumcos

    #rss = p * astra
    r = p * astra + rms

    del rms
    del p

    wa1 = 1.10363
    wa2 = -6.70122
    wx0 = 2.19777
    wdx = 0.51656
    bex = np.exp((g - wx0) / wdx)
    #sssss = (wa1 - wa2) / (1. + bex) + wa2

    # backscattering fraction
    arg = - 0.5 * one_g_tau / ((wa1 - wa2) / (1. + bex) + wa2)

    del bex
    del one_g_tau

    t1 = np.exp(arg / cos_sza)
    t2 = np.exp(arg / cos_vza)

    # SALBED
    #    ratm = salbed(tau, g)
    a_s = (.18016, -0.18229, 0.15535, -0.14223)
    bs = (.58331, -0.50662,  -0.09012, 0.0207)
    cs = (0.21475, -0.1, 0.13639, -0.21948)
    als = (0.16775, -0.06969, 0.08093, -0.08903)
    bets = (1.09188, 0.08994, 0.49647, -0.75218)

    a_cst = a_s[0] + a_s[1] * g
    b_cst = bs[0] + bs[1] * g
    c_cst = cs[0] + cs[1] * g
    al_cst = als[0] + als[1] * g
    bet_cst = bets[0] + bets[1] * g

    for n in range(2, 4):
        if n == 2:
            gg = g ** 2
        else:
            gg *= g
        a_cst += a_s[n] * gg
        b_cst += bs[n] * gg
        c_cst += cs[n] * gg
        al_cst += als[n] * gg
        bet_cst += bets[n] + gg
    del gg

    ratm = tau * (a_cst * np.exp(-tau / al_cst) + b_cst * np.exp(-tau / bet_cst) + c_cst)
    return t1, t2, ratm, r   # used to also return: , astra, rms

# %% snow_imputirities


def snow_impurities(alb_sph, bal):
    # analysis of snow impurities
    # ( the concentrations below 0.0001 are not reliable )
    # bf    normalized absorption coefficient of pollutants ay 1000nm ( in inverse mm)
    # bm    Angstroem absorption coefficient of pollutants ( around 1 - for soot, 3-7 for dust)

    ind_nonan = ~np.isnan(alb_sph.sel(band=0)) & ~np.isnan(alb_sph.sel(band=1))

    p1 = np.log(alb_sph.sel(band=0))**2
    p2 = np.log(alb_sph.sel(band=1))**2
    bm = np.log(p1 / p2) / np.log(wls.sel(band=1) / wls.sel(band=0))

    # type of pollutants
    ntype = 1 + (bm > 1.2)     # 1=soot, 2 = dust

    soda = (wls.sel(band=0)**bm).where(bm >= 0.1)
    bf = soda * p1 / bal

    # normalized absorption coefficient of pollutants at the wavelength  1000nm
    k_abs_0 = p1 / bal
    # bal   -effective absorption length in microns

    B_soot = 1.6        # enhancement factors for soot
    B_ice = 0.9       # enhancement factors for ice grains
    alfa_soot = 4. * np.pi * 0.47 / wls.sel(band=0)  # bulk soot absorption coefficient at 1000nm
    k_dust = 0.01       # volumetric absorption coefficient of dust

    # conc[ntype == 1] = B_soot*k_abs_0[ntype == 1] / B_ice / alfa_soot
    # conc[ntype == 2] = B_soot*k_abs_0[ntype == 2] / k_dust

    conc = (ntype == 1) * B_soot * k_abs_0 / B_ice / alfa_soot + \
           (ntype == 2) * B_soot * k_abs_0 / k_dust

    ntype = ntype.where(bm > 0.5, 3)  # type is other or mixture
    ntype = ntype.where(bm < 10., 4)  # type is other or mixture

    return ntype.where(ind_nonan), bf.where(ind_nonan), conc.where(ind_nonan)


# %% ===========================================================================
def alb2rtoa(a, t1, t2, r0, ak1, ak2, ratm, r):
    # Function that calculates the theoretical reflectance from a snow spherical albedo a
    # This function can then be solved to find optimal snow albedo
    # Inputs:
    # a                     Surface albedo
    # r0                    reflectance of a semi-infinite non-absorbing snow layer
    #
    # Outputs:
    # rs                  surface reflectance at specific channel
    surf = t1 * t2 * r0 * a**(ak1 * ak2 / r0) / (1 - a * ratm)
    rs = r + surf
    return rs

# %% ===========================================================================


def salbed(tau, g):
    # WARNING: NOT USED ANYMORE
    # SPHERICAL ALBEDO OF TERRESTRIAL ATMOSPHERE:
    # bav: replaced as by a_s
    # inputs:
    # tau               directional albedo ?
    # g                 asymetry coefficient
    # outputs:
    # salbed            spherical albedo
    a_s = (.18016,  -0.18229,  0.15535,     -0.14223)
    bs = (.58331,  -0.50662,  -0.09012,        0.0207)
    cs = (0.21475,   -0.1,  0.13639,            -0.21948)
    als = (0.16775, -0.06969,  0.08093,     -0.08903)
    bets = (1.09188,  0.08994,  0.49647,   -0.75218)

    a = a_s[0]*g**0 + a_s[1]*g**1 + a_s[2]*g**2 + a_s[3]*g**3
    b = bs[0]*g**0 + bs[1]*g**1 + bs[2]*g**2 + bs[3]*g**3
    c = cs[0]*g**0 + cs[1]*g**1 + cs[2]*g**2 + cs[3]*g**3
    al = als[0]*g**0 + als[1]*g**1 + als[2]*g**2 + als[3]*g**3
    bet = bets[0]*g**0 + bets[1]*g**1 + bets[2]*g**2 + bets[3]*g**3

    salbed = tau*(a*np.exp(-tau/al)+b*np.exp(-tau/bet)+c)
    return salbed


# %% =====================================================================

@numba.jit(nopython=True, cache=True)
def f(albedo, t1, t2, r0, ak1, ak2, ratm, r, toa_cor_o3):
    surf = t1 * t2 * r0 * albedo**(ak1 * ak2 / r0) / (1 - albedo * ratm)
    rs = r + surf
    return toa_cor_o3 - rs  # sl.alb2rtoa(albedo, t1, t2, r0, ak1, ak2, ratm, r)


@numba.jit(nopython=True, cache=True)
def zbrent(x0, x1, args=(), max_iter=100, tolerance=1e-6):
    # Equation solver using Brent's method
    # https://en.wikipedia.org/wiki/Brent%27s_method
    # Brent’s is essentially the Bisection method augmented with Inverse
    # Quadratic Interpolation whenever such a step is safe. At it’s worst case
    # it converges linearly and equal to Bisection, but in general it performs
    # superlinearly; it combines the robustness of Bisection with the speedy
    # convergence and inexpensive computation of Quasi-Newtonian methods.
    # Because of this, you’re likely to find Brent’s as a default root-finding
    # algorithm in popular libraries. For example, MATLAB’s fzero, used to find
    # the root of a nonlinear function, employs a variation of Brent’s.
    # Python script from https://nickcdryan.com/2017/09/13/root-finding-algorithms-in-python-line-search-bisection-secant-newton-raphson-boydens-inverse-quadratic-interpolation-brents/

    fx0 = f(x0, *args)
    fx1 = f(x1, *args)

#    print(str(fx0) + ", " + str(fx1))
    if ((fx0 * fx1) > 0):
        #        print("Root not bracketed "+str(fx0)+", "+str(fx1))
        #        assert ((fx0 * fx1) <= 0), ("-----Root not bracketed"+str(fx0)+", "+str(fx1))
        return -999

    if abs(fx0) < abs(fx1):
        x0, x1 = x1, x0
        fx0, fx1 = fx1, fx0

    x2, fx2 = x0, fx0

    d = x2  # any value is fine. It is not used at the first iteration

    mflag = True
    steps_taken = 0

    while steps_taken < max_iter and abs(x1-x0) > tolerance:
        fx0 = f(x0, *args)
        fx1 = f(x1, *args)
        fx2 = f(x2, *args)

        if fx0 != fx2 and fx1 != fx2:
            L0 = (x0 * fx1 * fx2) / ((fx0 - fx1) * (fx0 - fx2))
            L1 = (x1 * fx0 * fx2) / ((fx1 - fx0) * (fx1 - fx2))
            L2 = (x2 * fx1 * fx0) / ((fx2 - fx0) * (fx2 - fx1))
            new = L0 + L1 + L2

        else:
            new = x1 - ((fx1 * (x1 - x0)) / (fx1 - fx0))

        if ((new < ((3 * x0 + x1) / 4) or new > x1) or
            (mflag == True and (abs(new - x1)) >= (abs(x1 - x2) / 2)) or
            (mflag == False and (abs(new - x1)) >= (abs(x2 - d) / 2)) or
            (mflag == True and (abs(x1 - x2)) < tolerance) or
                (mflag == False and (abs(x2 - d)) < tolerance)):
            new = (x0 + x1) / 2
            mflag = True

        else:
            mflag = False

        fnew = f(new, *args)
        d, x2 = x2, x1

        if (fx0 * fnew) < 0:
            x1 = new
        else:
            x0 = new

        if abs(fx0) < abs(fx1):
            x0, x1 = x1, x0

        steps_taken += 1

    return x1

# %% =====================================================================


def funp(x, al, sph_calc, ak1):
    #     Spectral planar albedo
    # Inputs:
    # x                     input wavelength (should work with any)
    # ak1
    # al                    absorption length
    # sph_calc              sph_calc= 0 for planar =1 for spherical
    #
    # Constants:
    # xa(168),ya(168)       imaginary part (ya) of the refraction index at specified wavelength (xa)
    #
    # Outputs:
    # f1*funcs              ?
    #
    # bav 2020
    # using numpy interpolation

    y = np.interp(x, xa, ya)

    dega = 1000. * al * 4.*np.pi*y/x
    pow = np.sqrt(dega)

    if (pow >= 1.e-6):
        rsd = np.exp(-pow)
    else:
        rsd = 1.

    if (sph_calc == 0):
        rs = rsd**ak1
    elif (sph_calc == 1):
        rs = rsd

    if (x < 0.4):
        x = 0.4
    funcs = f0 + f1 * np.exp(-x * bet) + f2 * np.exp(-x * gam)

    return rs * funcs

# %% Approximation functions for BBA integration


def plane_albedo_sw_approx(D, cos_sza):
    anka = 0.7389 - 0.1783 * cos_sza + 0.0484 * cos_sza**2.
    banka = 0.0853 + 0.0414 * cos_sza - 0.0127 * cos_sza**2.
    canka = 0.1384 + 0.0762 * cos_sza - 0.0268 * cos_sza**2.
    diam1 = 187.89 - 69.2636 * cos_sza + 40.4821 * cos_sza**2.
    diam2 = 2687.25 - 405.09 * cos_sza + 94.5 * cos_sza**2.
    return anka + banka * np.exp(-1000 * D / diam1) + canka * np.exp(-1000 * D / diam2)


def spher_albedo_sw_approx(D):
    anka = 0.6420
    banka = 0.1044
    canka = 0.1773
    diam1 = 158.62
    diam2 = 2448.18
    return anka + banka * np.exp(-1000 * D / diam1) + canka * np.exp(-1000 * D / diam2)

# %%   CalCULATION OF BBA for clean pixels


def BBA_calc_clean(al, ak1):
    # for clean snow
    # plane albedo
    sph_calc = 0  # planar
    # visible(0.3-0.7micron)

    def func_integ(x):
        return funp(x, al, sph_calc, ak1)

    p1 = qsimp(func_integ, 0.3, 0.7)

    # near-infrared (0.7-2.4micron)
#        p2 = trapzd(func_integ,0.7,2.4, 20)
    p2 = qsimp(func_integ, 0.7, 2.4)

    # spherical albedo
    sph_calc = 1  # spherical calculation

    def func_integ(x):
        return funp(x, al, sph_calc, ak1)

    # visible(0.3-0.7micron)
#        s1 = trapzd(func_integ,0.3,0.7, 20)
    s1 = qsimp(func_integ, 0.3, 0.7)
    # near-infrared (0.7-2.4micron)
#        s2 = trapzd(func_integ,0.7,2.4, 20)
    s2 = qsimp(func_integ, 0.7, 2.4)
    # shortwave(0.3-2.4 micron)
    # END of clean snow bba calculation
    return p1, p2, s1, s2

# %% ===============================


@numba.jit(nopython=True)
def qsimp(func, a, b):
    # integrate function between a and b using simpson's method.
    # works as fast as scipy.integrate quad
    eps = 1.e-3
    jmax = 20
    ost = -1.e30
    os = -1.e30
    for j in range(jmax):
        if (j == 0):
            st = 0.5 * (b - a) * (func(a) + func(b))
        else:
            it = 2**(j - 1)
            tnm = it
            delta = (b - a) / tnm
            x = a+0.5 * delta
            sum = 0.
            for jj in range(it):
                sum = sum + func(x)
                x = x + delta
            st = 0.5*(st + (b - a) * sum / tnm)
        s = (4.*st - ost) / 3.
        if (j > 4):
            if (np.abs(s - os) < eps * np.abs(os)):
                return s
            if (s == 0) and (os == 0.):
                return s
        os = s
        ost = st
    print("Max iteration reached")
    return s

# %% Calculation f BBA for polluted snow


def BBA_calc_pol(alb, asol, sol1_pol, sol2, sol3_pol):
    # polluted snow
    # NEW CODE FOR BBA OF BARE ICE
    # alb is either the planar or spherical albedo

    # ANAlYTICal EQUATION FOR THE NOMINATOR
    # integration over 3 segments

    # segment 1
    # QUADRATIC POLYNOMIal for the range 400-709nm
    # input wavelength
    #    alam2=w[0]
    #    alam3=w[5]
    #    alam5=w[10]
    #    alam6=w[11]
    #    alam7=w[16]
    #    alam8=w[20]

    alam2 = 0.4
    alam3 = 0.56
    alam5 = 0.709
    alam6 = 0.753
    alam7 = 0.865
    alam8 = 1.02

    # input reflectances
    r2 = alb[0, :]
    r3 = alb[5, :]
    r5 = alb[10, :]
    r6 = alb[11, :]
    r7 = alb[16, :]
    r8 = alb[20, :]

    sa1, a1, b1, c1 = quad_func(alam2, alam3, alam5, r2, r3, r5)
    ajx1 = a1*sol1_pol
    ajx2 = b1*coef1
    ajx3 = c1*coef2

    aj1 = ajx1 + ajx2 + ajx3
    # segment 2.1
    # QUADRATIC POLYNOMIal for the range 709-865nm
    sa1, a2, b2, c2 = quad_func(alam5, alam6, alam7, r5, r6, r7)
    ajx1 = a2 * asol
    ajx2 = b2 * coef3
    ajx3 = c2 * coef4

    aj2 = ajx1 + ajx2 + ajx3    # segment 2.2
    # exponential approximation for the range 865- 2400 nm
    z1 = 0.865
    z2 = 2.4
    rati = r7 / r8
    alasta = (alam8 - alam7) / np.log(rati)
    an = 1. / alasta
    p = r7 * np.exp(alam7 / alasta)

    aj31 = (1. / an) * (np.exp(-an * z2) - np.exp(-an * z1))
    aj32 = (1. / (bet + an)) * (np.exp(-(bet + an) * z2) - np.exp(-(an + bet) * z1))
    aj33 = (1. / (gam + an)) * (np.exp(-(gam + an) * z2) - np.exp(-(an + gam) * z1))
    aj3 = (-f0 * aj31 - f1 * aj32 - f2 * aj33) * p

    BBA_vis = aj1 / sol1_pol
    BBA_nir = (aj2 + aj3) / sol2  # here segment 2.1 and 2.2 are summed
    BBA_sw = (aj1 + aj2 + aj3) / sol3_pol

    return BBA_vis, BBA_nir, BBA_sw

# %% ==========================================================================


def quad_func(x0, x1, x2, y0, y1, y2):
    # quadratic function used for the polluted snow BBA calculation
    # see BBA_calc_pol
    # compatible with arrays
    d1 = (x0 - x1) * (x0 - x2)
    d2 = (x1 - x0) * (x1 - x2)
    d3 = (x2 - x0) * (x2 - x1)

    a1 = x1 * x2 * y0 / d1 + x0 * x2 * y1 / d2 + x0 * x1 * y2 / d3
    b1 = -(x1 + x2) * y0 / d1 - (x0 + x2) * y1 / d2 - (x0 + x1) * y2 / d3
    c1 = y0 / d1 + y1 / d2 + y2 / d3
    x = x1
    sa = a1 + b1 * x + c1 * x * x
    return sa, a1, b1, c1
