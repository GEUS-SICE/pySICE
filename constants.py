# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:58:31 2019
Update 07032019
This script contains the constants needed by the pySICEv1.1 library.
@author: bav@geus.dk
"""

import numpy as np
import xarray as xr


bandcoord = ("band", np.arange(21))

# ICE REFRATIVE INDEX
xa = np.array(
    [
        2.010e-001,
        2.019e-001,
        2.100e-001,
        2.500e-001,
        3.00e-001,
        3.500e-001,
        3.900e-001,
        4.000e-001,
        4.100e-001,
        4.200e-001,
        4.300e-001,
        4.400e-001,
        4.500e-001,
        4.600e-001,
        4.700e-001,
        4.800e-001,
        4.900e-001,
        5.000e-001,
        5.100e-001,
        5.200e-001,
        5.300e-001,
        5.400e-001,
        5.500e-001,
        5.600e-001,
        5.700e-001,
        5.800e-001,
        5.900e-001,
        6.000e-001,
        6.100e-001,
        6.200e-001,
        6.300e-001,
        6.400e-001,
        6.500e-001,
        6.600e-001,
        6.700e-001,
        6.800e-001,
        6.900e-001,
        7.000e-001,
        7.100e-001,
        7.200e-001,
        7.300e-001,
        7.400e-001,
        7.500e-001,
        7.600e-001,
        7.700e-001,
        7.800e-001,
        7.900e-001,
        8.000e-001,
        8.100e-001,
        8.200e-001,
        8.300e-001,
        8.400e-001,
        8.500e-001,
        8.600e-001,
        8.700e-001,
        8.800e-001,
        8.900e-001,
        9.000e-001,
        9.100e-001,
        9.200e-001,
        9.300e-001,
        9.400e-001,
        9.500e-001,
        9.600e-001,
        9.700e-001,
        9.800e-001,
        9.900e-001,
        1.000e000,
        1.010e000,
        1.020e000,
        1.030e000,
        1.040e000,
        1.050e000,
        1.060e000,
        1.070e000,
        1.080e000,
        1.090e000,
        1.100e000,
        1.110e000,
        1.120e000,
        1.130e000,
        1.140e000,
        1.150e000,
        1.160e000,
        1.170e000,
        1.180e000,
        1.190e000,
        1.200e000,
        1.210e000,
        1.220e000,
        1.230e000,
        1.240e000,
        1.250e000,
        1.260e000,
        1.270e000,
        1.280e000,
        1.290e000,
        1.300e000,
        1.310e000,
        1.320e000,
        1.330e000,
        1.340e000,
        1.350e000,
        1.360e000,
        1.370e000,
        1.380e000,
        1.390e000,
        1.400e000,
        1.410e000,
        1.420e000,
        1.430e000,
        1.440e000,
        1.449e000,
        1.460e000,
        1.471e000,
        1.481e000,
        1.493e000,
        1.504e000,
        1.515e000,
        1.527e000,
        1.538e000,
        1.563e000,
        1.587e000,
        1.613e000,
        1.650e000,
        1.680e000,
        1.700e000,
        1.730e000,
        1.760e000,
        1.800e000,
        1.830e000,
        1.840e000,
        1.850e000,
        1.855e000,
        1.860e000,
        1.870e000,
        1.890e000,
        1.905e000,
        1.923e000,
        1.942e000,
        1.961e000,
        1.980e000,
        2.000e000,
        2.020e000,
        2.041e000,
        2.062e000,
        2.083e000,
        2.105e000,
        2.130e000,
        2.150e000,
        2.170e000,
        2.190e000,
        2.220e000,
        2.240e000,
        2.245e000,
        2.250e000,
        2.260e000,
        2.270e000,
        2.290e000,
        2.310e000,
        2.330e000,
        2.350e000,
        2.370e000,
        2.390e000,
        2.410e000,
        2.430e000,
        2.460e000,
        2.500e000,
    ]
)

ya = np.array(
    [
        3.249e-011,
        2.0e-011,
        2.0e-011,
        2.0e-011,
        2.0e-011,
        2.0e-011,
        2.0e-011,
        2.365e-011,
        2.669e-011,
        3.135e-011,
        4.140e-011,
        6.268e-011,
        9.239e-011,
        1.325e-010,
        1.956e-010,
        2.861e-010,
        4.172e-010,
        5.889e-010,
        8.036e-010,
        1.076e-009,
        1.409e-009,
        1.813e-009,
        2.289e-009,
        2.839e-009,
        3.461e-009,
        4.159e-009,
        4.930e-009,
        5.730e-009,
        6.890e-009,
        8.580e-009,
        1.040e-008,
        1.220e-008,
        1.430e-008,
        1.660e-008,
        1.890e-008,
        2.090e-008,
        2.400e-008,
        2.900e-008,
        3.440e-008,
        4.030e-008,
        4.300e-008,
        4.920e-008,
        5.870e-008,
        7.080e-008,
        8.580e-008,
        1.020e-007,
        1.180e-007,
        1.340e-007,
        1.400e-007,
        1.430e-007,
        1.450e-007,
        1.510e-007,
        1.830e-007,
        2.150e-007,
        2.650e-007,
        3.350e-007,
        3.920e-007,
        4.200e-007,
        4.440e-007,
        4.740e-007,
        5.110e-007,
        5.530e-007,
        6.020e-007,
        7.550e-007,
        9.260e-007,
        1.120e-006,
        1.330e-006,
        1.620e-006,
        2.000e-006,
        2.250e-006,
        2.330e-006,
        2.330e-006,
        2.170e-006,
        1.960e-006,
        1.810e-006,
        1.740e-006,
        1.730e-006,
        1.700e-006,
        1.760e-006,
        1.820e-006,
        2.040e-006,
        2.250e-006,
        2.290e-006,
        3.040e-006,
        3.840e-006,
        4.770e-006,
        5.760e-006,
        6.710e-006,
        8.660e-006,
        1.020e-005,
        1.130e-005,
        1.220e-005,
        1.290e-005,
        1.320e-005,
        1.350e-005,
        1.330e-005,
        1.320e-005,
        1.320e-005,
        1.310e-005,
        1.320e-005,
        1.320e-005,
        1.340e-005,
        1.390e-005,
        1.420e-005,
        1.480e-005,
        1.580e-005,
        1.740e-005,
        1.980e-005,
        3.442e-005,
        5.959e-005,
        1.028e-004,
        1.516e-004,
        2.030e-004,
        2.942e-004,
        3.987e-004,
        4.941e-004,
        5.532e-004,
        5.373e-004,
        5.143e-004,
        4.908e-004,
        4.594e-004,
        3.858e-004,
        3.105e-004,
        2.659e-004,
        2.361e-004,
        2.046e-004,
        1.875e-004,
        1.650e-004,
        1.522e-004,
        1.411e-004,
        1.302e-004,
        1.310e-004,
        1.339e-004,
        1.377e-004,
        1.432e-004,
        1.632e-004,
        2.566e-004,
        4.081e-004,
        7.060e-004,
        1.108e-003,
        1.442e-003,
        1.614e-003,
        1.640e-003,
        1.566e-003,
        1.458e-003,
        1.267e-003,
        1.023e-003,
        7.586e-004,
        5.255e-004,
        4.025e-004,
        3.235e-004,
        2.707e-004,
        2.228e-004,
        2.037e-004,
        2.026e-004,
        2.035e-004,
        2.078e-004,
        2.171e-004,
        2.538e-004,
        3.138e-004,
        3.858e-004,
        4.591e-004,
        5.187e-004,
        5.605e-004,
        5.956e-004,
        6.259e-004,
        6.820e-004,
        7.530e-004,
    ]
)


# OLCI channels
wls = xr.DataArray(
    [
        0.4000e00,
        0.4125e00,
        0.4425e00,
        0.4900e00,
        0.5100e00,
        0.5600e00,
        0.6200e00,
        0.6650e00,
        0.6737e00,
        0.6812e00,
        0.7088e00,
        0.7538e00,
        0.7613e00,
        0.7644e00,
        0.7675e00,
        0.7788e00,
        0.8650e00,
        0.8850e00,
        0.9000e00,
        0.9400e00,
        0.1020e01,
    ],
    coords=[bandcoord],
)

# Imaginary part of ice refrative index at OLCI channels
bai = xr.DataArray(
    [
        2.365e-11,
        2.7e-11,
        7.0e-11,
        4.17e-10,
        8.04e-10,
        2.84e-09,
        8.58e-09,
        1.78e-08,
        1.95e-08,
        2.1e-08,
        3.3e-08,
        6.23e-08,
        7.1e-08,
        7.68e-08,
        8.13e-08,
        9.88e-08,
        2.4e-07,
        3.64e-07,
        4.2e-07,
        5.53e-07,
        2.25e-06,
    ],
    coords=[bandcoord],
)
# in fortran code:
# c       imaginary part of ice refractive index at OLCI channels
#     DATA kappa/6.27E-10,5.78E-10,6.49E-10,
#  c 1.08E-9,1.46E-9,3.35E-09,
#  c 8.58E-09,1.78E-08,1.95E-08,2.1E-08,3.3E-08,6.23E-08,7.1E-08,
#  c 7.68E-08,8.13E-08,9.88E-08,2.4E-07,3.64E-07,4.2E-07,5.53e-07,
#  c       2.25E-06/

# BULK ice absorption coefficient at all OLCI channels  (1/nm)
alpha = 4 * np.pi * bai / wls


# ozone vertical optical density  at OLCI channels
# cabsoz = xr.DataArray([1.3782e-4, 3.0488e-4, 1.6457e-3, 8.9359e-3, 1.7505e-2,
#                        4.3471e-2, 4.4871e-2, 2.1016e-2, 1.7162e-2, 1.4663e-2,
#                        7.9830e-3, 3.8797e-3, 2.9338e-3, 2.7992e-3, 2.7297e-3,
#                        3.2560e-3, 8.9569e-4, 5.1888e-4, 6.7158e-4, 3.1278e-4,
#                        1.4088e-5], coords=[bandcoord])

# relative vertical optical density of ozone f at all OLCI channels
# (normalized to that at 620nm)
# f = cabsoz/cabsoz.isel(band=6)

# new constants defined by Alex
# aerosol properties(aot=aot500nm,ANNA=Angström exponent)
aot = 0.07
anna = 1.3
# MOLEC=0,1 - two versions of molecular optical thickness
molec = 0
jend = 100  # number of pixels to be printed (with spectral data)
# patchy snow threshold for the channel R(400nm)
thv0 = 0.5
# threshold for the measured and simulated spectra differences
thv1 = 10.0
# threshold for the grain diameter
thv2 = 0.07
# threshold for the total ozone column retrieved
thv3 = 1000.0


# %% Solar flux
# solar spectrum constants
f0 = 32.38
f1 = -160140.33
f2 = 7959.53
bet = 1.0 / 0.08534
gam = 1.0 / 0.40179


def sol(x):
    # SOLAR SPECTRUM at GROUND level
    # Inputs:
    # x         wave length in micrometer
    # Outputs:
    # sol       solar spectrum in W m-2 micrometer-1 (?)
    #    if (x < 0.4):
    #            x=0.4
    sol1a = f0 * x
    sol1b = -f1 * np.exp(-bet * x) / bet
    sol1c = -f2 * np.exp(-gam * x) / gam
    return sol1a + sol1b + sol1c


# solar flux calculation
# sol1      visible(0.3-0.7micron)
# Update 2022: same for clean and polluted
sol1 = sol(0.7) - sol(0.3)
# sol2      near-infrared (0.7-2.4micron)
# Update 2022: same for clean and polluted
sol2 = sol(2.4) - sol(0.7)
# sol3      shortwave(0.3-2.4 micron)
sol3 = sol1 + sol2

# asol specific band
asol = sol(0.865) - sol(0.7)

# analystical integration of the solar flux
def analyt_func(z1, z2):
    # see BBA_calc_pol
    # compatible with array
    ak1 = (z2 ** 2.0 - z1 ** 2.0) / 2.0
    ak2 = (z2 / bet + 1.0 / bet ** 2) * np.exp(-bet * z2) - (
        z1 / bet + 1.0 / bet ** 2
    ) * np.exp(-bet * z1)
    ak3 = (z2 / gam + 1.0 / gam ** 2) * np.exp(-gam * z2) - (
        z1 / gam + 1.0 / gam ** 2
    ) * np.exp(-gam * z1)

    am1 = (z2 ** 3.0 - z1 ** 3.0) / 3.0
    am2 = (z2 ** 2.0 / bet + 2.0 * z2 / bet ** 2 + 2.0 / bet ** 3) * np.exp(
        -bet * z2
    ) - (z1 ** 2.0 / bet + 2.0 * z1 / bet ** 2 + 2.0 / bet ** 3) * np.exp(-bet * z1)
    am3 = (z2 ** 2.0 / gam + 2.0 * z2 / gam ** 2 + 2.0 / gam ** 3.0) * np.exp(
        -gam * z2
    ) - (z1 ** 2.0 / gam + 2.0 * z1 / gam ** 2 + 2.0 / gam ** 3.0) * np.exp(-gam * z1)

    return (f0 * ak1 - f1 * ak2 - f2 * ak3), (f0 * am1 - f1 * am2 - f2 * am3)


# %% solar constant
# segment 1
coef1, coef2 = analyt_func(0.3, 0.7)
# segment 2
coef3, coef4 = analyt_func(0.7, 0.865)
