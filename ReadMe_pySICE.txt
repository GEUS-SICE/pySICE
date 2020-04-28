# pySICE
# 
# VERSION 2.3
# Oct. 25, 2019

# BAV 10-10-2019 (bav@geus.dk)
# Changes:
#   old variable:           Replaced by:
#   raa                     its formula: 180.-(vaa-saa)
#   AKOEF                   its formula: totadu/404.59
#   nv np and isk           i_channel
#   xxx                     ak1*ak2/r0
#   answer                  alb_sph
#   nsolo                   spher_calc
#   fun                     alb2rtoa
#   deltak                  removed (diagnostic variable)
#   sobthv                  specific case o alb2rtoa
#   funs                    removed (not used)
#   zbrent                  replaced by a python version
#   psi                     specific case of sol
#   wave                    removed (using vectorized w instead)
              
# This code retrieves snow/ice  albedo and related snow products for clean Arctic
# atmosphere. The errors increase with the load of pollutants in air.
# Alexander  KOKHANOVSKY
# a.kokhanovsky@vitrocisetbelgium.com
# **************************************************
               
# Inputs:
# sza                       solar zenith angle
# vza                       viewing zenith angle
# saa                       solar azimuthal angle
# vaa                       viewing azimuthal angle
# height                    height of underlying surface(meters)
# toa[i_channel]            spectral OLCI TOA reflectance at 21 channels (R=pi*I_reflec/cos(SZA)/E_0)
# tozon [i_channel]         spectral ozone vertical optical depth at the fixed ozone concentration 404.59DU ( wavelength, VOD)
# voda[i_channel]           spectral water vapour vertical optical depth at the fixed concentration 3.847e+22 molecules per square sm
# aot                       threshold value on aerosol optical thickness (aot) at 500nm

# Outputs: 
# Ozone retrieval:
# BXXX                      retrieved total ozone from OLCI measurements
# totadu                    ECMWF total column ozone in Dobson Unit
# toa                       ozone-corrected OLCI toa relfectances

# snow characteristics:
# isnow                     0 = clean snow, 1 = polluted snow
# ntype                     pollutant type: 1(soot), 2( dust), 3 and 4 (other or mixture)
# conc                      pollutant concentration is defined as the volumetric concentration 
#                           of pollutants devided by the volumetric concentration of ice grains
# bf                        normalized absorption coefficient of pollutants ay 1000nm ( in inverse mm)
# bm                        Angstroem absorption coefficient of pollutants ( around 1 - for soot, 3-7 for dust)
    
# alb_sph(i),i=1,21)        spherical albedo
# (rp(i),i=1,21)            planar albedo
# (refl(i),i=1,21)          relfectance (boar)

# D                         diamater of grains(mm)
# area                      specific surface area (kg/m/m)
# al                        effective absorption length(mm)
# r0                        reflectance of a semi-infinite non-absorbing snow layer
#
# plane BroadBand Albedo (BBA)
# rp1                       visible(0.3-0.7micron)
# rp2                       near-infrared (0.7-2.4micron)
# rp3                       shortwave(0.3-2.4 micron)shortwave(0.3-2.4 micron)

# spherical BBA
# rs1                       visible(0.3-0.7micron)
# rs2                       near-infrared (0.7-2.4micron)
# rs3                       shortwave(0.3-2.4 micron)shortwave(0.3-2.4 micron)

# Constants required:
# xa, ya                    ice refractive index ya at wavelength xa
# w                         OLCI channels
# bai                       Imaginary part of ice refrative index at OLCI channels

# Functions required:
# alb2rtoa                  calculates TOA reflectance from surface albedo
# salbed                    calculates ratm for albedo correction (?)
# zbrent                    equation solver
# sol                       solar spectrum
# analyt_func               calculation of surface radiance
# quad_func                 calculation of quadratic parameters
# trapzd                    trapezoidal rule for integral calculation
# funp                      snow spectral planar and spherical albedo function

#Removed parts compared to Fortran:
#--------------------------
## if Oa11 < ALR21 the surface is assumed to be dark ice and we move to the next pixel    
#    if(toa[21]<ALR21):
#        iice = 1  
#        # moving to next pixel
#        break
#    else:
#        iice = 0
#---------------------------
#    # calculation of NDSI, NDBI and flags
#    rr1=toa[17]   
#    rr2=toa[21]
#    arr1=toa[1]
#
#    # derivation of indices
#    indexs=0, indexi=0, indexd=0
#
#    andsi=(rr1-rr2)/(rr1+rr2)
#    andbi=(arr1-rr2)/(arr1+rr2)
#
#    if (andsi>0.03 and arr1>0.5): indexs=1
#    if (andbi>0.33):     indexi=1
#    if (andbi>0.66):    indexd=1
#
#	===============================0