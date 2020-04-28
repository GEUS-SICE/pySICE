sice.f


Version 2.1

A. KOKHANOVSKY
30.08.2019- the case with small values of R(1020nm) (below 0.5) is added

---------------------------
The FORTRAN77 code sice.f is aimed at the determination of snow/ice properties using spaceborne observations.
It has been prepared during the ESA SICE Project.

In this version the ozone and water vapour  corection has been performed.

Therefore, new input files tg_vod.dat for ozone and tg_water_vod.dat have appeared
The files contains ozone/water vapour  vertical depth at fixed ozone/water vapour  concentration for all
OLCI channels.





In addition limiting values of toa reflectance at 1020nm ( currently: 0.5) and limiting value of diameter of grains (currently: 0.1mm) are given in the file limits.dat
---------------------------



Please, use UNIX commands

gfortran -o sice.exe sice.f
./sice.exe

These commands will produce output files.
--------------------------------------------------------------------------------------------------------

INPUT files:

1)

** JEB edit *** olci_toa_newformat.dat (formerly olci_toa.dat))

structure:
line number ( e.g.,time of measurements),latitude,longitude,solar zenith angle, viewing zenith angle,solar azimuthal angle,viewing azimuthal angle,height of underlying surface in m,
21 OLCI TOA reflectances, ozone columne, water vapour column ( as provided in OLCI files - the same dimension)
** JEB edit *** 
here is output from s3_band_extract.py ... I have asked Alex 24 Sept if the dimensions look OK...
o3=0.0081 # according to sice.f o3 may be in g/m/m
wv=3.6513


note that ozone.dat is now not used
c               it is given now in the file olci_toa_newformat.dat                     
c                     open(1973,file='ozone.dat')

** JEB edit *** 


                    read(1,*) ns,alat,alon,sza,vza,saa,vaa,height,(toa(iks),iks=1,21),ozone,water
                    

2)
nlines.dat
1
the number of lines to be processed
1 000 000 pixels is processed during 5min

The following input files shall not be changed:
-----------------------------------------------------
3)
ice_index.dat
This files contains the following  3 columns: wavelength ( microns), real part of ice refractive index, imaginary part of ice refractive index

4)
thv.dat
This file contains the AOT(500nm) to be used in the automatic decision, if clean or polluted snow routine is used in the retrievals
Also this AOT is used in atmospheric correction routine

5)
tg_vod.dat
This file contains the spectral ozone vertical optical depth (VOD) at the fixed ozone concentration
404.59DU ( wavelength, VOD)

6)
tg_water_vod.dat
This file contains the spectral water vapour vertical optical depth (VOD) at the fixed concentration
3.847e+22 molecules per square sm ( wavelength, VOD)

The OLCI file gives the total water vapour column in kg/m**2
The transfer coefficient form molecules/sm**2 to kg/m/m is
2.99236e-22

7)
limits.dat
This file contains two numbers
1-limit for TOA reflectance at 1020nm ALR21
2-limit for the diameter of grains ALDI in mm
Default: 0.1 for ALR21
Deafault: 0.1mm for ALDI
If the respective values in the code are below these numbers the pixel is not processed
(bare ice/opean water or clouds)

8)retrieved_O3.dat
This file contains retrieved total ozone from OLCI measurements in DU among other data:
      write(1914,*) ns,alat,alon,BXXX,totadu,deltak,sza,vza,amf
      pixel number, latitude,longitude,retrieved total ozone, ECMWF total ozone from OLCI files, relative difference,sza,vza,airmass factor
      BXXX=999. means no retrieval / polluted snow

------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

OUTPUT files:
 1                open(157,  file= 'spherical_albedo.dat' )
 2                open(257,  file= 'planar_albedo.dat'     )
 3                open(357,  file= 'boar.dat'                    )
 4                open(1001,file= 'size.dat'                     )
 5                open(457,file=   'impurity.dat'              )
 6                open(557,file=   'bba.dat'                     )
 7                open(5570,file=   'bba_alex_reduced.dat')       
 8                open(55701,file=   'notsnow.dat')
 
               The file 'notsnow.dat' lists the lines which are not processed bacause they have clouds (first index=1) or bare ice (second index=1)
               The files *interm* give intermediate results needed only for the advanved users of the code. They do not provide essential output information.

1             write(157,*) ns,ndate(3),alat,alon,(answer(i),i=1,21),isnow

2             write(257,*) ns,ndate(3),alat,alon,(rp(i),i=1,21),isnow

3             write(357,*) ns,ndate(3),alat,alon,(refl(i),i=1,21),isnow

4             write(457,*) ns,ndate(3),alat,alon,ntype,conc,bf,bm,thv,
       c       toa(1),isnow

5             write(1001,*) ns,ndate(3),alat,alon,D,area,al,r0,
      c       andsi,andbi,indexs,indexi,indexd,isnow
 
6             write(557,33) ns,ndate(3),alat,alon,rp3,rp1,rp2,rs3,rs1,
     c      rs2,isnow
7              write(5570,333)ns,ndate(3),rp3,isnow
8             write(55701,*)ns,ndate(3),icloud,iice

   
 ---


The file 'bba_alex_reduced.dat' gives plane snow albedo ( 0.3-2.4 microns)  to be compared with PROMICE data:
pixel number, latitutde, longitude, BBA, polluted snow index ( if isnow=1- then we have the case II - polluted snow)

Meaning of parameters

isnow=0 clean snow

isnow=1 polluted snow

ntype- type of pollutants: 1(soot), 2( dust), 3 and 4 (other ot mixture)

conc-the concentration is defined as the vlolumetric concentration of pollutants devided by
the volumetric concentration of ice grains

bf - the normalized absorption coefficient of pollutants ay 1000nm ( in inverse mm)
bm-the Angstroem absorption coefficient of pollutants ( around 1 - for soot, 3-7 for dust)
thv-assumed minimal TOA  reflectance at 400nm  for a clean snow field
toa(1)-measured OLCI TOA reflectance at this channel (if  toa(1) is smaller than thv then the algorithm for polluted snow  is used)
D-diamater of grains(mm)
area-specific surface area ( kg/m/m)
al-effective absorption length(mm)

andsi-NDSI ( see ATBD)
andbi-NDBI ( see ATBD)
-------------------------------------------------
indexs=1  for snow
indexi=1   for bare ice
indexd=1  for polluted bare ice
icloud=1   for cloud
iice=1       for weakly reflecting ice ( no retrieval )
otherwise: 0
-------------------------------------------------

rp3,rp1,rp2-plane BBA, sequence: shortwave(0.3-2.4 micron), visible(0.3-0.7micron), near-infrared (0.7-2.4micron)
rs3,rs1,rs2-spherical BBA( the same sequence as above)






