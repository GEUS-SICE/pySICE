**Pre-operational Sentinel-3 snow and ice products**

Algorithm Theoretical Basis Document

Version 3.0

May 06, 2020

1. A. Kokhanovsky (1), J. Box (2), B. Vandecrux (2)

(1) VITROCISET Belgium SPRL, Bratustrasse 7, 64293 Darmstadt, Germany

(2) Geological Survey of Denmark and Greenland (GEUS)
 Øster Voldgade 10, 1350 Copenhagen, Denmark

[1.Introduction 3](#_Toc39758561)

[2.Overview 3](#_Toc39758562)

[2.1.Ocean and Land and Colour Instrument 4](#_Toc39758563)

[2.2.Generated Products 4](#_Toc39758564)

[2.3.Summary of assumptions 6](#_Toc39758565)

[3.Snow and ice property retrievals 7](#_Toc39758566)

[3.1.Definitions 7](#_Toc39758567)

[3.1.1.Geometry of the system 7](#_Toc39758568)

[3.1.2.Reflectance, spherical and plane albedos 7](#_Toc39758569)

[3.2.Overview 8](#_Toc39758570)

[3.3.Atmospheric correction 10](#_Toc39758571)

[3.3.1.Correction of the OLCI TOA reflectance for molecular absorption 10](#_Toc39758572)

[3.3.2.Molecular and aerosol scattering of light and effects on the top-of the atmosphere reflectance 11](#_Toc39758573)

[3.4.Retrieval of the surface characteristics 15](#_Toc39758574)

[3.4.1.Clean snow 15](#_Toc39758575)

[3.4.2.Polluted snow and ice 17](#_Toc39758576)

[3.5.Broadband albedo calculation 20](#_Toc39758577)

[3.5.1.General case 20](#_Toc39758578)

[3.5.1.Approximation used for clean snow 21](#_Toc39758579)

[3.5.2.Approximation used for polluted snow and ice 22](#_Toc39758580)

[4.NDSI and NDBI 24](#_Toc39758581)

[5.Appendix: Data tables used in the retrieval 25](#_Toc39758582)

[References 26](#_Toc39758583)

1.
# Introduction

This document is aimed at the description of the theoretical basis of the algorithms to determine properties of snow and ice from Sentinel-3 observations.

The following topics are covered:

- Atmospheric correction
- Snow extent determination/snow mask
- Bare ice extent determination
- Dirty ice extent determination
- Snow and ice albedo retrieval (spectral and broadband)
- Bottom-of-atmosphere snow reflectance retrieval
- Pollution load retrieval
- Size of snow grains determination
- Snow specific surface area determination

The code used for the SICE retrieval and its documentation can be found at [https://github.com/BaptisteVandecrux/pySICE](https://github.com/BaptisteVandecrux/pySICE).

1.
# Overview

Snow is composed of ice crystals in contact with each other and surrounded by air. Snow can include &#39;impurities&#39; such as dust, soot, algae (e.g., Skiles et al., 2018). Here, we refer to impurities as &#39;pollution&#39;. Snow can also contain liquid water. The volume concentration of snow grains is usually around 1/3 with 2/3 of the snow volume occupied by air (Proksch et al., 2016). The concentration of pollutants is often low, that is, below 100 ng/g especially in polar regions (Doherty et al., 2010).

The algorithms described here are dedicated to the retrieval of snow optical properties such as snow spectral and broadband albedo and also snow microstructure (snow specific surface area and effective optical grain size). We propose a snow mask based on the Normalized Difference Snow Index (NDSI) and a technique to retrieve the concentration of pollutants in snow, which is possible only for the cases with relatively heavy (above 1ppmv) pollution load (Warren, 2013).

  1.
## Ocean and Land and Colour Instrument

Ocean and Land and Colour Instrument (OLCI) is a 21 band spectrometer that measures solar radiation reflected by the Earth&#39;s atmosphere and surface with a ground spatial resolution of 300 m (see Table 1). The OLCI swath width is 1270 km. OLCI is installed on both Sentinel-3A and Sentinel-3B satellite platforms operated by the ESA in service to the EU Copernicus Programme. The Sentinel-3 A and B orbit at 802 km altitude, 98.6 orbital inclination and a 10:00 UTC sun-synchronous equatorial crossing time.

_**Table 1. Band characteristicsof the SENTINEL-3 Ocean and Land Colour Instrument (OLCI)
# 1**_

| Band | λ centre (nm) | Width (nm) |
| --- | --- | --- |
| 1 | 400 | 15 |
| 2 | 412.5 | 10 |
| 3 | 442.5 | 10 |
| 4 | 490 | 10 |
| 5 | 510 | 10 |
| 6 | 560 | 10 |
| 7 | 620 | 10 |
| 8 | 665 | 10 |
| 9 | 673.75 | 7.5 |
| 10 | 681.25 | 7.5 |
| 11 | 708.75 | 10 |
| 12 | 753.75 | 7.5 |
| 13 | 761.25 | 2.5 |
| 14 | 764.375 | 3.75 |
| 15 | 767.5 | 2.5 |
| 16 | 778.75 | 15 |
| 17 | 865 | 20 |
| 18 | 885 | 10 |
| 19 | 900 | 10 |
| 20 | 940 | 20 |
| 21 | 1 020 | 40 |

  1.
## Generated Products

We arrived at the list of planned products (Table 2) as being both based on a theory that is considered mature and comprising products of need by the global snow modelling and Earth Observation community.

In addition to the products listed in Table 2, we also provide cloud mask derived using an approach not discussed in this ATBD. Most of retrievals are based on the measurements at 865 and 1020nm_,_ where the influence of atmospheric light scattering and absorption processes on top-of-atmosphere signal as detected on a satellite over polar regions is weak.

_ **Table 2. SICE: Snow and ice products** _

|
 | **Snow product name** | **Related Section** | **Units** | **Expected range** | **Maximum acceptable uncertainty in modelling** | **Optimum uncertainty** |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Snow mask (based on NDSI) |
 | - | 0 - no snow, 1 -snow | \&lt;10% | 5 % |
| 2
3
4 | Spectral spherical snow albedoSpectral planar snow albedoSpectral surface reflectance | Retrieval of the surface characteristics | - | 0 -1.0 | \&lt;10% | 5 % |
| 3 | Broadband snow albedo(planar and spherical) | Broadband albedo calculation | - | 0-0.9 | \&lt;15%\*
 | \&lt;5%\* |
| 4 | Snow Specific Surface Area | Clean snow | m2 kg-1 | 20-200 | \&lt;15% | 5 % |
| 5 | Snow grain size (diameter) | Clean snow | mm | 0.02-0.2mm | \&lt;15% | 5 % |
| 6 | Concentration of pollutants | Pollutant characteristics | ppmv(10-6) | 0.1-10.0 | - | - |
| 7
8 | Normalized difference snow indexNormalized difference bare ice index | NDSI and NDBI | - | - | - | - |
| \* Source: GCOS (WMO, 2011) |

_Table 2 notes_

1. The snow mask is based on NDSI. We do not provide the fractional snow cover. However, we provide the NDSI in the output of the algorithm, which can be used to estimate the snow cover.
2. The retrieval of spectral snow albedo is based on processing of OLCI data.
3. SSA is separate from grain size because it has different accuracy associated with differing field measurement approach.
4. BBA (planar and spherical) is provided for three spectral broad bands.
5. To determine BOAR, the top-of-atmosphere reflectance is corrected for atmospheric aerosol/molecular scattering and absorption effects.
6. In addition to the products listed in in Table 2, the following products are given in output:

- reflectance from a semi-infinite non-absorbing snow layer
- absorption length
- parameters describing the absorption coefficient of pollutants in snow (see below).

1. We also derive the type of pollutants and provide flags for bare clean and polluted ice as discussed below.

  1.
## Summary of assumptions

A number of assumptions are made that produce an algorithm set that have been shown to be robust. The simplified snow layer model represents snow as a:

1. Horizontally homogeneous plane parallel turbid medium;
2. Vertically homogeneous layer;
3. Semi-infinite layer. Therefore, there is no need to account for the reflective properties of underlying surface.
4. Close packed effects are ignored (although ice fraction is roughly 30%).
5. Geometrical optics can be used to derive local optical snow characteristics.
6. Impurities (dust, soot, etc.) are located external to ice grains.
7. The single light scattering angular pattern is spectrally neutral in the spectral range given in Table 1.
8. Only pixels completely covered by snow are considered, i.e., pixels with ice and/or partially snow pixels are ignored.
9. The effects of slopes and snow roughness are not accounted for.

The output is provided if the OLCI reflectance at 1020nm is larger than 0.1 and the derived diameter of grains is larger than 0.1mm. These numbers are given in the configuration file and can be changed, if needed.

1.
# Snow and ice property retrievals

  1.
## Definitions

    1.
### Geometry of the system

![](RackMultipart20200513-4-1d51jl4_html_11ed9397f90fb802.jpg)

_**Figure 1. Definition of the solar zenith angle , azimuth angle , viewing zenith angle and relative azimuth angle . Illustration adapted from Hudson et al. (2006).**_

The angles describing the solar and satellite positions around the point observation are presented in Figure 1. From these we derive the cosine of the solar zenith angle , the cosine of the viewing zenith angle and the scattering angle :

|
 |


 | (3.1.1) |
| --- | --- | --- |

    1.
### Reflectance, spherical and plane albedos

The **top-of-atmosphere reflection function or reflectance** is defined as (Kokhanovsky, 2006):

|
 |
 | (3.1.2) |
| --- | --- | --- |

where, is the intensity of reflected light, is the solar flux at the top-of-atmosphere. Many satellite instruments simultaneously measure both and . Therefore, the reflection function can be easily found using Eq. (A1.1) as measured by satellite sensors. We shall consider only cloud free pixels in this work. The **bottom-of atmosphere reflectance or snow reflectance** is defined by Eq. (A1.1) when applied at the bottom of the atmosphere.

The reflectance depends on atmospheric effects due to molecular and aerosol scattering and absorption of solar radiation. For retrieval of surface optical properties, these effects must be removed. The _ **plane albedo** _ is defined as the integration of bottom-of atmosphere reflectance _R_ across all viewing azimuth and zenith angles:

|
 |
 | (3.1.3) |
| --- | --- | --- |

The _ **spherical albedo** _ is found by integration of _R_ over all incident angles :

|
 |
 | (3.1.4) |
| --- | --- | --- |

  1.
## Overview

After converting top of the atmosphere radiance to reflectance using the [SNAP](http://step.esa.int/main/toolboxes/snap/) Rad2Refl module. The top of the atmosphere reflectances are being affected by ozone and molecular scattering. Retrievals are approached in two ways, depending on a dynamic threshold in band 1. The threshold is derived from the synthetic radiative transfer calculations for the assumed aerosol optical thickness at 550nm.

_Clean snow retrieval approach_

If reflectance in OLCI band 1 is larger than the threshold value, it is assumed that the ground scene is covered by unpolluted snow (the majority of pixels in the terrestrial cryosphere). Then we provide an index for the clean snow in the output and derive snow spectral albedo in the spectral range 0.3-2.4 micrometres using the two-parameter analytical equation as described by Kokhanovsky et al. (2019) under assumption that atmospheric effects can be ignored at OLCI channels 865 and 1020nm (not in other channels). The unknown two parameters are derived from OLCI reflectances at 865 and 1020nm. This simple approach to atmospheric correction has appeared to produce highly accurate snow spectral albedo in the range 0.4-1.02micrometers with deviations from ground measurements below 0.01-0.02. The same is true for the broadband albedo.

_Polluted snow retrieval approach_

The atmospheric correction for the polluted snow case (low value of OLCI reflectance at 400nm) is treated in two ways depending on the OLCI reflectance at 1020nm.

_Case 1_

If OLCI reflectance at channel 21 is above 0.4, the retrievals for polluted snow are based

1. on the OLCI measurements at bands 16 and 21 (respectively 865 nm and 1020nm wavelengths) and extrapolation to the longer wavelengths using the analytical equation for the spectral surface albedo identical to that as used for a clean snow.
2. on the OLCI measurements at the wavelengths below 865nm corrected for gaseous absorption and light scattering by aerosol in the framework of the theory described below (see Eq. (10)). The albedo inside O2 absorption band is derived using the linear interpolation of results for neighbouring channels.

In this dark surface regime, we assume that scattering and absorption of light by surface impurities and atmosphere can be ignored at the wavelengths 865 and 1020nm. Such an assumption is similar to the approach used for the clean snow.

_Case 2_

In the case of dark snow and ice pixels and band 21 under 0.4 TOA reflectance, atmospheric correction of measurements _at all OLCI channels_ must be performed. We can no longer assume that pollutants and other effects have small influence on OLCI reflectance above 865nm. Here, we use the transcendent Eq. (10) to account for gaseous absorption and light scattering by aerosol for all OLCI channels. We have found that such an approach produces good results outside oxygen and water vapor absorption bands. Therefore, the albedo inside O2 and water vapor absorption bands is derived using the linear interpolation of results for neighbouring channels.

  1.
## Atmospheric correction

    1.
### Correction of the OLCI TOA reflectance for molecular absorption

The top-of-the-atmosphere reflectance is corrected for ozone absorption using the ozone transmittance function :

|
 |
 | (3.3.1) |
| --- | --- | --- |

Where the transmittance function is defined as in Rozanov and Rozanov (2010):

|
 |
 | (3.3.2) |
| --- | --- | --- |

Where is the air mass factor and is the ozone vertical optical depth (VOD) at the wavelength defined as:

|
 |
 | (3.3.3) |
| --- | --- | --- |

Here is the ozone absorption cross-section at the height _z_ and the wavelength is the concentration of the ozone molecules at the height _z_. Equation (3.3.3) depends on the vertical profile of and but in the absence of such information, we use the simplification by (Kokhanovsky et al., 2020c) where the ozone optical depth is expressed as the product the total column ozone and an reference ozone absorption cross-section :

|
 |

 | (3.3.4) |
| --- | --- | --- |

The value of was estimated by Kokhanovsky et al. ( n.d.) as follows. A reference vertical profile of pressure, temperature and ozone concentration was extracted from the 2D chemistry-climate model from Sinnhuber et al. (2009) at 75oN and for the month of August. For these reference profiles, was calculated using the parametrization of ozone cross-section from Serdyuchenko et al. (2014) and integrated vertically (Eq. (3.3.3) to calculate the reference optical depth (Figure 2). Eventually, the reference optical depth was normalized by the total ozone of the reference profile, = 405 DU = 8.6728e-3 kg m-2, to derive the reference ozone absorption cross-section:

|
 |
 | (3.3.5) |
| --- | --- | --- |

Eventually, the transmittance can be calculated for each pixel using Equation (3.3.2) and (3.3.4):

|
 |
 | (3.3.6) |
| --- | --- | --- |

Where is the ECMWF total column ozone concentration provided in the OLCI files.

We note that one should account for the instrument spectral response function because the measurements are usually performed not at a single wavelength but in narrow spectral range Therefore, the value of will differ for different instruments even if measured at the same central wavelength.

![](RackMultipart20200513-4-1d51jl4_html_3a2288af237fbabb.png)

_ **Figure 2. Optical depth of the total ozone column as a function of wavelength. The OLCI bands are highlighted in light blue.** _

    1.
### Molecular and aerosol scattering of light and effects on the top-of the atmosphere reflectance

The background atmospheric aerosol in Arctic is usually characterized by the low values of aerosol optical thickness and values of single scattering albedo close to one. Therefore, one can neglect light absorption by aerosol and assume that the atmosphere-underlying surface reflectance (due to molecular and aerosol scattering and reflectance from underlying surface) can be presented in the following way:

|
 |
 | (3.3.7) |
| --- | --- | --- |

where the surface spherical albedo and snow reflectance function are the quantities we want to quantify in this retrieval. But before that three characteristics of the atmosphere need to be quantified: the atmospheric reflectance , the spherical albedo of atmosphere , and the total atmospheric transmittance from the top-of-atmosphere to the surface and back to the satellite .

Before the atmospheric reflectance can be derived, several characteristics of the atmosphere should be described: the molecular and aerosol optical depth, which describe how opaque the atmosphere is at a given wavelength, the phase function and its derivatives: the asymmetry parameter and backscattering fraction. In this section we describe how we derive these characteristics and eventually present the atmospheric reflectance calculation in Section_ **:** __ **Calculation of the atmospheric reflectance, transmittance and spherical albedo** _ **.**

#### Molecular and aerosol optical depth

The atmospheric reflection function depends on the atmospheric optical thickness, which can be presented in the following form:

|
 |
 | (3.3.8) |
| --- | --- | --- |

The molecular optical depth can be approximated as

|
 |
 | (3.3.9) |
| --- | --- | --- |

where , _p_ is the site pressure, and the wavelength is in microns. We calculate the site pressure using the following equation: . Here the height of the underlying surface provided in OLCI files and H = 7.64 km is the scale height.

It follows for the aerosol optical depth:

|
 |
 | (3.3.10) |
| --- | --- | --- |

where . The pair represents the Angström parameters. Currently, we use the fixed values of in our retrievals. Due to low aerosol load in Arctic, this assumption does not lead to the substantial errors.

#### Phase function, asymmetry parameter and backscatter of the atmosphere

The phase function of a media define the light intensity scattered by the media at a given wavelength and towards the scattering angle . The phase function is normalized so that integrating for all gives one.The asymmetry parameter of the media is then defined as the intensity-weighted average cosine of the scattering angle (Hansen and Travis, 1974). It ranges from -1 for completely backscattered light to +1 for entirely forward scattered light and is defined as:

|
 |
 | (3.3.11) |
| --- | --- | --- |

In presence of both molecular scattering and aerosol, the phase function can be presented in the following form:

|
 |
 | (3.3.12) |
| --- | --- | --- |

where

|
 |
 | (3.3.13) |
| --- | --- | --- |

is the molecular scattering phase function and is the aerosol phase function. We shall represent this function as:

|
 |
 | (3.3.14) |
| --- | --- | --- |

Therefore, it follows for the asymmetry parameter:

|
 |
 | (3.3.15) |
| --- | --- | --- |

The parameter varies with the location, time, aerosol, type, etc. We shall assume that it can be approximated by the following equation:

|
 |
With | (3.3.16) |
| --- | --- | --- |

The coefficients in this equation (as derived from multiple year AERONET observations over Greenland).

Another useful quantity is the so-called backscattering fraction, meaning the fraction of light scattered in the backward direction defined as:

|
 |
 | (3.3.17) |
| --- | --- | --- |

Which, using Eq. (3.3.12), translated into:

|
 |

 | (3.3.18) |
| --- | --- | --- |

where =0.5 and (Kokhanovsky et al., 2020c).

It should be pointed that the system of equations given above enables the calculation of underlying snow-atmosphere reflectance as a function of the aerosol optical thickness for a known value of the snow spherical albedo.

#### Calculation of the atmospheric reflectance, transmittance and spherical albedo

The atmospheric reflectance due to coupled aerosol-molecular scattering can be presented within the framework of the Sobolev approximation (Sobolev, 1975) as the sum of the reflectance due to single scattering and the reflectance due to multiple scattering :

|
 |
 | (3.3.19) |
| --- | --- | --- |

Both contributions can be expressed as a function of the atmospheric optical depth and of the atmosphere&#39;s phase function , and its asymmetry parameter .

The single scattering contribution is expressed as:

|
 |
 | (3.3.20) |
| --- | --- | --- |

The multiple light scattering contribution is approximated as

|
 |
 | (3.3.21) |
| --- | --- | --- |

where

|
 |

 | (3.3.22) |
| --- | --- | --- |

The transmission function is approximated as follows:

|
 |
 | (3.3.23) |
| --- | --- | --- |

Where _B_ is the backscattering coefficient (Eq. and

The atmosphere&#39;s spherical albedo is found using the approximation proposed by Kokhanovsky et al. (2007):

|
 |
 | (3.3.24) |
| --- | --- | --- |

The coefficients of polynomial expansions of all coefficients (a, b, c, ) with respect to the value of _g_ are given by Kokhanovsky et al. (2005).

  1.
## Retrieval of the surface characteristics

    1.
### Clean snow

The asymptotic radiative transfer theory relates the reflectance of a medium to its spherical albedo and consequently allows for the determination of spherical albedo using reflectance observations for a given observation geometry. This approach was first developed and verified with airborne measurements of albedo and reflectance over a bright cloud field with the spherical albedo in the range 0.8-0.95 (Kokhanovsky et al., 2007). This relationship was later adapted to snow by Kokhanovsky et al. (2018, 2019a, 2020a) to retrieve the snow spherical albedo from the atmosphere-corrected OLCI reflectance as:

|
 |
 | (3.4.1) |
| --- | --- | --- |

where is the theoretical reflectance of snow in the absence of absorption (Kokhanovsky et al., 2019a, Appendix A), s the cosine of the solar zenith angle, is the cosine of the viewing zenith angle, and is the angular function (Kokhanovsky et al., 2019a) defined as:

|
 |
 | (3.4.2) |
| --- | --- | --- |

where is either or .

The spherical albedo for clean snow can be presented in the following form (Kokhanovsky et al, 2018):

|
 |
 | (3.4.3) |
| --- | --- | --- |

Where _l_ is the_effective absorption length_in snow and is the bulk absorption coefficient of ice calculated from the wavelength and, the imaginary parts of ice refractive index (Warren and Brandt, 2008) reported in Section _ **Appendix: Data tables used in the retrieval** _:

|
 |
 | (3.4.4) |
| --- | --- | --- |

The plane albedo can be derived eventually from spherical albedo. Namely, it follows (Kokhanovsky et al., 2019a):

|
 |
 | (3.4.5) |
| --- | --- | --- |

The effective absorption length _l_ _does not depend_ on the wavelength in the OLCI spectral range as demonstrated by Kokhanovsky et al. (2018). The same is truefor , the reflectance of non-absorbing snow layer. Therefore, we can derive both parameters from measurements of OLCI reflectance at two wavelengths and as in Kokhanovsky et al., (2018) that satisfy the following two criteria: i) the reflectance at these channels must be sensitive to the parameters of interest, and ii) these channels must be least influenced by atmospheric scattering and absorption processes. Consequently, the OLCI channels centred around 865 and 1020 nm are the best candidates for the retrieval.

The and parameters are subsequently defined as follows.

|
 |
 | (3.4.6) |
| --- | --- | --- |

And can then be used in Eq. (3.4.1) and (3.4.3) at the wavelength to determine _l_ :

|
 |
 | (3.4.7) |
| --- | --- | --- |

The derived value of _l_ can be used to determine the snow spherical/plane albedo and also snow reflection function (OLCI bottom of atmosphere reflectance) at any OLCI wavelength using Eqs. (3.1.1)-(3.3.5). The diameter _d_ of ice grains in snow is estimated using the effective absorption length (Kokhanovsky et al., 2019a):

|
 |
 | (3.4.8) |
| --- | --- | --- |

where the parameter depends on the type of snow/shape of grains. We assume that _A=0.06_ in the retrievalsas suggested by Kokhanovsky et al. (2019a). The snow specific surface area is derived as

|
 |
 | (3.4.9) |
| --- | --- | --- |

where g cm-3 is the bulk ice density.

    1.
### Polluted snow and ice

As in the case of clean snow, we assumed that scattering and absorption of light by atmosphere and impurities in snowpack can be ignored at the wavelengths 865 and 1020nm. This makes it possible to derive the parameters , _l, d,_ in the same way as for clean snow. However, for polluted snow, cannot be derived from Eq. (3.4.3) because the spectral reflectance in the visible is influenced not just by ice grains but also by various impurities (soot, dust, algae).

Nevertheless, Eq. (3.3.1) describing ozone absorption, Eq. (3.3.7) describing atmospheric correction and Eq. (3.4.1) that links observed reflectance to surface reflectance and spherical albedo can still be used at all channels except those subject to molecular absorption by oxygen: channels 13-15; and water vapor: channels 19 and 20. For the remaining channels, 1-12, 16-18 and 21, Equations (3.3.1) (3.3.7) and (3.4.1) give the following system:

|
 |
 | (3.4.10) |
| --- | --- | --- |

Where is the measured top-of-atmosphere reflectance function, is atmospheric contribution to the measured signal, is the spherical albedo of the atmosphere, is the bottom-of-atmosphere surface reflectance, is atmospheric transmittance from the top-of-atmosphere to the underlying surface and back to the satellite position, is the transmittance of purely gaseous atmosphere. Given that is measured and that, , , , can be calculated following the approach detailed above, the system presented in Eq. (3.4.10) has therefore only one unknown, , which cannot be presented in closed form. We consequently derive iteratively using Simpson&#39;s rule.

For the channels 13-15, affected by oxygen absorption, and 19-20, affected by water vapor absorption, the spherical albedos are linearly interpolated between the retrieved spherical albedo at channels 12 and 16 in the first case and channels 18 and 21 in the second.

If the underlying surface is not snow, the application of the equation relating the snow albedo to the snow grain size is not justified. In this case the spherical albedo is found for all OLCI channels (except 19 and 20) using Eq. (10). It is assumed that the value of reflectance for a non-absorbing surface can be approximated, as discussed by Kokhanovsky et al. (2019a), by the following expression:

|
 |
 | (3.4.11) |
| --- | --- | --- |

where _A_ = 1.247, _B_ = 1.186, _C_ = 5.157, and

|
 |
 | (3.4.12) |
| --- | --- | --- |

and θ is the scattering angle in degrees. This formulation of is then used when solving Eq. (3.4.10).

#### Pollutant characteristics

The concentration of pollutants in snow is estimated using the approach described below. It is assumed that the spherical albedo, solved in Eq. (3.4.10), can also be expressed as in Kokhanovsky et al. (2018):

|
 |
 | (3.4.13) |
| --- | --- | --- |

where is spectral absorption coefficient of impurities, is so-called absorption enhancement parameter for ice grains (Kokhanovsky et al., 2019), is the volumetric concentration of ice grains in snowpack and is the bulk absorption coefficient of ice defined in Eq. (3.4.4). We shall assume that _B __abs__ =1.6_ in the retrieval procedure.

In the visible spectrum, the absorption by ice particle can be neglected () and the polluted snow spherical albedo can be presented in the following form:

|
 |
or
 | (3.4.14) |
| --- | --- | --- |

Let us assume that the impurity absorption coefficient can also be expressed in following form:

|
 |
 | (3.4.15) |
| --- | --- | --- |

where The Angström absorption coefficient can then be derived from Eq. (3.4.14) and (3.4.15) evaluated at 400 nm and 412.5 nm where the spherical albedo was previously derived:

|
 |
 | (3.4.16) |
| --- | --- | --- |

The Angström coefficient can then be used to characterize the type of pollutant present in the snow. Since soot has a typical Angström coefficient around one while dust&#39;s Angström coefficient ranges from 3 to 7, we here assume that the snow is polluted by black carbon if _m_ \&lt; 1.2 and by dust pollution otherwise.

Eq. (3.4.15) evaluated for also gives:

|
 |
 | (3.4.17) |
| --- | --- | --- |

Eq. (3.4.13) can be used to derive , the normalized absorption coefficient of impurities:

|
 |
 | (3.4.18) |
| --- | --- | --- |

Once again, neglecting absorption by ice particle from the measurements at the wavelength 400 nm (), the relative volumetric concentration of pollutants in snow can be derived:

|
 |
 | (3.4.19) |
| --- | --- | --- |

where is the volumetric absorption coefficient of impurities.

In the case of soot, can be approximated as in Kokhanovsky et al., (2018):

|
 |
 | (3.4.20) |
| --- | --- | --- |

Here, is the enhancement is the bulk absorption coefficient of soot , is the imaginary part of soot refractive index, currently set at a constant .

In the case of dust pollution, we assume: to calculate the relative volumetric concentration of pollutants in snow .

  1.
## Broadband albedo calculation

    1.
### General case

The derived spectral albedo is used to integrate the planar and spherical broadband albedo (BBA) over any wavelength interval [:

|
 |
 | (3.5.1) |
| --- | --- | --- |

where is the incident solar flux at the snow surface, is plane (_p_) or spherical (_s_) albedo depending plane or spherical BBA is to be calculated. Currently, only shortwave spherical/plane BBA () is being retrieved but additional ranges may be added in the future depending on user demand.

Broadband albedo are only weakly sensitive to the variation of . The spectrum of incident solar flux at the snow surface is therefore assumed to be identical in all pixels and is approximated by the following analytical equation:

|
 |
 | (3.5.2) |
| --- | --- | --- |

where 3.238e+1, -1.6014033e+5, 7.95953e+3, 11.71, and 2.48. The coefficients have been derived using the code SBDART (Ricchiazzi et al., 1998) in the spectral range 300-2400 nm at the following assumptions.

_ **Table 3. Assumptions used in SBDART to derive the solar flux at the surface** _

| **Parameter** | **Value** |
| --- | --- |
| water vapor column | 2.085 g m-2 |
| --- | --- |
| ozone column | 0.35 atm-cm |
| tropospheric ozone | 0.0346atm-cm |
| aerosol model | rural (Shettle and Fenn, 1979) |
| vertical optical depth of boundary layer at 550nm | 0.1 |
| altitude | 825 m a.s.l. |
| solar zenith angle | 60 degrees |
| snow albedo at the surface | calculated using spherical grains of 0.25 mm diameter |

    1.
### Approximation used for clean snow

In the case of clean snow, the exact integration of BBA (Eq. (3.3.14)) is possible because the spectral reflectance is known for each of OLCI measurement wavelength. However, this integration was time consuming. To speed up the retrieval process, the shortwave spherical albedo can be directly expressed as a function of the retrieved grain diameter:

|
 |
 | (3.5.3) |
| --- | --- | --- |

where _d_ is expressed in microns and a = 0.642, b = 0.1044, c = 0.1773, 158.62 , and = 2448.18 .

For the shortwave broadband plane albedo, we use the same equation as for the spherical albedo. However, the second order polynomial is used to represent the dependence of the coefficients _a,b,c,_ , with respect to the cosine of the solar zenith angle : , … etc. The coefficients of parametrization are given in Table 4.

_ **Table 4. The coefficients of the parametrization for the shortwave plane albedo.** _

|
 |
 |
 |
 |
| --- | --- | --- | --- |
| a | 0.7389 | -0.1783 | 0.0484 |
| b | 0.0853 | 0.0414 | -0.0127 |
| c | 0.1384 | 0.0762 | -0.0268 |
| 𝛥, microns | 187.89 | -69.2636 | 40.4821 |
| , microns | 2687.25 | -405.09 | 94.5 |

    1.
### Approximation used for polluted snow and ice

For polluted snow, the spherical albedo and planar albedo cannot be expressed in a closed form as it is the solution of the system of equation (3.4.10). Nevertheless, ( and respectively) are known for each of the wavelength corresponding to OLCI channels. To circumvent this issue, we build functions of the wavelength that approximate the retrieved over three intervals:

1. Over 400-709 nm, we approximate spherical and planar albedo by a polynomial of the second order fitted to the retrieved 400 nm), 560 nm), and 709 nm).
2. Over 709-865 nm, we approximate spherical and planar albedo by a polynomial of the second order fitted to the retrieved 709 nm), 753 nm), and 865 nm).
3. Over 865-2400 nm, we approximate the spherical and planar albedo with an exponential function fitted to the retrieved 865 nm), and 1020nm).

These assumptions make it possible to derive the value of BBA analytically.

First, for all three intervals, the denominator of Eq. (3.5.1) can be calculated as:

|
 |
 ,
 | (3.5.4) |
| --- | --- | --- |

Over the intervals , either equal to 400-709 nm or 709-865 nm, the spherical albedo can be expressed using its polynomial approximation:

|
 |
 | (3.5.5) |
| --- | --- | --- |

Where a, b and c take different values whether they are fitted to derived or . With this formulation of , the numerator in Eq. (3.5.1) can be expressed in the following form:

|
 |
 | (3.5.6) |
| --- | --- | --- |

where

|
 |
, | (3.5.7) |
| --- | --- | --- |

And

|
 |


 | (3.5.8) |
| --- | --- | --- |

For wavelengths in the 865-2400 nm range, we use the exponential approximation of :

|
 |
 | (3.5.9) |
| --- | --- | --- |

Where and take different values whether they are fitted to derived or . And in that case, the numerator in Eq. (3.5.1) can be expressed in the following form:

|
 |
 | (3.5.10) |
| --- | --- | --- |

1.
# NDSI and NDBI

Several flags are introduced in the snow processor. They are explained in this section.

The snow flag is determined by the value of OLCI normalized difference snow index (NDSI):

|
 |
 | (3.5.1) |
| --- | --- | --- |

The snow flag is equal to one (100% snow – covered pixel), if NDSI is in the range 0.3-0.4 and R(400nm) is larger than 0.75.

The bare ice flag is determined by the value of OLCI normalized difference bare ice index (NDBI):

|
 |
 | (3.5.2) |
| --- | --- | --- |

The bare ice is classified in two steps. First, dark bare ice is identified where NDBI is less than 0.65 and R (400nm) is less than 0.75. Then for cases the dark bare ice flag is not set, the bare ice flag is equal to one (100% bare ice – covered pixel), if NDSI is larger than 0.33. Also is assumed that the dark dirty bare ice flag is equal to one (100% dark dirty bare ice – covered pixel), if NDBI is smaller than 0.65 and R (400nm) is smaller than 0.75 and that a land mask is used.

The values of NDSI and NDBI are provided in the output of the algorithm. In principle, the value of NDSI can be used for the estimation of snow fraction in the OLCI pixel.

1.
# Appendix: Data tables used in the retrieval

_**Table A 1. The ozone vertical optical thickness ( as function of wavelength in terrestrial atmosphere at the ozone concentration of 405 DU = 8.6728e-3**_ **kg m** -2_._

| (nm) | (-) |
| --- | --- |
| 400 | 1.38E-04 |
| --- | --- |
| 412.5 | 3.05E-04 |
| 442.5 | 1.65E-03 |
| 490 | 8.94E-03 |
| 510 | 1.75E-02 |
| 560 | 4.35E-02 |
| 620 | 4.49E-02 |
| 665 | 2.10E-02 |
| 673.75 | 1.72E-02 |
| 681.25 | 1.47E-02 |
| 708.75 | 7.98E-03 |
| 753.75 | 3.88E-03 |
| 761.25 | 2.92E-03 |
| 764.375 | 2.79E-03 |
| 767.5 | 2.73E-03 |
| 778.75 | 3.26E-03 |
| 865 | 8.96E-04 |
| 885 | 5.19E-04 |
| 900 | 6.72E-04 |
| 940 | 3.13E-04 |
| 1020 | 1.41E-05 |

_**Table A 2. The imaginary part of ice refractive index at OLCI bands and wavelength () as reported in Warren and Brandt (1994).**_

| Band | (nm) |
 |
| --- | --- | --- |
| Oa1 | 400 | 2.37E-11 |
| --- | --- | --- |
| Oa2 | 412 | 2.70E-11 |
| Oa3 | 442 | 7.00E-11 |
| Oa4 | 490 | 4.17E-10 |
| Oa5 | 510 | 8.04E-10 |
| Oa6 | 560 | 2.84E-09 |
| Oa7 | 620 | 8.58E-09 |
| Oa8 | 665 | 1.78E-08 |
| Oa9 | 673 | 1.95E-08 |
| Oa10 | 681 | 2.10E-08 |
| Oa11 | 708 | 3.30E-08 |
| Oa12 | 753 | 6.23E-08 |
| Oa13 | 761 | 7.10E-08 |
| Oa14 | 764 | 7.68E-08 |
| Oa15 | 767 | 8.13E-08 |
| Oa16 | 778 | 9.88E-08 |
| Oa17 | 865 | 2.40E-07 |
| Oa18 | 885 | 3.64E-07 |
| Oa19 | 900 | 4.20E-07 |
| Oa20 | 940 | 5.53E-07 |
| Oa21 | 1020 | 2.25E-06 |

# References

Doherty, S. J., Warren, S. G., Grenfell, T. C., Clarke, A. D. and Brandt, R. E.: Light-absorbing impurities in Arctic snow, Atmos. Chem. Phys., 10(23), 11647–11680, doi:10.5194/acp-10-11647-2010, 2010.

Hansen, J. E. and Travis, L. D.: Light scattering in planetary atmospheres, Space Sci. Rev., 16(4), 527–610, doi:10.1007/BF00168069, 1974.

Hudson, S. R., Warren, S. G., Brandt, R. E., Grenfell, T. C. and Six, D.: Spectral bidirectional reflectance of Antarctic snow: Measurements and parameterization, J. Geophys. Res. Atmos., 111(18), D18106, doi:10.1029/2006JD007290, 2006.

Kokhanovsky, A., Mayer, B., Von Hoyningen-Huene, W., Schmidt, S. and Pilewskie, P.: Retrieval of cloud spherical albedo from top-of-atmosphere reflectance measurements performed at a single observation angle, Atmos. Chem. Phys., 7(13), 3633–3637, doi:10.5194/acp-7-3633-2007, 2007.

Kokhanovsky, A., Lamare, M., Mauro, B. Di, Picard, G., Arnaud, L., Dumont, M., Tuzet, F., Brockmann, C. and Box, J. E.: On the reflectance spectroscopy of snow, , 2371–2382, 2018.

Kokhanovsky, A., Lamare, M., Danne, O., Brockmann, C., Dumont, M., Picard, G., Arnaud, L., Favier, V., Jourdain, B., Meur, E. Le, Di Mauro, B., Aoki, T., Niwano, M., Rozanov, V., Korkin, S., Kipfstuhl, S., Freitag, J., Hoerhold, M., Zuhr, A., Vladimirova, D., Faber, A. K., Steen-Larsen, H. C., Wahl, S., Andersen, J. K., Vandecrux, B., van As, D., Mankoff, K. D., Kern, M., Zege, E. and Box, J. E.: Retrieval of snow properties from the Sentinel-3 Ocean and Land Colour Instrument, Remote Sens., 11(19), 1–49, doi:10.3390/rs11192280, 2019.

Kokhanovsky, A., Box, J. E., Vandecrux, B., Mankoff, K. D., Lamare, M., Smirnov, A. and Kern, M.: The determination of snow albedo from satellite measurements using fast atmospheric correction technique, Remote Sens., 12(2), 1–18, doi:10.3390/rs12020234, 2020a.

Kokhanovsky, A., Box, J. E., Vandecrux, B., Mankoff, K. D., Lamare, M., Smirnov, A., Kern, M. and Darmstadt, D.: The Determination of Snow Albedo from Satellite Measurements Using Fast Atmospheric Correction Technique, , doi:10.3390/rs12020234, 2020b.

Kokhanovsky, A. A.: Scaling constant and its determination from simultaneous measurements of light reflection and methane adsorption by snow samples, Opt. Lett., 31(22), 3282, doi:10.1364/OL.31.003282, 2006.

Kokhanovsky, A. A., Lamare, M. and Rozanov, V.: Retrieval of the total ozone over snow fields using Sentinel-3 Ocean and Land Colour Instrument, 2020c.

Proksch, M., Rutter, N., Fierz, C. and Schneebeli, M.: Intercomparison of snow density measurements: Bias, precision, and vertical resolution, Cryosphere, 10(1), 371–384, doi:10.5194/tc-10-371-2016, 2016.

Ricchiazzi, P., Yang, S., Gautier, C. and Sowle, D.: SBDART: A Research and Teaching Software Tool for Plane-Parallel Radiative Transfer in the Earth&#39;s Atmosphere, Bull. Am. Meteorol. Soc., 79(10), 2101–2114, doi:10.1175/1520-0477(1998)079\&lt;2101:SARATS\&gt;2.0.CO;2, 1998.

Rozanov, V. V. and Rozanov, A. V.: Differential optical absorption spectroscopy (DOAS) and air mass factor concept for a multiply scattering vertically inhomogeneous medium: Theoretical consideration, Atmos. Meas. Tech., 3(3), 751–780, doi:10.5194/amt-3-751-2010, 2010.

Serdyuchenko, A., Gorshelev, V., Weber, M., Chehade, W. and Burrows, J. P.: High spectral resolution ozone absorption cross-sections &amp;ndash; Part 2: Temperature dependence, Atmos. Meas. Tech., 7(2), 625–636, doi:10.5194/amt-7-625-2014, 2014.

Shettle, E. P. and Fenn, R. W.: Models for the Aerosols of the Lower Atmosphere and the Effects of Humidity Variations on their Optical Properties, 1979.

Sinnhuber, B. M., Sheode, N., Sinnhuber, M., Chipperfield, M. P. and Feng, W.: The contribution of anthropogenic bromine emissions to past stratospheric ozone trends: A modelling study, Atmos. Chem. Phys., 9(8), 2863–2871, doi:10.5194/acp-9-2863-2009, 2009.

Skiles, S. M. K., Flanner, M., Cook, J. M., Dumont, M. and Painter, T. H.: Radiative forcing by light-absorbing particles in snow, Nat. Clim. Chang., 8(11), 964–971, doi:10.1038/s41558-018-0296-5, 2018.

Sobolev, V. V.: Light scattering in planetary atmospheres, Pergamon Press., 1975.

Warren, S. G.: Can black carbon in snow be detected by remote sensing?, J. Geophys. Res. Atmos., 118(2), 779–786, doi:10.1029/2012JD018476, 2013.

Warren, S. G. and Brandt, R. E.: Optical constants of ice from the ultraviolet to the microwave: A revised compilation, J. Geophys. Res. Atmos., 113(14), D14220, doi:10.1029/2007JD009744, 2008.

[1](#sdfootnote1anc)[_https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-olci/resolutions/radiometric_](https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-olci/resolutions/radiometric)

21