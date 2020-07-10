# pySICE
Python scripts behind the SICE toolchain for albedo retrieval.

by  B. Vandecrux (2), A. Kokhanovsky (1), J. Box (2)

(1) VITROCISET Belgium SPRL, Bratustrasse 7, 64293 Darmstadt, Germany
(2) Geological Survey of Denmark and Greenland (GEUS)
 Øster Voldgade 10, 1350 Copenhagen, Denmark


## Table of Contents  
* [Running environment](#running-environment)  
* [Theoretical background](#theoretical-background)  
* [Scripts description](#scripts-description)  
    * [Scripts overview](#scripts-overview)  
    * [Input preparation](#input-preparation)  
    * [Clean or polluted snow pixels](#clean-or-polluted-snow-pixels)  
    * [Bottom of the atmosphere reflectance and broadband albedo](#test)  
* [Installation](#installation)  
    * [Python](#python)
    * [Download pySICE](#download)
* [Examples](#examples)  
    * [pySICE](#pytsice)
    * [Python interface for the fortran script sice.f	](#sicef)


	
## Running environment
developped on Python 3.7.6
uses numpy, rasterio, time and sys packages

## Theoretical background

The snow surface characteristics retrieval is based on the following work:
[Kokhanovsky et al. (2018) On the reflectance spectroscopy of snow](https://tc.copernicus.org/articles/12/2371/2018/)
[Kokhanovsky et al. (2019) Retrieval of Snow Properties from the Sentinel-3
Ocean and Land Colour Instrument](http://dx.doi.org/10.3390/rs11192280)
[Kokhanovsky et al. (2018) The Determination of Snow Albedo from Satellite
Measurements Using Fast Atmospheric
Correction Technique](http://dx.doi.org/10.3390/rs12020234)
[Kokhanovsky et al. (2018) On the reflectance spectroscopy of snow](https://tc.copernicus.org/articles/12/2371/2018/)

The ozone total ozone retreival is described in 
[Kokhanovsky et al. (2020) Retrieval of the total ozone over Antarctica using Sentinel-3 ocean and land colour instrument](https://doi.org/10.1016/j.jqsrt.2020.107045)

The Algorithm Theoretical Basis Document is available [here](docs/atbd/FINAL_SICE_ATBD__v3.0_MAY06_2020.pdf)

## Scripts description

### Scripts overview
![](docs/atbd/ATBD_plots1.png)

### Input preparation
![](docs/atbd/SICE_overview1.png)

### Clean or polluted snow pixels
![](docs/atbd/SICE_overview2.png)

<a name="test"/>
### Bottom of the atmosphere reflectance and broadband albedo
![](docs/atbd/SICE_overview3.png)

## Installation

# Python

We recommend the use of [Anaconda](https://www.anaconda.com/products/individual) and recent version of Python (>3.7).


<a name="download"/>
# Download pySICE

Download the repository using the browser or typing in your command prompt
```
git clone https://github.com/BaptisteVandecrux/pySICE.git
```

## Examples

Test input files are available [here](https://www.dropbox.com/s/9pb9n0k54ev3yg4/S3_test_data.zip?dl=0).

### pySICE



<a name="sicef"/>
### Python interface for the fortran script sice.f



