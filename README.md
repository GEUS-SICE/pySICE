# pySICE
Python scripts behind the SICE toolchain for albedo retrieval.

by  B. Vandecrux (2), A. Kokhanovsky (1), J. Box (2)

(1) VITROCISET Belgium SPRL, Bratustrasse 7, 64293 Darmstadt, Germany
(2) Geological Survey of Denmark and Greenland (GEUS)
 Øster Voldgade 10, 1350 Copenhagen, Denmark


## Table of Contents  
*[Running environment](#running-environment)  
*[Theoretical background](#theoretical-background)  
*[Script description](#script-description)  
	*[Script overview](#script-overview)  
	*[Input preparation](#input-preparation)  
	*[Clean or polluted snow pixels](#clean-or-polluted-snow-pixels)  
	*[Bottom of the atmosphere reflectance and broadband albedo](#bottom-of-the-atmosphere-reflectance-and-broadband-albedo)  

## Running environment
developped on Python 3.7.6
uses numpy, rasterio, time and sys packages

## Theoretical background

The Algorithm Theoretical Basis Document is available [here](docs/atbd/FINAL_SICE_ATBD__v3.0_MAY06_2020.pdf)

## Script description

### Script overview
![](docs/atbd/ATBD_plots1.png)

### Input preparation
![](docs/atbd/SICE_overview1.png)

### Clean or polluted snow pixels
![](docs/atbd/SICE_overview2.png)

### Bottom of the atmosphere reflectance and broadband albedo
![](docs/atbd/SICE_overview3.png)

