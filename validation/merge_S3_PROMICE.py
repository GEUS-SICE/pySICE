# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 07:01:37 2020

@author: bav
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 08:58:21 2020
Data processing for the evaluation of SICE against PROMICE measurements
@author: bav
"""

import numpy as np
import pandas as pd
import os
import errno
from pathlib import Path
import wget # conda install: https://anaconda.org/anaconda/pywget

#%% preparing input file
path = './data/'
filename = 'S3_28072020'
data = pd.read_csv(path+filename+'.csv');

cols = data.columns
cols = [c for c in cols if c.find('radiance')<0]
cols = [c for c in cols if c.find('solar')<0]
cols = [c for c in cols if c.find('temperature')<0]
cols = [c for c in cols if c.find('spectral')<0]
cols = [c for c in cols if c.find('rBRR')<0]
cols = [c for c in cols if c.find('albedo')<0]
cols = [c for c in cols if c.find('variance')<0]
cols = [c for c in cols if c.find('grain')<0]
cols = [c for c in cols if c.find('snow')<0]
cols = [c for c in cols if c.find('ndsi')<0]
cols = [c for c in cols if c.find('ndbi')<0]
cols = [c for c in cols if c.find('wind')<0]
cols = [c for c in cols if c.find('OAA')<0]
cols = [c for c in cols if c.find('SAA')<0]
cols = [c for c in cols if c.find('OZA')<0]
cols = [c for c in cols if c.find('SZA')<0]
cols = [c for c in cols if c.find('plat')<0]
cols = [c for c in cols if c.find('aspect')<0]
cols = [c for c in cols if c.find('_y')<0]
cols = [c for c in cols if c.find('Unname')<0]
data=data[cols]

for c in cols:
    if c.find('_x')>0:
        data.rename(columns={c:c.replace('_x', '')}, inplace=True)


#%% Adding PROMICE observations
# Array information of stations available at PROMICE official site: https://promice.org/WeatherStations.html
PROMICE_stations = [('EGP',(75.6247,-35.9748), 2660), 
                   ('KAN_B',(67.1252,-50.1832), 350), 
                   ('KAN_L',(67.0955,-35.9748), 670), 
                   ('KAN_M',(67.0670,-48.8355), 1270), 
                   ('KAN_U',(67.0003,-47.0253), 1840), 
                   ('KPC_L',(79.9108,-24.0828), 370),
                   ('KPC_U',(79.8347,-25.1662), 870), 
                   ('MIT',(65.6922,-37.8280), 440), 
                   ('NUK_K',(64.1623,-51.3587), 710), 
                   ('NUK_L',(64.4822,-49.5358), 530),
                   ('NUK_U',(64.5108,-49.2692), 1120),
                   ('QAS_L',(61.0308,-46.8493), 280),
                   ('QAS_M',(61.0998,-46.8330), 630), 
                   ('QAS_U',(61.1753,-46.8195), 900), 
                   ('SCO_L',(72.2230,-26.8182), 460),
                   ('SCO_U',(72.3933,-27.2333), 970),
                   ('TAS_A',(65.7790,-38.8995), 890),
                   ('TAS_L',(65.6402,-38.8987), 250),
                   ('THU_L',(76.3998,-68.2665), 570),
                   ('THU_U',(76.4197,-68.1463), 760),
                   ('UPE_L',(72.8932,-54.2955), 220), 
                   ('UPE_U',(72.8878,-53.5783), 940)]

path_to_PROMICE = "./data/PROMICE"  

# Function for making directories if they do not exists. 
def mkdir_p(path):
    try:
        os.makedirs(path)
        return 'Path created.'
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            return 'Path already exists!'
        else:
            raise
            
mkdir_p(path_to_PROMICE)

# Goes through each station and fetch down data online. Necessary manipulations and sorting are made.
for ws in PROMICE_stations:
    if Path(f'{path_to_PROMICE}/{ws[0]}_hour_v03.txt').is_file():
        print('\nPROMICE.csv file already exists.')
        pass
    else:
        print(ws)
        url = f'https://promice.org/PromiceDataPortal/api/download/f24019f7-d586-4465-8181-d4965421e6eb/v03/hourly/csv/{ws[0]}_hour_v03.txt'
        filename = wget.download(url, out= path_to_PROMICE + f'/{ws[0]}_hour_v03.txt')
        
if Path(f'{path_to_PROMICE}/PROMICE.csv').is_file():
    print('\nPROMICE.csv file already exists.')
    pass
else: 
    # Create one big file "PROMICE.csv" with all the stations data. 
    # reading PROMICE file and joining them together
    for ws in PROMICE_stations:
        filepath = path_to_PROMICE + '/' + ws[0] + '_hour_v03.txt'
    
        df = pd.read_csv (filepath, delim_whitespace=True)
        df = df[['Year','MonthOfYear','DayOfYear','HourOfDay(UTC)','CloudCover',
                 'TiltToEast(d)','TiltToNorth(d)','Albedo_theta<70d']]
        df = df[df.Year > 2015]
        df['station_name'] = ws[0]
        df['latitude N'] = ws[1][0]
        df['longitude W'] = ws[1][1]
        df['elevation'] = float(ws[2])
        df['TiltToEast(d)'].replace(-999.0, np.nan, inplace=True)
        try :
            df['TiltToEast(d)'] = df['TiltToEast(d)'].interpolate(method='nearest', limit_direction='both')
        except:
            print('no tilt at '+ws[0])
        df['TiltToNorth(d)'].replace(-999.0, np.nan, inplace=True)
        try: 
            df['TiltToNorth(d)'] = df['TiltToNorth(d)'].interpolate(method='nearest', limit_direction='both')
        except:
            print('no tilt at '+ws[0])
    
        df.to_csv(path_to_PROMICE + '/' + ws[0] + '.csv', index=None)
    
    PROMICE = pd.DataFrame()
    filelist = [ f for f in os.listdir(path_to_PROMICE) if f.endswith(".csv") ]
    for f in filelist:
        print(f)
        PROMICE = PROMICE.append(pd.read_csv(f'{path_to_PROMICE}/{f}'))
    
    PROMICE['TiltToEast(d)'] = PROMICE['TiltToEast(d)'].abs()
    PROMICE['TiltToNorth(d)'] = PROMICE['TiltToNorth(d)'].abs()
    PROMICE['tilt'] = PROMICE[['TiltToEast(d)','TiltToNorth(d)']].values.max(1)
    PROMICE.drop(['TiltToNorth(d)','TiltToEast(d)'], axis=1, inplace=True)
    PROMICE.rename(columns={'Year': 'year', 
                            'MonthOfYear': 'month',
                            'DayOfYear':'dayofyear',
                            'HourOfDay(UTC)':'hour',
                            'CloudCover':'cloud',
                            'station_name':'station',
                            'Albedo_theta<70d':'PROMICE_alb'}, inplace=True)
    PROMICE.to_csv(f'{path_to_PROMICE}/PROMICE.csv', index=None)
    
#%% Merging
data_PROMICE = pd.read_csv(path_to_PROMICE+'/PROMICE.csv')
data.rename(columns={'site': 'station'}, inplace=True)

full_size = data.shape[0]
data = data.drop_duplicates(keep='first') #dropping duplicates except first sample, if any.
print(f'Removed duplicates: {full_size-data.shape[0]}')

# the PROMICE data is only given in hourly interval. Hence the hour have to be corrected in the s3_extract data.
#s3_extract['hour'] = s3_extract['hour'] - 1 if s3_extract['minute'].astype('int64') < 30 else (s3_extract['hour'] + 1) 

# Merge images with PROMICE data
data_all = pd.merge(left=data, right=data_PROMICE, how='inner',
                             left_on=['year','month','dayofyear','hour','station'], 
                             right_on=['year','month','dayofyear','hour', 'station'])
# data_all.drop(columns=['Unnamed: 0','Unnamed: 0.1'],inplace=True)
# keeping only good data
print(data_all.shape[0])
data_all = data_all.replace(to_replace=-999,value=np.nan)
data_all = data_all[data_all['PROMICE_alb'].notnull()]
data_all = data_all[data_all['PROMICE_alb']<1]
data_all = data_all[data_all['PROMICE_alb']>0]
data_all = data_all[data_all['sza']<70]
print(data_all.shape[0])
tmp =data_all.shape[0]

data_all['BBA_emp'] = 0.945*(data_all["Oa01_reflectance"] + data_all["Oa06_reflectance"] + data_all["Oa17_reflectance"] + data_all["Oa21_reflectance"]) / 4.0+ 0.055

# Create merged csv
data_all.to_csv('./data/S3_PROMICE.csv', index=None) 
