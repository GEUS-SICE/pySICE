# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 11:14:55 2020

@author: bav
"""

#%%
folder = "out/SICE_test_2020_v2/"
import time
from sice import sice

start_time = time.time()
sice(folder)
time_new = time.time() - start_time

#%%
from sice import sice_old

start_time = time.time()
sice_old(folder)
time_old = time.time() - start_time
