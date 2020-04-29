# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:09:38 2019

@author: bav
"""

from satpy.scene import Scene
from satpy import find_files_and_readers
from datetime import datetime
 
files = find_files_and_readers(sensor='olci',
                               start_time=datetime(2017, 7, 15, 0, 0, 0),
                               end_time=datetime(2017, 7, 16, 0, 0, 0),
                               base_dir=".\OLCI scenes",
                               reader='olci_l1b')
 
scn = Scene(filenames=files)
scn.load(['true_color'])

scn.save_dataset('true_color', filename='scene'+'.png')