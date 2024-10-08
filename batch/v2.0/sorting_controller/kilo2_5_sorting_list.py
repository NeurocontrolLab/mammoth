#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 17:05:20 2024

@author: lichenyang
"""

import subprocess
import os

path = '/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Data_recording'
for p in os.listdir(path):
    data_path = os.path.join(path,p)
    if 'formatted_data' in os.listdir(data_path):
        continue
    command = ''
