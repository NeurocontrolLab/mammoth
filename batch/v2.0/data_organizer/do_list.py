#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 18:46:27 2024

@author: lichenyang
"""

import subprocess
import os

path = '/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Brain_control'
for p in os.listdir(path):
    data_path = os.path.join(path,p)
    if 'formatted_data' in os.listdir(data_path):
        continue
    command = 'python /home/cuihe_lab/lichenyang/DATA_AMAX/code/data_organization.py -f {}'.format(data_path)
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)