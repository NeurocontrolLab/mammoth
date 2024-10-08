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
    command = 'python /home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Code/nwb_builder/data_integration_pipeline_traj/neural_data/br_share/pipeline.py -f {} -o {}'.format(data_path,data_path)
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    command = 'python /home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Code/nwb_builder/data_integration_pipeline_traj/trial_behavior/monkeylogic_share/pipeline.py -f {} -o {}'.format(data_path,data_path)
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    command = 'python /home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Code/nwb_builder/data_integration_pipeline_traj/continuous_behavior/aie_share/pipeline.py -f {} -o {}'.format(data_path,data_path)
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)