#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

import argparse
import os
import yaml
import copy
import pandas as pd
import numpy as np
import quantities as pq
import joblib
from SmartNeo.user_layer.dict_to_neo import templat_neo
from SmartNeo.interface_layer.nwb_interface import NWBInterface
from user_input_entry import AIEShare as As
import re
import shutil
import json
import hdf5storage
#%% parse the input arguments
parser = argparse.ArgumentParser(description='extract trial info')
parser.add_argument("-f", "--file", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Data_recording/20240403_interception_001', 
                    metavar='/the/path/your/data/located/in', help='input folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Data_recording/20240403_interception_001', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()
raw_dirname = args.file

#%% load template

FILEPATH = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(FILEPATH,'template.yml')) as f:
    Template = yaml.safe_load(f)

data_path = os.path.join(raw_dirname, 'continuous_data')
walk_file = [j for j in os.walk(data_path)]

for f_l in walk_file:
    rec_name = [f_n for f_n in f_l[2] if '.mat' in f_n]
    if len(rec_name) !=0:
        break
datafile = os.path.join(f_l[0], rec_name[0])

traj = hdf5storage.loadmat(datafile)

InputList = []

InputData = copy.deepcopy(Template)
InputData['IrregularSampledData'] = {}
InputData['IrregularSampledData']['irr'] = {}
InputData['IrregularSampledData']['irr'] = templat_neo['irr'].copy()
InputData['IrregularSampledData']['irr']['signal'] = traj['Unlabel'][1][0]*pq.mm
InputData['IrregularSampledData']['irr']['times'] = traj['Unlabel'][0][0].squeeze()/100*pq.s
InputData['name'] = 'Vicon motion'
InputList.append(InputData)

#%% operate the dicts
continuous_bhv_block = As.data_input(user_data = InputList,
                                        index = 3,
                                        name = 'continuous_behavior')

#%% save to appointed path
if not os.path.exists(os.path.join(args.output,'formatted_data')):
    os.mkdir(os.path.join(args.output,'formatted_data'))
    
nwb_saver = NWBInterface()
nwb_saver.save_nwb(blockdata = continuous_bhv_block, 
                    filename = os.path.join(args.output,'formatted_data','continuous_bhv_block.nwb'))
