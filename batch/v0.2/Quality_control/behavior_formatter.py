#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 15:58:49 2024

@author: lichenyang
"""

import argparse
import os
import pandas as pd
import numpy as np
import quantities as pq
from SmartNeo.interface_layer.nwb_interface import NWBInterface
from SmartNeo.analysis_layer.spike_preprocessor3 import SpikeStatistics
import matplotlib.pyplot as plt
from elephant.kernels import GaussianKernel

#%% read data
# data folder
unit_dict = {'s': pq.s}
parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-f", "--file", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Brain_control/20240307_BrUtahInterception120SemiBc_001', 
                    metavar='/the/path/your/data/located/in', help='input folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Brain_control/20240307_BrUtahInterception120SemiBc_001', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()
raw_dirname = args.file

formated_data_path = os.path.join(raw_dirname, 'formatted_data')
description_data_path = os.path.join(raw_dirname, 'description')

if not os.path.exists(description_data_path):
    os.mkdir(description_data_path)

# read formated SmartNeo data
nwb_saver = NWBInterface()
neural_data = nwb_saver.read_nwb(filename = os.path.join(formated_data_path,'neural_data.nwb'))

nwb_saver = NWBInterface()
trial_behavior = nwb_saver.read_nwb(filename = os.path.join(formated_data_path,'trial_behavior.nwb'))

index = [i.index for i in trial_behavior.segments]
sorted_index = np.argsort(index)
trials = [trial_behavior.segments[i] for i in sorted_index]

# index
trial_index = [i.index for i in trials]
# event
event = {}
event['Labels'] = [i.events[0].labels for i in trials]
event['Times'] = [i.events[0].times for i in trials]
df_event = pd.DataFrame(event)
df_event.columns=[len(df_event.columns)*['Event'], list(df_event.columns)]
df_event.index=trial_index

# description
description = {}
description['Block'] = [i.description['Block'] for i in trials]
description['AngularV'] = [float(i.description['UserVars']['angularV'].squeeze()) for i in trials]
description['AbsoluteTrialStartTime'] = [i.description['AbsoluteTrialStartTime'] for i in trials]
description['Condition'] = [i.description['Condition'] for i in trials]
description['TrialError'] = [i.description['TrialError'] for i in trials]
df_description = pd.DataFrame(description)
df_description.columns=[len(df_description.columns)*['Description'], list(df_description.columns)]
df_description.index=trial_index


# ObjectStatusRecord
ObjectStatusRecord = {}
obj_pos = [[j for j in i.irregularlysampledsignals if 'Position' in j.name][0] for i in trials]
obj_sta = [[j for j in i.irregularlysampledsignals if 'Status' in j.name][0] for i in trials]
ObjectStatusRecord['Position'] = [np.array(i) for i in obj_pos]
ObjectStatusRecord['Status'] = [np.array(i) for i in obj_sta]
df_ObjectStatusRecord = pd.DataFrame(ObjectStatusRecord)
df_ObjectStatusRecord.columns=[len(df_ObjectStatusRecord.columns)*['ObjectStatusRecord'], list(df_ObjectStatusRecord.columns)]
df_ObjectStatusRecord.index=trial_index

frames = [df_event, df_description, df_ObjectStatusRecord]
df_trials = pd.concat(frames,axis=1).rename_axis('TrialIndex')

df_trials.to_csv(os.path.join(description_data_path, 'df_trials.csv'))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


