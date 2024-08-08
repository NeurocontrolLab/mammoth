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

description_data_path = os.path.join(raw_dirname, 'description')
formated_data_path = os.path.join(raw_dirname, 'formatted_data')
continuous_saver = NWBInterface()
bhv_data = continuous_saver.read_nwb(filename = os.path.join(formated_data_path,'continuous_behavior.nwb'))

frame = [i for i in bhv_data.segments if i.name=='frame'][0]
frame_pd = {'frame':frame.events[0].labels}
frame_pd = pd.DataFrame(frame_pd)
frame_pd.index = frame.events[0].times.rescale(pq.s).magnitude
frame_pd.index.name = 'Times (s)'
frame_pd.to_csv(os.path.join(description_data_path, 'frame_pd.csv'))

coeff_pd = {}
for i in bhv_data.segments[1].irregularlysampledsignals:
    coeff_pd[i.name[1]] = i.as_array().squeeze()
coeff_pd = pd.DataFrame(coeff_pd)
# coeff_pd.rename_axis('Times')
coeff_pd.index = bhv_data.segments[1].irregularlysampledsignals[0].times.rescale(pq.s).magnitude
coeff_pd.index.name = 'Times (s)'
coeff_pd.to_csv(os.path.join(description_data_path, 'coeff_pd.csv'))

pos_pd = {}
for i in bhv_data.segments[3].irregularlysampledsignals:
    pos_pd[i.name[1]] = list(i.as_array().squeeze())

pos_pd = pd.DataFrame(pos_pd)
pos_pd.rename_axis('Times')
pos_pd.index = bhv_data.segments[3].irregularlysampledsignals[0].times.rescale(pq.s).magnitude
pos_pd.index.name = 'Times (s)'
pos_pd.to_csv(os.path.join(description_data_path, 'pos_pd.csv'))

trial_pd = {}
for i in bhv_data.segments[-1].irregularlysampledsignals:
    trial_pd[i.name[1]] = list(i.as_array().squeeze())

trial_pd = pd.DataFrame(trial_pd)
trial_pd.rename_axis('Times')
trial_pd.index = bhv_data.segments[-1].irregularlysampledsignals[0].times.rescale(pq.s).magnitude
trial_pd.index.name = 'Times (s)'
trial_pd.to_csv(os.path.join(description_data_path, 'trial_pd.csv'))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


