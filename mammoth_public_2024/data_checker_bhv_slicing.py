#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 15:58:49 2024

@author: cuilab
"""

import argparse
import os
import pandas as pd
import quantities as pq
from SmartNeo.interface_layer.nwb_interface import NWBInterface


#%% read data
# data folder

def run(data_dir, output_dir):
    
    continuous_saver = NWBInterface()
    bhv_data = continuous_saver.read_nwb(filename = os.path.join(data_dir,'continuous_behavior.nwb'))

    frame = [i for i in bhv_data.segments if i.name=='frame'][0]
    frame_df = {'frame':frame.events[0].labels}
    frame_df = pd.DataFrame(frame_df)
    frame_df.index = frame.events[0].times.rescale(pq.s).magnitude
    frame_df.index.name = 'Times (s)'
    frame_df.to_csv(os.path.join(output_dir, 'frame_df.csv'))

    coeff_df = {}
    for i in bhv_data.segments[1].irregularlysampledsignals:
        coeff_df[i.name[1]] = i.as_array().squeeze()
    coeff_df = pd.DataFrame(coeff_df)
    # coeff_df.rename_axis('Times')
    coeff_df.index = bhv_data.segments[1].irregularlysampledsignals[0].times.rescale(pq.s).magnitude
    coeff_df.index.name = 'Times (s)'
    coeff_df.to_csv(os.path.join(output_dir, 'coeff_df.csv'))

    pos_df = {}
    for i in bhv_data.segments[3].irregularlysampledsignals:
        pos_df[i.name[1]] = list(i.as_array().squeeze())
    pos_df = pd.DataFrame(pos_df)
    pos_df.rename_axis('Times')
    pos_df.index = bhv_data.segments[3].irregularlysampledsignals[0].times.rescale(pq.s).magnitude
    pos_df.index.name = 'Times (s)'
    pos_df.to_csv(os.path.join(output_dir, 'pos_df.csv'))

    trial_df = {}
    for i in bhv_data.segments[-1].irregularlysampledsignals:
        trial_df[i.name[1]] = list(i.as_array().squeeze())

    trial_df = pd.DataFrame(trial_df)
    trial_df.rename_axis('Times')
    trial_df.index = bhv_data.segments[-1].irregularlysampledsignals[0].times.rescale(pq.s).magnitude
    trial_df.index.name = 'Times (s)'
    trial_df.to_csv(os.path.join(output_dir, 'trial_df.csv'))


parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-d", "--data", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Brain_control/20240307_BrUtahInterception120SemiBc_001/formatted_data', 
                    metavar='/the/path/your/nwb/data/located/in', help='data folder')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Abel/Data_recording/20240705_centerOut_001/description', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()

    
# run(formatted_data_path, description_data_path)

run(args.data, args.output)
