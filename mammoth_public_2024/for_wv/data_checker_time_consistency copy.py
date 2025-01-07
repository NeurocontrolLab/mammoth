#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 15:09:19 2024

@author: cuilab
"""

import os
from dtw import accelerated_dtw
import json
import numpy as np
import argparse
import quantities as pq
import matplotlib.pyplot as plt
# from SmartNeo.interface_layer.nwb_interface import NWBInterface
from pynwb import NWBHDF5IO


def run(data_dir, output_dir):
    # load data
    filename = os.path.join(data_dir, 'neural_data_no_sort.nwb')
    neural_data = NWBHDF5IO(filename, mode='r').read()

    filename = os.path.join(data_dir, 'continuous_behavior.nwb')
    bhv_data = NWBHDF5IO(filename, mode='r').read()

    # get events from neural recording system
    neural_event_labels = neural_data.processing['ecephys']['RecordingSystemEvents'].data[:]
    neural_event_times = neural_data.processing['ecephys']['RecordingSystemEvents'].timestamps[:]*pq.s

    if neural_event_labels[0]>60000:
        neural_event_labels = np.array(neural_event_labels)-65280
    elif neural_event_labels[0]<-60000:
        neural_event_labels = np.array(neural_event_labels)+65280
    
    neural_event_times = neural_event_times[neural_event_labels!=0]
    neural_event_labels = neural_event_labels[neural_event_labels!=0]

    # get events from trial operating system
    ml_event_labels = bhv_data.processing['behavior']['MonkeyLogicEvents'].data[:]
    ml_event_times = bhv_data.processing['behavior']['MonkeyLogicEvents'].timestamps[:] * pq.s
    
    # compute distance between two events (time difference)
    li_distance = lambda x, y: 0 if np.abs(x - y)==0 else len(neural_event_labels)
        # manhattan_distance = lambda x, y: np.abs(x - y)
    
    _,_,_,path1 = accelerated_dtw(ml_event_labels, neural_event_labels, dist=li_distance)
    path1 = np.array(path1).T
    # zero_dir = (ml_event_labels[path1[:,0]]-neural_event_labels[path1[:,1]])==0
    # assert len(neural_event_labels)-sum(zero_dir)<(0.005*len(neural_event_labels)), 'wrong alignment'
    # a = (len(neural_event_labels)-sum(zero_dir))
    diff_time = ml_event_times[path1[:,0]]-neural_event_times[path1[:,1]]

    def remove_outliers(data, factor=1.5):
        q1 = np.percentile(data, 25)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        lower_bound = q1 - (factor * iqr)
        upper_bound = q3 + (factor * iqr)
        return data[(data > lower_bound) & (data < upper_bound)]
    
    diff_time = remove_outliers(diff_time, factor=1.5)
    plt.boxplot(diff_time-np.mean(diff_time))
    # len(diff_time)/len(neural_event_labels)

    #%% save to appointed path
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    plt.savefig(os.path.join(output_dir, 'Time_consistency_check.png'), dpi=300, bbox_inches = 'tight')

    diff_time_mean = np.mean(diff_time)

    with open(os.path.join(output_dir, 'diff_time_mean.txt'), 'w') as file:
        file.write(f'{diff_time_mean}\n')

    with open(os.path.join(output_dir, 'qc_summary.json'), 'w') as file:
        summary = {}
        summary['diff_time [q1, median, q3]'] = [np.percentile(diff_time-diff_time_mean, 25),
                                                 np.percentile(diff_time-diff_time_mean, 50),
                                                 np.percentile(diff_time-diff_time_mean, 75)]
        json.dump(summary, file)



parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-d", "--data", type=str,
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/interception/Abel/data_recording/20240501_Interception_001/formatted_data', 
                    metavar='/the/path/your/data/located/in', help='data folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/interception/Abel/data_recording/20240501_Interception_001/description', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()

run(args.data, args.output)
