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
import matplotlib.pyplot as plt
from SmartNeo.interface_layer.nwb_interface import NWBInterface


def run(data_dir, output_dir):
    # load data
    nwb_saver = NWBInterface()
    neural_data = nwb_saver.read_nwb(filename = os.path.join(data_dir, 'neural_data_no_sort.nwb'))

    nwb_saver = NWBInterface()
    bhv_data = nwb_saver.read_nwb(filename = os.path.join(data_dir, 'continuous_behavior.nwb'))

    # get events from neural recording system
    neural_event = [i for i in neural_data.segments if i.name=='RecordingSystemEvent'][0]
    neural_event_labels = neural_event.events[0].labels
    neural_event_times = neural_event.events[0].times

    # get events from trial operating system
    ml_event_labels0 = bhv_data.segments[0].events[0].labels
    ml_event_labels = np.array(
        [json.loads(i)['Marker'] for i in ml_event_labels0 if isinstance(json.loads(i), dict)])
    ml_event_times = bhv_data.segments[0].events[0].times

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
    len(diff_time)/len(neural_event_labels)

    #%% save to appointed path
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    plt.savefig(os.path.join(output_dir, 'Time_consistency_check.png'), dpi=300, bbox_inches = 'tight')

    diff_time_mean = np.mean(diff_time)

    with open(os.path.join(output_dir, 'diff_time_mean.txt'), 'w') as file:
        file.write(f'{diff_time_mean}\n')


parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-d", "--data", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20240920_interception_003/formatted_data', 
                    metavar='/the/path/your/data/located/in', help='data folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20240920_interception_003/description', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()

run(args.data, args.output)
