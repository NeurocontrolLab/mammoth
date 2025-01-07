#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


import scipy
import numpy as np
import pandas as pd
import brpylib
import argparse
from datetime import datetime
from zoneinfo import ZoneInfo 
from probeinterface import read_probeinterface

from pynwb import NWBHDF5IO, NWBFile
from pynwb.ecephys import ElectricalSeries
from pynwb.file import Subject
from pynwb.core import DynamicTable


#%% define functions
def get_timestamp(root_dir):
    walk_file = [j for j in os.walk(root_dir)]

    try:
        for f_l in walk_file:
            rec_name = [f_n for f_n in f_l[2] if '.ns6' in f_n]
            if len(rec_name) !=0:
                break
        datafile = os.path.join(f_l[0], rec_name[0])
    except:
        sys.exit()

    nsx_file = brpylib.NsxFile(str(datafile))

    # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
    cont_data = nsx_file.getdata(full_timestamps=True)

    # Close the nsx file now that all data is out
    nsx_file.close()

    data = np.concatenate(cont_data['data'],1).T
    timestamp = np.concatenate([i['Timestamp'] for i in cont_data['data_headers']])

    return data, timestamp


def format_file(root_dir, map_path, output_dir):
   
    #%% Set basic info and create NWB file
    root_directory = os.path.join(root_dir, 'bhv')
    file_pattern = [i for i in os.listdir(root_directory) if ('csv' in i) and ('meta' in i)][0]
    assert os.path.exists(os.path.join(root_directory, file_pattern)), "Please define the metadata first!"
    meta_dict_pd = pd.read_csv(os.path.join(root_directory, file_pattern)).to_dict()
    meta_dict = {}
    for i,j in zip(meta_dict_pd['Key'], meta_dict_pd['Value']):
        meta_dict[meta_dict_pd['Key'][i]] = meta_dict_pd['Value'][j]

    subject_dict = {}
    nwb_dict = {}
    for i in meta_dict:
        key = i.split('/')[1]
        if 'Subject' in i:
            subject_dict[key] = meta_dict[i]
        else:
            nwb_dict[key] = meta_dict[i]

    if not 'description' in subject_dict:
        subject_dict['description'] = 'None'

    subject = Subject(
        **subject_dict
        )

    session_id = os.path.basename(root_dir)
    nwb_dict['identifier'] = session_id
    nwb_dict['session_id'] = nwb_dict['session_number']
    del nwb_dict['session_number']
    nwb_dict['lab'] = nwb_dict['lab_name']
    del nwb_dict['lab_name']
    nwb_dict['session_start_time'] = \
        datetime.strptime(nwb_dict['session_start_time'], 
                        '%Y/%m/%d %H:%M').replace(tzinfo=ZoneInfo('Asia/Shanghai'))
    nwb_dict['subject'] = subject
    nwb_dict['keywords'] = ["ecephys", "monkey", "motor control"]
    del nwb_dict['experiment_name']
    if 'task' in nwb_dict.keys():
        del nwb_dict['task']

    for k, v in nwb_dict.items():
        if isinstance(v, float):
            if np.isnan(v):
                nwb_dict[k] = ''

    nwbfile = NWBFile(
        **nwb_dict
    )
    
    #%% get raw data path
    try:
        walk_file = [j for j in os.walk(root_dir)]
        for f_l in walk_file:
            rec_name = [f_n for f_n in f_l[2] if '.ns6' in f_n]
            if len(rec_name) !=0:
                break
        rec_dir = f_l[0]
    except:
        sys.exit()
    
    raw_dir = os.path.join(root_dir, rec_dir)

    # get timestamps
    data, timestamp = get_timestamp(raw_dir)


    # %% Add ecephys data
    # Load probe
    probegroup = read_probeinterface(map_path)

    probe_info_str = os.path.basename(map_path).split(".")[0]
    probe_info = probe_info_str.split("_")
    device_dict = {'name': '%s array' % probe_info[1], 'description': '%s' % probe_info[2],
                   'manufacturer': '%s Microsystem' % probe_info[4]} 

    # set device and electrode table
    device = nwbfile.create_device(name=device_dict['name'],
                                description=device_dict['description'],
                                manufacturer=device_dict['manufacturer'])

    # get electrode table
    # nshanks = len(probegroup.probes)
    shank_location = probe_info[3].split("-")
    ETR_list = []

    for ishank, p in enumerate(probegroup.probes):

        electrode_group = nwbfile.create_electrode_group(
            name="shank{}".format(ishank),
            description="electrode group for shank {}".format(ishank),
            device=device,
            location=shank_location[ishank]
        )
        ETR_list.append(electrode_group)

        # add electrodes to the electrode table
        for ielec in range(p.contact_positions.shape[0]):
            nwbfile.add_electrode(
                x=float(p.contact_positions[ielec, 0]),
                y=float(p.contact_positions[ielec, 1]),
                group=electrode_group,
                id = p.device_channel_indices[ielec],
                location=shank_location[ishank])

    # view the electrodes table in a pandas DataFrame
    # nwbfile.electrodes.to_dataframe()


    #%% add raw data
    electrode_index = nwbfile.electrodes.id.data
    all_table_region = nwbfile.create_electrode_table_region(
        region=list(range(len(electrode_index))), 
        description="all electrodes",)

    elec_data = data[:, electrode_index]

    raw_electrical_series = ElectricalSeries(
        name="raw",
        description="raw data",
        data=elec_data,
        electrodes=all_table_region,
        timestamps=timestamp
    )

    ecephys_module = nwbfile.create_processing_module(
        name="ecephys", 
        description="original extracellular electrophysiology data"
    )
    ecephys_module.add(raw_electrical_series)


    # %% Add time difference data
    # load saved time difference data
    filename_ = os.path.join(root_dir, 'description', 'diff_time_mean.txt')
    assert os.path.exists(filename_), "Please run time consistency check first!"
    with open(filename_, 'r') as file:
        dt_str = file.read()  

    diff_time_mean = float(dt_str.split(' ')[0])

    # add time difference
    behavior_ecephys_module = nwbfile.create_processing_module(
        name='behavior_ecephys_analysis', 
        description='Quality control and pre-analysis results')

    time_diff_table = DynamicTable(
        name='TimeDifference',
        description='Time difference between behavior and ecephys')

    time_diff_table.add_column(
        name='time_difference',
        description='Behavior - Ecephys (s)'
    )

    time_diff_table.add_row(id=0, time_difference=diff_time_mean)
    behavior_ecephys_module.add(time_diff_table)


    #%% Writing standard NWB file
    save_path = os.path.join(output_dir, "raw_data.nwb")
    if os.path.exists(save_path):
        os.remove(save_path)

    with NWBHDF5IO(save_path, "w") as io:
        io.write(nwbfile)


#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20240925_interception_004', 
                    metavar='/the/path/your/data/located/in', help='input folder')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20240925_interception_004/formatted_data', 
                    metavar='/the/path/you/want/to/save', help='output folder')

# parser.add_argument('-mp', '--map_path', 
#                     default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Abel_Utah_64x4_PMd-M1-S1-A7_BlackRock.json')

parser.add_argument('-mp', '--map_path', 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Bohr_Utah_96x2_PMd-M1_BlackRock.json')


args = parser.parse_args()


format_file(args.root, args.map_path, args.output)        
        