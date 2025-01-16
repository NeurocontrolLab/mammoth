#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 10:02:11 2024

@author: cuilab
"""

#reference: https://pynwb.readthedocs.io/en/stable/tutorials/domain/ecephys.html#sphx-glr-tutorials-domain-ecephys-py

import os
# import copy
# import json
import numpy as np
# import pandas as pd
import argparse
# import quantities as pq

from pynwb import NWBHDF5IO, NWBFile, TimeSeries
from pynwb.behavior import BehavioralEvents
from pynwb.core import DynamicTable
from pynwb.ecephys import LFP, ElectricalSeries

# from SmartNeo.interface_layer.nwb_interface import NWBInterface
from probeinterface import read_probeinterface


#%%

def run(root_dir, map_path, output_dir):
   
    # %% Create a new NWB file
    # load ecephys data
    filename = os.path.join(root_dir, 'formatted_data', 'neural_data.nwb')
    neural_data = NWBHDF5IO(filename, mode='r').read()
    
    # synchronize session and subject information
    session_id = os.path.basename(root_dir)
    nwbfile = NWBFile(    
       experiment_description=neural_data.experiment_description+'ecephys + bhv',
       experimenter=neural_data.experimenter,
    #    file_create_data = (),
       identifier=session_id + 'collection',
       institution=neural_data.institution,
       keywords=neural_data.keywords,
       lab=neural_data.lab,
       session_id=neural_data.session_id,
       session_start_time=neural_data.session_start_time,
       session_description=neural_data.session_description,
       timestamps_reference_time=neural_data.timestamps_reference_time,
    )

    if neural_data.subject:
        source_subject = neural_data.subject
        new_subject = source_subject.__class__(
            subject_id=source_subject.subject_id,
            age=source_subject.age,
            description=source_subject.description,
            genotype=source_subject.genotype,
            sex=source_subject.sex,
            species=source_subject.species,
            weight=source_subject.weight,
            date_of_birth=source_subject.date_of_birth,
            strain=source_subject.strain,
        )
        nwbfile.subject = new_subject


    # %% Add ecephys data
      # load probe
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


    # # view the electrodes table in a pandas DataFrame
    # nwbfile.electrodes.to_dataframe()
  

    # define unit columns
    neural_units = neural_data.units.to_dataframe()

    for col in neural_data.units.colnames:
            if col not in ['spike_times', 'waveform_mean', 'electrode_group']:
                nwbfile.add_unit_column(name=col, 
                                        description=neural_data.units[col].description)
    
    for i, row in neural_units.iterrows():
        new_dict = {k:v for k, v in row.to_dict().items() if k!='electrode_group'}
        # new_dict['electrode_group'] = ETR_list[int([k for k, n in enumerate(ETR_list) if n.name==neural_units['electrode_group'][i].name][0])],
        nwbfile.add_unit(**new_dict, 
                         electrode_group=ETR_list[int([k for k, n in enumerate(ETR_list) if n.name==neural_units['electrode_group'][i].name][0])])
    

    ## validate 
    neural_units_v = neural_units.copy()
    neural_units_v['electrode_group'] = neural_units_v['electrode_group'].apply(lambda x: x.name)
    neural_units_v['spike_times'] = neural_units_v['spike_times'].apply(lambda x: x.tolist())

    nwbfile_units_v = nwbfile.units.to_dataframe().copy()
    nwbfile_units_v['electrode_group'] = nwbfile_units_v['electrode_group'].apply(lambda x: x.name)
    assert neural_units_v['chn_id'].equals(nwbfile_units_v['chn_id'])


    # add LFP data 
    ecephys_module = nwbfile.create_processing_module(
            name="ecephys", 
            description="processed extracellular electrophysiology data"
        )
    
    electrode_index = nwbfile.electrodes.id.data   

    all_table_region = nwbfile.create_electrode_table_region(
        region=list(range(len(electrode_index))), 
        description="all electrodes",
    )

    lfp_electrical_series = ElectricalSeries(
        name="LFP",
        description="LFP data",
        data=np.expand_dims(neural_data.processing['ecephys']['LFP']['LFP'].data[:], axis=0),
        electrodes=all_table_region,
        timestamps=neural_data.processing['ecephys']['LFP']['LFP'].timestamps[:]
    )

    lfp = LFP(electrical_series=lfp_electrical_series)

    ecephys_module.add(lfp)


    # add RecordingSystemEvent    
    event_marker_series = TimeSeries(
        name='RecordingSystemEvents',
        data=neural_data.processing['ecephys']['RecordingSystemEvents'].data,
        unit='NA', 
        timestamps=neural_data.processing['ecephys']['RecordingSystemEvents'].timestamps,
        description='Event markers recorded by recording system'
    )

    ecephys_module.add(event_marker_series)


    # %% Add time difference data
    # load saved time difference data
    filename_ = os.path.join(root_dir,'description', 'diff_time_mean.txt')
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

    global time_difference
    time_difference = diff_time_mean


    # %% Add behavioral data 
    # load data
    filename = os.path.join(root_dir, 'formatted_data', 'continuous_behavior.nwb')
    bhv_data = NWBHDF5IO(filename, mode='r').read()

    behavior_module = nwbfile.create_processing_module(
            name="behavior", 
            description="Behavioral data. "
    )


    event_marker = bhv_data.processing['behavior']['MonkeyLogicEvents'].data[:]
    event_time = bhv_data.processing['behavior']['MonkeyLogicEvents'].timestamps[:] 


    '''---Align behavioral data time to neural data time---''' 
    event_time = event_time - time_difference
    # global marker_time
    # if 24 in event_marker:
    #     marker_time = event_time[np.array(event_marker) == 24][0] # first vicon recorded frame is 24
    # marker_time = event_time[np.array(event_marker) == 1][0] # for bohr_data_recording_0910
    '''----------------------------------------------------'''

    # write into
    behavioral_events = BehavioralEvents(name='MonkeyLogicEvents')

    behavioral_events.create_timeseries(
        name='BehaviorMarkers',
        data=event_marker,
        unit='NA',
        timestamps = event_time,
        description='Event markers recorded by PsychoPy'
    )

    behavior_module.add(behavioral_events)


    # %% Collect trial info
    bhvdf = bhv_data.trials.to_dataframe()

    '''---Align behavioral data time to neural data time---'''
    # info_time = info_time - time_difference
    bhvdf['AbsoluteTrialStartTime'] = bhvdf['AbsoluteTrialStartTime'].apply(lambda x: x - time_difference)
    '''----------------------------------------------------'''

    skip_list = ['start_time', 'stop_time'] 

    for key in bhvdf.columns:
        if key not in skip_list:
            nwbfile.add_trial_column(
                name=key,
                description=''
            )

    s = "nwbfile.add_trial(start_time=bhvdf.loc[i, 'start_time'], stop_time=bhvdf.loc[i, 'stop_time'], " + \
        ','.join(["%s=bhvdf.loc[i, '%s']" % (key, key) for key in bhvdf.columns if key not in skip_list]) + ')'
    for i in bhvdf.index: 
       exec(s)


    #%% Writing standard NWB file
    save_path = os.path.join(output_dir, "standard_data.nwb")
    if os.path.exists(save_path):
        os.remove(save_path)

    with NWBHDF5IO(os.path.join(output_dir, "standard_data.nwb"), "w") as io:
        io.write(nwbfile)


# %% Parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/interception/Abel/data_recording/20240528_Interception_001', 
                    metavar='/the/path/your/data/located/in', help='root folder')

parser.add_argument('-mp', '--map_path', 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Abel_Utah_64x4_PMd-M1-S1-A7_BlackRock.json')

# parser.add_argument('-mp', '--map_path', 
#                     default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Bohr_Utah_96x2_PMd-M1_BlackRock.json')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/interception/Abel/data_recording/20240528_Interception_001/formatted_data', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()

run(args.root, args.map_path, args.output)