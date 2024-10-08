#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 10:02:11 2024

@author: cuilab
"""

#reference: https://pynwb.readthedocs.io/en/stable/tutorials/domain/ecephys.html#sphx-glr-tutorials-domain-ecephys-py

import os
import json
import numpy as np
import pandas as pd
import argparse
from datetime import datetime
from zoneinfo import ZoneInfo 
import quantities as pq

from pynwb import NWBHDF5IO, NWBFile, TimeSeries
from pynwb.ecephys import LFP, ElectricalSeries
from pynwb.behavior import BehavioralEvents,BehavioralTimeSeries,SpatialSeries
from pynwb.file import Subject
from pynwb.core import DynamicTable

from SmartNeo.interface_layer.nwb_interface import NWBInterface
from probeinterface import read_probeinterface


# %% Parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument("-f", "--file", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20240925_interception_004', 
                    metavar='/the/path/your/data/located/in', help='input folder')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20240925_interception_004', 
                    metavar='/the/path/you/want/to/save', help='output folder')

parser.add_argument('-mp', '--map_path', 
                    default='/home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Spike_sorting/SN+11386-000049.cmp')

args = parser.parse_args()
raw_dirname = args.file
map_path = args.map_path
formated_data_path = os.path.join(raw_dirname,'formatted_data')


# unit_dict = {'s': pq.s}
# unit = unit_dict[dt_str.split(' ')[1][0:-1]]


# %% Set basic info and create NWB file
root_directory = os.path.join(raw_dirname,'bhv')
file_pattern = [i for i in os.listdir(root_directory) if ('csv' in i) and ('meta' in i)][0]
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

session_id = os.path.basename(raw_dirname)
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


nwbfile = NWBFile(
    **nwb_dict
)


# %% Add ecephys data
# Load probe
probegroup = read_probeinterface(map_path)

device_dict = {'name': 'Utah array', 'description': '96 x 2',
               'manufacturer': 'BlackRock Microsystem'} #TODO

# set device and electrode table
device = nwbfile.create_device(name=device_dict['name'],
                               description=device_dict['description'],
                               manufacturer=device_dict['manufacturer'])

# get electrode table
nshanks = len(probegroup.probes)
shank_location = ['PMd', 'M1']  #TODO
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

# load data
nwb_saver = NWBInterface()
neural_data = nwb_saver.read_nwb(filename = os.path.join(formated_data_path,'neural_data.nwb'))


# define unit columns
nwbfile.add_unit_column(name="chn_id", 
                        description="source channel of unit")
nwbfile.add_unit_column(name="sorter", 
                        description="sorting method")
nwbfile.add_unit_column(name="sorting_info", 
                        description="metadata of sorting results")
nwbfile.add_unit_column(name="time_unit", 
                        description="time scale")

# add spike data sorted by kilosort into nwbfile.units
kilo = [i for i in neural_data.segments if i.name=='kilosort2.5'][0]
nunits = len(kilo.spiketrains)

# b = kilo.spiketrains[0].description["chn_meta"]
for iunit in range(nunits):
    spike_times = kilo.spiketrains[iunit]
    nwbfile.add_unit(spike_times=np.array(spike_times.times.rescale(pq.s).magnitude), #should be array rather than quantities 
                     id=iunit,
                     chn_id=int(spike_times.description["chn"]),
                     waveform_mean=spike_times.description["mean_waveform"],
                     electrode_group=ETR_list[int(spike_times.description["electrode"])],
                     electrode=int(spike_times.description["chn_meta"]),
                     time_unit='seconds',
                     sorting_info=json.dumps(spike_times.description["chn_meta"]),
                     # sorting_info=spike_times.description["chn_meta"],
                     sorter = 'kilosort2.5')

# add TCR into nwbfile.units
TCR = [i for i in neural_data.segments if i.name=='TCR'][0]
for iunit in range(nunits):
    spike_times = TCR.spiketrains[iunit]
    nwbfile.add_unit(spike_times=np.array(spike_times.times.rescale(pq.s).magnitude), #should be array rather than quantities 
                     id=int(iunit+nunits),
                     chn_id=int(spike_times.description["chn"]),
                     waveform_mean=spike_times.description["mean_waveform"],
                     electrode_group=ETR_list[int(spike_times.description["electrode"])],
                     electrode=int(spike_times.description["chn_meta"]),
                     time_unit='seconds',
                     sorting_info = 'NA',
                     sorter = 'TCR')
    
# view the unit table in a pandas DataFrame
# nwbfile.units.to_dataframe()

# add LFP data
electrode_index = nwbfile.electrodes.id.data
lfp_data = [i for i in neural_data.segments if i.name=='LFP'][0]
lfp_array = np.array(lfp_data.irregularlysampledsignals[0])[:,electrode_index]

all_table_region = nwbfile.create_electrode_table_region(
    region=list(range(len(electrode_index))), 
    description="all electrodes",
)

lfp_electrical_series = ElectricalSeries(
    name="LFP",
    description="LFP data",
    data=lfp_array,
    electrodes=all_table_region,
    timestamps=lfp_data.irregularlysampledsignals[0].times.rescale(pq.s).magnitude
)

lfp = LFP(electrical_series=lfp_electrical_series)

ecephys_module = nwbfile.create_processing_module(
    name="ecephys", 
    description="processed extracellular electrophysiology data"
)
ecephys_module.add(lfp)

# add Recording System Event
event = [i for i in neural_data.segments if i.name=='RecordingSystemEvent'][0].events[0]
event_label = np.array(event.labels)-65280
event_times = event.times.rescale(pq.s).magnitude
event_times = event_times[event_label!=0]
event_label = event_label[event_label!=0]
event_marker_series = TimeSeries(
    name='RecordingSystemEvents',
    data=event_label,
    unit='NA', 
    timestamps=event_times,
    description='Event markers recorded by recording system'
)

behavioral_events = BehavioralEvents(name='Recording system Events')
behavioral_events.add_timeseries(event_marker_series)
ecephys_module.add(behavioral_events)

# maybe add sbp and raw data
del neural_data


#%% Add time difference data
# load saved time difference data
filename_ = os.path.join(raw_dirname,'description', 'diff_time_mean.txt')
with open(filename_, 'r') as file:
    dt_str = file.read()  

diff_time_mean = float(dt_str.split(' ')[0])

# # add time difference
# behavior_ecephys_module = nwbfile.create_processing_module(
#     name='behavior_ecephys_analysis', 
#     description='Quality control and pre-analysis results')

# time_diff_table = DynamicTable(
#     name='TimeDifference',
#     description='Time difference between behavior and ecephys')

# time_diff_table.add_column(
#     name='time_difference',
#     description='Behavior - Ecephys (s)'
# )

# time_diff_table.add_row(id=0, time_difference=diff_time_mean)

# behavior_ecephys_module.add(time_diff_table)
global time_difference
time_difference=diff_time_mean

# %% Add behavioral data 
# load data
nwb_saver = NWBInterface()
bhv_data = nwb_saver.read_nwb(filename=os.path.join(formated_data_path, 
                                                    'continuous_behavior.nwb'))

# build behavior module for saving behavior data
bhv_par = [i for i in bhv_data.segments 
            if i.name=='Task parameters'][0].events[0] # get behavior parameters
behavior_module = nwbfile.create_processing_module(
        name="behavior", 
        description="Processed behavioral data. The task parameters are {}".format(bhv_par.labels[0])
)

# add behavior event marker
bhv_event = [i for i in bhv_data.segments 
             if i.name=='Event marker'][0].events[0]
# marker = [json.loads(i) for i in bhv_event.labels]

event_marker = []
event_label = []
event_time = []
for i,j in zip(bhv_event.labels, bhv_event.times):
    event_recorded = json.loads(i)
    if isinstance(event_recorded, dict):
        event_marker.append(event_recorded['Marker'])
        event_label.append(event_recorded['Event'])
        event_time.append(j.rescale(pq.s).magnitude)
event_time = np.array(event_time).squeeze()

'''---Align behavioral data time to neural data time---''' 
event_time = event_time - time_difference
global marker_time
marker_time = event_time[np.array(event_marker) == 24][0] # first vicon recorded frame is 24
'''----------------------------------------------------'''

# write into
behavioral_events = BehavioralEvents(name='PsychoPy Events')

behavioral_events.create_timeseries(
    name='BehaviorMarkers',
    data=event_marker,
    unit='NA',
    timestamps = event_time,
    description='Event markers recorded by PsychoPy'
)

behavioral_events.create_timeseries(
    name='BehaviorLabels',
    data=event_label,
    unit='NA',
    timestamps = event_time,
    description='Event labels recorded by PsychoPy'
)
behavior_module.add(behavioral_events)


# add Vicon data
vm = [i for i in bhv_data.segments 
      if i.name=='Vicon motion'][0].irregularlysampledsignals[0]

vm_time = vm.times.rescale(pq.s).magnitude

'''---Align Vicon Motion data time to aligned behavioral data time---'''
vicon_start_time = vm_time[0]
marker_vicon_difference = marker_time - vicon_start_time
vm_time = vm_time + marker_vicon_difference
'''------------------------------------------------------------------'''

vicon_pos_series = SpatialSeries(
    name='ViconMotion',
    data=np.array(vm),
    unit='mm',
    reference_frame='Zero-position was set by Vicon before every experiment.\
        The position was closed to the monkey sitting position',
    timestamps=vm_time, #vm.times.rescale(pq.s).magnitude,
    description='Finger position recorded by Vicon Motion System.\
        The Marker 24 represents the start of recording. \
        The three columns are x, y, and z positions'
)
behavior_module.add(vicon_pos_series)


# add object info in each frame
def unpack_objects(d):
    objects_list = d.get('objects', [])
    for idx, obj in enumerate(objects_list, 1):
        d[f'object{idx}'] = obj
    del d['objects']
    return d

bhv_frame = [i for i in bhv_data.segments if i.name=='Frame'][0].events[0]
bhv_frame_ = [unpack_objects(json.loads(i)) for i in bhv_frame.labels]
frame_time = np.array(bhv_frame.times.rescale(pq.s).magnitude).squeeze()

'''---Align behavioral data time to neural data time---'''
frame_time = frame_time - time_difference
'''----------------------------------------------------'''

behavioral_times = BehavioralTimeSeries(name='FrameInfo')

behavioral_times.create_timeseries(
    name='FrameInfoTable',
    data=[json.dumps(i) for i in bhv_frame_],
    unit='NA',
    timestamps = frame_time,
    description='JSON encoded data.\
        There are three objects: \
            center is the object that appears in the center; \
            when displayed, the monkey needs to place its hand at the central position.\
            Target is the goal that the monkey needs to touch. \
            Feedback represents the position on the screen where the macaque \
                clicks and will display "wrong" or "right" to indicate whether \
                the action was incorrect or correct. Other attributes are \
                measured in centimeters. If the content is empty, it means \
                that the cursor did not appear.'
)
behavior_module.add(behavioral_times)


# %% Collect trial info
bhv_trial_info = [i for i in bhv_data.segments if i.name=='Trial info'][0].events[0]
info = [json.loads(i)['Trial info'].copy() for i in bhv_trial_info.labels]
info_time = np.array(bhv_trial_info.times.rescale(pq.s).magnitude).squeeze()

'''---Align behavioral data time to neural data time---'''
info_time = info_time - time_difference
'''----------------------------------------------------'''

# add as time series
behavioral_events = BehavioralEvents(name='TrialInfo')

behavioral_events.create_timeseries(
    name='TrialInfoTable',
    data=[json.dumps(i) for i in info],
    unit='NA',
    timestamps = info_time,
    description='JSON encoded data. Timestamps of trial start and trial end. Summary info of each trial'
)

behavior_module.add(behavioral_events)

trial_ts_view = pd.DataFrame.from_dict(info)

# add as trials
trial_view = pd.DataFrame(columns=['start_time', 'stop_time', 'target_on_time',
                                   'go_cue_time', 'movement_onset_time',
                                   'touch_time', 'init_pos', 'target_speed',
                                   'result_marker','result_label',
                                   ])

key_label_dict = {'target_on_time': 'Target on', 'go_cue_on_time': 'Gosignal on',
                  'movement_onset_time': 'Movement onset', 
                  'touch_time': 'Target touched'}
for itrial in range(max(trial_ts_view['trial_number'])):
    
    start = info_time[
        trial_ts_view[(trial_ts_view['status']=='start') &
                      (trial_ts_view['trial_number']==itrial)].index][0]
    end = info_time[
        trial_ts_view[(trial_ts_view['status']=='end') &
                      (trial_ts_view['trial_number']==itrial)].index][0]
    
    trial_view.loc[itrial, 'start_time'] = start
    trial_view.loc[itrial, 'stop_time'] = end 
    
    for k, v in key_label_dict.items():
        ts_ls = [event_time[ind] for ind, i in enumerate(event_label) 
                 if (i==v) and 
                 ((event_time[ind]>start)&(event_time[ind]<end))]
    
        if len(ts_ls)==1:
            trial_view.loc[itrial, k] = ts_ls[0] 

    trial_view.loc[itrial, 'init_pos_x'] = \
        trial_ts_view.loc[
        (trial_ts_view['trial_number']==itrial)&
        (trial_ts_view['status']=='end'), 'init_pos'].tolist()[0][0]
        
    trial_view.loc[itrial, 'init_pos_y'] = \
        trial_ts_view.loc[
        (trial_ts_view['trial_number']==itrial)&
        (trial_ts_view['status']=='end'), 'init_pos'].tolist()[0][1]
    
    trial_view.loc[itrial, 'result_label'] = trial_ts_view.loc[
        (trial_ts_view['trial_number']==itrial)&
        (trial_ts_view['status']=='end'), 'wrong'].values[0]['type']
    
    trial_view.loc[itrial, 'result_marker'] = trial_ts_view.loc[
        (trial_ts_view['trial_number']==itrial)&
        (trial_ts_view['status']=='end'), 'wrong'].values[0]['number']

trial_view['init_pos'] = [[x, y] for x, y in 
                          zip(trial_view['init_pos_x'], 
                              trial_view['init_pos_y'])]

assert trial_view['start_time']==[event_time[i] for i in range(len(event_time))
                                    if event_label[i]=='Center on'][0]
aaa = trial_view['start_time'].tolist()
bbb = [event_time[i] for i in range(len(event_time))
                                    if event_label[i]=='Center on']
[i for i in range(len(bbb)) if (aaa[i]-bbb[i])>1e-6]

# set trials
nwbfile.add_trial_column(
    name="target_on_time",
    description="when target appeared",
)
nwbfile.add_trial_column(
    name="go_cue_time",
    description="when GO-cue appeared",
)
nwbfile.add_trial_column(
    name="movement_onset_time",
    description="when hands-off",
)
nwbfile.add_trial_column(
    name="touch_time",
    description="when touch screen",
)
nwbfile.add_trial_column(
    name="result_marker",
    description="Marker for behavioral result",
)
nwbfile.add_trial_column(
    name="result_label",
    description="Label for behavioral result",
)
nwbfile.add_trial_column(
    name="init_pos",
    description="Initial position for the moving target",
)

ntrial = len(trial_view)
for itrial in range(ntrial): 
    nwbfile.add_trial(
        start_time=trial_view.loc[itrial, 'start_time'],
        stop_time=trial_view.loc[itrial, 'stop_time'],
        target_on_time=trial_view.loc[itrial, 'target_on_time'],
        go_cue_time=trial_view.loc[itrial, 'go_cue_time'],
        movement_onset_time=trial_view.loc[itrial, 'movement_onset_time'],
        touch_time=trial_view.loc[itrial, 'touch_time'],
        result_marker=trial_view.loc[itrial, 'result_marker'],
        result_label=trial_view.loc[itrial, 'result_label'],
        init_pos=trial_view.loc[itrial, 'init_pos'],
        )
nwbfile.trials.to_dataframe()


del bhv_data


#%% Writing standard NWB file
save_path = os.path.join(formated_data_path, "standard_data.nwb")
if os.path.exists(save_path):
    os.remove(save_path)

with NWBHDF5IO(os.path.join(formated_data_path, "standard_data.nwb"), "w") as io:
    io.write(nwbfile)