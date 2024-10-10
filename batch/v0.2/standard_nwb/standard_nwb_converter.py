#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 10:02:11 2024

@author: cuilab
"""

#reference: https://pynwb.readthedocs.io/en/stable/tutorials/domain/ecephys.html#sphx-glr-tutorials-domain-ecephys-py

from datetime import datetime
from zoneinfo import ZoneInfo 
import numpy as np
from pynwb import NWBHDF5IO, NWBFile, TimeSeries
from pynwb.ecephys import LFP, ElectricalSeries
from pynwb.behavior import BehavioralEvents,BehavioralTimeSeries,SpatialSeries
from pynwb.file import Subject
from SmartNeo.interface_layer.nwb_interface import NWBInterface
import os
from probeinterface import Probe, ProbeGroup
import argparse
import quantities as pq
import pandas as pd
import json
from pynwb.core import DynamicTable

#%% parse the input arguments
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

#%% construct probe
map_path = args.map_path
with open(map_path,'r') as f:
    probe_info = f.readlines()[14::]

probe_info = [i.split() for i in probe_info][0:-1]

array_name = ['elec'] if '-' not in probe_info[0][-1]\
    else list(set([i[-1].split('-')[0] for i in probe_info]))

array_name.sort()

nchannels_per_array = len(probe_info)
electrode_counter = []
probegroup = ProbeGroup()

for array in array_name:
    # create an electrode group for this shank
    
    probe_2d = Probe(ndim=2, si_units='um')
    probe_2d.annotate(
        name = array, 
        manufacturer="blackrock microsystem",
        escription = 'one 96 Utah array'
        )
    positions = []
    device_channel = []
    # test2 = []
    # add electrodes to the electrode table
    for ielec in probe_info:
        if array not in ielec[-1]:
            continue
    
        positions.append([float(ielec[0])*400, float(ielec[1])*400])
        device_channel.append((ord(ielec[2])-65)*32+int(ielec[3])-1)
        # test2.append(ielec[2]+ielec[3])
        
    
    probe_2d.set_contacts(positions=np.array(positions), 
                          shapes='circle', 
                          shape_params={'radius': 20})
    probe_2d.set_device_channel_indices(device_channel)
    probe_2d.create_auto_shape(probe_type='tip')
    probegroup.add_probe(probe_2d)


#%% Set basic info
root_directory = os.path.join(raw_dirname,'bhv')
file_pattern = [i for i in os.listdir(root_directory) if ('csv' in i) and ('meta' in i)][0]
meta_dict_pd = pd.read_csv(os.path.join(root_directory, file_pattern)).to_dict()
meta_dict = {}
for i,j in zip(meta_dict_pd['Key'],meta_dict_pd['Value']):
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
    datetime.strptime(nwb_dict['session_start_time'], '%Y/%m/%d %H:%M').replace(tzinfo=ZoneInfo('Asia/Shanghai'))
nwb_dict['subject'] = subject
nwb_dict['keywords'] = ["ecephys", "monkey", "motor control"]
del nwb_dict['experiment_name']
#%% Create NWB files
nwbfile = NWBFile(
    **nwb_dict
)

#%% add electrode
device_dict = {'name': 'Utah array', 'description': '2 x 96',
               'manufacturer': 'blackrock microsystem'}

# set device and electrode table
device = nwbfile.create_device(name=device_dict['name'],
                               description=device_dict['description'],
                               manufacturer=device_dict['manufacturer'])

# nwbfile.add_electrode_column(name="2 x 96 Utah array", description="multi-channel electrode array")

# get electrode table
# nshanks = len(pd.unique(elec_info['sh']))
nshanks = len(probegroup.probes)
# p.device_channel_indices[0] #TODO
shank_location = ['PMd', 'M1']
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
    for ielec in range(p.contact_positions.shape[0]): #p.device_channel_indices:
        nwbfile.add_electrode(
            x=float(p.contact_positions[ielec, 0]),
            y=float(p.contact_positions[ielec, 1]),
            group=electrode_group,
            id = p.device_channel_indices[ielec],
            location=shank_location[ishank])

# view the electrodes table in a pandas DataFrame
nwbfile.electrodes.to_dataframe()
unit_dict = {'s': pq.s}
filename_ = os.path.join(raw_dirname,'description', 'diff_time_mean.txt')
with open(filename_, 'r') as file:
    dt_str = file.read()  # 读取整个文件内容
unit = unit_dict[dt_str.split(' ')[1][0:-1]]
diff_time_mean = float(dt_str.split(' ')[0])

behavior_ecephys_module = nwbfile.create_processing_module(name='behavior_ecephys_analysis', 
                                           description='Quality control and pre-analysis results')

# 创建 DynamicTable 保存时间差
time_diff_table = DynamicTable(name='TimeDifference',
                               description='Time difference between behavior and ecephys. (behavior - ecephys)')

# 添加一行，将时间差保存到表中

time_diff_table.add_column(
    name='time_difference',
    description='Time difference between behavior and ecephys (behavior - ecephys). The unit is second.'
)

time_diff_table.add_row(id=0, time_difference=diff_time_mean)

# 将 DynamicTable 添加到 ProcessingModule
behavior_ecephys_module.add(time_diff_table)

#%% add ecephys data
formated_data_path = os.path.join(raw_dirname,'formatted_data')

nwb_saver = NWBInterface()
neural_data = nwb_saver.read_nwb(filename = os.path.join(formated_data_path,'neural_data.nwb'))

nwb_saver = NWBInterface()
bhv_data = nwb_saver.read_nwb(filename = os.path.join(formated_data_path, 'continuous_behavior.nwb'))

# add spike data sorted by kilosort into nwbfile.units
kilo = [i for i in neural_data.segments if i.name=='kilosort2.5'][0]
nunits = len(kilo.spiketrains)

nwbfile.add_unit_column(name="chn_id", description="source channel of unit")
nwbfile.add_unit_column(name="sorter", description="sorting method")
nwbfile.add_unit_column(name="sorting_info", description="metadata of sorting results")
nwbfile.add_unit_column(name="time_unit", description="time scale")


b = kilo.spiketrains[0].description["chn_meta"]
for iunit in range(nunits):
    spike_times = kilo.spiketrains[iunit]
    nwbfile.add_unit(spike_times=np.array(spike_times.times.rescale(pq.s).magnitude), #should be array rather than quantities 
                     # quality='good',
                     id = iunit,
                     chn_id = int(spike_times.description["chn"]),
                     waveform_mean=spike_times.description["mean_waveform"],
                     electrode_group = ETR_list[int(spike_times.description["electrode"])],
                     time_unit='seconds',
                     sorting_info = json.dumps(spike_times.description["chn_meta"]),
                     # sorting_info = spike_times.description["chn_meta"],
                     sorter = 'kilosort2.5')

# add TCR into nwbfile.units
TCR = [i for i in neural_data.segments if i.name=='TCR'][0]
for iunit in range(nunits):
    spike_times = TCR.spiketrains[iunit]
    nwbfile.add_unit(spike_times=np.array(spike_times.times.rescale(pq.s).magnitude), #should be array rather than quantities 
                     # quality='good',
                     id = int(iunit+nunits),
                     chn_id = int(spike_times.description["chn"]),
                     waveform_mean=spike_times.description["mean_waveform"],
                     electrode_group = ETR_list[int(spike_times.description["electrode"])],
                     time_unit='seconds',
                     sorting_info = 'NA',
                     sorter = 'TCR')
    
# view the unit table in a pandas DataFrame
nwbfile.units.to_dataframe()

electrode_index = nwbfile.electrodes.id.data
# add LFP data
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
    name="ecephys", description="processed extracellular electrophysiology data"
)
ecephys_module.add(lfp)

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
    description='Event markers recording by recording system'
)

behavioral_events = BehavioralEvents(name='Recording system Events')
behavioral_events.add_timeseries(event_marker_series)
ecephys_module.add(behavioral_events)

# maybe add sbp and raw data
del neural_data

#%% For Behavior Data
bhv_par = [i for i in bhv_data.segments if i.name=='Task parameters'][0].events[0] # get behavior parameters

# build behavior module for saving behavior data
behavior_module = nwbfile.create_processing_module(
        name="behavior", description="Processed behavioral data. The task parameters are {}".format(bhv_par.labels[0])
) # save task parameters

# convert vicon motion data
vm = [i for i in bhv_data.segments if i.name=='Vicon motion'][0].irregularlysampledsignals[0]

vicon_pos_series = SpatialSeries(
    name='ViconMotion',
    data=np.array(vm),
    unit='mm',
    reference_frame = 'zero-position was set by vicon before every experiment. The position was closed to the monkey sitting position',
    timestamps=vm.times.rescale(pq.s).magnitude,
    description='finger position recorded by vicon motion system. 24 marker represents recording start. Three columns are x, y, z position'
)

behavior_module.add(vicon_pos_series)

# add behavior event marker
bhv_event = [i for i in bhv_data.segments if i.name=='Event marker'][0].events[0]
marker = [json.loads(i) for i in bhv_event.labels]
num_marker = []
event_marker = []
marker_time = []

for i,j in zip(bhv_event.labels,bhv_event.times):
    event_recorded = json.loads(i)
    if isinstance(event_recorded,dict):
        num_marker.append(event_recorded['Marker'])
        event_marker.append(event_recorded['Event'])
        marker_time.append(j.rescale(pq.s).magnitude)
marker_time = np.array(marker_time).squeeze()

behavioral_events = BehavioralEvents(name='Psychopy Events')

behavioral_events.create_timeseries(
    name='BehaviorMarkers',
    data=num_marker,
    unit='NA',
    timestamps = marker_time,
    description='marker numbers recorded by psychopy'
)

behavioral_events.create_timeseries(
    name='BehaviorLabels',
    data=event_marker,
    unit='NA',
    timestamps = marker_time,
    description='marker labels recorded by psychopy'
)

# behavioral_events.events = events_table
behavior_module.add(behavioral_events)

# add trial info
bhv_trial_info = [i for i in bhv_data.segments if i.name=='Trial info'][0].events[0]
data = [json.loads(i)['Trial info'].copy() for i in bhv_trial_info.labels]
info_time = np.array(bhv_trial_info.times.rescale(pq.s).magnitude).squeeze()

behavioral_events = BehavioralEvents(name='TrialInfo')

behavioral_events.create_timeseries(
    name='TrialInfoTable',
    data=[json.dumps(i) for i in data],
    unit='NA',
    timestamps = info_time,
    description='JSON encoded data. Timestamps of trial start and trial end. Summary info of each trial'
)

behavior_module.add(behavioral_events)


# add object info in each frame
def unpack_objects(d):
    objects_list = d.get('objects', [])
    for idx, obj in enumerate(objects_list, 1):
        d[f'object{idx}'] = obj
    del d['objects']
    return d

bhv_frame = [i for i in bhv_data.segments if i.name=='Frame'][0].events[0]
bhv_frame_ = [unpack_objects(json.loads(i)) for i in bhv_frame.labels]
data = bhv_frame_
info_time = np.array(bhv_frame.times.rescale(pq.s).magnitude).squeeze()

# nwbfile.add_acquisition(frame_table)
behavioral_times = BehavioralTimeSeries(name='FrameInfo')

behavioral_times.create_timeseries(
    name='FrameInfoTable',
    data=[json.dumps(i) for i in data],
    unit='NA',
    timestamps = info_time,
    description='JSON encoded data. There are three objects: center is the object that appears in the center; when displayed, the monkey needs to place its hand at the central position. Target is the goal that the monkey needs to touch. Feedback represents the position on the screen where the macaque clicks and will display "wrong" or "right" to indicate whether the action was incorrect or correct. Other attributes are measured in centimeters. If the content is empty, it means that the cursor did not appear.'
)

behavior_module.add(behavioral_times)

del bhv_data

#%% Writing electrophysiology data
save_path = os.path.join(formated_data_path, "standard_data.nwb")
if os.path.exists(save_path):
    os.remove(save_path)

with NWBHDF5IO(os.path.join(formated_data_path, "standard_data.nwb"), "w") as io:
    io.write(nwbfile)