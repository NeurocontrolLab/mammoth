#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 10:02:11 2024

@author: cuilab
"""

#reference: https://pynwb.readthedocs.io/en/stable/tutorials/domain/ecephys.html#sphx-glr-tutorials-domain-ecephys-py

from datetime import datetime
from uuid import uuid4
from zoneinfo import ZoneInfo 

import numpy as np

from pynwb import NWBHDF5IO, NWBFile
from pynwb.ecephys import LFP, ElectricalSeries
from pynwb.behavior import SpatialSeries
from pynwb.file import Subject

from SmartNeo.interface_layer.nwb_interface import NWBInterface
import os

from probeinterface import read_probeinterface


#%% Define path

# session_id = "20231205_interception_001"
# data_dir = "/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Qianqian/Data_recording/{}/formatted_data".format(session_id)

session_id = "20240401_centerOut_001"
data_dir = "/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Data_recording/{}/formatted_data".format(session_id)


#%% Set basic info
device_dict = {'name': 'array', 'description': 'for PMd and M1',
               'manufacturer': 'Utah'}

shank_location = ['M1']

map_path = '/AMAX/cuihe_lab/lichenyang/Neucyber-NC-2023-A-01/Spike_sorting/nezha-gai.json'
# map_path = '/AMAX/cuihe_lab/lichenyang/Neucyber-NC-2023-A-01/Qianqian/Spike_sorting/qianqian_map_zong.json'
probegroup = read_probeinterface(map_path)

subject = Subject(
    subject_id="Nezha",
    age="P6Y",
    description="Monkey",
    species="Macaca mulatta",
    sex="M",
    )

#%% Create NWB files
nwbfile = NWBFile(
    session_description="Monkey interception",
    identifier=str(uuid4()),
    session_start_time=datetime(int(session_id[:4]), 
                                int(session_id[4:6]), 
                                int(session_id[6:8]),
                                tzinfo=ZoneInfo("Asia/Shanghai")),
    experimenter=[
        "Chenyang Li",
    ],
    lab="Cui Lab",
    institution="Chinese Institute for Brain Research, Beijing",
    experiment_description="Flexible Interception",
    keywords=["ecephys", "monkey", "arm movement"],
    related_publications="",
    subject=subject
)

#%% ecephys data
# load
nwb_saver = NWBInterface()
neural_data = nwb_saver.read_nwb(filename = os.path.join(data_dir,
                                                         'neural_data.nwb'))

# set device and electrode table
device = nwbfile.create_device(name=device_dict['name'],
                               description=device_dict['description'],
                               manufacturer=device_dict['manufacturer'])

nwbfile.add_electrode_column(name="label", description="label of electrode")

# get electrode table
# nshanks = len(pd.unique(elec_info['sh']))
nshanks = len(probegroup.probes)
electrode_counter = 0
# p.device_channel_indices[0] #TODO

for ishank in range(nshanks):        
    # create an electrode group for this shank
    p = probegroup.probes[ishank]
    
    electrode_group = nwbfile.create_electrode_group(
        name="shank{}".format(ishank),
        description="electrode group for shank {}".format(ishank),
        device=device,
        location=shank_location[0]
    )

    # add electrodes to the electrode table
    for ielec in range(p.contact_positions.shape[0]): #p.device_channel_indices:
        nwbfile.add_electrode(
            x=float(p.contact_positions[ielec, 0]),
            y=float(p.contact_positions[ielec, 1]),
            group=electrode_group,
            label="shank{}elec{}".format(ishank, ielec),
            location="M1")
        
        electrode_counter += 1

# view the electrodes table in a pandas DataFrame
nwbfile.electrodes.to_dataframe()


# add spike data into nwbfile.units
kilo = [i for i in neural_data.segments if i.name=='kilosort2.5'][0]
nunits = len(kilo.spiketrains)

# nwbfile.add_unit_column(name="quality", description="sorting quality")

for iunit in [0]: #range(nunits):
    spike_times = kilo.spiketrains[iunit]
    nwbfile.add_unit(spike_times=np.array(spike_times.times), #should be array rather than quantities 
                     # quality='good',
                      waveform_mean=spike_times.description["mean_waveform"],
                     electrodes=[int(kilo.spiketrains[iunit].description["chn"])]) 

# view the unit table in a pandas DataFrame
nwbfile.units.to_dataframe()


# add LFP data
lfp_data = [i for i in neural_data.segments if i.name=='LFP'][0]
assert electrode_counter == lfp_data.irregularlysampledsignals[0].shape[1]

all_table_region = nwbfile.create_electrode_table_region(
    region=list(range(electrode_counter)), 
    description="all electrodes",
)

lfp_electrical_series = ElectricalSeries(
    name="ElectricalSeries",
    description="LFP data",
    data=lfp_data.irregularlysampledsignals[0],
    electrodes=all_table_region,
    starting_time=0.0, 
    rate=30000.0,
)

lfp = LFP(electrical_series=lfp_electrical_series)

ecephys_module = nwbfile.create_processing_module(
    name="ecephys", description="processed extracellular electrophysiology data"
)
ecephys_module.add(lfp)

del neural_data


#%% For Behavior Data
# load
trial_saver = NWBInterface()
bhv_data = trial_saver.read_nwb(filename = os.path.join(data_dir,
                                                        'trial_behavior.nwb'))

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
    name="TrialError",
    description="whether the trial was correct",
)
nwbfile.add_trial_column(
    name="Condition",
    description="task condition",
)
nwbfile.add_trial_column(
    name="TargetAngularSpeed",
    description="target angular speed",
)


ntrial = len(bhv_data.segments)
for itrial in range(ntrial):
    
    event_labels = bhv_data.segments[itrial].events[0].labels
    event_times = bhv_data.segments[itrial].events[0].times
    
    nwbfile.add_trial(
        start_time=float(event_times[event_labels==9][0]),
        stop_time=float(event_times[event_labels==18][0]),
        target_on_time=float(event_times[event_labels==3][0]) if 3 in event_labels else -1.0,
        go_cue_time=float(event_times[event_labels==4][0]) if 4 in event_labels else -1.0,
        movement_onset_time=float(event_times[event_labels==5][0]) if 5 in event_labels else -1.0,
        touch_time=float(event_times[event_labels==6][0]) if 6 in event_labels else -1.0,
        TrialError=float(bhv_data.segments[itrial].description['TrialError']),
        Condition=float(bhv_data.segments[itrial].description['Condition']),
        TargetAngularSpeed=float(bhv_data.segments[itrial].description['UserVars']['angularV'][0, 0]),
        )
nwbfile.trials.to_dataframe()


# set other behavior data
behavior_module = nwbfile.create_processing_module(
    name="behavior", description="Processed behavioral data"
)

# target position
# Row: status x Columns: objects
# [Start, TO, GO, Touch, End] x [center dot, Go-cue, touch dot, -, target dot]

position_data = np.array([pos[4, :] 
                     for i in bhv_data.segments 
                     for pos in i.irregularlysampledsignals[0]])
ts = np.array([ts for i in bhv_data.segments for ts in i.irregularlysampledsignals[0].times])
ts2 = np.array([ts for i in bhv_data.segments for ts in i.irregularlysampledsignals[1].times]) 
assert np.unique(ts==ts2)[0]==True 

target_position_spatial_series = SpatialSeries(
    name="TargetPosition",
    description="Target Position (x, y) in the touch screen.",
    data=position_data,
    timestamps=ts,
    reference_frame="(0,0) is center",
    unit='cm',
)

behavior_module.add(target_position_spatial_series)

del target_position_spatial_series

#call: nwbfile.processing["behavior"]["TargetPosition"].data[:]
#call: nwbfile.processing["behavior"]["TargetPosition"].timestamps


# touch position
position_data = np.array([pos[2, :] 
                     for i in bhv_data.segments 
                     for pos in i.irregularlysampledsignals[0]])
ts = np.array([ts for i in bhv_data.segments for ts in i.irregularlysampledsignals[0].times])
ts2 = np.array([ts for i in bhv_data.segments for ts in i.irregularlysampledsignals[1].times]) 
assert np.unique(ts==ts2)[0]==True 

touch_position_spatial_series = SpatialSeries(
    name="TouchPosition",
    description="Touch Position (x, y) in the touch screen.",
    data=position_data,
    timestamps=ts,
    reference_frame="(0,0) is center",
    unit='cm',
)

behavior_module.add(touch_position_spatial_series)

del touch_position_spatial_series

#call: nwbfile.processing["behavior"]["TouchPosition"].data[:]
#call: nwbfile.processing["behavior"]["TouchPosition"].timestamps

del bhv_data 


#%% Writing electrophysiology data
with NWBHDF5IO(os.path.join(data_dir, "data.nwb"), "w") as io:
    io.write(nwbfile)


#%% Reading electrophysiology data
io = NWBHDF5IO(os.path.join(data_dir, "data.nwb"), mode="r")
read_nwbfile = io.read()
    
# read electrodes
read_nwbfile.electrodes.to_dataframe()

# read spike time
read_nwbfile.units.to_dataframe()
   
# read LFP
# read_nwbfile.processing["ecephys"]["LFP"]["ElectricalSeries"].data[:]

import pynapple as nap

nwb = nap.NWBFile(read_nwbfile)
print(nwb)

all_trials = nwb["trials"]
corr_trials = all_trials[(all_trials["TrialError"]==0)]
print(corr_trials)
