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

from SmartNeo.interface_layer.nwb_interface import NWBInterface
from probeinterface import read_probeinterface


#%%

def run(root_dir, map_path, output_dir):
   
    # unit_dict = {'s': pq.s}
    # unit = unit_dict[dt_str.split(' ')[1][0:-1]]

    # %% Set basic info and create NWB file
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


    nwbfile = NWBFile(
        **nwb_dict
    )


    # %% Add ecephys data
    # Load probe
    probegroup = read_probeinterface(map_path)

    device_dict = {'name': 'Utah array', 'description': '64 x 4',
                   'manufacturer': 'BlackRock Microsystem'} #TODO

    # set device and electrode table
    device = nwbfile.create_device(name=device_dict['name'],
                                description=device_dict['description'],
                                manufacturer=device_dict['manufacturer'])

    # get electrode table
    # nshanks = len(probegroup.probes)
    shank_location = ['PMd', 'M1', 'S1', 'A7']  #TODO
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
    neural_data = nwb_saver.read_nwb(filename = os.path.join(root_dir, 'formatted_data', 'neural_data.nwb'))


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
    kilo = [i for i in neural_data.segments if 'kilosort2.5' in i.name][0]
    nunits = len(kilo.spiketrains)

    # b = kilo.spiketrains[0].description["chn_meta"]
    for iunit in range(nunits):
        spike_times = kilo.spiketrains[iunit]
        nwbfile.add_unit(spike_times=np.array(spike_times.times.rescale(pq.s).magnitude), #should be array rather than quantities 
                        id=iunit,
                        chn_id=int(spike_times.description["chn"]),
                        waveform_mean=spike_times.description["mean_waveform"],
                        electrode_group=ETR_list[int(spike_times.description["electrode"])],
                        time_unit='seconds',
                        sorting_info=json.dumps(spike_times.description["chn_meta"]),
                        # sorting_info=spike_times.description["chn_meta"],
                        sorter = 'kilosort2.5')

    # add TCR into nwbfile.units
    TCR = [i for i in neural_data.segments if i.name=='TCR'][0]
    for iunit in range(len(TCR.spiketrains)):
        spike_times = TCR.spiketrains[iunit]
        nwbfile.add_unit(spike_times=np.array(spike_times.times.rescale(pq.s).magnitude), #should be array rather than quantities 
                        id=int(iunit+nunits),
                        chn_id=int(spike_times.description["chn"]),
                        waveform_mean=spike_times.description["mean_waveform"],
                        electrode_group=ETR_list[int(spike_times.description["electrode"])],
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


    # %% Add time difference data
    # load saved time difference data
    filename_ = os.path.join(root_dir,'description', 'diff_time_mean.txt')
    assert os.path.exists(filename_), "Please run time consistency check first!"
    with open(filename_, 'r') as file:
        dt_str = file.read()  

    diff_time_mean = float(dt_str.split(' ')[0])

    # add time difference
    from pynwb.core import DynamicTable
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
    time_difference=diff_time_mean

    # %% Add behavioral data 
    # load data
    nwb_saver = NWBInterface()
    bhv_data = nwb_saver.read_nwb(filename=os.path.join(root_dir, 'formatted_data', 
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

    if len(vm.times)>0:
        
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
        del vm, vm_time

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
    del bhv_frame

    '''---Align behavioral data time to neural data time---'''
    frame_time = frame_time - time_difference
    '''----------------------------------------------------'''

    behavioral_times = BehavioralTimeSeries(name='FrameInfo')

    behavioral_times.create_timeseries(
        name='FrameInfoTable',
        data=[json.dumps(i) for i in bhv_frame_],
        unit='NA',
        timestamps=frame_time,
        description='JSON encoded data. There are three objects: "center" is the object that appears in the center; when displayed, the monkey needs to place its hand at the central position. "Target" is the goal that the monkey needs to touch. "Feedback" represents the position on the screen where the macaque clicks and will display "wrong" or "right" to indicate whether the action was incorrect or correct. Other attributes are measured in centimeters. If the content is empty, it means that the cursor did not appear.'
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
                                    'go_cue_on_time', 'movement_onset_time',
                                    'touch_time', 'feedback_on_time',
                                    'target_pos_init', 'target_speed',
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
        
        ts_ls = [event_time[ind] for ind, i in enumerate(event_label) 
                if (i in ['Right feedback on', 'Wrong feedback on']) and 
                ((event_time[ind]>start)&(event_time[ind]<end))]

        if len(ts_ls)==1:
            trial_view.loc[itrial, 'feedback_on_time'] = ts_ls[0] 

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

    trial_view['target_pos_init'] = [[x, y] for x, y in 
                            zip(trial_view['init_pos_x'], 
                                trial_view['init_pos_y'])]

    # assert trial_view['start_time']==[event_time[i] for i in range(len(event_time))
    #                                     if event_label[i]=='Center on'][0]
    # aaa = trial_view['start_time'].tolist()
    # bbb = [event_time[i] for i in range(len(event_time))
    #                                     if event_label[i]=='Center on']
    # [i for i in range(len(bbb)) if (aaa[i]-bbb[i])>1e-6]

    frame_view = pd.DataFrame.from_dict(bhv_frame_)
    frame_view['time'] = frame_time

    pos_df = pd.DataFrame(columns=['target_pos_go_x', 'target_pos_go_y'
                                'target_pos_touch_x', 'target_pos_touch_y',
                                'feedback_pos_touch_x' 'feedback_pos_touch_y'])

    for itrial in range(len(trial_view.index)):
        if not pd.isna(trial_view.loc[itrial, 'go_cue_on_time']):        
            
            pls1 = frame_view.loc[
                abs(frame_view['time']-trial_view.loc[itrial, 'go_cue_on_time']).nsmallest(1).index,
                'object2'].values[0]['pos']
            pos_df.loc[itrial, 'target_pos_go_x'] = pls1[0]
            pos_df.loc[itrial, 'target_pos_go_y'] = pls1[1]
        else:
            pos_df.loc[itrial, 'target_pos_go_x'] = np.nan
            pos_df.loc[itrial, 'target_pos_go_y'] = np.nan

        if not pd.isna(trial_view.loc[itrial, 'touch_time']):
            
            pls2 =  frame_view.loc[
                abs(frame_view['time']-trial_view.loc[itrial, 'touch_time']).nsmallest(1).index,
                'object2'].values[0]['pos']
            
            pos_df.loc[itrial, 'target_pos_touch_x'] = pls2[0]
            pos_df.loc[itrial, 'target_pos_touch_y'] = pls2[1]
        else:
            pos_df.loc[itrial, 'target_pos_touch_x'] = np.nan
            pos_df.loc[itrial, 'target_pos_touch_y'] = np.nan
        
        if not pd.isna(trial_view.loc[itrial, 'feedback_on_time']):
            pls3_candidate =  frame_view.loc[
                abs(frame_view['time']-trial_view.loc[itrial, 'feedback_on_time']).nsmallest(5).index,
                'object3']
            pls3 = [i for i in pls3_candidate if not pd.isna(i)][0]['pos']
            
            pos_df.loc[itrial, 'feedback_pos_touch_x'] = pls3[0]
            pos_df.loc[itrial, 'feedback_pos_touch_y'] = pls3[1]
        else:
            pos_df.loc[itrial, 'feedback_pos_touch_x'] = np.nan
            pos_df.loc[itrial, 'feedback_pos_touch_y'] = np.nan
            
    trial_view['target_pos_go'] = [[x, y] 
                                for x, y in zip(pos_df['target_pos_go_x'], 
                                                pos_df['target_pos_go_y'])]

    trial_view['target_pos_touch'] = [[x, y] 
                                    for x, y in zip(pos_df['target_pos_touch_x'], 
                                                    pos_df['target_pos_touch_y'])]

    trial_view['feedback_pos'] = [[x, y] 
                                for x, y in zip(pos_df['feedback_pos_touch_x'], 
                                                pos_df['feedback_pos_touch_y'])]


    all_w = eval(bhv_par.labels[0])['all_w']
    for itrial in range(len(trial_view.index)):
        
        if not pd.isna(trial_view.loc[itrial, 'target_on_time']):
            trial_frame = \
                frame_view[(frame_view['time']>trial_view.loc[itrial, 'target_on_time'])&
                    (frame_view['time']<trial_view.loc[itrial, 'stop_time'])]
            center = trial_frame['object1'].values[0]['pos']
            
            if len(trial_frame)>1:
                target_frame1 = trial_frame['object2'].values[0]['pos']
                target_frame2 = trial_frame['object2'].values[1]['pos']
            
                a1 = np.mod(np.arctan2(target_frame1[1]-center[1], target_frame1[0]-center[0]),
                            np.pi*2)
                a2 = np.mod(np.arctan2(target_frame2[1]-center[1], target_frame2[0]-center[0]),
                            np.pi*2)
                
                w = abs((a2-a1)/(trial_frame['time'].values[1] - trial_frame['time'].values[0])
                        /np.pi*180)*np.sign(np.cross(target_frame1, target_frame2))
            else:
                w = np.array([0])
            
            w_cond = all_w[np.argmin(abs(all_w-w))]
        
            trial_view.loc[itrial, 'target_speed'] = w_cond
        
        else:
            trial_view.loc[itrial, 'target_speed'] = 0


    # set trials
    nwbfile.add_trial_column(
        name="target_on_time",
        description="when target appeared",
    )
    nwbfile.add_trial_column(
        name="go_cue_on_time",
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
        name="target_pos_init",
        description="Initial position for the moving target",
    )
    nwbfile.add_trial_column(
        name="target_pos_go",
        description="The position of the moving target at Go-signal on",
    )
    nwbfile.add_trial_column(
        name="target_pos_touch",
        description="The position of the moving target at Target touched",
    )
    nwbfile.add_trial_column(
        name="feedback_pos",
        description="The feedback (touch) position of the hand",
    )
    nwbfile.add_trial_column(
        name="target_speed",
        description="Angular speed of the moving target",
    )

    ntrial = len(trial_view)
    for itrial in range(ntrial): 
        nwbfile.add_trial(
            start_time=trial_view.loc[itrial, 'start_time'],
            stop_time=trial_view.loc[itrial, 'stop_time'],
            target_on_time=trial_view.loc[itrial, 'target_on_time'],
            go_cue_on_time=trial_view.loc[itrial, 'go_cue_on_time'],
            movement_onset_time=trial_view.loc[itrial, 'movement_onset_time'],
            touch_time=trial_view.loc[itrial, 'touch_time'],
            result_marker=trial_view.loc[itrial, 'result_marker'],
            result_label=trial_view.loc[itrial, 'result_label'],
            target_pos_init=trial_view.loc[itrial, 'target_pos_init'],
            target_pos_go=trial_view.loc[itrial, 'target_pos_go'],
            target_pos_touch=trial_view.loc[itrial, 'target_pos_touch'],
            feedback_pos=trial_view.loc[itrial, 'feedback_pos'],
            target_speed=trial_view.loc[itrial, 'target_speed'],
            )
    nwbfile.trials.to_dataframe()


    del bhv_data


    #%% Writing standard NWB file
    save_path = os.path.join(output_dir, "standard_data.nwb")
    if os.path.exists(save_path):
        os.remove(save_path)

    with NWBHDF5IO(os.path.join(output_dir, "standard_data.nwb"), "w") as io:
        io.write(nwbfile)


# %% Parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Data_recording/20241008_Interception_001', 
                    metavar='/the/path/your/data/located/in', help='root folder')

parser.add_argument('-mp', '--map_path', 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Abel_Utah_64x4_PMd-M1-S1-A7_BlackRock.json')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Data_recording/20241008_Interception_001/formatted_data', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()

run(args.root, args.map_path, args.output)