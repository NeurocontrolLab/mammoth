#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

import os
import glob
import json
import numpy as np
import pandas as pd
import argparse
import datetime
from dateutil.tz import tzlocal
import quantities as pq
import hdf5storage
from pynwb import NWBFile, NWBHDF5IO, TimeSeries

TrialError_annotation = {
0: 'Correct',
1: 'no touch',
2: 'break touch center before target shows',
3: 'break touch center before go-signal shows',
4: 'wrong place of touch center',
5: 'no target touch',
6: 'break touch target',
7: 'wrong place of touch target',
9: 'false start, reaction time after go-signal too short',
}

Event_Label = {
9: 'trial start',
1: 'center on', 
2: 'touch center',
3: 'target on',
11: 'memory on',
4: 'go signal',
5: 'movement onset',
6: 'touch target',
7: 'touch wrong',
8: 'touch right',
10: 'leave target',
18: 'trial end',
24: 'TTL, trigger EMG and Vicon',
13:	'system frame skip',
14: 'tbc',
}

def format_file(root_dir, output_dir):
    # FILEPATH = os.path.dirname(os.path.abspath(__file__))
   
    #%% extract info
    # load converted .mat monkeylogic data
    bhv_name = glob.glob(os.path.join(os.path.join(root_dir,'bhv'),'*.bhv2'))[0].split('/')[-1].split('.')[0]
    bhv_mat = [i for i in glob.glob(os.path.join(root_dir,'bhv','*.mat')) if bhv_name in i][0]
    bhvsave = hdf5storage.loadmat(bhv_mat)
    bhvkey = [i for i in bhvsave if '_' not in i][0]
    bhvsave = bhvsave[bhvkey].squeeze()

    #%% prepare trial dataframe
    # convert to dataframe
    bhvdf = pd.DataFrame.from_dict(bhvsave)


    # flatten nested columns
    def flatten_df(df):
        def flatten_columns(df, key):
            new_columns = {}
            for subkey in df[key][0].dtype.names:
                new_key = key + '_' + subkey 
                new_value = [trial[subkey].squeeze() for trial in df[key]] 
                new_columns[new_key] = new_value
            return new_columns
        
        while any(isinstance(df[key][0], np.ndarray) and df[key][0].dtype.names is not None for key in df.columns):
            for key in df.columns:
                if isinstance(df[key][0], np.ndarray) and df[key][0].dtype.names is not None:
                    flattened_columns = flatten_columns(df, key)
                    for nk, nv in flattened_columns.items():
                        df[nk] = nv 
    
                    df = df.drop(columns=key)

        return df


    bhvdf = flatten_df(bhvdf)

    # flatten elemental arrays and convert ndarray to dict
    def flatten(arr):
        if isinstance(arr, np.ndarray):
            if arr.shape==(1, 1):
                return arr[0, 0]
            else:
                return arr.squeeze()
        return arr
    
    def arr_to_dict(arr):
        if isinstance(arr, np.ndarray):
            if not bhvdf[key][0].shape and bhvdf[key][0].item().dtype.names:
                da = {}
                narr = arr.item()
                for k in narr.dtype.names:
                    v = narr[k]
                    da[k] = flatten(v)
                return da
            return arr
        return arr


    dict_columns = []
    for key in bhvdf.columns:
        if bhvdf[key][0].dtype==object:
            if not bhvdf[key][0].shape and bhvdf[key][0].item().dtype.names:
                dict_columns.append(key)

    while any(isinstance(bhvdf[key][0], np.ndarray) and bhvdf[key][0].shape==(1, 1) for key in bhvdf.columns):
        for key in bhvdf.columns:
            if key not in dict_columns:
                bhvdf[key] = bhvdf[key].apply(lambda x: flatten(x))
            else:
                bhvdf[key] = bhvdf[key].apply(lambda x: arr_to_dict(x))


    #%% set time unit
    for i in bhvdf.columns:
        if 'time' in i.lower():
            bhvdf[i] = bhvdf[i].apply(lambda x: x*pq.ms.rescale(pq.s).magnitude)


    for i in bhvdf.index:
        bhvdf.loc[i, 'start_time'] = (bhvdf.loc[i, 'BehavioralCodes_CodeTimes'][0] + bhvdf.loc[i, 'AbsoluteTrialStartTime'])
        bhvdf.loc[i, 'stop_time'] = (bhvdf.loc[i, 'BehavioralCodes_CodeTimes'][-1] + bhvdf.loc[i, 'AbsoluteTrialStartTime'])

    bhvdf['TrialError_annotation'] = bhvdf['TrialError'].apply(lambda x: TrialError_annotation[x])
    

    # add events
    event_markers = [i for _, trial in bhvdf.iterrows() for i in trial['BehavioralCodes_CodeNumbers']]
    event_times = [i+trial['AbsoluteTrialStartTime'] for _, trial in bhvdf.iterrows() for i in trial['BehavioralCodes_CodeTimes']]
    event_times = np.array(event_times) # *pq.ms.rescale(pq.s).magnitude


    #%% convert elemental array to list and then string
    def convert_to_list(obj):
        if isinstance(obj, np.ndarray):
            return [convert_to_list(x) for x in obj.tolist()]
        elif isinstance(obj, list):
            return [convert_to_list(x) for x in obj]
        elif isinstance(obj, dict):
            return {k: convert_to_list(v) for k, v in obj.items()}
        else:
            return obj
    
    # convert array to string
    for key in bhvdf.columns:
        if isinstance(bhvdf[key][0], (np.ndarray, dict)):  
            bhvdf[key] = bhvdf[key].apply(lambda x: json.dumps(convert_to_list(x))) 
    

    # unify shape
    ucol = []
    for col in bhvdf.columns:
        if bhvdf[col].apply(lambda x: isinstance(x, (list, np.ndarray))).any():
            ucol.append(col)
    
    def unify_column_to_list(x):
        if isinstance(x, list):
            return x
        elif isinstance(x, (float, int)):
            return [x]
        elif isinstance(x, np.ndarray):
            return x.tolist()
        else:
            return []

    for col in ucol:
        bhvdf[col] = bhvdf[col].apply(lambda x: json.dumps(unify_column_to_list(x)))


    #%% create an NWB file
    nwbfile = NWBFile(
        session_description="behavioral recording",
        identifier=root_dir,
        session_start_time=datetime.datetime.combine(
            datetime.date(*bhvsave[0]['TrialDateTime'].astype(int)[0, 0:3]),
            datetime.time(*bhvsave[0]['TrialDateTime'].astype(int)[0, 3:]), tzlocal()),
        lab="Cui Lab",
        institution="CIBR",
        experiment_description="interception",
        session_id=root_dir.split(os.sep)[-1],
    )

    # create a behavioral module
    behavior_module = nwbfile.create_processing_module(
        name="behavior", 
        description="Raw behavioral data."
    )

    behavioral_events = TimeSeries(
        name='MonkeyLogicEvents',
        data=event_markers,
        unit='NA',
        timestamps = event_times,
        description='Event markers recorded by MokeyLogic'
    )
    behavior_module.add(behavioral_events)

    behavioral_events2 = TimeSeries(
        name='MonkeyLogicEvents_Labels',
        data=[Event_Label[i] for i in event_markers],
        unit='NA',
        timestamps = event_times,
        description='Event labels correspond to event markers'
    )
    behavior_module.add(behavioral_events2)


    # add trial info          
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
    
    # view trial info
    # nwbfile.trials.to_dataframe()


    #%% Writing standard NWB file
    if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    save_path = os.path.join(output_dir, "continuous_behavior.nwb")
    if os.path.exists(save_path):
        os.remove(save_path)

    with NWBHDF5IO(os.path.join(output_dir, "continuous_behavior.nwb"), "w") as io:
        io.write(nwbfile)

# df = nwbfile.trials.to_dataframe() 
# for col in df.columns:
#     print(f"Column: {col}")
#     print(df[col].apply(type).value_counts())
#     print("-" * 50)

# for col in df.columns:
#     if df[col].apply(lambda x: isinstance(x, (list, np.ndarray))).any():
#         print(f"Column: {col}")
#         print(df[col].apply(lambda x: np.shape(x) if isinstance(x, (list, np.ndarray)) else None).value_counts())
#         print("-" * 50)

# for col in df.columns:
#     try:
#         df[[col]].to_hdf(os.path.join(output_dir, 'test.h5'), key='test', mode='w')
#         print(f"Column {col}: OK")
#     except Exception as e:
#         print(f"Column {col}: Error - {e}")


#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/interception/Abel/data_recording/20240602_Interception_001', 
                    metavar='/the/root/path/your/data/located/in', help='root folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/interception/Abel/data_recording/20240602_Interception_001/formatted_data', 
                    metavar='/the/path/you/want/to/save', help='output folder')
args = parser.parse_args()

format_file(args.root, args.output)
