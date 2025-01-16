#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 15:18:10 2024

@author: cuilab
"""

import argparse
import os
import json
import pandas as pd
import numpy as np
import quantities as pq
from SmartNeo.interface_layer.nwb_interface import NWBInterface
from SmartNeo.analysis_layer.spike_preprocessor3 import SpikeStatistics
import matplotlib.pyplot as plt
from elephant.kernels import GaussianKernel
import numpy as np
from sklearn import linear_model
from sklearn.model_selection import cross_val_score
from pynwb import NWBHDF5IO

unit_dict = {'s': pq.s}

def run(data_dir, output_dir, description_dir):

    # load preprocessed data
    filename = os.path.join(data_dir, 'neural_data_no_sort.nwb')
    neural_data = NWBHDF5IO(filename, mode='r').read()

    filename = os.path.join(data_dir, 'continuous_behavior.nwb')
    bhv_data = NWBHDF5IO(filename, mode='r').read()

    # nwb_saver = NWBInterface()
    # neural_data = nwb_saver.read_nwb(filename = os.path.join(data_dir, 'neural_data_no_sort.nwb'))

    # nwb_saver = NWBInterface()
    # bhv_data = nwb_saver.read_nwb(filename = os.path.join(data_dir, 'continuous_behavior.nwb'))

    filename_ = os.path.join(description_dir, 'diff_time_mean.txt')
    assert os.path.exists(filename_), "Please check time consistency first!"
    with open(filename_, 'r') as file:
        dt_str = file.read()  # 读取整个文件内容
    unit = unit_dict[dt_str.split(' ')[1][0:-1]]
    diff_time_mean = float(dt_str.split(' ')[0])*unit

    # index
    # index = [i.index for i in trial_behavior.segments]
    # sorted_index = np.argsort(index)
    # trials = [trial_behavior.segments[i] for i in sorted_index]
    # trial_index = [i.index for i in trials]

    #%% load bhv data
    # event
    event_marker = bhv_data.processing['behavior']['MonkeyLogicEvents'].data[:]
    event_time = bhv_data.processing['behavior']['MonkeyLogicEvents'].timestamps[:] * pq.s.magnitude

    event = {}
    event['Markers'] = event_marker
    event['Times'] = event_time

    df_event = pd.DataFrame(event)
    df_event.columns=[len(df_event.columns)*['Event'], list(df_event.columns)]
    # df_event.index=trial_index

    # description
    description = {}
    bhv_trial_info = bhv_data.trials.to_dataframe()
    
    description['AbsoluteTrialStartTime'] = \
        bhv_trial_info['AbsoluteTrialStartTime']
    description['TrialError'] = \
        bhv_trial_info['TrialError']
   
    df_description = pd.DataFrame(description)
    df_description.columns=[len(df_description.columns)*['Description'], list(df_description.columns)]
    # df_description.index=trial_index


    # ObjectStatusRecord
    ObjectStatusRecord = {}
    ObjectStatusRecord['Status'] = [json.loads(i) for i in bhv_trial_info['ObjectStatusRecord_Status']]
    ObjectStatusRecord['Position'] = [json.loads(i) for i in bhv_trial_info['ObjectStatusRecord_Position']]

    df_ObjectStatusRecord = pd.DataFrame(ObjectStatusRecord)
    df_ObjectStatusRecord.columns=[len(df_ObjectStatusRecord.columns)*['ObjectStatusRecord'], list(df_ObjectStatusRecord.columns)]
    # df_ObjectStatusRecord.index=trial_index

    trials = [df_description[:len(df_ObjectStatusRecord)], df_ObjectStatusRecord]
    df_trials = pd.concat(trials, axis=1).rename_axis('TrialIndex')

    ## load neural data
    st = neural_data.units.to_dataframe()['spike_times'].tolist()
    # st.sort()
    st_shift = [(i+diff_time_mean.magnitude)*pq.s for i in st]


    #%% data slicing
    kwargs_list = []
    trial_ind = []
    move_direction = []
    # speed = []

    for ind, i in enumerate(df_trials.iloc):
        Description = i['Description']
        # Event = i['Event']
        ObjectStatusRecord = i['ObjectStatusRecord']
        if Description['TrialError']!=0 :
            continue
        
        start_time = bhv_trial_info.loc[ind, 'start_time']
        end_time = bhv_trial_info.loc[ind, 'stop_time']
        aligned_time = event_time[(np.array(event_marker)==5)*(event_time>start_time)&(event_time<end_time)]

        kwargs = {
                # 't_start': Event['Times'][Event['Labels']==5]-1*pq.s+Description['AbsoluteTrialStartTime'],
                # 't_stop': Event['Times'][Event['Labels']==5]+1*pq.s+Description['AbsoluteTrialStartTime'],
                't_start': aligned_time*pq.s-1*pq.s,
                't_stop': aligned_time*pq.s+1*pq.s,
                'aligned_marker':[5],
                'aligned_marker2':[5],
                'trial_index': ind}
        kwargs_list.append(kwargs)
        trial_ind.append(ind)

        md = np.array(ObjectStatusRecord['Position'])[np.array(ObjectStatusRecord['Status'])[:, 2]==1, 2, :]
        move_direction.append(md)
        # speed.append(Description['AngularV'])

    trial_ind = np.array(trial_ind)
    move_direction = np.array(move_direction).squeeze()
    kwargs = {
        'kernel' : GaussianKernel(50*pq.ms),
        'sampling_period' : 50*pq.ms
    }

    # sliced_st = SpikeStatistics.preprocessing('spike_time', kwargs_list, st_shift)
    sliced_is = SpikeStatistics.preprocessing('instantaneous_rate', kwargs_list, st_shift, **kwargs)
    trial_ind = trial_ind[trial_ind<sliced_is.shape[0]]
    move_direction = move_direction[0:len(trial_ind)]
    sliced_is = sliced_is[trial_ind,:,5:-5]

    #%% data regression
    time_suc = []
    for t_i in range(sliced_is.shape[-1]):
        X = sliced_is[:,:,t_i].reshape((sliced_is.shape[0],-1))
        Y = move_direction
        reg = linear_model.MultiTaskLasso()
        scores = np.mean(cross_val_score(reg, X, Y, cv=5, scoring='r2',n_jobs=-1))
        time_suc.append(scores)

    print(time_suc)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    with open(os.path.join(output_dir, 'correlation_test_score.txt'), 'a') as file:
        file.write(f'{scores}\n')


    #%% data shuffiled slicing
    shuffled_st_shift = [(i+diff_time_mean.magnitude+2)*pq.s for i in st]
    shuffled_sliced_is = SpikeStatistics.preprocessing('instantaneous_rate', kwargs_list, shuffled_st_shift, **kwargs)
    shuffled_sliced_is = shuffled_sliced_is[trial_ind,:,5:-5]

    #%% data shuffiled regression
    shuffled_time_suc = []
    for t_i in range(sliced_is.shape[-1]):
        X = shuffled_sliced_is[:,:,t_i].reshape((sliced_is.shape[0],-1))
        Y = move_direction
        reg = linear_model.MultiTaskLasso()
        scores = np.mean(cross_val_score(reg, X, Y, cv=5, scoring='r2',n_jobs=-1))
        shuffled_time_suc.append(scores)

    import seaborn as sns
    fig, ax1 = plt.subplots()

    data = {}
    data["Times"] = list(np.arange(0,30,1))+list(np.arange(0,30,1))
    data['R^2'] = time_suc + shuffled_time_suc
    data['Mismatch'] = ['Aligned']*len(time_suc)+['2 sec lag']*len(shuffled_time_suc)
    sns.lineplot(x="Times", y='R^2',data=data, hue='Mismatch', ax=ax1)
    sta, mid, end = data["Times"][4], data["Times"][14], data["Times"][24]
    ax1.set_xticks([4,14,24],['-0.5','0','0.5'],fontsize=15,rotation=45)
    ax1.set_yticks([0,0.5,1],[0,0.5,1],fontsize=15)
    ax1.set_ylabel('$R^2$',fontsize=15)
    ax1.set_xlabel('Times (s)',fontsize=15)
    ax1.set_xlim([0,29])
    # ax1.set_ylim([-0.2,1])
        
    plt.savefig(os.path.join(output_dir, 'sliding_R2.png'), dpi=300, bbox_inches = 'tight')

    with open(os.path.join(description_dir, 'qc_summary.json'), "r") as file:
        summary = json.load(file)

        summary['neural correlation R^2 (MO+/-1s)'] = time_suc
    
    with open(os.path.join(output_dir, 'qc_summary.json'), 'w') as file:
        json.dump(summary, file)


parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-d", "--data", type=str,
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/interception/Abel/data_recording/20240528_Interception_001/formatted_data', 
                    metavar='/the/path/your/nwb/data/located/in', help='data folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/interception/Abel/data_recording/20240528_Interception_001/description', 
                    metavar='/the/path/you/want/to/save', help='output folder')
parser.add_argument("-s", "--description", type=str,
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/interception/Abel/data_recording/20240528_Interception_001/description', 
                    metavar='/the/path/your/descriptive/data/located/in', help='root folder')

args = parser.parse_args()

run(args.data, args.output, args.description)
