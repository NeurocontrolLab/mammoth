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

unit_dict = {'s': pq.s}

def run(data_dir, output_dir, description_dir):

    # load preprocessed data
    nwb_saver = NWBInterface()
    neural_data = nwb_saver.read_nwb(filename = os.path.join(data_dir, 'neural_data_no_sort.nwb'))

    nwb_saver = NWBInterface()
    bhv_data = nwb_saver.read_nwb(filename = os.path.join(data_dir, 'continuous_behavior.nwb'))

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

    # event
    bhv_event = [i for i in bhv_data.segments 
                if i.name=='Event marker'][0].events[0]

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

    event = {}
    event['Markers'] = event_marker
    event['Times'] = event_time

    df_event = pd.DataFrame(event)
    df_event.columns=[len(df_event.columns)*['Event'], list(df_event.columns)]
    # df_event.index=trial_index

    # description
    description = {}
    bhv_trial_info = [i for i in bhv_data.segments if i.name=='Trial info'][0].events[0]
    info = [json.loads(i)['Trial info'].copy() for i in bhv_trial_info.labels]
    info_time = np.array(bhv_trial_info.times.rescale(pq.s).magnitude).squeeze()  
    info_df = pd.DataFrame.from_dict(info)  
    
    if info_df.iloc[-1, 1] == 'start':
        info_df = info_df.iloc[:-1]

    description['AbsoluteTrialStartTime'] = \
        info_time[info_df[info_df['status']=='start'].index.tolist()]
    description['TrialError'] = \
        [i['number'] for i in info_df.loc[info_df['status']=='end', 'wrong']]
   
    df_description = pd.DataFrame(description)
    df_description.columns=[len(df_description.columns)*['Description'], list(df_description.columns)]
    # df_description.index=trial_index


    # ObjectStatusRecord
    ObjectStatusRecord = {}
    # obj_pos = [[j for j in i.irregularlysampledsignals if 'Position' in j.name][0] for i in trials]
    # obj_sta = [[j for j in i.irregularlysampledsignals if 'Status' in j.name][0] for i in trials]
    # ObjectStatusRecord['Position'] = [np.array(i) for i in obj_pos]
    # ObjectStatusRecord['Status'] = [np.array(i) for i in obj_sta]

    def unpack_objects(d):
        objects_list = d.get('objects', [])
        for idx, obj in enumerate(objects_list, 1):
            d[f'object{idx}'] = obj
        del d['objects']
        return d
    
    bhv_frame = [i for i in bhv_data.segments if i.name=='Frame'][0].events[0]
    bhv_frame_ = [unpack_objects(json.loads(i)) for i in bhv_frame.labels]
    frame_view = pd.DataFrame.from_dict(bhv_frame_)
    frame_time = np.array(bhv_frame.times.rescale(pq.s).magnitude).squeeze()
    frame_view['time'] = frame_time

    pos_ls = []
    for itrial in range(info_df['trial_number'].values[-1]):
        start_time = info_time[info_df[(info_df['trial_number']==itrial)&(info_df['status']=='start')].index]
        end_time = info_time[info_df[(info_df['trial_number']==itrial)&(info_df['status']=='end')].index]
        
        if 'object3' in frame_view.loc[(frame_time>start_time)&(frame_time<end_time), :].columns:
            pls3_candidate =  frame_view.loc[(frame_time>start_time)&(frame_time<end_time), 'object3']

            pos0 = [[i['name'], i['pos']] for i in pls3_candidate if (not pd.isna(i))]
            if len(pos0)>0:
                pos = [i[1] for i in pos0 if 'feedback' in i[0]][0]
                pos_ls.append(pos)
            
            else:
                pos_ls.append(np.nan)
        
        else:
            pos_ls.append(np.nan)

    ObjectStatusRecord['Position'] = pos_ls
    df_ObjectStatusRecord = pd.DataFrame(ObjectStatusRecord)
    df_ObjectStatusRecord.columns=[len(df_ObjectStatusRecord.columns)*['ObjectStatusRecord'], list(df_ObjectStatusRecord.columns)]
    # df_ObjectStatusRecord.index=trial_index

    frames = [df_description[:len(df_ObjectStatusRecord)], df_ObjectStatusRecord]
    df_trials = pd.concat(frames,axis=1).rename_axis('TrialIndex')

    st = [i for i in neural_data.segments if i.name=='TCR'][0].spiketrains
    st_shift = [i.time_shift(diff_time_mean) for i in st]

    #%% data slicing
    kwargs_list = []
    trial_ind = []
    move_direction = []
    # speed = []

    for ind, i in enumerate(df_trials.iloc):
        Description = i['Description']
        # Event = i['Event']
        ObjectStatusRecord = i['ObjectStatusRecord']
        if 0!= Description['TrialError']:
            continue
        
        start_time = info_time[info_df[(info_df['trial_number']==ind)&(info_df['status']=='start')].index][0]
        end_time = info_time[info_df[(info_df['trial_number']==ind)&(info_df['status']=='end')].index][0]
        aligned_time = event_time[(np.array(event_marker)==5)*(event_time>start_time)&(event_time<end_time)][0]

        kwargs = {
                # 't_start': Event['Times'][Event['Labels']==5]-1*pq.s+Description['AbsoluteTrialStartTime'],
                # 't_stop': Event['Times'][Event['Labels']==5]+1*pq.s+Description['AbsoluteTrialStartTime'],
                't_start': aligned_time*pq.s-0.05*pq.s,
                't_stop': aligned_time*pq.s+0*pq.s,
                'aligned_marker':[5],
                'aligned_marker2':[5],
                'trial_index': ind}
        kwargs_list.append(kwargs)
        trial_ind.append(ind)
        move_direction.append(ObjectStatusRecord['Position'])
        # speed.append(Description['AngularV'])

    trial_ind = np.array(trial_ind)
    move_direction = np.array(move_direction)
    kwargs = {
        'kernel' : GaussianKernel(50*pq.ms),
        'sampling_period' : 10*pq.ms
    }

    # sliced_st = SpikeStatistics.preprocessing('spike_time', kwargs_list, st_shift)
    sliced_th = SpikeStatistics.preprocessing('time_histogram', kwargs_list, st_shift, binsize=10*pq.ms)
    sliced_is = SpikeStatistics.preprocessing('instantaneous_rate', kwargs_list, st_shift, **kwargs)
    print(sliced_is.shape)
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
    shuffled_st_shift = [i.time_shift(diff_time_mean+2*pq.s) for i in st]
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
    ax1.set_ylim([-0.2,1])
        
    plt.savefig(os.path.join(output_dir, 'sliding_R2.png'), dpi=300, bbox_inches = 'tight')

    with open(os.path.join(description_dir, 'qc_summary.json'), "r") as file:
        summary = json.load(file)

        summary['neural correlation R^2 (MO+/-1s)'] = time_suc
    
    with open(os.path.join(output_dir, 'qc_summary.json'), 'w') as file:
        json.dump(summary, file)


parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-d", "--data", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Data_recording/20240924_Interception_001/formatted_data', 
                    metavar='/the/path/your/nwb/data/located/in', help='data folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Data_recording/20240924_Interception_001/description', 
                    metavar='/the/path/you/want/to/save', help='output folder')
parser.add_argument("-s", "--description", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Data_recording/20240924_Interception_001/description', 
                    metavar='/the/path/your/descriptive/data/located/in', help='root folder')

args = parser.parse_args()

run(args.data, args.output, args.description)
