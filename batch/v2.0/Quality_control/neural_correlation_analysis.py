#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 15:58:49 2024

@author: lichenyang
"""

import argparse
import os
import pandas as pd
import numpy as np
import quantities as pq
from SmartNeo.interface_layer.nwb_interface import NWBInterface
from SmartNeo.analysis_layer.spike_preprocessor3 import SpikeStatistics
import matplotlib.pyplot as plt
from elephant.kernels import GaussianKernel

#%% read data
# data folder
unit_dict = {'s': pq.s}
parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-f", "--file", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Abel/Data_recording/20240705_centerOut_001', 
                    metavar='/the/path/your/data/located/in', help='input folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Abel/Data_recording/20240705_centerOut_001', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()
raw_dirname = args.file

formated_data_path = os.path.join(raw_dirname, 'formatted_data')

filename_ = os.path.join(raw_dirname,'description', 'diff_time_mean.txt')
with open(filename_, 'r') as file:
    dt_str = file.read()  # 读取整个文件内容
unit = unit_dict[dt_str.split(' ')[1][0:-1]]
diff_time_mean = float(dt_str.split(' ')[0])*unit
# read formated SmartNeo data
nwb_saver = NWBInterface()
neural_data = nwb_saver.read_nwb(filename = os.path.join(formated_data_path,'neural_data_no_sort.nwb'))

continuous_saver = NWBInterface()
bhv_data = nwb_saver.read_nwb(filename = os.path.join(formated_data_path,'continuous_behavior.nwb'))

coeff_pd = {}
for i in bhv_data.segments[1].irregularlysampledsignals:
    coeff_pd[i.name[1]] = i.as_array().squeeze()
coeff_pd = pd.DataFrame(coeff_pd)
# coeff_pd.rename_axis('Times')
coeff_pd.index = bhv_data.segments[1].irregularlysampledsignals[0].times
coeff_pd.index.name = 'Times'
coeff_pd.head()

pos_pd = {}
for i in bhv_data.segments[3].irregularlysampledsignals:
    pos_pd[i.name[1]] = list(i.as_array().squeeze())

pos_pd = pd.DataFrame(pos_pd)
pos_pd.rename_axis('Times')
pos_pd.index = bhv_data.segments[3].irregularlysampledsignals[0].times
pos_pd.index.name = 'Times'
pos_unit = bhv_data.segments[3].irregularlysampledsignals[0].times.units

pos_pd.head()

trial_pd = {}
for i in bhv_data.segments[-1].irregularlysampledsignals:
    trial_pd[i.name[1]] = list(i.as_array().squeeze())

trial_pd = pd.DataFrame(trial_pd)
trial_pd.rename_axis('Times')
trial_pd.index = bhv_data.segments[-1].irregularlysampledsignals[0].times
trial_pd.index.name = 'Times'
trial_pd.head()

events = bhv_data.segments[0].events[0]
marker_24_pos = np.where(events.labels==24)[0]
marker_5_pos = np.where(events.labels==5)[0]

trial_event = [events[i:j+1] for i,j in zip(marker_24_pos,marker_5_pos)]

marker_24_time = events.times[np.where(events.labels==24)[0]]
marker_5_time = events.times[np.where(events.labels==5)[0]]
trial_pos = [pos_pd.iloc[np.where(pos_pd.index<j * (pos_pd.index>i))[0]] \
                 for i,j in zip(marker_24_time,marker_5_time)]
    
trial_coeff = [coeff_pd.iloc[np.where(coeff_pd.index<j * (coeff_pd.index>i))[0]] \
                   for i,j in zip(marker_24_time,marker_5_time)]

neural_event = [i for i in neural_data.segments if i.name=='RecordingSystemEvent'][0]


events = bhv_data.segments[0].events[0]
marker_24_pos = np.where(events.labels==3)[0]
marker_5_pos = np.where(events.labels==5)[0]
trial_event = [events[i:j+1] for i,j in zip(marker_24_pos,marker_5_pos)]

marker_24_time = events.times[np.where(events.labels==3)[0]]
marker_5_time = events.times[np.where(events.labels==5)[0]]
trial_pos = [pos_pd.iloc[np.where(pos_pd.index<j * (pos_pd.index>i))[0]] \
             for i,j in zip(marker_24_time,marker_5_time)]
trial_coeff = [coeff_pd.iloc[np.where(coeff_pd.index<j * (coeff_pd.index>i))[0]] for i,j in zip(marker_24_time,marker_5_time)]

trial = [events[i:j+1] for i,j in zip(marker_24_pos,marker_5_pos)]
# duration = [i.times[-1]-i.times[3] for i in trial]
trial_success = [20 in i.labels for i in trial]
trial_time = [(i.times[-1]-trial[0].times[0]).magnitude/60 for i in trial]

flag_times = bhv_data.segments[1].irregularlysampledsignals[-1].times
training_flag = np.array(bhv_data.segments[1].irregularlysampledsignals[-1])
trial_flag = [training_flag[np.where(j[0]<flag_times * (j[-1]>flag_times))[0]] for j in trial]
trial_flag_times = [flag_times[np.where(j[0]<flag_times * (j[-1]>flag_times))[0]] for j in trial]
trial_flag_times = [(i-trial[0].times[0]).magnitude.squeeze()/60 for i in trial_flag_times]

st = [i for i in neural_data.segments if i.name=='TCR'][0].spiketrains
st_shift = [i.time_shift(diff_time_mean) for i in st]

from elephant.kernels import GaussianKernel
#%% data slicing
kwargs_list = []
trial_ind = []
move_direction = []
pos_input = []
is_kwargs = {
    'kernel' : GaussianKernel(20*pq.ms),
    'sampling_period' : 20*pq.ms
}

flag_list = []
for ind, i in enumerate(trial_pos):
    if not sum(trial_coeff[ind]['training flag'].tolist())==0:
        continue
    if not 20 in trial[ind].labels:
        continue
    pos_time = np.array(i.index)*pos_unit
    kwargs = {'t_start': pos_time[0]-5*is_kwargs['sampling_period'],
              't_stop': pos_time[-1]+5*is_kwargs['sampling_period'],
              'aligned_marker':[3],
              'aligned_marker2':[5],
              'trial_index': ind}
    pos_input.append(i)
    kwargs_list.append(kwargs)
    trial_ind.append(ind)
    move_direction.append(np.array(i['pos'].to_list()).squeeze()[-1])
    flag_list.append(trial_flag[ind])

trial_ind = np.array(trial_ind)
sliced_is = SpikeStatistics.preprocessing('instantaneous_rate', 
                                          kwargs_list, st_shift, **is_kwargs)
trial_ind = trial_ind[np.where(trial_ind<sliced_is.shape[0])[0]]
sliced_is = sliced_is[trial_ind]
pos_input = pos_input[0:len(trial_ind)]
kwargs_list = kwargs_list[0:len(trial_ind)]
flag_list = flag_list[0:len(trial_ind)]

bin_num = [(i['t_stop'] - i['t_start']).rescale(is_kwargs['sampling_period'].units)/is_kwargs['sampling_period'] for i in kwargs_list]
select_is = [sliced_is[ind,:,0:int(i.magnitude)] for ind, i in enumerate(bin_num)]
select_is = [i[:,5:-5] for i in select_is]
fr_sum = [np.sum(i.magnitude) for i in select_is]
if len(np.where(np.array(fr_sum)==0)[0])!=0:
    select_ind = max(np.where(np.array(fr_sum)==0)[0])+1

    select_is = select_is[select_ind::]
    pos_input = pos_input[select_ind::]
    
pos = []
speed = []
from scipy import signal
from scipy.ndimage import gaussian_filter

for i,j in zip(pos_input,select_is):
    traj = np.array(i['pos'].tolist())  
    vel = np.array(i['vel'].tolist()) 
    traj_x = gaussian_filter(signal.resample(traj[:,0],j.shape[1]+1)[1::],sigma=1)
    traj_y = gaussian_filter(signal.resample(traj[:,1],j.shape[1]+1)[1::],sigma=1)
    vel_x = gaussian_filter(signal.resample(vel[:,0],j.shape[1]+1)[1::],sigma=1)
    vel_y = gaussian_filter(signal.resample(vel[:,1],j.shape[1]+1)[1::],sigma=1)
    pos.append(np.array([traj_x,traj_y]).T)
    speed.append(np.array([vel_x,vel_y]).T)
    plt.plot(traj_x, traj_y)

X_pos = np.concatenate(pos,axis=0)
X_speed = np.concatenate(speed,axis=0)
Y_neu = np.concatenate(select_is,axis=1).T

from sklearn.linear_model import MultiTaskLasso
from sklearn.model_selection import cross_val_score
lasso_speed = MultiTaskLasso(fit_intercept=False)
scores = cross_val_score(lasso_speed, Y_neu, X_pos, scoring='r2',cv=5,n_jobs=32)
print(scores)

with open(os.path.join(args.output, 'description', 'correlation_test_score.txt'), 'a') as file:
    file.write(f'{scores}\n')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


