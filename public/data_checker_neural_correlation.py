#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 15:18:10 2024

@author: cuilab
"""

import argparse
import os
import pandas as pd
import numpy as np
import quantities as pq
from SmartNeo.interface_layer.nwb_interface import NWBInterface
from SmartNeo.analysis_layer.spike_preprocessor3 import SpikeStatistics
from elephant.kernels import GaussianKernel
from sklearn.linear_model import MultiTaskLasso
from sklearn.model_selection import cross_val_score
from elephant.kernels import GaussianKernel
from scipy import signal
from scipy.ndimage import gaussian_filter


def run(data_dir, output_dir, description_dir):

    #%% read data
    # read formated NWB data (TCR & BCI behaviour)
    nwb_saver = NWBInterface()
    neural_data = nwb_saver.read_nwb(filename = os.path.join(data_dir,'neural_data_no_sort.nwb'))

    continuous_saver = NWBInterface()
    bhv_data = continuous_saver.read_nwb(filename = os.path.join(data_dir,'continuous_behavior.nwb'))

    # read diff time from description
    unit_dict = {'s': pq.s} # bad convertor of quantities
    filename_ = os.path.join(description_dir, 'diff_time_mean.txt')
    with open(filename_, 'r') as file:
        dt_str = file.read()
    unit = unit_dict[dt_str.split(' ')[1][0:-1]]
    diff_time_mean = float(dt_str.split(' ')[0])*unit

    # read pandas.DataFrame data
    pos_pd = pd.read_csv(os.path.join(description_dir, 'pos_pd.csv'), header=0, index_col=0)
    coeff_pd = pd.read_csv(os.path.join(description_dir, 'coeff_pd.csv'), header=0, index_col=0)

    #%% slice behavioral data
    # get trial time markers
    events = bhv_data.segments[0].events[0]
    marker_start_pos = np.where(events.labels==3)[0]
    marker_end_pos = np.where(events.labels==5)[0]
    marker_start_time = events.times[marker_start_pos].rescale(pq.s).magnitude
    marker_end_time = events.times[marker_end_pos].rescale(pq.s).magnitude

    # slice pos_pd and coeff_pd by start and end time
    trial_pos = [pos_pd.iloc[np.where(pos_pd.index<j * (pos_pd.index>i))[0]] \
                for i,j in zip(marker_start_time, marker_end_time)]
    trial_coeff = [coeff_pd.iloc[np.where(coeff_pd.index<j * (coeff_pd.index>i))[0]] \
                    for i,j in zip(marker_start_time,marker_end_time)]
    trial = [events[i:j+1] for i,j in zip(marker_start_pos, marker_end_pos)]

    #%% slice neural data
    # get TCR
    st = [i for i in neural_data.segments if i.name=='TCR'][0].spiketrains # state spike trains (TCR)
    st_shift = [i.time_shift(diff_time_mean) for i in st]

    # data slicing
    kwargs_list = [] # slicing parameters container
    trial_ind = [] # success trials container
    move_direction = [] # trial related direction (end direction)
    pos_input = [] # sliced instantaneous postion

    is_kwargs = {
        'kernel' : GaussianKernel(20*pq.ms),
        'sampling_period' : 20*pq.ms
    } # spike statistic parameters

    def str2list(s):
        s = s.strip('[]').split(' ')
        return [float(i) for i in s if len(i)!=0]

    for ind, i in enumerate(trial_pos):
        # this branch means the 'training flag' needs to be zero
        if not sum(trial_coeff[ind]['training flag'].tolist())==0:
            continue
        # '20' represents the monkey complete the trial successfully
        if not 20 in trial[ind].labels:
            continue
        pos_time = np.array(i.index)*pq.s
        # parameters for slicing
        kwargs = {'t_start': pos_time[0]-5*is_kwargs['sampling_period'],
                  't_stop': pos_time[-1]+5*is_kwargs['sampling_period'],
                  'aligned_marker':[3],
                  'aligned_marker2':[5],
                  'trial_index': ind}
        # trial info saver
        pos_input.append(i)
        kwargs_list.append(kwargs)
        trial_ind.append(ind)
        move_direction.append(np.array([str2list(ele) for ele in i['pos'].to_list()]).squeeze()[-1]) # end pos

    trial_ind = np.array(trial_ind) # for indexing
    sliced_is = SpikeStatistics.preprocessing('instantaneous_rate', 
                                              kwargs_list, st_shift, **is_kwargs) # spike slicing and preprocessing
    trial_ind = trial_ind[np.where(trial_ind<sliced_is.shape[0])[0]]
    sliced_is = sliced_is[trial_ind]
    pos_input = pos_input[0:len(trial_ind)]
    kwargs_list = kwargs_list[0:len(trial_ind)]

    # left align
    bin_num = [(i['t_stop'] - i['t_start']).rescale(is_kwargs['sampling_period'].units)/is_kwargs['sampling_period'] for i in kwargs_list]
    select_is = [sliced_is[ind, :, 0:int(i.magnitude)] for ind, i in enumerate(bin_num)] # slice neural data, remove zero padding
    select_is = [i[:, 5:-5] for i in select_is] # remove left and right border of convolution
    # remove zero firing trial (unnecessary)
    fr_sum = [np.sum(i.magnitude) for i in select_is]
    if len(np.where(np.array(fr_sum)==0)[0])!=0:
        select_ind = max(np.where(np.array(fr_sum)==0)[0])+1
        select_is = select_is[select_ind::]
        pos_input = pos_input[select_ind::]

    #%% neural correlation analysis
    # container of processing result
    pos = []

    for i,j in zip(pos_input, select_is):
        traj = np.array([str2list(ele) for ele in i['pos'].tolist()]) # not sure whethe scipy can handle pandas
        vel = np.array([str2list(ele) for ele in i['vel'].tolist()])
        # x and y axis of pos and vel, performing resampling and smoothing
        traj_x = gaussian_filter(signal.resample(traj[:,0], j.shape[1]+1)[1::], sigma=1)
        traj_y = gaussian_filter(signal.resample(traj[:,1], j.shape[1]+1)[1::], sigma=1)
        pos.append(np.array([traj_x,traj_y]).T)
        
    # build train & test set
    X_pos = np.concatenate(pos, axis=0)
    Y_neu = np.concatenate(select_is, axis=1).T

    lasso_speed = MultiTaskLasso(fit_intercept=False)
    scores = cross_val_score(lasso_speed, Y_neu, X_pos, scoring='r2', cv=5, n_jobs=32)
    print(scores)

    with open(os.path.join(output_dir, 'correlation_test_score.txt'), 'w') as file:
        file.write(f'{scores}\n')


parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Brain_control/20240307_BrUtahInterception120SemiBc_001', 
                    metavar='/the/root/path/your/data/located/in', help='root folder')
parser.add_argument("-d", "--data", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Brain_control/20240307_BrUtahInterception120SemiBc_001/formatted_data', 
                    metavar='/the/path/your/data/located/in', help='data folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Brain_control/20240307_BrUtahInterception120SemiBc_001/description', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()

run(args.data, args.output, args.output)