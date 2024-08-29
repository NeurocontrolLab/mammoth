#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:35:16 2024

@author: cuilab
"""
import os
import pandas as pd
import argparse
import numpy as np

sf = 30000
def converted(spike_times, sampling_frequency=sf):
    binlen = 1/1000/(1/sampling_frequency)  #ms
    spike_times_ms = np.ceil(spike_times/binlen)
    
    return spike_times_ms


parser = argparse.ArgumentParser()
parser.add_argument('-dp', '--data_path',
                    default=r'/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001')
args = parser.parse_args()

data_path = args.data_path

# data_path = r'/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001'
sorted_path = os.path.join(data_path, 'sorted_data')

ks_path = os.path.join(sorted_path, 'kilosort2_5_output')
skcs_path = os.path.join(sorted_path, 'spykingcircus_output')
assert os.path.exists(ks_path) & os.path.exists(skcs_path)
kshank_list = next(os.walk(ks_path))[1]

# for phy
# kilosort
# cluster_KSLabel.tsv gives suggesting label by kilosort2.5
# cluster_group.tsv gives extra/external labels
# cluster_info.tsv gives sorting summary but only after phy save
# Therefor, we use 'KSLabel' to mark 'group' automatically
# And use 'info' to check if there is manual operation involved

# mountainsort
# cluster_si_unit_id.tsv gives suggesting label?
# cluster_group.tsv gives extr/external labels
# cluster_info.tsv gives sorting summary but only after phy save

# spykingcircus
# cluster_channel_group.tsv
# cluster_group.tsv
# cluster_si_unit_ids.tsv

for kshank in kshank_list:
    
    ks_shank_path = os.path.join(ks_path, '%s/sorter_output'%kshank)
    # ks_label = pd.read_csv(os.path.join(ks_shank_path, 'cluster_KSLabel.tsv'),
    #                         sep='\t')
    # ks_group = pd.read_csv(os.path.join(ks_shank_path, 'cluster_group.tsv'),
    #                       sep='\t')
    
    skcshank = 'spykingcircus'+kshank.split('kilo')[1]
    skcs_shank_path = os.path.join(skcs_path, '%s/output-phy'%skcshank)

    k_spike_times = np.load(os.path.join(ks_shank_path, 'spike_times.npy')).squeeze()
    k_spike_clusters = np.load(os.path.join(ks_shank_path, 'spike_clusters.npy')).squeeze()
    
    skc_spike_times = np.load(os.path.join(skcs_shank_path, 'spike_times.npy')).squeeze()
    skc_spike_clusters = np.load(os.path.join(skcs_shank_path, 'spike_clusters.npy')).squeeze()
     
    skc_sc_list = [converted(skc_spike_times)[skc_spike_clusters==i]
                   for i in np.unique(skc_spike_clusters)]
    
    good = []
    for cluster in np.unique(k_spike_clusters):
        k_st = converted(k_spike_times)[k_spike_clusters==cluster]
        
        # overlap_rate = max(
        #     [len(np.intersect1d(k_st, skc_st_i)) for skc_st_i in skc_sc_list])/len(k_st)
        # print(overlap_rate)
        overlap_rate = max(
            [len(np.intersect1d(k_st, skc_st_i)) for skc_st_i in skc_sc_list])/len(k_st)
        # print(overlap_rate)
        if overlap_rate>0.7:
            good.append(cluster)
    
    ks_group = pd.read_csv(os.path.join(ks_shank_path, 'cluster_KSLabel.tsv'),
                          sep='\t')
    
    ks_good = ks_group[ks_group['KSLabel']=='good'].index.tolist()
    both_good = np.intersect1d(good, ks_good)
    
    ks_group_update = ks_group
    ks_group_update.loc[both_good.tolist(), 'KSLabel'] = 'best'
    ks_group_update.loc[[i for i in good if i not in both_good], 'KSLabel'] = 'good'
    os.rename(os.path.join(ks_shank_path, 'cluster_KSLabel.tsv'),
              os.path.join(ks_shank_path, 'cluster_KSLabel.tsv.origin')) 
    ks_group_update.to_csv(os.path.join(ks_shank_path, 'cluster_KSLabel.tsv'),
                           sep='\t')

marktxt = next(os.walk(ks_path))[2]
os.rename(os.path.join(ks_path, marktxt[0]),
          os.path.join(ks_path, 'autokilo_voted'+marktxt[0].split('autokilo')[1])) 
