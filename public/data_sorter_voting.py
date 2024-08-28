#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 9:44:09 2024

@author: cuilab
"""


#%% 
import os
from tqdm import tqdm
import pandas as pd
import argparse 
import spikeinterface.extractors as se
import spikeinterface.comparison as sc 


def compare(ks_dir, skc_dir):
    
    shank_num = len(next(os.walk(ks_dir))[1])
    # slice recorded data to shank for sorting
    for ind in tqdm(range(shank_num)):
        sp_ks = os.path.join(ks_dir, "kilo_shank_{}/sorter_output".format(ind))
        sp_sc = os.path.join(skc_dir, "spykingcircus_shank_{}/sorter_output".format(ind))
    
        # Then run two spike sorters and compare their outputs.
        sorting_ks = se.read_kilosort(folder_path=sp_ks)
        sorting_sc = se.read_spykingcircus(folder_path=sp_sc)
        
        # Run the comparison
        # Let's see how to inspect and access this matching.
        cmp_KS_SC = sc.compare_multiple_sorters(sorting_list=[sorting_sc, sorting_ks],
                                                name_list=['SpykingCircus', 'Kilosort_2_5'])
        agr = cmp_KS_SC.get_agreement_sorting(minimum_agreement_count=2)
        # agr_id = agr.get_unit_ids()
        agr_ks_id = [int(agr.get_property('unit_ids')[i]['Kilosort_2_5'])
                     for i in range(len(agr.get_property('unit_ids')))]
        
        ks_group = pd.read_csv(os.path.join(sp_ks, 'cluster_KSLabel.tsv'), sep='\t')
        
        ks_group_update = ks_group
        ks_group_update.loc[agr_ks_id, 'KSLabel'] = 'best'
        os.rename(os.path.join(sp_ks, 'cluster_KSLabel.tsv'),
                  os.path.join(sp_ks, 'cluster_KSLabel.tsv.origin')) 
        ks_group_update.to_csv(os.path.join(sp_ks, 'cluster_KSLabel.tsv'),
                               sep='\t')

    marktxt = next(os.walk(data_path_ks))[2]
    os.rename(os.path.join(data_path_ks, marktxt[0]),
              os.path.join(data_path_ks, 'autokilo_voted'+marktxt[0].split('autokilo')[1])) 


# os.chdir('/AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller')
parser = argparse.ArgumentParser(
    prog = 'Vote for sorting results', 
    description = 'Compare sorting results from Kilosort 2.5 and SpykingCircus',
)

parser.add_argument('-dp', '--data_path',
                    default='/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001')

args = parser.parse_args()

data_path = args.data_path

# make folder for kilo data saving
data_path_ks = os.path.join(data_path,'sorted_data','kilosort2_5_output')
data_path_sc = os.path.join(data_path,'sorted_data','spykingcircus_output')

compare(data_path_ks, data_path_sc)
