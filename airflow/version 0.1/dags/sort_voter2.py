#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 13:35:05 2024

@author: chenyun
"""


import spikeinterface.extractors as se
import spikeinterface.comparison as sc 
import numpy as np
from tqdm import tqdm
import pandas as pd 

#%% perform sorting
import os
import argparse
os.chdir('/AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller')

parser = argparse.ArgumentParser(
    prog = 'Nezha kilosort 2.5 pipeline', 
    description = 'Using 1024 flexible array with spikegadgets recording system',
    epilog = 'Neucyber-NC-2023-A-01'
)

parser.add_argument('-dp', '--data_path',
                    default=r'/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001')

args = parser.parse_args()

data_path = args.data_path


# make folder for kilo data saving
data_path_ks = os.path.join(data_path,'sorted_data','kilosort2_5_output')
data_path_sc = os.path.join(data_path,'sorted_data','spykingcircus_output')

# slice recorded data to shank for sorting
for ind, border in tqdm(enumerate(range(0,1024,128))):
    sp_ks = os.path.join(data_path_ks, "kilo_shank_{}/sorter_output".format(ind))
    sp_sc = os.path.join(data_path_sc, 
                         "spykingcircus_shank_{}/sorter_output".format(ind))
   
    # Then run two spike sorters and compare their outputs.
    sorting_ks = se.read_kilosort(folder_path=sp_ks)
    sorting_sc = se.read_spykingcircus(folder_path=sp_sc)
    
    # Run the comparison
    # Let's see how to inspect and access this matching.
    cmp_KS_SC = sc.compare_multiple_sorters(
        sorting_list=[sorting_sc, sorting_ks],
        name_list=['SpykingCircus', 'Kilosort_2_5']
    )
    agr = cmp_KS_SC.get_agreement_sorting(minimum_agreement_count=2)
    # agr_id = agr.get_unit_ids()
    agr_ks_id = [int(agr.get_property('unit_ids')[i]['Kilosort_2_5'])
             for i in range(len(agr.get_property('unit_ids')))]
    
    ks_group = pd.read_csv(os.path.join(sp_ks, 'cluster_KSLabel.tsv'),
                           sep='\t')
    
    ks_good = ks_group[ks_group['KSLabel']=='good'].index.tolist()
    
    ks_group_update = ks_group
    ks_group_update.loc[agr_ks_id, 'KSLabel'] = 'best'
    os.rename(os.path.join(sp_ks, 'cluster_KSLabel.tsv'),
              os.path.join(sp_ks, 'cluster_KSLabel.tsv.origin')) 
    ks_group_update.to_csv(os.path.join(sp_ks, 'cluster_KSLabel.tsv'),
                           sep='\t')

marktxt = next(os.walk(data_path_ks))[2]
os.rename(os.path.join(data_path_ks, marktxt[0]),
          os.path.join(data_path_ks, 'autokilo_voted'+marktxt[0].split('autokilo')[1])) 

    