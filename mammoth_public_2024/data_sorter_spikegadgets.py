#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 15:49:09 2024

@author: cuilab
"""

import os
import argparse
import numpy as np
from tqdm import tqdm
import spikeinterface.extractors as se
import spikeinterface.sorters as ss
from probeinterface import read_probeinterface


def sorting(sorter, root_dir, map_path, output_dir, container_dir):
    #%% load probe
    probegroup = read_probeinterface(map_path) # read probe json
    
    #%% perform sorting
    # find .rec file
    walk_file = [j for j in os.walk(root_dir)]
    for f_l in walk_file:
        rec_name = [f_n for f_n in f_l[2] if '.rec' in f_n]
        if len(rec_name) !=0:
            break
   
    # find kilosort file
    for kilo_data in os.listdir(f_l[0].replace(' ','\ ')):
        if 'kilosort' in kilo_data:
            break
        
    kilo_path = os.path.join(f_l[0].replace(' ','\ '), kilo_data)

    for kilo_data in os.listdir(kilo_path):
        if 'group0.dat' in kilo_data:
            break
        
    kilo_path = os.path.join(kilo_path, kilo_data)

    # build recording instance for sorting
    recording = se.BinaryRecordingExtractor(kilo_path,
                                            sampling_frequency=30000,
                                            dtype=np.int16,
                                            num_channels=int(len(probegroup.probes[-1].shank_ids)*len(probegroup.probes)))
    recording = recording.set_probegroup(probegroup)

    # make folder for kilo data saving
    data_path_ = os.path.join(output_dir, '{}_output'.format(sorter))

    if not os.path.exists(data_path_):
        os.makedirs(data_path_)

    os.chdir(container_dir)
    wrong_flag = 0
    # slice recorded data to shank for sorting
    for ind, p in tqdm(enumerate(probegroup.probes)):
        sp = os.path.join(data_path_, "{}_shank_{}".format(sorter[:4], ind))
        sliced_recording = recording.channel_slice(list(p.device_channel_indices))
        
        if os.path.exists(sp):
            continue
        
        try:
            sorting = ss.run_sorter(sorter_name=sorter,
                                    recording=sliced_recording,     
                                    output_folder=sp,
                                    singularity_image=True,
                                    verbose = True,
                                    delete_tmp_files=False,
                                    delete_recording_dat = True,
                                    n_jobs = 12)
        except:
            os.rename(sp,sp+'_wrong') # if something wrong happened, change the folder name
            wrong_flag += 1

    if wrong_flag == 0:
        with open(os.path.join(data_path_, 'autokilo.txt'), 'w') as file:
            file.write("")    
    else:
        with open(os.path.join(data_path_, 'wrong_shanks_%d.txt' % wrong_flag), 'w') as file:
            file.write("")  


#%%
parser = argparse.ArgumentParser()

parser.add_argument('-sorter', '--sorter_name', 
                    default='kilosort2_5')

parser.add_argument('-r', '--root', 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Data_recording/20240319_centerOut_001')

parser.add_argument('-mp', '--map_path', 
                    default='/AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/nezha-gai.json')

parser.add_argument('-o', '--output',
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Data_recording/20240319_centerOut_001/sorted_data')

parser.add_argument('-cp', '--container', 
                    default='/AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller')

args = parser.parse_args()

sorting(args.sorter_name, args.root, args.map_path, args.output, args.container)
