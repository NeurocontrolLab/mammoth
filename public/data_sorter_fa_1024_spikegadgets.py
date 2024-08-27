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
from probeinterface import Probe, ProbeGroup
from probeinterface import read_probeinterface


def sorting(root_dir, map_path, container_dir):
    #%% load probe
    probegroup2 = read_probeinterface(map_path) # read probe json

    pos_list = []
    shapes_list = []
    shape_params_list = []
    probe_2d = Probe(ndim=2, si_units='um')

    # reset channel postion
    for ind, i in enumerate(probegroup2.probes):
        pos = i.contact_positions[:,[0,2]].copy()
        pos[:,-1] = pos[:,-1]*-1+10*400*ind
        sha_par = i.contact_shape_params
        pos_list.extend(list(pos.copy()))
        shapes_list.extend(i.contact_shapes)
        shape_params_list.extend(sha_par)
        
    # build probe instance
    probe_2d.set_contacts(positions=pos_list, 
                          shapes=shapes_list, 
                          shape_params=shape_params_list)

    # set device
    probe_2d.set_device_channel_indices(np.arange(0, 1024))
    
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
                                            dtype=np.int16,num_channels=1024)
    recording = recording.set_probe(probe_2d)

    # make folder for kilo data saving
    data_path_ = os.path.join(root_dir, 'sorted_data', 'kilosort2_5_output')

    if not os.path.exists(data_path_):
        os.makedirs(data_path_)

    os.chdir(container_dir)
    wrong_flag = 0
    # slice recorded data to shank for sorting
    for ind, border in tqdm(enumerate(range(0, 1024, 128))):
        sp = os.path.join(data_path_, "kilo_shank_{}".format(ind))
        sliced_recording = recording.channel_slice(list(range(border, border+128)))
        
        if os.path.exists(sp):
            continue
        
        try:
            sorting = ss.run_sorter(sorter_name='kilosort2_5',
                                    recording=sliced_recording,     
                                    output_folder=sp,
                                    singularity_image=True,
                                    verbose = True,
                                    delete_tmp_files=False,
                                    delete_recording_dat = True,
                                    n_jobs = 12)
                                    #NT=64*256+64)
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
parser = argparse.ArgumentParser(
    prog = 'Nezha kilosort 2.5 pipeline', 
    description = 'Using 1024 flexible array with spikegadgets recording system',
    epilog = 'Neucyber-NC-2023-A-01'
)

parser.add_argument('-dp', '--data_path', 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Data_recording/20240319_centerOut_001')

parser.add_argument('-mp', '--map_path', 
                    default='/AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/nezha-gai.json')

parser.add_argument('-cp', '--container_path', 
                    default='/AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller')

args = parser.parse_args()

sorting(args.data_path, args.map_path, args.container_path)
