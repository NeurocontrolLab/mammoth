#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 15:49:09 2024

@author: cuilab
"""

import os
import brpylib
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
        rec_name = [f_n for f_n in f_l[2] if '.ns6' in f_n]
        if len(rec_name) !=0:
            break
    datafile = os.path.join(f_l[0], rec_name[0])

    # find temp.bin 
    if not os.path.exists(os.path.join(root_dir, 'raw_data', 'temp.bin')):
        print('make bin file')
        # Open file and extract headers
        nsx_file = brpylib.NsxFile(str(datafile))
        
        # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
        cont_data = nsx_file.getdata(full_timestamps=True)
        
        # Close the nsx file now that all data is out
        nsx_file.close()
        
        data = np.concatenate(cont_data['data'], 1).T
        
        fp = np.memmap(os.path.join(root_dir, 'raw_data', 'temp.bin'), 
                       dtype = data.dtype, 
                       mode='w+', shape = data.shape) # creat memmap
        
        batch_chn = 32
        nchannels_per_array = data.shape[-1]
        
        for j in range(0,nchannels_per_array, batch_chn):
                
            for k in tqdm(range(0,len(data), 20000)):
                if (k+20000)<len(data):
                    fp[k:k+20000,j:j+batch_chn] = data[k:k+20000,j:j+batch_chn]
                else:
                    fp[k::,j:j+batch_chn] = data[k::,j:j+batch_chn]
    else:
        # Open file and extract headers
        nsx_file = brpylib.NsxFile(str(datafile))
        
        # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
        cont_data = nsx_file.getdata(full_timestamps=True)
        
        # Close the nsx file now that all data is out
        nsx_file.close()
        data = cont_data['data'][0]

    recording = se.BinaryRecordingExtractor(os.path.join(root_dir, 'raw_data', 'temp.bin'),
                                            sampling_frequency=30000,
                                            dtype=data.dtype,
                                            num_channels=cont_data['data'][0].shape[0])

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

    # if wrong_flag == 0:
    #     with open(os.path.join(data_path_, 'autokilo.txt'), 'w') as file:
    #         file.write("")    
    # else:
    #     with open(os.path.join(data_path_, 'wrong_shanks_%d.txt' % wrong_flag), 'w') as file:
    #         file.write("")  


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