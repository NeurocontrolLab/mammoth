# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import spikeinterface.extractors as se
import spikeinterface.sorters as ss
import spikeinterface.core as sc
from spikeinterface.exporters import export_to_phy
from probeinterface import Probe, ProbeGroup
from probeinterface import read_probeinterface
import numpy as np
from tqdm import tqdm
import shutil

#%% perform sorting
import os
import argparse
os.chdir('/AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller')

parser = argparse.ArgumentParser(
    prog = 'Nezha mountainsort 5 pipeline', 
    description = 'Using 1024 flexible array with spikegadgets recording system',
    epilog = 'Neucyber-NC-2023-A-01'
)

parser.add_argument('-dp', '--data_path', 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Data_recording/20240319_centerOut_001')

parser.add_argument('-mp', '--map_path', 
                    default='/AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/nezha-gai.json')

args = parser.parse_args()

# data_path = r'/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001'
data_path = args.data_path
map_path = args.map_path
# target_path = args.target_path

probegroup2 = read_probeinterface(map_path) # read probe json

probegroup = ProbeGroup()
pos_list = []
shapes_list = []
shape_params_list = []
p = []
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
probe_2d.set_device_channel_indices(np.arange(0,1024))
    

# find .rec file

walk_file = [j for j in os.walk(data_path)]
for f_l in walk_file:
    rec_name = [f_n for f_n in f_l[2] if '.rec' in f_n]
    if len(rec_name) !=0:
        break
rec_file = os.path.join(f_l[0].replace(' ','\ '), rec_name[0])

'''
command = 'ulimit -n 10240 && module load trodes/2.4.2 && trodesexport -rec {} -lfp -lfplowpass 300 -dio -kilosort -spikeband -spikehighpass 300 -spikelowpass 6000 -spikes -thresh 50'.format(rec_file)
# command = r'~/Trodes_2-4-1_Ubuntu2004/trodesexport -rec {} -dio -kilosort'.format(rec_file)
# call bash command
return_info = subprocess.Popen(command, shell=True,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
return_info.wait()
'''

# find kilosort file
for kilo_data in os.listdir(f_l[0].replace(' ','\ ')):
    if 'kilosort' in kilo_data:
        break
    
kilo_path = os.path.join(f_l[0].replace(' ','\ '),kilo_data)

for kilo_data in os.listdir(kilo_path):
    if 'group0.dat' in kilo_data:
        break
    
kilo_path = os.path.join(kilo_path,kilo_data)

# build recording instance for sorting
recording = se.BinaryRecordingExtractor(kilo_path,
                                 sampling_frequency=30000,
                                 dtype=np.int16,num_channels=1024)


recording = recording.set_probe(probe_2d)

# make folder for kilo data saving
data_path_ = os.path.join(data_path,'sorted_data','spykingcircus_output')

if not os.path.exists(data_path_):
    os.makedirs(data_path_)


wrong_flag = 0
# slice recorded data to shank for sorting
for ind, border in tqdm(enumerate(range(0,1024,128))):
    sp = os.path.join(data_path_,"spykingcircus_shank_{}".format(ind))
    sliced_recording = recording.channel_slice(list(range(border,border+128)))
    # if os.path.exists(sp):
    #     continue

    try:
       
        sorting = ss.run_sorter(sorter_name='spykingcircus',
                                singularity_image="spyking-circus-base",
                                recording=sliced_recording,     
                                output_folder=sp)
      
        sliced_recording.annotate(is_filtered=True)
        
        wep = os.path.join(sp, 'waveforms')
        phyp = os.path.join(sp, 'output-phy')
        if os.path.exists(wep):
            shutil.rmtree(wep)
        
        if os.path.exists(phyp):
            shutil.rmtree(phyp)
            
        we = sc.extract_waveforms(sliced_recording, sorting, wep, sparse=True)
        export_to_phy(we, output_folder=phyp)
        
    except:
        os.rename(sp,sp+'_wrong') # if something wrong happened, change the folder name
        wrong_flag += 1

if wrong_flag == 0:
    with open(os.path.join(data_path_, 'autoskc.txt'), 'w') as file:
        file.write("")      
else:
    with open(os.path.join(data_path_, 'wrong_shanks_%d.txt' % wrong_flag), 'w') as file:
        file.write("")  
