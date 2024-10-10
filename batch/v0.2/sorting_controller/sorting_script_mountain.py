# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import spikeinterface.extractors as se
import spikeinterface.sorters as ss
from probeinterface import Probe, ProbeGroup
import spikeinterface.preprocessing as spre
import numpy as np
from tqdm import tqdm
import os
import brpylib
import argparse
import spikeinterface as si

parser = argparse.ArgumentParser(
    prog = 'Bohr kilosort 2.5 pipeline', 
    description = 'Using 2*96 Utah array with blackrock recording system',
    epilog = 'Neucyber-NC-2023-A-01'
)

parser.add_argument('-dp', '--data_path', 
                    default='/home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Data_recording/20231213')

parser.add_argument('-mp', '--map_path', 
                    default='/home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Spike_sorting/SN+11386-000049.cmp')

parser.add_argument('-tp', '--target_path', 
                    default='/home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Data_recording/20231213')

args = parser.parse_args()

data_path = args.data_path
map_path = args.map_path
target_path = args.target_path

walk_file = [j for j in os.walk(data_path)]

for f_l in walk_file:
    rec_name = [f_n for f_n in f_l[2] if '.ns6' in f_n]
    if len(rec_name) !=0:
        break
datafile = os.path.join(f_l[0], rec_name[0])

with open(map_path,'r') as f:
    probe_info = f.readlines()[14::]

probe_info = [i.split() for i in probe_info][0:-1]
array_name = ['elec'] if '-' not in probe_info[0][-1]\
    else list(set([i[-1].split('-')[0] for i in probe_info]))

nchannels_per_array = len(probe_info)
electrode_counter = []
probegroup = ProbeGroup()

for array in array_name:
    # create an electrode group for this shank
    
    probe_2d = Probe(ndim=2, si_units='um')
    probe_2d.annotate(
        name = array, 
        manufacturer="blackrock microsystem",
        escription = 'one 96 Utah array'
        )
    positions = []
    device_channel = []
    # add electrodes to the electrode table
    for ielec in probe_info:
        if array not in ielec[-1]:
            continue
    
        positions.append([float(ielec[0])*400, float(ielec[1])*400])
        device_channel.append((ord(ielec[2])-65)*32+int(ielec[3])-1)
        
    
    probe_2d.set_contacts(positions=np.array(positions), 
                          shapes='circle', 
                          shape_params={'radius': 20})
    probe_2d.set_device_channel_indices(device_channel)
    probe_2d.create_auto_shape(probe_type='tip')
    probegroup.add_probe(probe_2d)


if not os.path.exists(os.path.join(data_path, 'temp.bin')):
  print('make bin file')
  # Open file and extract headers
  nsx_file = brpylib.NsxFile(str(datafile))
  
  # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
  cont_data = nsx_file.getdata(full_timestamps=True)
  
  # Close the nsx file now that all data is out
  nsx_file.close()
  
  data = np.concatenate(cont_data['data'],1).T
  
  fp = np.memmap(
      os.path.join(data_path, 'temp.bin'), 
      dtype = data.dtype, 
      mode='w+', shape = data.shape
      ) # creat memmap
  
  batch_chn = 32
  nchannels_per_array = data.shape[-1]
  
  for j in range(0,nchannels_per_array,batch_chn):
          
      for k in tqdm(range(0,len(data),20000)):
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

recording = se.BinaryRecordingExtractor(os.path.join(data_path, 'temp.bin'),
                                 sampling_frequency=30000,
                                 dtype=data.dtype,num_channels=cont_data['data'][0].shape[0])

recording = recording.set_probegroup(probegroup)
# location = np.array([nwbfile.electrodes['x'][:],nwbfile.electrodes['y'][:]]).T
# recording.set_property('location',location)

# nsfilePath = glob(os.path.join(data_path, nd, '*.ns3'))[0]
# shutil.copy(nsfilePath, os.path.split(nsfilePath)[-1])
if not os.path.exists(os.path.join(data_path, 'mountain5_output')):
  os.mkdir(os.path.join(data_path, 'mountain5_output'))

sp1 = os.path.join(data_path, 'mountain5_output')

recording_filtered = spre.bandpass_filter(recording, freq_min=300, freq_max=6000, dtype=data.dtype)
recording_preprocessed: si.BaseRecording = spre.whiten(recording_filtered,int_scale=200)

# slice recorded data to shank for sorting
for ind, p in tqdm(enumerate(probegroup.probes)):
    sp = os.path.join(sp1,"ms_shank_{}".format(ind))
    sliced_recording = recording.channel_slice(list(p.device_channel_indices))
    if os.path.exists(sp):
        continue
    try:
        sorting = ss.run_sorter(sorter_name='mountainsort5',
            recording=sliced_recording,     
            output_folder=sp,
            singularity_image=True,
            verbose = True)
            #NT=64*256+64)
    except:
        os.rename(sp,sp+'_wrong') # if something wrong happened, change the folder name

# shutil.move("neural_data.nwb",os.path.join(target_folder, "neural_data.nwb"))
# shutil.move("kilo3",os.path.join(target_folder, "kilo3"))
    
    