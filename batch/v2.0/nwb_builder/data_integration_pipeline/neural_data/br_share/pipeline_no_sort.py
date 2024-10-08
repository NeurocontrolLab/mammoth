#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

from user_input_entry import BRShare as bs
from SmartNeo.user_layer.dict_to_neo import templat_neo
from SmartNeo.interface_layer.nwb_interface import NWBInterface
import readTrodesExtractedDataFile3 as trodesReader
from probeinterface import read_probeinterface
from probeinterface import Probe, ProbeGroup

import neo
import joblib
import argparse
import quantities as pq
import glob
import os
from neo import io
import yaml
import copy
import numpy as np
import sys
from tqdm import tqdm
import pandas as pd
import brpylib

os.chdir('/AMAX/cuihe_lab/cuilab_share/Abel/Code/nwb_builder/data_integration_pipeline/neural_data/br_share')

#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument("-f", "--file", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Abel/Data_recording/20240705_centerOut_001', 
                    metavar='/the/path/your/data/located/in', help='input folder')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Abel/Data_recording/20240705_centerOut_001', 
                    metavar='/the/path/you/want/to/save', help='output folder')

parser.add_argument('-mp', '--map_path', 
                    default='/home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Spike_sorting/SN+11386-000049.cmp')

args = parser.parse_args()
raw_dirname = args.file

#%% load template
FILEPATH = os.path.dirname(os.path.abspath(__file__))
Template = yaml.safe_load(open(os.path.join(FILEPATH,'template.yml')))
Template['LFP'] = templat_neo['ana']
Template['RecordingSystemEvent'] = templat_neo['event']
Template['Spike'] = {'SorterName' : {},'kargs' : {}}

map_path = args.map_path
with open(map_path,'r') as f:
    probe_info = f.readlines()[14::]

probe_info = [i.split() for i in probe_info][0:-1]

array_name = ['elec'] if '-' not in probe_info[0][-1]\
    else list(set([i[-1].split('-')[0] for i in probe_info]))

array_name.sort()

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
    # test2 = []
    # add electrodes to the electrode table
    for ielec in probe_info:
        if array not in ielec[-1]:
            continue
    
        positions.append([float(ielec[0])*400, float(ielec[1])*400])
        device_channel.append((ord(ielec[2])-65)*32+int(ielec[3])-1)
        # test2.append(ielec[2]+ielec[3])
        
    
    probe_2d.set_contacts(positions=np.array(positions), 
                          shapes='circle', 
                          shape_params={'radius': 20})
    probe_2d.set_device_channel_indices(device_channel)
    probe_2d.create_auto_shape(probe_type='tip')
    probegroup.add_probe(probe_2d)
    
rs_pos = np.concatenate([i.contact_positions for i in probegroup.probes])

InputList = []

#%% convert Spike
InputData = copy.deepcopy(Template)
InputData['RecordingSystemEvent'] = 'null'
InputData['TCR'] = 'null'
InputData['SBP'] = 'null'
InputData['LFP'] = 'null'

InputData['name'] = 'kilosort2.5'
InputData['Spike'] = {}

kilo_path = os.path.join(raw_dirname,'sorted_data','kilosort2_5_output')

walk_file = [j for j in os.walk(raw_dirname)]

for f_l in walk_file:
    rec_name = [f_n for f_n in f_l[2] if '.ns6' in f_n]
    if len(rec_name) !=0:
        break
datafile = os.path.join(f_l[0], rec_name[0])

nsx_file = brpylib.NsxFile(str(datafile))

# Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
cont_data = nsx_file.getdata(full_timestamps=True)

# Close the nsx file now that all data is out
nsx_file.close()

data = np.concatenate(cont_data['data'],1).T
timestamp = np.concatenate([i['Timestamp'] for i in cont_data['data_headers']])

import numpy as np

fs = 30000.0
lowcut = 300.0
highcut = 6000.0

from scipy.signal import butter, lfilter
import scipy

def butter_bandpass(lowcut, highcut, fs, order=5):
    return butter(order, [lowcut, highcut], fs=fs, btype='band')

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

#%% convert TCR
InputData = copy.deepcopy(Template)
InputData['RecordingSystemEvent'] = 'null'
InputData['Spike'] = 'null'
InputData['SBP'] = 'null'
InputData['LFP'] = 'null'

InputData['name'] = 'TCR'
InputData['TCR'] = {}

# a = np.memmap(os.path.join(raw_dirname, rec_file, tcr_path,data_name), dtype=np.int16, mode='r+').reshape((-1,1024))

fs = 30000.0
lowcut = 300.0
highcut = 6000.0


from scipy.signal import find_peaks


for i, p in tqdm(enumerate(probegroup.probes)):
    
    shank_ind = str(i)
    p = probegroup.probes[int(shank_ind)]
    
    for ch_ind, chn in tqdm(enumerate(p.device_channel_indices)):
        
        f_ch = butter_bandpass_filter(data[:,chn], lowcut, highcut, fs, order=5)
        ch1_mad = scipy.stats.median_abs_deviation(f_ch)
        thred = np.median(f_ch)-(3.5/0.6745)*ch1_mad
        peaks, _ = find_peaks(-f_ch, height=-thred)
        ind=np.array([peaks-24+i for i in range(64)]).T
        ch_spike = f_ch[ind[1:-10]]
        
        spike_description = {'clu':0.0,
                             'chn':float(chn),
                             'mean_waveform':ch_spike.mean(0).astype(float),
                             'pos' : p.contact_positions[ch_ind],
                             'electrode' : float(shank_ind),
                             'annotations' : '',
                             'chn_meta' : ''}
        sampling_rate = 30000.0
        ptp_t = timestamp[peaks]/1e9
        
        if not str(chn) in InputData['TCR']:
            InputData['TCR'][str(chn)] = {}
        InputData['TCR'][str(chn)]['0'] = {}
        InputData['TCR'][str(chn)]['0']['spk'] = templat_neo['spk'].copy()
        InputData['TCR'][str(chn)]['0']['spk']['times'] = ptp_t * pq.s
        InputData['TCR'][str(chn)]['0']['spk']['t_stop'] = timestamp[-1]/1e9 * pq.s
        InputData['TCR'][str(chn)]['0']['spk']['t_start'] = timestamp[0]/1e9 * pq.s
        InputData['TCR'][str(chn)]['0']['spk']['sampling_rate'] = float(sampling_rate) * pq.Hz
        InputData['TCR'][str(chn)]['0']['spk']['description'] = spike_description
    
InputList.append(InputData)

#%% convert RecordingSystemEvent
walk_file = [j for j in os.walk(raw_dirname)]

for f_l in walk_file:
    rec_name = [f_n for f_n in f_l[2] if ('.nev' in f_n) and ('NSP' in f_n)]
    if len(rec_name) !=0:
        break
datafile = os.path.join(f_l[0], rec_name[0])

nev_file = brpylib.NevFile(str(datafile))

# Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
nsp_data = nev_file.getdata()

# Close the nsx file now that all data is out
nev_file.close()

# data = np.concatenate(lfp_data['data'],1).T
# timestamp = np.concatenate([i['Timestamp'] for i in lfp_data['data_headers']])/1e9
    
InputData = copy.deepcopy(Template)   

InputData['LFP'] = 'null'
InputData['Spike'] = 'null'
InputData['TCR'] = 'null'
InputData['SBP'] = 'null'

ptp_t = np.array(nsp_data['digital_events']['TimeStamps'])/1e9
InputData['RecordingSystemEvent']['times'] = ptp_t*pq.sec
InputData['RecordingSystemEvent']['labels'] = list(np.array(nsp_data['digital_events']['UnparsedData']).astype('float'))
InputData['RecordingSystemEvent']['description'] = ''
InputData['name'] = 'RecordingSystemEvent'

InputList.append(InputData)


#%% operate the dicts
neuralblock = bs.data_input(user_data = InputList, index = 1, name = 'neural_data')

del InputList
import gc
gc.collect()

#%% save to appointed path
if not os.path.exists(os.path.join(args.output,'formatted_data')):
    os.mkdir(os.path.join(args.output,'formatted_data'))
    
nwb_saver = NWBInterface()
nwb_saver.save_nwb(blockdata = neuralblock, 
                   filename = os.path.join(args.output, 'formatted_data', 'neural_data_no_sort.nwb'))
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
