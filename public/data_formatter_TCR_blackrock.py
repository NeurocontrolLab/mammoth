#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

import os
import yaml
import copy
from tqdm import tqdm
import numpy as np
import scipy
from scipy.signal import butter, lfilter, find_peaks
import brpylib
import argparse
from get_probe_bohr_utah96 import get as get_probe
from user_input_entry_collection import BRShare as bs
from SmartNeo.user_layer.dict_to_neo import templat_neo
from SmartNeo.interface_layer.nwb_interface import NWBInterface


def format_file(root_dir, map_path, output_dir):
    #%% load template
    FILEPATH = os.path.dirname(os.path.abspath(__file__))
    Template = yaml.safe_load(open(os.path.join(FILEPATH,'template_neural_data.yml')))
    Template['LFP'] = templat_neo['ana']
    Template['RecordingSystemEvent'] = templat_neo['event']
    Template['Spike'] = {'SorterName' : {},'kargs' : {}}

    probegroup = get_probe(map_path)

    #%% read ns6 file
    walk_file = [j for j in os.walk(root_dir)]
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

    data = np.concatenate(cont_data['data'], 1).T
    timestamp = np.concatenate([i['Timestamp'] for i in cont_data['data_headers']])

    #%% filter
    def butter_bandpass(lowcut, highcut, fs, order=5):
        return butter(order, [lowcut, highcut], fs=fs, btype='band')

    def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        y = lfilter(b, a, data)
        return y
    

    #%% convert TCR
    InputList = []
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
    walk_file = [j for j in os.walk(root_dir)]

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
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        
    nwb_saver = NWBInterface()
    nwb_saver.save_nwb(blockdata = neuralblock, 
                       filename = os.path.join(output_dir, 'neural_data_no_sort.nwb'))


#%% 
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Abel/Data_recording/20240705_centerOut_001', 
                    metavar='/the/path/your/data/located/in', help='input folder')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Abel/Data_recording/20240705_centerOut_001/formatted_data', 
                    metavar='/the/path/you/want/to/save', help='output folder')

parser.add_argument('-mp', '--map_path', 
                    default='/home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Spike_sorting/SN+11386-000049.cmp')

args = parser.parse_args()

format_file(args.root, args.map_path, args.output)
