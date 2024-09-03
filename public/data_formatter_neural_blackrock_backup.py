#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""


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
from scipy.signal import butter, lfilter, find_peaks
import scipy
from user_input_entry import BRShare as bs
from SmartNeo.user_layer.dict_to_neo import templat_neo
from SmartNeo.interface_layer.nwb_interface import NWBInterface
from get_probe_bohr_utah96 import get as get_probe


def format_file(root_dir, map_path, output_dir):
    #%% load template
    FILEPATH = os.path.dirname(os.path.abspath(__file__))
    Template = yaml.safe_load(open(os.path.join(FILEPATH,'template_neural_data.yml')))
    Template['LFP'] = templat_neo['ana']
    Template['RecordingSystemEvent'] = templat_neo['event']
    Template['Spike'] = {'SorterName' : {}, 'kargs' : {}}


    #%% load probe
    probegroup = get_probe(map_path)


    #%% define bandpass functions
    def butter_bandpass(lowcut, highcut, fs, order=5):
        return butter(order, [lowcut, highcut], fs=fs, btype='band')

    def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        y = lfilter(b, a, data)
        return y


    #%% Load timestamps
    InputList = []

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

    data = np.concatenate(cont_data['data'],1).T
    timestamp = np.concatenate([i['Timestamp'] for i in cont_data['data_headers']])

    #%% convert Spike
    InputData = copy.deepcopy(Template)
    InputData['RecordingSystemEvent'] = 'null'
    InputData['TCR'] = 'null'
    InputData['SBP'] = 'null'
    InputData['LFP'] = 'null'

    InputData['name'] = 'kilosort2.5'
    InputData['Spike'] = {}

    fs = 30000.0
    lowcut = 300.0
    highcut = 6000.0

    kilo_path = os.path.join(root_dir,'sorted_data','kilosort2_5_output')

    for i in os.listdir(kilo_path):
        data_path = os.path.join(kilo_path, i, 'sorter_output')
        assert os.path.exists(os.path.join(data_path, 'cluster_info.tsv')), "cannot find cluster info"
        
    for i in os.listdir(kilo_path):
        
        shank_ind = i.split("_")[2]
        data_path = os.path.join(kilo_path, i, 'sorter_output')
        cluster_info = pd.read_csv(os.path.join(data_path, 'cluster_info.tsv'),sep="\t")
        spike_times = np.load(os.path.join(data_path, 'spike_times.npy')).squeeze()
        spike_clusters = np.load(os.path.join(data_path, 'spike_clusters.npy')).squeeze()
        #whitened_data = np.array(np.memmap(os.path.join(data_path, 'temp_wh.dat'),mode='r',dtype=np.int16).reshape((-1,96)))
        
        p = probegroup.probes[int(shank_ind)]
        
        for clu in tqdm(np.unique(spike_clusters)):
            ind = np.where(cluster_info.cluster_id==clu)[0][0]
            ci = cluster_info.iloc[ind].to_dict()
            st = spike_times[spike_clusters==clu][0:-1]
            swi = [st-24+i for i in range(64)]
            bch = p.device_channel_indices[ci['ch']]
            
            # mean_waveform = whitened_data[:,ci['ch']][swi].squeeze().mean(1).astype(float)
                
            if not 'good' == ci['group']:
                continue
            
            f_ch = butter_bandpass_filter(data[:,bch], lowcut, highcut, fs, order=5)
            mean_waveform = f_ch[swi].squeeze().mean(1).astype(float)
            # mean_waveform = data[swi,bch].squeeze().mean(1).astype(float)
            
            spike_description = {'clu':float(ci['cluster_id']),
                                'chn':float(bch),
                                'mean_waveform':mean_waveform,
                                'pos' : p.contact_positions[ci['ch']],
                                'electrode' : float(shank_ind),
                                'annotations' : '',
                                'chn_meta' : ci}
            sampling_rate = 30000
            
            ptp_t = timestamp[st]/1e9
            if not 'shank ' + shank_ind in InputData['Spike']:
                InputData['Spike']['shank ' + shank_ind] = {}
            if not 'chn ' + str(ci['ch']) in InputData['Spike']['shank ' + shank_ind]:
                InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])] = {}
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])] = {}
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk'] = templat_neo['spk'].copy()
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk']['times'] = ptp_t * pq.s
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk']['t_stop'] = timestamp[-1]/1e9 * pq.s
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk']['t_start'] = timestamp[0]/1e9 * pq.s
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk']['sampling_rate'] = float(sampling_rate) * pq.Hz
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk']['description'] = spike_description

    InputList.append(InputData)

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

    for i in os.listdir(kilo_path):
        
        shank_ind = i.split("_")[2]
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

    #%% convert LFP
    walk_file = [j for j in os.walk(root_dir)]

    for f_l in walk_file:
        rec_name = [f_n for f_n in f_l[2] if '.ns2' in f_n]
        if len(rec_name) !=0:
            break
    datafile = os.path.join(f_l[0], rec_name[0])

    lfp_file = brpylib.NsxFile(str(datafile))

    # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
    lfp_data = lfp_file.getdata(full_timestamps=True)

    # Close the nsx file now that all data is out
    lfp_file.close()

    data = np.concatenate(lfp_data['data'],1).T
    timestamp = np.concatenate([i['Timestamp'] for i in lfp_data['data_headers']])/1e9

    # LFPRecordingParam = blackrock_data.extended_headers
    InputData = copy.deepcopy(Template)
    InputData['RecordingSystemEvent'] = 'null'
    InputData['Spike'] = 'null'
    InputData['SBP'] = 'null'
    InputData['TCR'] = 'null'

    InputData['LFP'] = templat_neo['irr'].copy()
    InputData['LFP']['signal'] = data*pq.uV
    InputData['LFP']['times'] = timestamp*pq.s
    InputData['LFP']['t_start'] = timestamp[0]*pq.s
    InputData['LFP']['description'] = ''
    InputData['name'] = 'LFP'
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
                       filename = os.path.join(output_dir, 'neural_data.nwb'))
        

#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Data_recording/20240402_interception_001', 
                    metavar='/the/path/your/data/located/in', help='input folder')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Data_recording/20240402_interception_001/formatted_data', 
                    metavar='/the/path/you/want/to/save', help='output folder')

parser.add_argument('-mp', '--map_path', 
                    default='/home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Spike_sorting/SN+11386-000049.cmp')

args = parser.parse_args()

format_file(args.root, args.map_path, args.output)        
        