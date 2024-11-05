#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import yaml
import copy
from tqdm import tqdm
import scipy
from scipy.signal import butter, lfilter, find_peaks
import numpy as np
import pandas as pd
import brpylib
import argparse
import quantities as pq
from SmartNeo.user_layer.dict_to_neo import templat_neo
from SmartNeo.interface_layer.nwb_interface import NWBInterface
from probeinterface import read_probeinterface
from dependencies.user_input_entry_collection import BRShare as bs

#%% define functions
def butter_bandpass(lowcut, highcut, fs, order=5):
    return butter(order, [lowcut, highcut], fs=fs, btype='band')


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def get_timestamp(root_dir):
    walk_file = [j for j in os.walk(root_dir)]

    try:
        for f_l in walk_file:
            rec_name = [f_n for f_n in f_l[2] if '.ns6' in f_n]
            if len(rec_name) !=0:
                break
        datafile = os.path.join(f_l[0], rec_name[0])
    except:
        sys.exit()

    nsx_file = brpylib.NsxFile(str(datafile))

    # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
    cont_data = nsx_file.getdata(full_timestamps=True)

    # Close the nsx file now that all data is out
    nsx_file.close()

    data = np.concatenate(cont_data['data'],1).T
    timestamp = np.concatenate([i['Timestamp'] for i in cont_data['data_headers']])

    return data, timestamp


def convert_spike(data, timestamp, sorter_output_path, data_template, probegroup, bandpass_params):
    InputData = copy.deepcopy(data_template)
    InputData['RecordingSystemEvent'] = 'null'
    InputData['TCR'] = 'null'
    InputData['SBP'] = 'null'
    InputData['LFP'] = 'null'

    InputData['name'] = sorter_output_path.split(os.sep)[-1].lower().replace('_', '.')
    InputData['Spike'] = {}

    fs, lowcut, highcut = bandpass_params
    # fs = 30000.0
    # lowcut = 300.0
    # highcut = 6000.0

    # check if cluster_info.tsv exists
    for i in os.listdir(sorter_output_path):
        data_path = os.path.join(sorter_output_path, i, 'sorter_output')
        assert os.path.exists(os.path.join(data_path, 'cluster_info.tsv')), "cannot find cluster info"

    # assemble data    
    for i in os.listdir(sorter_output_path):
        shank_ind = i.split("_")[2]
        data_path = os.path.join(sorter_output_path, i, 'sorter_output')
        cluster_info = pd.read_csv(os.path.join(data_path, 'cluster_info.tsv'), sep="\t")
        spike_times = np.load(os.path.join(data_path, 'spike_times.npy')).squeeze()
        spike_clusters = np.load(os.path.join(data_path, 'spike_clusters.npy')).squeeze()
        #whitened_data = np.array(np.memmap(os.path.join(data_path, 'temp_wh.dat'),mode='r',dtype=np.int16).reshape((-1,96)))
        
        p = probegroup.probes[int(shank_ind)]
        
        for clu in tqdm(np.unique(spike_clusters)):
            ind = np.where(cluster_info.cluster_id==clu)[0][0]
            ci = cluster_info.iloc[ind].to_dict()
            st = spike_times[spike_clusters==clu]
            # choose 24 sampling points before and 40, sampling rate = 30k Hz
            # => [-0.5ms, +1 ms], referring to WaveClus
            swi = [st-24+i for i in range(64)]  
            bch = p.device_channel_indices[ci['ch']]
            
            # mean_waveform = whitened_data[:,ci['ch']][swi].squeeze().mean(1).astype(float)
                
            if not 'good' == ci['group']:
                continue
            
            f_ch = butter_bandpass_filter(data[:, bch], lowcut, highcut, fs, order=5)
            mean_waveform = f_ch[swi].squeeze().mean(1).astype(float)
            # mean_waveform = data[swi,bch].squeeze().mean(1).astype(float)
            
            spike_description = {'clu': float(ci['cluster_id']),
                                 'chn': float(bch),
                                 'mean_waveform': mean_waveform,
                                 'pos': p.contact_positions[ci['ch']],
                                 'electrode': float(shank_ind),
                                 'annotations': '',
                                 'chn_meta': ci}
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
        
    return InputData
    

def convert_TCR(data, timestamp, data_template, probegroup, bandpass_params):    
    InputData = copy.deepcopy(data_template)
    InputData['RecordingSystemEvent'] = 'null'
    InputData['Spike'] = 'null'
    InputData['SBP'] = 'null'
    InputData['LFP'] = 'null'

    InputData['name'] = 'TCR'
    InputData['TCR'] = {}

    # a = np.memmap(os.path.join(raw_dirname, rec_file, tcr_path,data_name), dtype=np.int16, mode='r+').reshape((-1,1024))

    fs, lowcut, highcut = bandpass_params
    # fs = 30000.0
    # lowcut = 300.0
    # highcut = 6000.0

    for shank_ind, p in enumerate(probegroup.probes):
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
    return InputData


def convert_LFP(root_dir, data_template):
    walk_file = [j for j in os.walk(root_dir)]

    for f_l in walk_file:
        rec_name = [f_n for f_n in f_l[2] if '.ns2' in f_n]
        if len(rec_name) !=0:
            break
    
    if len(rec_name) >0:
        
        datafile = os.path.join(f_l[0], rec_name[0])

        lfp_file = brpylib.NsxFile(str(datafile))

        # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
        lfp_data = lfp_file.getdata(full_timestamps=True)

        # Close the nsx file now that all data is out
        lfp_file.close()

        data = np.concatenate(lfp_data['data'], 1).T
        timestamp = np.concatenate([i['Timestamp'] for i in lfp_data['data_headers']])/1e9
    
    else:
        for f_l in walk_file:
            rec_name = [f_n for f_n in f_l[2] if '.ns6' in f_n]
            if len(rec_name) !=0:
                break

        datafile = os.path.join(f_l[0], rec_name[0])

        raw_file = brpylib.NsxFile(str(datafile))

        lfp_data = raw_file.getdata(downsample=15)
        # Note: .ns3 file only from 0.3-250 Hz
        # ButterWorthFilter lowpass and highpass

        raw_file.close()

        data = np.concatenate(lfp_data['data'], 1).T
        timestamp = np.concatenate([i['Timestamp'] for i in lfp_data['data_headers']])/1e9


    # LFPRecordingParam = blackrock_data.extended_headers
    InputData = copy.deepcopy(data_template)
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

    return InputData


def convert_RSE(root_dir, data_template):
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
        
    InputData = copy.deepcopy(data_template)   
    InputData['LFP'] = 'null'
    InputData['Spike'] = 'null'
    InputData['TCR'] = 'null'
    InputData['SBP'] = 'null'

    ptp_t = np.array(nsp_data['digital_events']['TimeStamps'])/1e9
    InputData['RecordingSystemEvent']['times'] = ptp_t*pq.sec
    InputData['RecordingSystemEvent']['labels'] = list(np.array(nsp_data['digital_events']['UnparsedData']).astype('float'))
    InputData['RecordingSystemEvent']['description'] = ''
    InputData['name'] = 'RecordingSystemEvent'

    return InputData


def format_file(root_dir, map_path, output_dir, content_list, sorter=None):
    #%% prepare
    # load template
    FILEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    Template = yaml.safe_load(open(os.path.join(FILEPATH, 'dependencies', 'template_neural_data.yml')))
    
    # load probe
    probegroup = read_probeinterface(map_path)
    
    # get raw data path
    try:
        walk_file = [j for j in os.walk(root_dir)]
        for f_l in walk_file:
            rec_name = [f_n for f_n in f_l[2] if '.ns6' in f_n]
            if len(rec_name) !=0:
                break
        rec_dir = f_l[0]
    except:
        sys.exit()
    
    raw_dir = os.path.join(root_dir, rec_dir)

    # get timestamps
    data, timestamp = get_timestamp(raw_dir)

    # get sorter_output path
    if len([i for i in content_list if i.lower()=='spike'])>0:
        if not sorter:
            print("The data have not yet sorted.")
            sys.exit()

    if sorter:
        sorter = sorter.lower()
        if '.' in sorter:
            sorter = sorter.replace(".", "_")
        sorter_output_path = os.path.join(root_dir,'sorted_data','{}_output'.format(sorter))

    # set bandpass parameters
    fs = 30000.0
    lowcut = 300.0
    highcut = 6000.0
    bandpass_params = (fs, lowcut, highcut)
    
    # initialize data list
    InputList = []

    #%% convert Spike
    if len([i for i in content_list if i.lower()=='spike'])>0:
        Template['Spike'] = {'SorterName' : {}, 'kargs' : {}}
        InputData = convert_spike(data, timestamp, sorter_output_path, Template, probegroup, bandpass_params)
        InputList.append(InputData)
        print('spike included')

    #%% convert TCR
    if len([i for i in content_list if i.upper()=='TCR'])>0:
        InputData = convert_TCR(data, timestamp, Template, probegroup, bandpass_params)    
        InputList.append(InputData)
        print('TCR included')

    #%% convert LFP
    if len([i for i in content_list if i.upper()=='LFP'])>0:
        Template['LFP'] = templat_neo['ana']
        InputData = convert_LFP(raw_dir, Template)
        InputList.append(InputData)
        print('LFP included')

    #%% convert RecordingSystemEvent
    Template['RecordingSystemEvent'] = templat_neo['event']
    InputData = convert_RSE(raw_dir, Template)
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

    if len(content_list) == 1:
        nwb_saver.save_nwb(blockdata = neuralblock, 
                        filename = os.path.join(output_dir, 'neural_data_no_sort.nwb'))
    elif 'spike' in content_list:
         nwb_saver.save_nwb(blockdata = neuralblock, 
                        filename = os.path.join(output_dir, 'neural_data.nwb'))

#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20241011_interception_002', 
                    metavar='/the/path/your/data/located/in', help='input folder')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20241011_interception_002/formatted_data', 
                    metavar='/the/path/you/want/to/save', help='output folder')

# parser.add_argument('-mp', '--map_path', 
#                     default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Abel_Utah_64x4.json')

parser.add_argument('-mp', '--map_path', 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Bohr_Utah_96x2_PMd-M1_BlackRock.json')

parser.add_argument('-flag', '--sort_flag', 
                    default='0')

parser.add_argument('-sorter', '--sorter_name', 
                    default='kilosort2_5')


args = parser.parse_args()

if args.sort_flag == '1':
    content_list = ['spike', 'TCR', 'LFP']
elif args.sort_flag == '0':
    content_list = ['TCR']

format_file(args.root, args.map_path, args.output, content_list, args.sorter_name)        
        