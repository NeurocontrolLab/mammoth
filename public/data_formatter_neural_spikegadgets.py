#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

import os
import sys
import yaml
import copy
from tqdm import tqdm
import scipy
from scipy.signal import butter, lfilter, find_peaks
import numpy as np
import pandas as pd
import argparse
import quantities as pq
import readTrodesExtractedDataFile3 as trodesReader
from probeinterface import read_probeinterface
from user_input_entry import BRShare as bs
from SmartNeo.user_layer.dict_to_neo import templat_neo
from SmartNeo.interface_layer.nwb_interface import NWBInterface



#%% define functions
def butter_bandpass(lowcut, highcut, fs, order=5):
    return butter(order, [lowcut, highcut], fs=fs, btype='band')


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def get_timestamp(root_dir):
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
        cluster_info = pd.read_csv(os.path.join(data_path, 'cluster_info.tsv'),sep="\t")
        spike_times = np.load(os.path.join(data_path, 'spike_times.npy')).squeeze()
        spike_clusters = np.load(os.path.join(data_path, 'spike_clusters.npy')).squeeze()
        # whitened_data = np.array(np.memmap(os.path.join(data_path, 'temp_wh.dat'),mode='r',dtype=np.int16).reshape((-1,256)))
        
        p = probegroup.probes[int(shank_ind)]

        for clu in tqdm(np.unique(spike_clusters)):
            ind = np.where(cluster_info.cluster_id==clu)[0][0]
            ci = cluster_info.iloc[ind].to_dict()
            st = spike_times[spike_clusters==clu]
            # choose 24 sampling points before and 40, sampling rate = 30k Hz
            # => [-0.5ms, +1 ms], referring to WaveClus
            swi = list(np.array([st-24+i for i in range(64)])[:, 0:-10])
            if not 'good' == ci['group']:
                continue
            
            dev_ch = p.device_channel_indices[ci['ch']]
            spk_data_name = [ch_name for ch_name in spikeband_files if str(dev_ch+1) == ch_name.split('_')[-1][2:-7]][0]
            data = trodesReader.readTrodesExtractedDataFile(os.path.join(rec_dir, spikeband_path, spk_data_name))
            mean_waveform = data['data']['voltage'][swi].squeeze().mean(1).astype(float)
            
            spike_description = {'clu':float(ci['cluster_id']),
                                'chn':float(p.device_channel_indices[ci['ch']+int(shank_ind)*128]),
                                'mean_waveform':mean_waveform,
                                'pos' : p.contact_positions[ci['ch']+int(shank_ind)*128],
                                'electrode' : float(shank_ind),
                                'annotations' : '',
                                'chn_meta' : ci}
            
            sampling_rate = int(timestamp['clockrate'])
            if not 'shank ' + shank_ind in InputData['Spike']:
                InputData['Spike']['shank ' + shank_ind] = {}
            if not 'chn ' + str(ci['ch']) in InputData['Spike']['shank ' + shank_ind]:
                InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])] = {}
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])] = {}
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk'] = templat_neo['spk'].copy()
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk']['times'] = timestamp['data']['time'][st]/sampling_rate * pq.s
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk']['t_stop'] = timestamp['data']['time'][-1]/sampling_rate * pq.s
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk']['t_start'] = timestamp['data']['time'][0]/sampling_rate * pq.s
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk']['sampling_rate'] = float(sampling_rate) * pq.Hz
            InputData['Spike']['shank ' + shank_ind]['chn ' + str(ci['ch'])]\
                ['clu ' + str(ci['cluster_id'])]['spk']['description'] = spike_description
        
        return InputData
    

def convert_TCR(data, timestamp, sorter_output_path, data_template, probegroup, bandpass_params):    
    tcr_path = [i for i in os.listdir(rec_file) if 'kilosort' in i][0]
    tcr_name_list = os.listdir(os.path.join(raw_dirname, rec_file, tcr_path))
    data_name = [i for i in tcr_name_list if 'group0.dat' in i][0]
    InputData = copy.deepcopy(Template)
    InputData['RecordingSystemEvent'] = 'null'
    InputData['Spike'] = 'null'
    InputData['SBP'] = 'null'
    InputData['LFP'] = 'null'

    InputData['name'] = 'TCR'
    InputData['TCR'] = {}

    # a = np.memmap(os.path.join(raw_dirname, rec_file, tcr_path,data_name), dtype=np.int16, mode='r+').reshape((-1,1024))

    spk_path = [i for i in os.listdir(rec_file) if 'spikeband' in i][0]
    spk_name_list = os.listdir(os.path.join(raw_dirname, rec_file,spk_path))

    for i in tqdm(spk_name_list):
        data = trodesReader.readTrodesExtractedDataFile(os.path.join(raw_dirname, rec_file,spk_path,i))
        # sbp_time = data
        des = {}
        for des_key in data:
            if 'data' == des_key:
                continue
            des[des_key] = data[des_key]
        
        if 'timestamps' in i:
            continue
            
        f_ch = data['data']['voltage']
        chn = float(i.split('_')[-1][2:-7])
        ch1_mad = scipy.stats.median_abs_deviation(f_ch)
        thred = np.median(f_ch)-(3.5/0.6745)*ch1_mad
        peaks, _ = find_peaks(-f_ch, height=-thred)
        ind=np.array([peaks-24+i for i in range(64)]).T
        ch_spike = f_ch[ind[1:-100]]
        
        spike_description = {'clu':0.0,
                            'chn':chn,
                            'mean_waveform':ch_spike.mean(0).astype(float),
                            'pos' : probe_2d.contact_positions[int(chn)-1],
                            'electrode' : float(1),
                            'annotations' : '',
                            'chn_meta' : des}
        
        sampling_rate = int(sbp_time['clockrate'])
        
        if not str(chn) in InputData['TCR']:
            InputData['TCR'][str(chn)] = {}
            
        InputData['TCR'][str(chn)]['0'] = {}
        InputData['TCR'][str(chn)]['0']['spk'] = templat_neo['spk'].copy()
        InputData['TCR'][str(chn)]['0']['spk']['times'] = sbp_time['data']['time'][peaks]/sampling_rate * pq.s
        InputData['TCR'][str(chn)]['0']['spk']['t_stop'] = sbp_time['data']['time'][-1]/sampling_rate * pq.s
        InputData['TCR'][str(chn)]['0']['spk']['t_start'] = sbp_time['data']['time'][0]/sampling_rate * pq.s
        InputData['TCR'][str(chn)]['0']['spk']['sampling_rate'] = float(sampling_rate) * pq.Hz
        InputData['TCR'][str(chn)]['0']['spk']['description'] = spike_description

    return InputData


def convert_LFP(root_dir, data_template):
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


def format_file(root_dir, map_path, output_dir, content_list, sorter):
    #%% prepare
    # load template
    FILEPATH = os.path.dirname(os.path.abspath(__file__))
    Template = yaml.safe_load(open(os.path.join(FILEPATH,'template_neural_data.yml')))
    
    # load probe
    probegroup = read_probeinterface('/AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/nezha-gai.json')
    
    # get data path
    global rec_dir, spikeband_path, spikeband_files, time_file, timestamp
    try:
        walk_file = [j for j in os.walk(root_dir)]
        for f_l in walk_file:
            rec_name = [f_n for f_n in f_l[2] if '.rec' in f_n]
            if len(rec_name) !=0:
                break
        rec_dir = f_l[0]
    except:
        sys.exit()
    
    spikeband_path = [i for i in os.listdir(rec_dir) if 'spikeband' in i][0]
    spikeband_files = os.listdir(os.path.join(rec_dir, spikeband_path))
    time_file = [t for t in spikeband_files if 'timestamps' in t][0]
    timestamp = trodesReader.readTrodesExtractedDataFile(os.path.join(rec_dir, spikeband_path, time_file))
    
    

    # get sorter_output path
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

    #%% convert TCR
    if len([i for i in content_list if i.upper()=='TCR'])>0:
        InputData = convert_TCR(data, timestamp, sorter_output_path, Template, probegroup, bandpass_params)    
        InputList.append(InputData)

    #%% convert LFP
    if len([i for i in content_list if i.upper()=='LFP'])>0:
        Template['LFP'] = templat_neo['ana']
        InputData = convert_LFP(os.path.join(root_dir, 'raw_data'), Template)
        InputList.append(InputData)

    #%% convert RecordingSystemEvent
    Template['RecordingSystemEvent'] = templat_neo['event']
    InputData = convert_RSE(os.path.join(root_dir, 'raw_data'), Template)
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