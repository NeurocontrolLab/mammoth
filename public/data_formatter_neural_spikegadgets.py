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
def convert_spike(raw_dir, sorter_output_path, data_template, probegroup):
    InputData = copy.deepcopy(data_template)
    InputData['RecordingSystemEvent'] = 'null'
    InputData['TCR'] = 'null'
    InputData['SBP'] = 'null'
    InputData['LFP'] = 'null'

    InputData['name'] = sorter_output_path.split(os.sep)[-1].lower().replace('_', '.') 
    InputData['Spike'] = {}

    spikeband_path = [i for i in os.listdir(raw_dir) if 'spikeband' in i][0]
    spikeband_files = os.listdir(os.path.join(raw_dir, spikeband_path))
    time_file = [t for t in spikeband_files if 'timestamps' in t][0]
    timestamp = trodesReader.readTrodesExtractedDataFile(os.path.join(raw_dir, spikeband_path, time_file))

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
            data = trodesReader.readTrodesExtractedDataFile(os.path.join(raw_dir, spikeband_path, spk_data_name))
            mean_waveform = data['data']['voltage'][swi].squeeze().mean(1).astype(float)
            
            spike_description = {'clu':float(ci['cluster_id']),
                                 'chn':float(p.device_channel_indices[ci['ch']]),
                                 'mean_waveform':mean_waveform,
                                 'pos' : p.contact_positions[ci['ch']],
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
    

def convert_TCR(raw_dir, data_template, probegroup):    
    InputData = copy.deepcopy(data_template)
    InputData['RecordingSystemEvent'] = 'null'
    InputData['Spike'] = 'null'
    InputData['SBP'] = 'null'
    InputData['LFP'] = 'null'

    InputData['name'] = 'TCR'
    InputData['TCR'] = {}

    spikeband_path = [i for i in os.listdir(raw_dir) if 'spikeband' in i][0]
    spikeband_files = os.listdir(os.path.join(raw_dir, spikeband_path))
    time_file = [t for t in spikeband_files if 'timestamps' in t][0]
    timestamp = trodesReader.readTrodesExtractedDataFile(os.path.join(raw_dir, spikeband_path, time_file))

    for shank_ind, p in enumerate(probegroup.probes):
        for dev_ch in p.device_channel_indices:
            spk_data_name = [ch_name for ch_name in spikeband_files if str(dev_ch+1) == ch_name.split('_')[-1][2:-7]][0]
            data = trodesReader.readTrodesExtractedDataFile(os.path.join(raw_dir, spikeband_path, spk_data_name))
            # sbp_time = data
            des = {}
            for des_key in data:
                if 'data' == des_key:
                    continue
                des[des_key] = data[des_key]
                
            f_ch = data['data']['voltage']
            chn = dev_ch+1
            ch1_mad = scipy.stats.median_abs_deviation(f_ch)
            thred = np.median(f_ch)-(3.5/0.6745)*ch1_mad
            peaks, _ = find_peaks(-f_ch, height=-thred)
            ind=np.array([peaks-24+i for i in range(64)]).T
            ch_spike = f_ch[ind[1:-100]]
            
            spike_description = {'clu':0.0,
                                 'chn':dev_ch,
                                 'mean_waveform':ch_spike.mean(0).astype(float),
                                 'pos' : p.contact_positions[int(chn)-1],
                                 'electrode' : float(1),
                                 'annotations' : '',
                                 'chn_meta' : des}
            
            sampling_rate = int(timestamp['clockrate'])
            
            if not str(chn) in InputData['TCR']:
                InputData['TCR'][str(chn)] = {}
                
            InputData['TCR'][str(chn)]['0'] = {}
            InputData['TCR'][str(chn)]['0']['spk'] = templat_neo['spk'].copy()
            InputData['TCR'][str(chn)]['0']['spk']['times'] = timestamp['data']['time'][peaks]/sampling_rate * pq.s
            InputData['TCR'][str(chn)]['0']['spk']['t_stop'] = timestamp['data']['time'][-1]/sampling_rate * pq.s
            InputData['TCR'][str(chn)]['0']['spk']['t_start'] = timestamp['data']['time'][0]/sampling_rate * pq.s
            InputData['TCR'][str(chn)]['0']['spk']['sampling_rate'] = float(sampling_rate) * pq.Hz
            InputData['TCR'][str(chn)]['0']['spk']['description'] = spike_description

    return InputData


def convert_LFP(raw_dir, data_template):
    # LFPRecordingParam = blackrock_data.extended_headers
    InputData = copy.deepcopy(data_template)
    InputData['RecordingSystemEvent'] = 'null'
    InputData['Spike'] = 'null'
    InputData['SBP'] = 'null'
    InputData['TCR'] = 'null'

    lfp_path = [i for i in os.listdir(raw_dir) if 'LFP' in i][0]
    lfp_name_list = os.listdir(os.path.join(raw_dir, lfp_path))
    
    des_list = []
    data_list = []
    for i in tqdm(lfp_name_list):
        data = trodesReader.readTrodesExtractedDataFile(os.path.join(raw_dir, lfp_path, i))
        des = {}
        for des_key in data:
            if 'data' == des_key:
                continue
            des[des_key] = data[des_key]
        
        if 'timestamps' in i:
            data_time = data
            continue
        des_list.append(des)
        key = data['data'].dtype.names[0]
        data_list.append(data['data'][key])
        assert len(data['data'].dtype.names)==1, 'too many key'

    cutoff = min([len(i) for i in data_list])
    for i,j in zip(des_list, data_list):
        i['cutoff'] = len(j)-cutoff
    data_list = [i[0:cutoff] for i in data_list]

    InputData['LFP'] = templat_neo['irr'].copy()
    InputData['LFP']['signal'] = np.array(data_list).T*pq.uV
    InputData['LFP']['times'] = data_time['data']['time'][0:cutoff]/int(data_time['clockrate'])*pq.s
    InputData['LFP']['t_start'] = 0*pq.s
    InputData['LFP']['description'] = des_list
    InputData['name'] = 'LFP'

    return InputData


def convert_RSE(raw_dir, data_template):
    event_path = [i for i in os.listdir(raw_dir) if 'DIO' in i][0]
    event_name_list = os.listdir(os.path.join(raw_dir, event_path))
    
    des_list = []
    data_list = []
    time_list = []
    for i in sorted(event_name_list):
        data = trodesReader.readTrodesExtractedDataFile(os.path.join(raw_dir, event_path, i))
        des = {}
        for des_key in data:
            if 'data' == des_key:
                continue
            des['des_key'] = data[des_key]
        des_list.append(des)
        
        if 'timestamps' in i:
            continue

        data_list.append(data['data']['state'])
        time_list.append(data['data']['time'])
        assert len(data['data'].dtype.names)==2, 'too many key'
        
    InputData = copy.deepcopy(data_template)   
    InputData['LFP'] = 'null'
    InputData['Spike'] = 'null'
    InputData['TCR'] = 'null'
    InputData['SBP'] = 'null'

    timeq=np.unique(np.concatenate(time_list))
    marker=np.zeros((len(timeq),6))
    marker[0]=[i[0] for i in data_list]

    for i in range(1,len(timeq)):
        marker[i,:]=marker[i-1,:] 
        for j in range(6):
            if timeq[i] in time_list[j]:
                marker[i,j]=data_list[j][time_list[j]==timeq[i]]   

    marker = marker[1::]
    # markertime=timeq-timeq[0]
    markertime=timeq
    markertime = markertime[1::]
    marker=np.fliplr(marker)

    bin2dec = lambda x: x.dot(2**np.arange(x.shape[1])[::-1])
    decmarker=bin2dec(marker)
    sampling_rate = int(data['clockrate'])

    InputData['RecordingSystemEvent']['times'] = markertime/sampling_rate*pq.sec
    InputData['RecordingSystemEvent']['labels'] = list(decmarker.astype('float'))
    InputData['RecordingSystemEvent']['description'] = des_list
    InputData['name'] = 'RecordingSystemEvent'

    return InputData


def format_file(root_dir, map_path, output_dir, content_list, sorter=None):
    #%% prepare
    # load template
    FILEPATH = os.path.dirname(os.path.abspath(__file__))
    Template = yaml.safe_load(open(os.path.join(FILEPATH,'template_neural_data.yml')))
    
    # load probe
    probegroup = read_probeinterface(map_path)
    
    # get data path
    try:
        walk_file = [j for j in os.walk(root_dir)]
        for f_l in walk_file:
            rec_name = [f_n for f_n in f_l[2] if '.rec' in f_n]
            if len(rec_name) !=0:
                break
        rec_dir = f_l[0]
    except:
        sys.exit()
    
    raw_dir = os.path.join(root_dir, rec_dir)
     
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

    # initialize data list
    InputList = []

    #%% convert Spike
    if len([i for i in content_list if i.lower()=='spike'])>0:
        Template['Spike'] = {'SorterName' : {}, 'kargs' : {}}
        InputData = convert_spike(raw_dir, sorter_output_path, Template, probegroup)
        InputList.append(InputData)

    #%% convert TCR
    if len([i for i in content_list if i.upper()=='TCR'])>0:
        InputData = convert_TCR(raw_dir, sorter_output_path, Template, probegroup)    
        InputList.append(InputData)

    #%% convert LFP
    if len([i for i in content_list if i.upper()=='LFP'])>0:
        Template['LFP'] = templat_neo['ana']
        InputData = convert_LFP(raw_dir, Template)
        InputList.append(InputData)

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
    nwb_saver.save_nwb(blockdata = neuralblock, 
                       filename = os.path.join(output_dir, '{}.nwb'.format('_'.join(content_list))))
        

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

parser.add_argument('-cl', '--content_list', 
                    default="['spike', 'TCR', 'LFP']")

parser.add_argument('-sorter', '--sorter_name', 
                    default='kilosort2_5')


args = parser.parse_args()

format_file(args.root, args.map_path, args.output, args.content_list, args.sorter_name)        