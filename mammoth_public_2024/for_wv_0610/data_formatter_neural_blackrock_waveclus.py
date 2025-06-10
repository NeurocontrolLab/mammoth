#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

#%%
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import json
from tqdm import tqdm
import scipy
import scipy.io as sio
from scipy.signal import butter, lfilter, find_peaks
import numpy as np
import pandas as pd
import brpylib
import argparse
from datetime import datetime
from zoneinfo import ZoneInfo 
import quantities as pq
from probeinterface import read_probeinterface

from pynwb import NWBHDF5IO, NWBFile, TimeSeries
from pynwb.ecephys import LFP, ElectricalSeries
# from pynwb.behavior import BehavioralEvents
from pynwb.file import Subject


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


def format_file(root_dir, map_path, output_dir, content_list):
    #%% prepare
    # load template
    # FILEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # Template = yaml.safe_load(open(os.path.join(FILEPATH, 'dependencies', 'template_neural_data.yml')))
    
     # %% Set basic info and create NWB file
    root_directory = os.path.join(root_dir, 'bhv')
    # file_pattern = [i for i in os.listdir(root_directory) if ('csv' in i) and ('meta' in i)][0]
    # assert os.path.exists(os.path.join(root_directory, file_pattern)), "Please define the metadata first!"
    # meta_dict_pd = pd.read_csv(os.path.join(root_directory, file_pattern)).to_dict()
    # meta_dict = {}
    # for i,j in zip(meta_dict_pd['Key'], meta_dict_pd['Value']):
    #     meta_dict[meta_dict_pd['Key'][i]] = meta_dict_pd['Value'][j]

    # metadata default
    meta_dict = {
        'Subject/subject_id': 'Caesar',
        'Subject/age': '',
        'Subject/species': 'Rehsus monkey',
        'Subject/sex': 'male',
        'Subject/description': 'None',
        'NWB/session_description': '',
        'NWB/session_start_time': '2020/10/09 09:00',  
        'NWB/session_number': '',
        'NWB/lab_name': 'cuilab',
        'NWB/experimenter': '',
        'NWB/institution': '',
        'NWB/related_publications': ''
    }
    
    subject_dict = {}
    nwb_dict = {}
    for i in meta_dict:
        key = i.split('/')[1]
        if 'Subject' in i:
            subject_dict[key] = meta_dict[i]
        else:
            nwb_dict[key] = meta_dict[i]

    if not 'description' in subject_dict:
        subject_dict['description'] = 'None'

    subject = Subject(
        **subject_dict
        )

    session_id = os.path.basename(root_dir)
    nwb_dict['identifier'] = session_id
    nwb_dict['session_id'] = nwb_dict['session_number']
    del nwb_dict['session_number']
    nwb_dict['lab'] = nwb_dict['lab_name']
    del nwb_dict['lab_name']
    nwb_dict['session_start_time'] = \
        datetime.strptime(nwb_dict['session_start_time'], 
                        '%Y/%m/%d %H:%M').replace(tzinfo=ZoneInfo('Asia/Shanghai'))
    nwb_dict['subject'] = subject
    nwb_dict['keywords'] = ["ecephys", "monkey", "motor control"]
    # del nwb_dict['experiment_name']
    if 'task' in nwb_dict.keys():
        del nwb_dict['task']

    for k, v in nwb_dict.items():
        if isinstance(v, float):
            if np.isnan(v):
                nwb_dict[k] = ''
    nwbfile = NWBFile( **nwb_dict)

    # load probe
    probegroup = read_probeinterface(map_path)

    probe_info_str = os.path.basename(map_path).split(".")[0]
    probe_info = probe_info_str.split("_")
    device_dict = {'name': '%s array' % probe_info[1], 'description': '%s' % probe_info[2],
                   'manufacturer': '%s Microsystem' % probe_info[4]} 

    # set device and electrode table
    device = nwbfile.create_device(name=device_dict['name'],
                                description=device_dict['description'],
                                manufacturer=device_dict['manufacturer'])

    # get electrode table
    # nshanks = len(probegroup.probes)
    shank_location = probe_info[3].split("-")
    ETR_list = []

    for ishank, p in enumerate(probegroup.probes):

        electrode_group = nwbfile.create_electrode_group(
            name="shank{}".format(ishank),
            description="electrode group for shank {}".format(ishank),
            device=device,
            location=shank_location[ishank]
        )
        ETR_list.append(electrode_group)

        # add electrodes to the electrode table
        for ielec in range(p.contact_positions.shape[0]):
            nwbfile.add_electrode(
                x=float(p.contact_positions[ielec, 0]),
                y=float(p.contact_positions[ielec, 1]),
                group=electrode_group,
                id = p.device_channel_indices[ielec],
                location=shank_location[ishank])

    # view the electrodes table in a pandas DataFrame
    # nwbfile.electrodes.to_dataframe()

    # define unit columns
    nwbfile.add_unit_column(name="chn_id", 
                            description="source channel of unit")
    nwbfile.add_unit_column(name="sorter", 
                            description="sorting method")
    nwbfile.add_unit_column(name="sorting_info", 
                            description="metadata of sorting results")
    nwbfile.add_unit_column(name="time_unit", 
                            description="time scale")
    
    ecephys_module = nwbfile.create_processing_module(
            name="ecephys", 
            description="processed extracellular electrophysiology data"
        )
    
    #%% convert RecordingSystemEvent
    walk_file = [j for j in os.walk(root_dir)]

    for f_l in walk_file:
        rec_name = [f_n for f_n in f_l[2] if ('.nev' in f_n) or ('NSP' in f_n)]
        if len(rec_name) !=0:
            break
    datafile = os.path.join(f_l[0], rec_name[0])

    nev_file = brpylib.NevFile(str(datafile))

    # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
    nsp_data = nev_file.getdata()

    # Close the nsx file now that all data is out
    nev_file.close()

    ptp_t = np.array(nsp_data['digital_events']['TimeStamps'])/1e9
    global sys_start_time
    sys_start_time = ptp_t[0]
    
    # add Recording System Event
    event_label = list(np.array(nsp_data['digital_events']['UnparsedData']).astype('float'))
    event_label = np.array(event_label)
    if event_label[0]>60000:
        event_label = event_label-65280

    event_times = ptp_t*pq.s
    event_times = event_times.rescale(pq.s).magnitude
    
    event_times = event_times[event_label!=0]
    event_label = event_label[event_label!=0]

    event_marker_series = TimeSeries(
        name='RecordingSystemEvents',
        data=event_label,
        unit='NA', 
        timestamps=event_times.flatten(),
        description='Event markers recorded by recording system'
    )

    ecephys_module.add(event_marker_series)
    
    #%% load .ns6 data
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
    sorter_output_path = os.path.join(root_dir, 'sorted_data', 'WaveClus', 'selected')

    # set bandpass parameters
    fs = 30000.0
    lowcut = 300.0
    highcut = 6000.0
    bandpass_params = (fs, lowcut, highcut)
    
    # # initialize data list
    # InputList = []
    # Template = {}

    #%% convert Spike
    if len([i for i in content_list if i.lower()=='spike'])>0:
       
        # def convert_spike(data, timestamp, sorter_output_path, data_template, probegroup, bandpass_params):
        InputData_Spike = {}

        sorter_output_path = os.path.join(root_dir, 'sorted_data', 'WaveClus', 'selected')
        sorting_files = [i for i in os.listdir(sorter_output_path) if 'times' in i]

        # fs, lowcut, highcut = bandpass_params
        # fs = 30000.0
        # lowcut = 300.0
        # highcut = 6000.0


        # assemble data
        cluster_id = 0     
        for i in tqdm(sorting_files):
            shank_ind = int(i.split('_elec')[-1].split('-')[0]) - 1

            chn = i.split('_ele_')[0].split('ch_')[-1]
            elec_chn = i.split('_elec')[-1].split('-')[-1].split('.mat')[0]
           
            # data_path = os.path.join(sorter_output_path, i, 'sorter_output')
            # cluster_info = pd.read_csv(os.path.join(data_path, 'cluster_info.tsv'), sep="\t")
            
            sorting_data = sio.loadmat(os.path.join(sorter_output_path, i))
            cluster_data = sorting_data['cluster_class']
            cluster_data = cluster_data[cluster_data[:, 0]>0, :] # exclude noise

            spike_times = cluster_data[:, 1].squeeze()
            spike_clusters = cluster_data[:, 0].squeeze()
            
            p = probegroup.probes[shank_ind]
            
            if len(np.unique(spike_clusters))<1:
                continue

            for clu in np.unique(spike_clusters):
                # ind = np.where(cluster_info.cluster_id==clu)[0][0]
                # ci = cluster_info.iloc[ind].to_dict()
                st_in_sec = spike_times[spike_clusters==clu] # spiketrains in ms

                sampling_rate = 30000
                st = (st_in_sec/1e3*sampling_rate).astype(int) # frame

                # choose 24 sampling points before and 40 after, sampling rate = 30k Hz
                # => [-0.5ms, +1 ms], referring to WaveClus
                swi = [st-24+i for i in range(64)]  
                bch = int(chn)-1
                # bch = int(elec_chn)
                # assert len(np.where(p.device_channel_indices==int(chn)))==1
                
                # mean_waveform = whitened_data[:,ci['ch']][swi].squeeze().mean(1).astype(float)
                
                f_ch = butter_bandpass_filter(data[:, bch], lowcut, highcut, fs, order=5)
                try:
                    mean_waveform = f_ch[swi].squeeze().mean(1).astype(float)
                except:
                    swi = [i[:-1] for i in swi]
                    # maxind = min([np.argmax(i>=len(data)) for i in swi])
                    # maxind = maxind if maxind>0 else -1
                    # swi = [i[:maxind] for i in swi]
                    mean_waveform = f_ch[swi].squeeze().mean(1).astype(float)

                # mean_waveform = data[swi,bch].squeeze().mean(1).astype(float)
                
                spike_description = {'clu': float(cluster_id),
                                    'chn': float(chn),
                                    'mean_waveform': mean_waveform,
                                    'pos': p.contact_positions[p.device_channel_indices==int(chn), :],
                                    'electrode': float(shank_ind),
                                    'annotations': '',
                                    'chn_meta': 'elec_chn %s' % elec_chn}
                
                # maxind = np.argmax(st>=len(timestamp))
                # if maxind>0:
                #     ptp_t = timestamp[st[:maxind]]/1e9
                # else:
                #     ptp_t = timestamp[st]/1e9
                ptp_t = timestamp[st]/1e9 # sec

                InputData_Spike['clu ' + str(cluster_id)] = {}
                InputData_Spike['clu ' + str(cluster_id)]['shankid'] = shank_ind
                InputData_Spike['clu '  + str(cluster_id)]['chn'] = str(chn)

                InputData_Spike['clu ' + str(cluster_id)]['times'] = ptp_t * pq.s
                InputData_Spike['clu ' + str(cluster_id)]['t_stop'] = timestamp[-1]/1e9 * pq.s
                InputData_Spike['clu ' + str(cluster_id)]['t_start'] = timestamp[0]/1e9 * pq.s
                InputData_Spike['clu ' + str(cluster_id)]['sampling_rate'] = float(sampling_rate) * pq.Hz
                InputData_Spike['clu ' + str(cluster_id)]['description'] = spike_description
                
                cluster_id += 1
                    
        nunits = len(InputData_Spike)

        # b = kilo.spiketrains[0].description["chn_meta"]
        for iunit, unit in enumerate(InputData_Spike.keys()):
            spiketimes = InputData_Spike[unit]
            nwbfile.add_unit(spike_times=np.array(spiketimes['times'].rescale(pq.s).magnitude), #should be array rather than quantities 
                            id=iunit,
                            chn_id=int(spiketimes['description']["chn"]),
                            waveform_mean=spiketimes['description']["mean_waveform"],
                            electrode_group=ETR_list[int(spiketimes['description']["electrode"])],
                            time_unit='seconds',
                            sorting_info=json.dumps(spiketimes['description']["chn_meta"]),
                            # sorting_info=spiketimes['description']["chn_meta"],
                            sorter = 'WaveClus')


    #%% convert TCR
    if len([i for i in content_list if i.upper()=='TCR'])>0:
    
        InputData_TCR = {}

        for shank_ind, p in enumerate(probegroup.probes):
            for ch_ind, chn in tqdm(enumerate(p.device_channel_indices)):
                
                f_ch = butter_bandpass_filter(data[:, chn], lowcut, highcut, fs, order=5)
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
                
                InputData_TCR[str(chn)] = {}
                InputData_TCR[str(chn)]['times'] = ptp_t * pq.s
                InputData_TCR[str(chn)]['t_stop'] = timestamp[-1]/1e9 * pq.s
                InputData_TCR[str(chn)]['t_start'] = timestamp[0]/1e9 * pq.s
                InputData_TCR[str(chn)]['sampling_rate'] = float(sampling_rate) * pq.Hz
                InputData_TCR[str(chn)]['description'] = spike_description

        if len([i for i in content_list if i.lower()=='spike'])>0:
            nunits = len(InputData_Spike)
        else:
            nunits = 0

        # add TCR into nwbfile.units
        for iunit, unit in enumerate(InputData_TCR.keys()):
            spiketimes = InputData_TCR[unit]
            nwbfile.add_unit(spike_times=np.array(spiketimes['times'].rescale(pq.s).magnitude), #should be array rather than quantities 
                            id=int(iunit+nunits),
                            chn_id=int(spiketimes['description']["chn"]),
                            waveform_mean=spiketimes['description']["mean_waveform"],
                            electrode_group=ETR_list[int(spiketimes['description']["electrode"])],
                            time_unit='seconds',
                            sorting_info = 'NA',
                            sorter = 'TCR')
            
        # view the unit table in a pandas DataFrame
        # nwbfile.units.to_dataframe()
    

    #%% convert LFP
    if len([i for i in content_list if i.upper()=='LFP'])>0:
       
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

            lfp_data = np.concatenate(lfp_data['data'], 1).T
            # lfp_timestamp = np.concatenate([i['Timestamp'] for i in lfp_data['data_headers']])/1e9
            lfp_timestamp = lfp_data['data_headers'][0]['Timestamp'][0]/1e9
        
        else:
            for f_l in walk_file:
                rec_name = [f_n for f_n in f_l[2] if '.ns6' in f_n]
                if len(rec_name) !=0:
                    break

            datafile = os.path.join(f_l[0], rec_name[0])

            raw_file = brpylib.NsxFile(str(datafile))

            raw_lfp_data = raw_file.getdata(downsample=15)
            # Note: .ns3 file only from 0.3-250 Hz
            # ButterWorthFilter lowpass and highpass

            raw_file.close()

            lfp_data = np.concatenate(raw_lfp_data['data'], 1).T
            # lfp_timestamp = np.concatenate([i['Timestamp'] for i in raw_lfp_data['data_headers']])/1e9
            if not isinstance(raw_lfp_data['data_headers'][0]['Timestamp'], list):
                lfp_timestamp = raw_lfp_data['data_headers'][0]['Timestamp']/1e9
            else:
                lfp_timestamp = raw_lfp_data['data_headers'][0]['Timestamp'][0]/1e9
        
        lfp_timestamp = [lfp_timestamp]

        # add LFP data
        electrode_index = nwbfile.electrodes.id.data   
        lfp_array = lfp_data[:, electrode_index]

        all_table_region = nwbfile.create_electrode_table_region(
            region=list(range(len(electrode_index))), 
            description="all electrodes",
        )

        lfp_electrical_series = ElectricalSeries(
            name="LFP",
            description="LFP data",
            data=lfp_array,
            electrodes=all_table_region,
            timestamps=lfp_timestamp
        )

        lfp = LFP(electrical_series=lfp_electrical_series)

        ecephys_module.add(lfp)
    

    #%% Writing standard NWB file
    if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    if len(content_list) == 1:
        save_path = os.path.join(output_dir, "neural_data_no_sort.nwb")
        if os.path.exists(save_path):
            os.remove(save_path)

        with NWBHDF5IO(os.path.join(output_dir, "neural_data_no_sort.nwb"), "w") as io:
            io.write(nwbfile)
            
    elif 'spike' in content_list:
        save_path = os.path.join(output_dir, "neural_data.nwb")
        if os.path.exists(save_path):
            os.remove(save_path)

        with NWBHDF5IO(os.path.join(output_dir, "neural_data.nwb"), "w") as io:
            io.write(nwbfile)


#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/double_reach/Caesar/data_recording/20201009_DoubleReach_001_TestUDP', 
                    metavar='/the/path/your/data/located/in', help='input folder')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/double_reach/Caesar/data_recording/20201009_DoubleReach_001_TestUDP/formatted_data', 
                    metavar='/the/path/you/want/to/save', help='output folder')

parser.add_argument('-mp', '--map_path', 
                    default='/AMAX/cuihe_lab/share_rw/CuiLab-Database/interception/Caesar/Caesar_Utah_128x2_PMd-M1_BlackRock.json')

parser.add_argument('-flag', '--sort_flag', 
                    default='0')


args = parser.parse_args()

if args.sort_flag == '1':
    content_list = ['spike', 'TCR', 'LFP']
elif args.sort_flag == '0':
    content_list = ['TCR']

format_file(args.root, args.map_path, args.output, content_list)        
        