#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:28:49 2024

@author: cuilab
"""

import os
import json
from tqdm import tqdm
import numpy as np
import shutil
import pandas as pd
import seaborn as sns
import brpylib
import argparse
import quantities as pq
import matplotlib.pyplot as plt
from SmartNeo.interface_layer.nwb_interface import NWBInterface
from scipy.signal import butter, lfilter


def run(data_dir, output_dir, root_dir):
    """
    data_dir: directory for nwb file
    output_dir: directory to save results
    root_dir: directory to search for ns6 file
    """
    #%% read data

    # read formated SmartNeo data
    nwb_saver = NWBInterface()
    neural_data = nwb_saver.read_nwb(filename = os.path.join(data_dir,'neural_data.nwb'))

    # get TCR and sorted spike data
    TCR = [i for i in neural_data.segments if 'TCR' in i.name][0]
    kilo = [i for i in neural_data.segments if 'kilosort2.5' in i.name][0]
    
    #%% check consistency
    # compare time (pairwise)
    data_list = []
    for i in TCR.spiketrains:
        for j in kilo.spiketrains:
            data = {}
            if i.description['chn']==(j.description['chn']):
                tcr_time = (i.times.rescale(pq.s) -i.t_start.rescale(pq.s)).magnitude*1000
                kilo_time = (j.times.rescale(pq.s) -j.t_start.rescale(pq.s)).magnitude*1000
                rescale_time = lambda x: np.where(np.histogram(x,range=[0,np.floor(x[-1])],bins = int(np.floor(x[-1])/1))[0]!=0)[0]
                chn_per = len(set(rescale_time(tcr_time)) & set(rescale_time(kilo_time)))/len(j)
                data['chn'] = 'unshuffled'
                data['chn_per'] = chn_per
                data_list.append(data)

    # compare time (shuffled)
    ch_num = len(TCR.spiketrains)
    for j,ind in zip(kilo.spiketrains, np.random.choice(range(ch_num),ch_num)):
        i  = TCR.spiketrains[ind]
        data = {}
        tcr_time = (i.times.rescale(pq.s) -i.t_start.rescale(pq.s)).magnitude*1000
        kilo_time = (j.times.rescale(pq.s) -j.t_start.rescale(pq.s)).magnitude*1000
        rescale_time = lambda x: np.where(np.histogram(x,range=[0,np.floor(x[-1])],bins = int(np.floor(x[-1])/1))[0]!=0)[0]  
        chn_per = len(set(rescale_time(tcr_time)) & set(rescale_time(kilo_time)))/len(j)
        data['chn'] = 'shuffled'
        data['chn_per'] = chn_per
        data_list.append(data)
    
    # save results
    df = pd.DataFrame(data_list)
    sns.barplot(x="chn",y="chn_per",data=df)
    save_path = f'{output_dir}/chn_consis_summary.png'
    plt.savefig(save_path)

    with open(os.path.join(output_dir, 'qc_summary.json'), "r") as file:
        summary = json.load(file)
        summary['channel consistency'] = {"unshuffled": list(df[df['chn']=='unshuffled']['chn_per'].values),
                                          "shuffled": list(df[df['chn']=='shuffled']['chn_per'].values)}
    
    with open(os.path.join(output_dir, 'qc_summary.json'), 'w') as file:
        json.dump(summary, file)


    #%% single-channel waveform

    # load raw data
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
  
    # bandpass
    def butter_bandpass(lowcut, highcut, fs, order=5):
        return butter(order, [lowcut, highcut], fs=fs, btype='band')

    def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        y = lfilter(b, a, data)
        return y
    
    fs = 30000.0
    lowcut = 300.0
    highcut = 6000.0

    # get new folder
    folder_name = os.path.join(output_dir, 'chn_waveform')
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)

    os.mkdir(folder_name)
    plt.ioff()

    # save results
    kilo = [i for i in neural_data.segments if 'kilosort2.5' in i.name][0]
    for j in tqdm(kilo.spiketrains):
        #data = {}
        dev_ch = int(j.description['chn'])
        times = ((j.times-j.t_start).rescale(pq.s).magnitude*30000).astype(int)
        ind=np.array([times-24+i for i in range(64)]).T
        f_ch = butter_bandpass_filter(data[:,dev_ch], lowcut, highcut, fs, order=5)
        ch_spike = f_ch[ind[0:-10]]
        
        for i in range(100):
            plt.plot(ch_spike[i])
        plt.title(f'chn_{dev_ch}') 
        
        save_path = f'{folder_name}/chn_{dev_ch}.png'
        plt.savefig(save_path)
        plt.close()

    plt.ion()


parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Data_recording/20241008_Interception_001', 
                    metavar='/the/root/path/your/data/located/in', help='root folder')
parser.add_argument("-d", "--data", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Data_recording/20241008_Interception_001/formatted_data', 
                    metavar='/the/path/your/data/located/in', help='data folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Data_recording/20241008_Interception_001/description', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()

run(args.data, args.output, args.root)

