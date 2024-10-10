import argparse
import os
import json
import copy
import pandas as pd
import numpy as np
import quantities as pq
import joblib
import re
import shutil
import seaborn as sns
import scipy
from dtw import accelerated_dtw
from SmartNeo.analysis_layer.tools.dataset_maker.dataset_maker import DatasetMaker
from SmartNeo.user_layer.dict_to_neo import templat_neo
from SmartNeo.interface_layer.nwb_interface import NWBInterface
from SmartNeo.analysis_layer.spike_preprocessor3 import SpikeStatistics
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from elephant.kernels import GaussianKernel
from probeinterface import read_probeinterface,Probe, ProbeGroup
from probeinterface.plotting import plot_probe
from tqdm import tqdm

#%% read data
# data folder
#%% read data
# data folder
unit_dict = {'s': pq.s}
parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-f", "--file", type=str,
                    default='/media/lenovo/Extreme Pro/Bohr data/20240201_br_utah_interception_120_semi_brain_control_4_001', 
                    metavar='/the/path/your/data/located/in', help='input folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Abel/Data_recording/20240705_centerOut_001', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()
raw_dirname = args.file

formated_data_path = os.path.join(raw_dirname, 'formatted_data')

# read formated SmartNeo data
nwb_saver = NWBInterface()
neural_data = nwb_saver.read_nwb(filename = os.path.join(formated_data_path,'neural_data.nwb'))

TCR = [i for i in neural_data.segments if i.name=='TCR'][0]
kilo = [i for i in neural_data.segments if i.name=='kilosort2.5'][0]
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
description_path = os.path.join(raw_dirname, 'description')
df = pd.DataFrame(data_list)
sns.barplot(x="chn",y="chn_per",data=df)
save_path = f'{description_path}/chn_consis_summary.png'
plt.savefig(save_path)

import brpylib

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

import matplotlib.pyplot as plt
import shutil

folder_name = os.path.join(description_path, 'chn_waveform')
kilo = [i for i in neural_data.segments if i.name=='kilosort2.5'][0]

if os.path.exists(folder_name):
    shutil.rmtree(folder_name)

os.mkdir(folder_name)
plt.ioff()

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