
import argparse
import os
import numpy as np
from dtw import accelerated_dtw
from SmartNeo.interface_layer.nwb_interface import NWBInterface
import matplotlib.pyplot as plt
import quantities as pq
import pandas as pd
import json

def string_to_float_list(string):

    string = string.strip('[]').strip()
    str_list = string.split()
    float_list = [float(item) for item in str_list]
    
    return float_list

def scientific_string_to_float_list(string):
    string = string.strip('[]').replace('\n', '').strip().rstrip('s').replace(']', '')
    str_list = string.split()
    float_list = [float(item) for item in str_list]
    return float_list

parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-f", "--file", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20240925_interception_004', 
                    metavar='/the/path/your/data/located/in', help='input folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20240925_interception_004', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()
raw_dirname = args.file

description_data_path = os.path.join(raw_dirname,'description')
nwb_saver = NWBInterface()
bhv_data = nwb_saver.read_nwb(filename = os.path.join(raw_dirname,'formatted_data', 'continuous_behavior.nwb'))

nwb_saver = NWBInterface()
neural_data = nwb_saver.read_nwb(filename = os.path.join(raw_dirname,'formatted_data', 'neural_data.nwb'))

bhv_event = [i for i in bhv_data.segments if i.name=='Event marker'][0].events[0]
marker = [json.loads(i) for i in bhv_event.labels]
num_marker = []
event_marker = []
marker_time = []

for i,j in zip(bhv_event.labels,bhv_event.times):
    event_recorded = json.loads(i)
    if isinstance(event_recorded,dict):
        num_marker.append(event_recorded['Marker'])
        event_marker.append(event_recorded['Event'])
        marker_time.append(j.rescale(pq.s).magnitude)


neural_event = [i for i in neural_data.segments if i.name=='RecordingSystemEvent'][0]
neural_event_labels = neural_event.events[0].labels
neural_event_times = neural_event.events[0].times

event = [i for i in neural_data.segments if i.name=='RecordingSystemEvent'][0].events[0]
if event.labels[0]>60000:
    event_label = np.array(event.labels)-65280
event_times = event.times.rescale(pq.s).magnitude
neural_event_times = event_times[event_label!=0] * pq.s
neural_event_labels = event_label[event_label!=0]

ml_event_labels = np.array(num_marker)
ml_event_times = np.array(marker_time) * pq.s



li_distance = lambda x, y: 0 if np.abs(x - y)==0 else len(neural_event_labels)
    # manhattan_distance = lambda x, y: np.abs(x - y)

_,_,_,path1 = accelerated_dtw(ml_event_labels, neural_event_labels, dist=li_distance)
path1 = np.array(path1).T

zero_dir = (ml_event_labels[path1[:,0]]-neural_event_labels[path1[:,1]])==0
# assert len(neural_event_labels)-sum(zero_dir)<(0.005*len(neural_event_labels)), 'wrong alignment'
# a = (len(neural_event_labels)-sum(zero_dir))
diff_time = ml_event_times[path1[:,0]]-neural_event_times[path1[:,1]]
def remove_outliers(data, factor=1.5):
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    lower_bound = q1 - (factor * iqr)
    upper_bound = q3 + (factor * iqr)
    return data[(data > lower_bound) & (data < upper_bound)]
diff_time = remove_outliers(diff_time, factor=1.5)
plt.boxplot(diff_time-np.mean(diff_time))
len(diff_time)/len(neural_event_labels)

#%% save to appointed path
if not os.path.exists(os.path.join(args.output,'description')):
    os.mkdir(os.path.join(args.output,'description'))
plt.savefig(os.path.join(args.output, 'description', 'Time_consistency_check.png'), dpi=300, bbox_inches = 'tight')

diff_time_mean = np.mean(diff_time)

with open(os.path.join(args.output, 'description', 'diff_time_mean.txt'), 'w') as file:
    file.write(f'{diff_time_mean}\n')