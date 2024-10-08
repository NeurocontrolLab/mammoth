from pynwb import NWBHDF5IO
import os
import argparse
import pynapple as nap
import json
import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter1d
from SmartNeo.analysis_layer.spike_preprocessor3 import SpikeStatistics
import quantities as pq
from elephant.kernels import GaussianKernel
import neo

def tsd_json_to_tsdframe(tsd):

    timestamps = tsd.index
    json_data = tsd.values

    parsed_data = [json.loads(entry) for entry in json_data]

    df = pd.json_normalize(parsed_data, sep='_')
    df.index = timestamps

    return df

#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument("-f", "--file", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20240925_interception_004', 
                    metavar='/the/path/your/data/located/in', help='input folder')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/20240925_interception_004', 
                    metavar='/the/path/you/want/to/save', help='output folder')

parser.add_argument('-mp', '--map_path', 
                    default='/home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Spike_sorting/SN+11386-000049.cmp')

#%% read nwb data
args = parser.parse_args()
raw_dirname = args.file
formated_data_path = os.path.join(raw_dirname,'formatted_data')
data_path = os.path.join(formated_data_path, "standard_data.nwb")

io = NWBHDF5IO(data_path, mode='r')
nwbfile = io.read()
behavior = nwbfile.processing['behavior']
ecephys = nwbfile.processing['ecephys']
behavior_ecephys_analysis = nwbfile.processing['behavior_ecephys_analysis']
time_series_data = nap.load_file(data_path)
time_difference = behavior_ecephys_analysis.data_interfaces['TimeDifference']['time_difference'][0] # behavior - ecephys

#%% convert nwb data to pynapple/pandas
units = time_series_data['units']
markers = time_series_data['BehaviorMarkers']
labels = time_series_data['BehaviorLabels']
vicon_motion = time_series_data['ViconMotion']
frame_df = tsd_json_to_tsdframe(time_series_data['FrameInfoTable'])
trial_info_df = tsd_json_to_tsdframe(time_series_data['TrialInfoTable'])

#%% multi-time stream alignment
markers_shift = nap.Tsd(t=markers.index - time_difference, d=markers.values)
labels_shift = nap.Tsd(t=labels.index - time_difference, d=labels.values)
f = labels_shift.d
marker_time = markers_shift.index[markers_shift.values == 24][0] # first vicon recorded frame is 24
vicon_start_time = vicon_motion.index[0]
marker_vicon_difference = marker_time - vicon_start_time
vicon_motion_shift = nap.TsdFrame(
    t=vicon_motion.index + marker_vicon_difference,
    d=vicon_motion.values
)
frame_df.index = frame_df.index - time_difference
trial_info_df.index = trial_info_df.index - time_difference

#%% slice data for analysis

def extract_trial_info(df):
    trial_data = []

    # 遍历每个 unique 的 trial_number
    for trial in df['trial_number'].unique():
        # 获取当前 trial 的所有数据
        trial_df = df[df['trial_number'] == trial]

        # 找到 start 和 end 对应的时间
        start_time = trial_df[trial_df['status'] == 'start'].index.min()
        end_time = trial_df[trial_df['status'] == 'end'].index.max()

        # 获取结束时的其他信息
        end_info = trial_df.loc[end_time].drop(['status', 'trial_number'])

        # 将结果存储
        trial_data.append({
            'trial_number': trial,
            'start_time': start_time,
            'end_time': end_time,
            **end_info.to_dict()
        })

    # 转换为 DataFrame，并设置 trial_number 作为索引
    result_df = pd.DataFrame(trial_data).set_index('trial_number')

    return result_df

def add_slices_to_trial_info(result_df, labels, marker):

    labels_slices = []
    marker_slices = []

    for trial in result_df.index:
        start_time = result_df.at[trial, 'start_time']
        end_time = result_df.at[trial, 'end_time']

        labels_slice = labels.get(start_time, end_time)
        marker_slice = marker.get(start_time, end_time)

        labels_slices.append(labels_slice)
        marker_slices.append(marker_slice)

    result_df['labels_slice'] = labels_slices
    result_df['markers_slice'] = marker_slices

    return result_df

def extract_MO_data(df):
    result_list = []

    for idx, row in df[df['wrong_number'] == 0].iterrows():
        label_slice = row['labels_slice']
        time_marker = float(label_slice.t[label_slice.d==b'Movement onset'])
        
        result_list.append({
            'trial_index': idx,
            'init_pos': row['init_pos'],
            't_left': time_marker - 0.5,
            't_right': time_marker + 0.5,
            'marker_time':time_marker,
            'marker_label':'Movement onset'
        })
    
    return pd.DataFrame(result_list)

def compute_instantaneous_rate(MO_pd, tsd_group, bin_size=50, sigma_seconds=50):

    spike_trains = []
    for unit in tsd_group:
        spike_times = unit.t
        # 创建 neo SpikeTrain 对象
        spiketrain = neo.SpikeTrain(
            times=spike_times, 
            units='s', 
            t_start=spike_times.min(), 
            t_stop=spike_times.max()
        )
        spike_trains.append(spiketrain)

    trials = nap.IntervalSet(
        start=MO_pd['t_left'].values,
        end=MO_pd['t_right'].values
    )

    is_kwargs = {
        'kernel' : GaussianKernel(sigma_seconds*pq.ms),
        'sampling_period' : bin_size*pq.ms
    }

    trial_ind = []
    kwargs_list = []
    for trial_idx, trial in enumerate(trials):
        kwargs = {'t_start': trial.start*pq.s-5*is_kwargs['sampling_period'],
                't_stop': trial.end*pq.s+5*is_kwargs['sampling_period'],
                'trial_index': trial_idx}
        kwargs_list.append(kwargs)
        trial_ind.append(trial_idx)


    sliced_is = SpikeStatistics.preprocessing('instantaneous_rate', 
                                            kwargs_list, spike_trains, **is_kwargs)[np.array(trial_ind),:,5:-5]

    return sliced_is

trial_info_df_se = extract_trial_info(trial_info_df)
trial = add_slices_to_trial_info(trial_info_df_se, labels_shift, markers_shift)
MO_pd = extract_MO_data(trial)
kilosort_unit = units[units['sorter'] == "kilosort2.5"]
input_array = compute_instantaneous_rate(MO_pd, kilosort_unit)
print(kilosort_unit)