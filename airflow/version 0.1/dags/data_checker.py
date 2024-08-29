#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 16:26:50 2024

@author: chenyun
"""

import os
from string import Template

code_dict = {}
code_dict['brain_control'] = Template(r'''{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "import argparse\n",
    "import os\n",
    "import json\n",
    "import copy\n",
    "import pandas as pd\n",
    "import numpy as np\n",
    "import quantities as pq\n",
    "import joblib\n",
    "import re\n",
    "import shutil\n",
    "import seaborn as sns\n",
    "from dtw import accelerated_dtw\n",
    "from SmartNeo.analysis_layer.tools.dataset_maker.dataset_maker import DatasetMaker\n",
    "from SmartNeo.user_layer.dict_to_neo import templat_neo\n",
    "from SmartNeo.interface_layer.nwb_interface import NWBInterface\n",
    "from SmartNeo.analysis_layer.spike_preprocessor3 import SpikeStatistics\n",
    "import matplotlib.pyplot as plt\n",
    "from matplotlib.gridspec import GridSpec\n",
    "from probeinterface import read_probeinterface,Probe, ProbeGroup\n",
    "from probeinterface.plotting import plot_probe"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "raw_dirname = '${session0}'\n",
    "\n",
    "nwb_saver = NWBInterface()\n",
    "bhv_data = nwb_saver.read_nwb(filename = os.path.join(raw_dirname,'formatted_data', 'continuous_behavior.nwb'))\n",
    "\n",
    "nwb_saver = NWBInterface()\n",
    "neural_data = nwb_saver.read_nwb(filename = os.path.join(raw_dirname,'formatted_data', 'neural_data_no_sort.nwb'))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Channel Map"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "description_path = os.path.join(raw_dirname, 'description')\n",
    "cmp_name = [i for i in os.listdir(description_path) if 'cmp' in i][0]\n",
    "\n",
    "with open(os.path.join(raw_dirname, 'description', cmp_name),'r') as f:\n",
    "    probe_info = f.readlines()[14::]\n",
    "\n",
    "probe_info = [i.split() for i in probe_info][0:-1]\n",
    "array_name = ['elec'] if '-' not in probe_info[0][-1]\\\n",
    "    else list(set([i[-1].split('-')[0] for i in probe_info]))\n",
    "\n",
    "nchannels_per_array = len(probe_info)\n",
    "electrode_counter = []\n",
    "probegroup = ProbeGroup()\n",
    "\n",
    "for array in array_name:\n",
    "    # create an electrode group for this shank\n",
    "    \n",
    "    probe_2d = Probe(ndim=2, si_units='um')\n",
    "    probe_2d.annotate(\n",
    "        name = array, \n",
    "        manufacturer=\"blackrock microsystem\",\n",
    "        escription = 'one 96 Utah array'\n",
    "        )\n",
    "    positions = []\n",
    "    device_channel = []\n",
    "    # test2 = []\n",
    "    # add electrodes to the electrode table\n",
    "    for ielec in probe_info:\n",
    "        if array not in ielec[-1]:\n",
    "            continue\n",
    "    \n",
    "        positions.append([float(ielec[0])*400, float(ielec[1])*400])\n",
    "        device_channel.append((ord(ielec[2])-65)*32+int(ielec[3])-1)\n",
    "        # test2.append(ielec[2]+ielec[3])\n",
    "        \n",
    "    \n",
    "    probe_2d.set_contacts(positions=np.array(positions), \n",
    "                          shapes='circle', \n",
    "                          shape_params={'radius': 20})\n",
    "    probe_2d.set_device_channel_indices(device_channel)\n",
    "    probe_2d.create_auto_shape(probe_type='tip')\n",
    "    probegroup.add_probe(probe_2d)\n",
    "\n",
    "fig, (ax1,ax2) = plt.subplots(2,1)\n",
    "plot_probe(probegroup.probes[0], ax=ax1)\n",
    "plot_probe(probegroup.probes[1], ax=ax2)\n",
    "ax2.title.set_text('')"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Trial behavior"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "coeff_pd = {}\n",
    "for i in bhv_data.segments[1].irregularlysampledsignals:\n",
    "    coeff_pd[i.name[1]] = i.as_array().squeeze()\n",
    "coeff_pd = pd.DataFrame(coeff_pd)\n",
    "# coeff_pd.rename_axis('Times')\n",
    "coeff_pd.index = bhv_data.segments[1].irregularlysampledsignals[0].times\n",
    "coeff_pd.index.name = 'Times'\n",
    "coeff_pd.head()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "pos_pd = {}\n",
    "for i in bhv_data.segments[3].irregularlysampledsignals:\n",
    "    pos_pd[i.name[1]] = list(i.as_array().squeeze())\n",
    "\n",
    "pos_pd = pd.DataFrame(pos_pd)\n",
    "pos_pd.rename_axis('Times')\n",
    "pos_pd.index = bhv_data.segments[3].irregularlysampledsignals[0].times\n",
    "pos_pd.index.name = 'Times'\n",
    "pos_unit = bhv_data.segments[3].irregularlysampledsignals[0].times.units\n",
    "\n",
    "pos_pd.head()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "trial_pd = {}\n",
    "for i in bhv_data.segments[-1].irregularlysampledsignals:\n",
    "    trial_pd[i.name[1]] = list(i.as_array().squeeze())\n",
    "\n",
    "trial_pd = pd.DataFrame(trial_pd)\n",
    "trial_pd.rename_axis('Times')\n",
    "trial_pd.index = bhv_data.segments[-1].irregularlysampledsignals[0].times\n",
    "trial_pd.index.name = 'Times'\n",
    "trial_pd.head()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "events = bhv_data.segments[0].events[0]\n",
    "marker_24_pos = np.where(events.labels==24)[0]\n",
    "marker_5_pos = np.where(events.labels==5)[0]\n",
    "\n",
    "trial_event = [events[i:j+1] for i,j in zip(marker_24_pos,marker_5_pos)]\n",
    "\n",
    "marker_24_time = events.times[np.where(events.labels==24)[0]]\n",
    "marker_5_time = events.times[np.where(events.labels==5)[0]]\n",
    "trial_pos = [pos_pd.iloc[np.where(pos_pd.index<j * (pos_pd.index>i))[0]] for i,j in zip(marker_24_time,marker_5_time)]\n",
    "trial_coeff = [coeff_pd.iloc[np.where(coeff_pd.index<j * (coeff_pd.index>i))[0]] for i,j in zip(marker_24_time,marker_5_time)]"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Time consistency check"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "neural_event = [i for i in neural_data.segments if i.name=='RecordingSystemEvent'][0]\n",
    "neural_event_labels = neural_event.events[0].labels\n",
    "neural_event_times = neural_event.events[0].times\n",
    "neural_event_labels[neural_event_labels==18]=20\n",
    "\n",
    "ml_event_labels = bhv_data.segments[0].events[0].labels\n",
    "ml_event_times = bhv_data.segments[0].events[0].times\n",
    "\n",
    "d_f = [ml_event_labels[i:i+len(neural_event_labels)]-neural_event_labels for i in range(len(ml_event_labels)-len(neural_event_labels))]\n",
    "start_ind = np.argmin([sum(np.abs(i)) for i in d_f])\n",
    "\n",
    "ml_event_labels = bhv_data.segments[0].events[0].labels[start_ind::]\n",
    "ml_event_times = bhv_data.segments[0].events[0].times[start_ind::]\n",
    "\n",
    "li_distance = lambda x, y: 0 if np.abs(x - y)==0 else len(neural_event_labels)\n",
    "    # manhattan_distance = lambda x, y: np.abs(x - y)\n",
    "\n",
    "_,_,_,path1 = accelerated_dtw(ml_event_labels, neural_event_labels, dist=li_distance)\n",
    "path1 = np.array(path1).T\n",
    "zero_dir = (ml_event_labels[path1[:,0]]-neural_event_labels[path1[:,1]])==0\n",
    "diff_time = ml_event_times[path1[:,0]]-neural_event_times[path1[:,1]]\n",
    "def remove_outliers(data, factor=1.5):\n",
    "    q1 = np.percentile(data, 25)\n",
    "    q3 = np.percentile(data, 75)\n",
    "    iqr = q3 - q1\n",
    "    lower_bound = q1 - (factor * iqr)\n",
    "    upper_bound = q3 + (factor * iqr)\n",
    "    return data[(data > lower_bound) & (data < upper_bound)]\n",
    "diff_time = remove_outliers(diff_time, factor=1.5)\n",
    "plt.boxplot(diff_time-np.mean(diff_time))\n",
    "len(diff_time)/len(neural_event_labels)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Behavior analysis"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "events = bhv_data.segments[0].events[0]\n",
    "marker_24_pos = np.where(events.labels==24)[0]\n",
    "marker_5_pos = np.where(events.labels==5)[0]\n",
    "\n",
    "trial = [events[i:j+1] for i,j in zip(marker_24_pos,marker_5_pos)]\n",
    "duration = [i.times[-1]-i.times[3] for i in trial]\n",
    "trial_success = [20 in i.labels for i in trial]\n",
    "trial_time = [(i.times[-1]-trial[0].times[0]).magnitude/60 for i in trial]\n",
    "\n",
    "flag_times = bhv_data.segments[1].irregularlysampledsignals[-1].times\n",
    "training_flag = np.array(bhv_data.segments[1].irregularlysampledsignals[-1])\n",
    "trial_flag = [training_flag[np.where(j[0]<flag_times * (j[-1]>flag_times))[0]] for j in trial]\n",
    "trial_flag_times = [flag_times[np.where(j[0]<flag_times * (j[-1]>flag_times))[0]] for j in trial]\n",
    "trial_flag_times = [(i-trial[0].times[0]).magnitude.squeeze()/60 for i in trial_flag_times]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "mean_success = []\n",
    "sliding_time = []\n",
    "mean_duration = []\n",
    "flag = []\n",
    "epoch = []\n",
    "paradigm_start = trial[0].times[0]\n",
    "\n",
    "for i in range(0,len(trial)-75):\n",
    "    mean_success.append(trial_success[i:i+75].count(True)/75)\n",
    "    mean_duration.append(np.mean(duration[i:i+75]))\n",
    "    sliding_time.append(trial_time[i+75])\n",
    "    epoch.append(i)\n",
    "    \n",
    "labels = ['Success rate']*len(mean_success) + ['Duration time']*len(mean_duration)\n",
    "value = mean_success + list(np.array(mean_duration)/max(mean_duration))\n",
    "time = sliding_time + sliding_time\n",
    "df = pd.DataFrame(zip(value, time, labels), \n",
    "                  columns = ['Percentage', 'Time (min)', 'labels'])\n",
    "df.head()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "plt.figure(figsize = (8,4))\n",
    "\n",
    "ax = sns.lineplot(data = df, x = 'Time (min)', y = 'Percentage', \n",
    "                  hue = 'labels', lw = 2.5,\n",
    "                 style = 'labels', dashes = False, markersize = 8 ,\n",
    "                 palette = ['firebrick', 'gray'])\n",
    "ax.set_xlim([df['Time (min)'][0], 38])\n",
    "ax.set_ylim([0,1.02])\n",
    "\n",
    "flag = []\n",
    "for i,j in zip(trial_flag,trial_flag_times):\n",
    "    if sum(i.squeeze()) != 0:\n",
    "        flag = flag + list(j.squeeze())\n",
    "    if sum(i.squeeze()) == 0 and len(flag) != 0:    \n",
    "        x=flag\n",
    "        y1=[0]*len(x)\n",
    "        y2=[1.02]*len(x)\n",
    "\n",
    "        ax.fill_between(x, y1, y2, #上限，下限\n",
    "                facecolor='gray', #填充颜色\n",
    "                alpha=0.3) #透明度\n",
    "        flag=[]\n",
    "\n",
    "for axis in ['bottom', 'left']:\n",
    "    ax.spines[axis].set_linewidth(2.5)\n",
    "    ax.spines[axis].set_color('0.2')\n",
    "    \n",
    "ax.spines['top'].set_visible(False)\n",
    "ax.spines['right'].set_visible(False)\n",
    "\n",
    "ax.tick_params(width = 2.5, color = '0.2')\n",
    "\n",
    "plt.xticks(size = 14, weight = 'bold', color = '0.2')\n",
    "plt.yticks(size = 14, weight = 'bold', color = '0.2')\n",
    "\n",
    "ax.set_xlabel(ax.get_xlabel(), fontsize = 18, weight = 'bold', color = '0.2')\n",
    "ax.set_ylabel(ax.get_ylabel(), fontsize = 18, weight = 'bold', color = '0.2')\n",
    "\n",
    "plt.legend(frameon = False, prop = {'weight':'bold', 'size':14}, \n",
    "           labelcolor = '0.2',\n",
    "           loc=4)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# trial_pos\n",
    "succ = []\n",
    "\n",
    "for i,j,k in zip(trial_pos, trial_coeff, trial):\n",
    "    \n",
    "    if (i.index[0]-trial_pos[0].index[0])/60 > 30:\n",
    "        continue\n",
    "        \n",
    "    if not sum(j['training flag'].tolist())==0:\n",
    "        continue\n",
    "        \n",
    "    succ.append(0)\n",
    "    if not 20 in k.labels:\n",
    "        continue\n",
    "    succ[-1]=1\n",
    "    traj = np.array(i['pos'].tolist())\n",
    "    plt.plot(traj[:,0], traj[:,1])\n",
    "\n",
    "print('success rate: {}, trial num: {}, sucess trial num: {} (in 30 min)'.format(sum(succ)/len(succ),len(succ), sum(succ)))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Neural correlation analysis"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "diff_time_mean = np.mean(diff_time)\n",
    "st = [i for i in neural_data.segments if i.name=='TCR'][0].spiketrains\n",
    "st_shift = [i.time_shift(diff_time_mean) for i in st if len(i.times)>3000]\n",
    "# len(st_shift)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "events = bhv_data.segments[0].events[0]\n",
    "marker_24_pos = np.where(events.labels==3)[0]\n",
    "marker_5_pos = np.where(events.labels==5)[0]\n",
    "trial_event = [events[i:j+1] for i,j in zip(marker_24_pos,marker_5_pos)]\n",
    "\n",
    "marker_24_time = events.times[np.where(events.labels==24)[0]]\n",
    "marker_5_time = events.times[np.where(events.labels==5)[0]]\n",
    "trial_pos = [pos_pd.iloc[np.where(pos_pd.index<j * (pos_pd.index>i))[0]] for i,j in zip(marker_24_time,marker_5_time)]\n",
    "trial_coeff = [coeff_pd.iloc[np.where(coeff_pd.index<j * (coeff_pd.index>i))[0]] for i,j in zip(marker_24_time,marker_5_time)]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "from elephant.kernels import GaussianKernel\n",
    "#%% data slicing\n",
    "kwargs_list = []\n",
    "trial_ind = []\n",
    "move_direction = []\n",
    "pos_input = []\n",
    "is_kwargs = {\n",
    "    'kernel' : GaussianKernel(20*pq.ms),\n",
    "    'sampling_period' : 20*pq.ms\n",
    "}\n",
    "\n",
    "for ind, i in enumerate(trial_pos):\n",
    "    if not sum(trial_coeff[ind]['training flag'].tolist())==0:\n",
    "        continue\n",
    "    if not 20 in trial[ind].labels:\n",
    "        continue\n",
    "    pos_time = np.array(i.index)*pos_unit\n",
    "    kwargs = {'t_start': pos_time[0]-5*is_kwargs['sampling_period'],\n",
    "              't_stop': pos_time[-1]+5*is_kwargs['sampling_period'],\n",
    "              'aligned_marker':[3],\n",
    "              'aligned_marker2':[5],\n",
    "              'trial_index': ind}\n",
    "    pos_input.append(i)\n",
    "    kwargs_list.append(kwargs)\n",
    "    trial_ind.append(ind)\n",
    "    move_direction.append(np.array(i['pos'].to_list()).squeeze()[-1])\n",
    "\n",
    "sliced_is = SpikeStatistics.preprocessing('instantaneous_rate', \n",
    "                                          kwargs_list, st_shift, **is_kwargs)[np.array(trial_ind)]\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "bin_num = [(i['t_stop'] - i['t_start']).rescale(is_kwargs['sampling_period'].units)/is_kwargs['sampling_period'] for i in kwargs_list]\n",
    "select_is = [sliced_is[ind,:,0:int(i.magnitude)] for ind, i in enumerate(bin_num)]\n",
    "select_is = [i[:,5:-5] for i in select_is]\n",
    "fr_sum = [np.sum(i.magnitude) for i in select_is]\n",
    "select_ind = max(np.where(np.array(fr_sum)==0)[0])+1\n",
    "\n",
    "select_is = select_is[select_ind::]\n",
    "pos_input = pos_input[select_ind::]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "pos = []\n",
    "speed = []\n",
    "from scipy import signal\n",
    "from scipy.ndimage import gaussian_filter\n",
    "\n",
    "for i,j in zip(pos_input,select_is):\n",
    "    traj = np.array(i['pos'].tolist())  \n",
    "    vel = np.array(i['vel'].tolist()) \n",
    "    traj_x = gaussian_filter(signal.resample(traj[:,0],j.shape[1]+1)[1::],sigma=1)\n",
    "    traj_y = gaussian_filter(signal.resample(traj[:,1],j.shape[1]+1)[1::],sigma=1)\n",
    "    vel_x = gaussian_filter(signal.resample(vel[:,0],j.shape[1]+1)[1::],sigma=1)\n",
    "    vel_y = gaussian_filter(signal.resample(vel[:,1],j.shape[1]+1)[1::],sigma=1)\n",
    "    pos.append(np.array([traj_x,traj_y]).T)\n",
    "    speed.append(np.array([vel_x,vel_y]).T)\n",
    "    plt.plot(traj_x, traj_y)\n",
    "X_pos = np.concatenate(pos,axis=0)\n",
    "X_speed = np.concatenate(speed,axis=0)\n",
    "Y_neu = np.concatenate(select_is,axis=1).T\n",
    "\n",
    "mat_data = {}\n",
    "mat_data['neu_data'] = [i.magnitude for i in select_is]\n",
    "mat_data['pos'] = pos\n",
    "mat_data['speed'] = speed\n",
    "\n",
    "import hdf5storage\n",
    "hdf5storage.savemat(os.path.join(description_path,'dataset_not_sort.mat'), mat_data)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "scrolled": false
   },
   "outputs": [],
   "source": [
    "from sklearn.linear_model import MultiTaskLasso\n",
    "from sklearn.model_selection import cross_val_score\n",
    "lasso_speed = MultiTaskLasso(fit_intercept=False)\n",
    "scores = cross_val_score(lasso_speed, Y_neu, X_speed, scoring='r2',cv=5,n_jobs=24)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "print(scores)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "lasso_speed.fit(Y_neu, X_speed)\n",
    "hist, bin_edges = np.histogram(lasso_speed.coef_, bins=50)\n",
    "# plt.plot(lasso_speed.coef_)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "import seaborn as sns\n",
    "import matplotlib.pyplot as plt\n",
    "import numpy as np\n",
    "import warnings\n",
    "\n",
    "warnings.filterwarnings(\"ignore\", category=UserWarning)\n",
    "# data = np.random.normal(size=1000)\n",
    "\n",
    "sns.distplot(lasso_speed.coef_, hist=True, kde=False, rug=False)\n",
    "\n",
    "plt.show()"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
''')

code_dict['data_recording'] = ''' '''

def check_file(session):
    
    dpath = os.path.join(session, 'description')
    if os.path.exists(os.path.join(dpath, 'quality_control_no_sort.ipynb')):
        print('Check code already exists')
        return
        
    paradigm = ''
    paradigm = 'brain_control' if 'brain_control' in session.lower() else paradigm
    paradigm = 'data_recording' if 'data_recording' in session.lower() else paradigm
    
    quality_control = code_dict[paradigm].substitute(session0=session)
    
    if (os.path.exists(os.path.join(os.path.join(session, 'formatted_data'),
                                    'neural_data_no_sort.nwb')) & 
        os.path.exists(os.path.join(os.path.join(session, 'formatted_data'),
                                    'continuous_behavior.nwb'))):
         
        if not os.path.exists(dpath):
            os.mkdir(dpath)
            
        with open(os.path.join(dpath, 'quality_control_no_sort.ipynb'), 'w') as file:
            file.write(quality_control)

# jupyter nbconvert --execute --to notebook --inplace --allow-errors quality_control.ipynb
