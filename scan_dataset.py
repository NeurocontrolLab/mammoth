#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 3 12:49:11 2024

@author: cuilab
"""

import os
import time
import pandas as pd
import argparse
import json


#%% define functions
def get_sessions(root_dir):
    data_menu = pd.DataFrame(columns=['subject', 'type', 'session', 'path'])
    for dir_path, dir_names, _ in os.walk(root_dir, topdown=True):
        for dir_name_i in dir_names:
            dir_path_i = os.path.join(dir_path, dir_name_i)
            relative_path = dir_path_i.split(root_dir)[1].split(os.sep)
            relative_depth = len(relative_path)-1
            if relative_depth == 3:          
                data_menu = data_menu._append({'subject': relative_path[1],
                                               'type': relative_path[2],
                                               'session': relative_path[3],
                                               'path': dir_path_i},
                                               ignore_index=True)
                                               
            else:
                break                                 
    
    return data_menu
    # data_menu.to_csv(os.path.join(meta_path, '%s_all_sessions_list.csv'
    #                  % time.strftime("%Y%m%d", time.localtime(time.time()))))
    # print('Already got all session folders.')
    

def scan_sessions(root_dir, output_dir):
    f_session_df = get_sessions(root_dir)
    f_session_df = f_session_df[f_session_df['subject'].isin(
        ['Abel', 'Bohr', 'Darwin', 'Galileo', 'Leibniz', 'Maxwell'])]
     
    for i in f_session_df.index:
        spath = f_session_df.loc[i, 'path']
        
        dir_list = next(os.walk(spath))[1]
        
        # organized
        f_session_df.loc[i, 'organized'] = int(
            len([f.name for f in os.scandir(spath) if not f.is_dir()])==0)
        
        # metadata
        f_session_df.loc[i, 'metadata'] = (
            1 if os.path.exists(os.path.join(spath, 'bhv', 'metadata.csv')) else 0)
       
        # formatted
        f_session_df.loc[i, 'TCR'] = (
            1 if os.path.exists(os.path.join(spath, 'formatted_data', 'neural_data_no_sort.nwb')) else 0)
        
        f_session_df.loc[i, 'bhv'] = (
            1 if os.path.exists(os.path.join(spath, 'formatted_data', 'continuous_behavior.nwb')) else 0)

        # checked
        f_session_df.loc[i, 'TCR_checked'] = (
             0 if not os.path.exists(os.path.join(spath, 'description')) 
             else int('qc_summary.json' in 
                      os.listdir(os.path.join(spath, 'description'))))
        
        # sorted
        if not os.path.exists(os.path.join(spath, 'sorted_data')):
            f_session_df.loc[i, 'sorted'] = '0%'
            continue
        
        kspath = os.path.join(spath, 'sorted_data/kilosort2_5_output') 
               
        subfolders = [f.name for f in os.scandir(kspath) if f.is_dir()]

        if len(subfolders)==0:
            f_session_df.loc[i, 'sorted'] = 'wrong'
            continue

        wrong_shanks = [i for i in next(os.walk(kspath))[1] if 'wrong' in i]
        
        unsort_shanks = []
        for shank in sorted(next(os.walk(kspath))[1]):
            if shank not in wrong_shanks:
                if 'cluster_info.tsv' not in os.listdir(
                        os.path.join(kspath, shank+'/sorter_output')):
                    unsort_shanks.append('shank_%s' % shank.split("_")[-1])

        wrong_shanks = ["shank_%s" % i.split("_")[-2] for i in wrong_shanks]

        f_session_df.loc[i, 'wrong_shanks'] = (
            ','.join(wrong_shanks) if len(wrong_shanks)>0 else 'None'
        )

        f_session_df.loc[i, 'unsort_shanks'] = (
            ','.join(unsort_shanks) if len(unsort_shanks)>0 else 'None'
        )

        f_session_df.loc[i, 'sorted'] = '%d%%' % int((1 - (len(unsort_shanks)+len(wrong_shanks))/len(next(os.walk(kspath))[1]))*100)

        # formatted
        f_session_df.loc[i, 'Spike+TCR+LFP'] = (
            1 if os.path.exists(os.path.join(spath, 'formatted_data', 'neural_data.nwb')) else 0)
        
        # quality control
        f_session_df.loc[i, 'Spike_checked'] = (
             0 if not os.path.exists(os.path.join(spath, 'description')) 
             else int('chn_consis_summary.png' in 
                      os.listdir(os.path.join(spath, 'description'))))
        
        # standardize             
        f_session_df.loc[i, 'standardNWB'] = (
            1 if os.path.exists(os.path.join(spath, 'formatted_data', 'standard_data.nwb')) else 0)
        
        if os.path.exists(os.path.join(spath, 'description', 'qc_summary.json')):
            with open(os.path.join(spath, 'description', 'qc_summary.json'), "r") as file:
                summary = json.load(file)
                
                if len([i for i in summary.keys() if 'diff_time' in i])>0:
                    diff_time = summary[[i for i in summary.keys() if 'diff_time' in i][0]]
                    f_session_df.loc[i, 'time_consistent'] = 1 if abs(diff_time[2])<30/1000 else 0
                
                if len([i for i in summary.keys() if 'neural correlation' in i])>0:
                    neural_corr = summary[[i for i in summary.keys() if 'neural correlation' in i][0]]
                    f_session_df.loc[i, 'neural_correlation'] = 1 if len([i for i in neural_corr if i>0.5])/len(neural_corr)>0.5 else 0

                if len([i for i in summary.keys() if 'channel consistency' in i])>0:
                    ch_consist = summary[[i for i in summary.keys() if 'channel consistency' in i][0]]
                    f_session_df.loc[i, 'channel_consistency'] = 1 if (
                        sum(ch_consist['unshuffled'])/len(ch_consist['unshuffled'])>0.5 and 
                        sum(ch_consist['shuffled'])/len(ch_consist['shuffled'])<0.2) else 0
        
        if not pd.isna(f_session_df.loc[i, 'time_consistent']):
            if (f_session_df.loc[i, 'time_consistent'] + f_session_df.loc[i, 'neural_correlation'] <2):
                os.rename(spath, spath+"_check")

    print('Already scan fit sessions.')  

    older = [f for f in os.listdir(output_dir) if ('.csv' in f) and ('dataset_overview' in f)]
    if len(older)>0:
        for o in older:
            os.remove(os.path.join(output_dir, o))

    f_session_df.to_csv(os.path.join(output_dir, 'dataset_overview_%s.csv' % time.strftime("%Y%m%d", time.localtime(time.time()))))
    # return f_session_df  

    
#%% 
parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01', 
                    metavar='/the/root/path/your/data/located/in', help='root folder')
parser.add_argument("-o", "--output", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01', 
                    metavar='/the/output/path', help='output folder')

args = parser.parse_args()

scan_sessions(args.root, args.output)

