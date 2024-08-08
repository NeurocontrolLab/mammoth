#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 16:49:17 2024

@author: lichenyang
"""
import argparse
import os
import shutil

def handle_file(target_file_name, suffix_name_list):
    rec_file = [os.path.join(i[0], j) for j in i[2]\
                if os.path.splitext(j)[-1] in suffix_name_list]
    if len(rec_file)>0:
        raw_data_name = os.path.join(raw_dirname, target_file_name)
        if not os.path.exists(raw_data_name):
            os.mkdir(raw_data_name)
        for s_f in rec_file:
            shutil.move(s_f, raw_data_name)

#%% give interface
parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-f", "--file", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Brain_control/20240416_sg_fa_2d_centerout_brain_control', 
                    metavar='/the/path/your/data/located/in', help='input folder')

args = parser.parse_args()
raw_dirname = args.file
target_file_name_list = ['raw_data','sorted_data','formatted_data','bhv',
                    'video','hand_trajectory','eye_tracking','emg','description',
                    'metadata.json','quality_control.ipynb','']

rec_name = ['.ns6','.ns3','.ns2','.ccf','.nev','.rec']
video_name = ['.mkv', '.mp4', 'webm']
bhv_name = ['.png', '.mat', '.log', '.bhv2', '.psydat']
#%% orgnize file
walk_file = [j for j in os.walk(raw_dirname)]
for i in walk_file:
    if len([j for j in target_file_name_list if j in i[0].split(raw_dirname)[1]])>1:
        continue
    
    #%% orgize recording file
    handle_file(target_file_name_list[0], rec_name)
    
    #%% orgize video file
    handle_file(target_file_name_list[4], video_name)
    
    #%% bhv
    handle_file(target_file_name_list[3], bhv_name)

for i in os.listdir(raw_dirname):
    if (not i in target_file_name_list) and (not '.' in i):
        shutil.rmtree(os.path.join(raw_dirname, i))
    
        
