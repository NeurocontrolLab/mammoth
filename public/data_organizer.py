#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 19:34:42 2024

@author: cuilab
"""

import os
import shutil
import argparse

#%%
global target_folder_name_list, rec_name, video_name, bhv_name
target_folder_name_list = ['raw_data','sorted_data','formatted_data','bhv',
                           'video','hand_trajectory','eye_tracking','emg',
                           'description','metadata.json','quality_control.ipynb',
                           '']

rec_name = ['.ns6','.ns3','.ns2','.ccf','.nev','.rec']
video_name = ['.mkv', '.mp4', 'webm']
bhv_name = ['.png', '.mat', '.log', '.bhv2', '.psydat']


#%%
def handle_file(target_folder_name, suffix_name_list, file_walk, root_dir):
    '''
    This function aims to distribute certain files into specified folders. 

    Parameters
    ----------
    target_folder_name : str
        The name of the target folder.
    suffix_name_list : list
        The list of the suffix names which should be included in the target folder.
    file_walk: tuple
        This results from os.walk, as (current directory, sub directory, files) 
    root_path: str
        The parent directory of target folder

    Returns
    -------
    None.

    '''
    
    # if files' formats in suffix_name_list, their paths are collected
    rec_file = [os.path.join(file_walk[0], j) for j in file_walk[2]\
                if os.path.splitext(j)[-1] in suffix_name_list]
    
    if len(rec_file)>0:
        target_folder_dir = os.path.join(root_dir, target_folder_name)
        if not os.path.exists(target_folder_dir):
            os.mkdir(target_folder_dir)
        for s_f in rec_file:
            shutil.move(s_f, target_folder_dir)


#%%
def organize_file(root_dir):
    file_walk_list = [j for j in os.walk(root_dir)]
    for i in file_walk_list:
        # if i is one of target folders, continue
        if len([j for j in target_folder_name_list if j in i[0].split(root_dir)[1]])>1:
            continue
                
        # organize recording files
        handle_file(target_folder_name_list[0], rec_name, i, root_dir)
        
        # organize video files
        handle_file(target_folder_name_list[4], video_name, i, root_dir)
        
        # organize bhv files
        handle_file(target_folder_name_list[3], bhv_name, i, root_dir)
    
    # delete directory with unclassified names
    for i in os.listdir(root_dir):
        if (not i in target_folder_name_list) and (not '.' in i):
            shutil.rmtree(os.path.join(root_dir, i))
            
    print('File organization in %s completed.' % root_dir)


#%%
if __name__ == "__main()__":
    parser = argparse.ArgumentParser(argument_default=None)
    parser.add_argument("-f", "--file", type=str,
                        default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Brain_control/20240416_sg_fa_2d_centerout_brain_control', 
                        metavar='/the/path/your/data/located/in', help='input folder')

    args = parser.parse_args()
    raw_dirname = args.file
    organize_file(raw_dirname)
