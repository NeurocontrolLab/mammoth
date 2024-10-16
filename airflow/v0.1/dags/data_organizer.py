#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 19:34:42 2024

@author: lichenyang
@modifier: chenyun
"""

import os
import shutil


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
        The name of target folder.
    suffix_name_list : list
        The list of suffix names which should be included.
    file_walk: tuple
        This results from os.walk, as (current directory, sub directory, files) 
    root_path: str
        The session directory, i.e. the parent directory of target folder
    Returns
    -------
    None.

    '''
    
    i = file_walk
    
    # if files' formats in suffix_name_list, their paths are collected
    rec_file = [os.path.join(i[0], j) for j in i[2]\
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
    

