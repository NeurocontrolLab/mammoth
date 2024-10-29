#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

import argparse
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import yaml
import copy
import pandas as pd
import numpy as np
import quantities as pq
import joblib
from SmartNeo.user_layer.dict_to_neo import templat_neo
from SmartNeo.interface_layer.nwb_interface import NWBInterface
from dependencies.user_input_entry_collection import AIEShare as As
import re
import shutil
import json
import os
import fnmatch
from datetime import datetime

def find_files_with_name(root_dir, pattern):
    """查找指定目录及其子目录中符合给定模式的文件"""
    matches = []
    for root, dirs, files in os.walk(root_dir):
        for filename in fnmatch.filter(files, pattern):
            matches.append(os.path.join(root, filename))
    return matches

def datetime_to_seconds(time_str):
    # 定义时间字符串的格式
    time_format = "%Y-%m-%d %H:%M:%S.%f"
    
    # 将字符串解析为 datetime 对象
    dt_obj = datetime.strptime(time_str, time_format)
    
    # 计算当天00:00:00到这个时间的总秒数
    total_seconds = dt_obj.hour * 3600 + dt_obj.minute * 60 + dt_obj.second + dt_obj.microsecond / 1_000_000
    
    return total_seconds


def format_file(root_dir, output_dir):
    #%% load template

    FILEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    with open(os.path.join(FILEPATH, 'dependencies', 'template_bhv_data.yml')) as f:
        Template = yaml.safe_load(f)
        
    #%% convert AnalogData
    root_directory = os.path.join(root_dir,'bhv')
    file_pattern = '*.log'
    matching_files = find_files_with_name(root_directory, file_pattern)
    with open(matching_files[0],'r') as f:
        bhv_data = f.readlines()

    times = []
    bhv_info = []
    name_info = []

    for i in bhv_data:

        match = re.match(r'^([\d-]+\s[\d:.]+)\s(.+)$', i)
        timestamp = match.group(1)  # 时间戳部分
        content = match.group(2)    # 后面的内容部分
        times.append(datetime_to_seconds(timestamp))
        content_dict = json.loads(content)

        if len(list(content_dict.keys()))==1:
            name_info.append(list(content_dict.keys())[0])
            name = list(content_dict.keys())[0]
            info = content_dict[name]
        else:
            name_info.append('info')
            info = content_dict

        bhv_info.append(json.dumps(info))


    bhv_df = pd.DataFrame({
        'times': times,
        'name_info': name_info,
        'bhv_info': bhv_info
    })

    grouped_dfs = {name: group for name, group in bhv_df.groupby('name_info')}
    InputList = []

    for i in grouped_dfs:
        InputData = copy.deepcopy(Template)
        InputData['Event'] = {}
        InputData['Event']['event'] = {}
        InputData['Event']['event'] = templat_neo['event'].copy()
        InputData['Event']['event']['labels'] = grouped_dfs[i]['bhv_info'].to_list()
        InputData['Event']['event']['times'] = np.array(grouped_dfs[i]['times'].to_list())*pq.s
        InputData['name'] = i
        InputList.append(InputData)


    root_directory = os.path.join(root_dir,'bhv')
    file_pattern = [i for i in os.listdir(root_directory) if ('csv' in i) and (not 'meta' in i)]
    if len(file_pattern)!=0:
        file_path = os.path.join(root_directory, file_pattern[0])  # 这里替换成你的C3D文件路径
        df = pd.read_csv(file_path, skiprows=3, index_col=0)
        df.dropna(how='all', inplace=True)

        InputData = copy.deepcopy(Template)
        InputData['IrregularSampledData'] = {}
        InputData['IrregularSampledData']['irr'] = {}
        InputData['IrregularSampledData']['irr'] = templat_neo['irr'].copy()
        InputData['IrregularSampledData']['irr']['signal'] = df[['X', 'Y', 'Z']].iloc[1:].to_numpy().astype(float)*pq.mm
        InputData['IrregularSampledData']['irr']['times'] = np.array(df.iloc[1:].index.to_list())/100*pq.s
        InputData['name'] = 'Vicon motion'
        InputList.append(InputData)

    #%% operate the dicts
    continuous_bhv_block = As.data_input(user_data = InputList,
                                            index = 3,
                                            name = 'continuous_behavior')

    #%% save to appointed path
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        
    nwb_saver = NWBInterface()
    nwb_saver.save_nwb(blockdata = continuous_bhv_block, 
                        filename = os.path.join(output_dir,'continuous_behavior.nwb'))


#%% parse the input arguments
parser = argparse.ArgumentParser(description='extract trial info')
parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Data_recording/20240925_Interception_002', 
                    metavar='/the/path/your/data/located/in', help='root folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Data_recording/20240925_Interception_002/formatted_data', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()

format_file(args.root, args.output)
