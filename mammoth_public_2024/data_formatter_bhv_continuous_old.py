#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""


import os
import re
import yaml
import copy
import pandas as pd
import numpy as np
import numba
from numba.typed import List
import fnmatch
import argparse
import quantities as pq
from user_input_entry_collection import AIEShare as As
from SmartNeo.user_layer.dict_to_neo import templat_neo
from SmartNeo.interface_layer.nwb_interface import NWBInterface


def find_files_with_name(root_dir, pattern):
    """find files with name fitting given pattern in given root and sub directory"""
    matches = []
    for root, dirs, files in os.walk(root_dir):
        for filename in fnmatch.filter(files, pattern):
            matches.append(os.path.join(root, filename))
    return matches


def format_file(root_dir, output_dir):
    #%% load template
    FILEPATH = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(FILEPATH,'template_bhv_data.yml')) as f:
        Template = yaml.safe_load(f)
        
    #%% convert AnalogData
    root_directory = os.path.join(root_dir, 'bhv')
    file_pattern = '*.log'
    matching_files = find_files_with_name(root_directory, file_pattern)
    with open(matching_files[0], 'r') as f:
        bhv_data = f.readlines()

    @numba.jit
    def get_list(l1, l2, word): 
        for i in l1:
            if word in i:
                l2.append(i)
                
    bhv_numba_data = List(bhv_data)
    frame_numba_data = List(['a'])
    get_list(bhv_numba_data, frame_numba_data, 'frame info')

    pos_numba_data = List(['a'])
    get_list(bhv_numba_data, pos_numba_data, 'pos:')

    state_numba_data = List(['a'])
    get_list(bhv_numba_data, state_numba_data, 'state:')

    coeff_numba_data = List(['a'])
    get_list(bhv_numba_data, coeff_numba_data, 'coefficient:')

    event_numba_data = List(['a'])
    get_list(bhv_numba_data, event_numba_data, 'Event marker')

    trial_numba_data = List(['a'])
    get_list(bhv_numba_data, trial_numba_data, 'Trial info')

    frame_data = list(frame_numba_data[1::])
    pos_data = list(pos_numba_data[1::])
    state_data = list(state_numba_data[1::])
    coeff_data = list(coeff_numba_data[1::])
    event_data = list(event_numba_data[1::])
    trial_data = list(trial_numba_data[1::])

    var_name_extract = lambda x: [[k for k in j if '  ' not in k and not 'array' in k] \
                                for j in [re.findall(r"[A-Z_\ a-z][A-Z_\ a-z]+", i) \
                                        for i in x.split(':')] \
                                            if not 'INFO' in j and len(j)!=0]
    filter_var_name = lambda x: [i[0] for i in var_name_extract(x) if len(i)!=0]
    make_re_word = lambda x: r'\[INFO\] \[(.*)\] \[.*'+''.join(['('+i+'.*).*' for i in filter_var_name(x) \
                                                                if not 'Publishing' in i])+'\n' 

    frame_re_word = r'\[INFO\] \[(.*)\] \[.*frame info: (.*).*'
    pos_re_word = make_re_word(pos_data[1])
    state_re_word = make_re_word(state_data[1])
    coeff_re_word = make_re_word(coeff_data[1])
    event_re_word = make_re_word(event_data[1])
    trial_re_word = make_re_word(trial_data[1])

    def parsor(data, re_word):
        parse_data = []
        for i in data:
            gps = list(re.search(re_word,i).groups())
            gps[0] = 'time: '+gps[0]
            ele = {}
            for j in gps:
                input_j = j[1::] if j[0]==' ' else j
                #input_j = input_j.split(',')[0]
                key, val = input_j.split(':')
                if len(val)==0:
                    continue
                val = val[1::] if val[0]==' ' else val
                val = val.replace(',','')
                val = val.replace('\n','')
                try:
                    val = float(val) if not '[' in val else [float(k) for k in re.search(r'\[(.*)\]',val).groups()[0].split()]
                except:
                    val = [float(k) for k in re.search(r'\[(.*)\]',val).groups()[0].split()]
                ele[key] = val
            parse_data.append(ele)
        return parse_data

    # frame = [json.loads(i) for i in frame_data]
    frame = [list(re.search(frame_re_word,i).groups()) for i in frame_data]
    frame = [{'time':float(i[0]),'frame info':i[1]} for i in frame]
    frame = pd.DataFrame(frame)
    pos = pd.DataFrame(parsor(pos_data, pos_re_word))
    state = pd.DataFrame(parsor(state_data, state_re_word))
    coeff = pd.DataFrame(parsor(coeff_data, coeff_re_word))
    event = pd.DataFrame(parsor(event_data, event_re_word))
    trial = pd.DataFrame(parsor(trial_data, trial_re_word))


    InputList = []

    InputData = copy.deepcopy(Template)
    InputData['Event'] = {}
    InputData['Event']['event'] = {}
    InputData['Event']['event'] = templat_neo['event'].copy()
    InputData['Event']['event']['labels'] = event['Event marker'].tolist()
    InputData['Event']['event']['times'] = np.array(event['time'].tolist())*pq.s
    InputData['name'] = 'Event'
    InputList.append(InputData)

    InputData = copy.deepcopy(Template)
    InputData['Event'] = {}
    InputData['Event']['event'] = {}
    InputData['Event']['event'] = templat_neo['event'].copy()
    InputData['Event']['event']['labels'] = frame['frame info'].tolist()
    InputData['Event']['event']['times'] = np.array(frame['time'].tolist())*pq.s
    InputData['name'] = 'frame'
    InputList.append(InputData)

    #%% convert IrregularSampledData
    irr_list = [pos,state,coeff,trial]
    irr_name = ['pos','state','coeff','trial']
    for irr,irr_name in zip(irr_list,irr_name):
        InputData = copy.deepcopy(Template)

        InputData['AnalogData'] = 'null'
        InputData['Event'] = 'null'

        for ind, i in enumerate(irr.columns):
            if i=='time':
                continue
            InputData['IrregularSampledData'][i] = {}
            InputData['IrregularSampledData'][i]['irr'] = templat_neo['irr'].copy()
            a = irr[i]
            try:
                InputData['IrregularSampledData'][i]['irr']['signal'] = \
                                                            np.array(irr[i].tolist())*pq.dimensionless
            except:
                max_len = max([len(d) for d in irr[i].tolist()])
                buffer = np.zeros((len(irr[i]),max_len)).astype(float)
                for irr_ind, d in enumerate(irr[i].tolist()):
                    buffer[irr_ind,0:len(d)] = d.copy()
                    
                InputData['IrregularSampledData'][i]['irr']['signal'] = buffer*pq.dimensionless
                
            InputData['IrregularSampledData'][i]['irr']['times'] = \
                                                        np.array(irr['time'].tolist())*pq.s

        InputData['name'] = irr_name
        InputList.append(InputData)

    #%% operate the dicts
    continuous_bhv_block = As.data_input(user_data = InputList,
                                         index = 3,
                                         name = 'continuous_behavior')

    #%% save to appointed path
    nwb_saver = NWBInterface()
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    nwb_saver = NWBInterface()
    nwb_saver.save_nwb(blockdata = continuous_bhv_block, 
                       filename = os.path.join(output_dir, 'continuous_behavior.nwb'))


#%% parse the input arguments
parser = argparse.ArgumentParser(description='extract trial info')
parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/lichenyang/Neucyber-NC-2023-A-01/Bohr/Brain_control/20240306_br_utah_interception_120_semi_brain_control_001', 
                    metavar='/the/path/your/data/located/in', help='root folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/lichenyang/Neucyber-NC-2023-A-01/Bohr/Brain_control/20240306_br_utah_interception_120_semi_brain_control_001/formatted_data', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()

format_file(args.root, args.output)
