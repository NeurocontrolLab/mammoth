#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

import os
import glob
import argparse
import hdf5storage
from dependencies.user_input_entry_collection import Ml2NeoTrial as mt
from SmartNeo.interface_layer.nwb_interface import NWBInterface

def format_file(root_dir, output_dir):
    FILEPATH = os.path.dirname(os.path.abspath(__file__))
    #%% load and convert .mat monkeylogic data
    bhv_name = glob.glob(os.path.join(os.path.join(root_dir,'bhv'),'*.bhv2'))[0].split('/')[-1].split('.')[0]
    bhv_mat = [i for i in glob.glob(os.path.join(root_dir,'bhv','*.mat')) if bhv_name in i][0]
    bhvsave = hdf5storage.loadmat(bhv_mat)
    bhvkey = [i for i in bhvsave if '_' not in i][0]
    bhvsave = bhvsave[bhvkey]

    #%% operate the dicts
    mlblock = mt.data_input(user_data = bhvsave.squeeze(), index = 2, name = 'trial_behavior')

    #%% save to appointed path
    nwb_saver = NWBInterface()
    nwb_saver.save_nwb(blockdata = mlblock, 
                    filename = os.path.join(output_dir,'trial_behavior.nwb'))


#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-r", "--root", type=str,
                    default='/AMAX/cuihe_lab/lichenyang/Neucyber-NC-2023-A-01/Qianqian/Data_recording/20231207_kilo/', 
                    metavar='/the/root/path/your/data/located/in', help='root folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/lichenyang/Neucyber-NC-2023-A-01/Qianqian/Data_recording/20231207_kilo/formatted_data', 
                    metavar='/the/path/you/want/to/save', help='output folder')
args = parser.parse_args()

format_file(args.root, args.output)
