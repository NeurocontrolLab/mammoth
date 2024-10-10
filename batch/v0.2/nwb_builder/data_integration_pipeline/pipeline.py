#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""
import argparse
import subprocess
import os

#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument("-f", "--file", type=str,
                    default='/media/lenovo/Extreme Pro/Bohr data/20240201_br_utah_interception_120_semi_brain_control_4_001', 
                    metavar='/the/path/your/data/located/in', help='input folder')

parser.add_argument('-o', '--output', type=str, 
                    default='/media/lenovo/Extreme Pro/Bohr data/20240201_br_utah_interception_120_semi_brain_control_4_001', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()
raw_dirname = os.path.abspath(args.file)
if not os.path.exists(args.output):
    os.makedirs(args.output)

#%% call pipeline of three type of neo data
FILEPATH = os.path.dirname(os.path.abspath(__file__))

#------------------------------------------------------------------------------
# convert neural data (spiketrain, lfp, and event marker)
#------------------------------------------------------------------------------

command = 'python ' + os.path.join(FILEPATH,'neural_data',
                                   'br_share','pipeline.py')+\
                                    ' -f '+ raw_dirname\
                                    +' -o '+os.path.join(args.output,
                                                 'neural_data.db')

ret = subprocess.run(command,shell=True,
                     stdout=subprocess.PIPE,stderr=subprocess.PIPE)

#------------------------------------------------------------------------------
# convert continuous behavior data (EMG and jrajectory)
#------------------------------------------------------------------------------

command = 'python ' + os.path.join(FILEPATH,'continuous_behavior','aie_share',
                                   'pipeline.py')+' -f '+raw_dirname\
                                    +' -o '+os.path.join(args.output,
                                                    'continuous_behavior.db')

ret = subprocess.run(command,shell=True,
                     stdout=subprocess.PIPE,stderr=subprocess.PIPE)