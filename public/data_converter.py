#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 18:46:27 2024

@author: cuilab
"""

import argparse
import subprocess
# trodes and trodesexport are necessary

def convert_file(data_dir):
    command = """dir='{}'
    if ! find $dir -type d -name "*kilosort" | read; then
      echo "convert .rec files"
      rec_files=$(find $dir -type f -name "*.rec" -exec printf "%s\n" {} +)
      echo $rec_files
      ulimit -n 10240 && module load trodes/2.4.2 && trodesexport -rec $rec_files -lfp -lfplowpass 300 -dio -kilosort -spikeband -spikehighpass 300 -spikelowpass 6000 -spikes -thresh 50
    else
      echo "already converted"
    fi""".format(data_dir)
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


#%%
parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-d", "--data", type=str,
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Brain_control/20240416_sg_fa_2d_centerout_brain_control/raw_data', 
                    metavar='/the/path/your/data/located/in', help='raw_data folder')

args = parser.parse_args()

convert_file(args.data)
