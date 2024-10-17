#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 13:21:20 2024

@author: cuilab
"""
import os
import numpy as np
import argparse
from probeinterface import Probe, ProbeGroup, write_probeinterface

def get_probe_json(map_path, output_path, file_name):
    with open(map_path,'r') as f:
        probe_info = f.readlines()[14::]
    
    probe_info = [i.split() for i in probe_info][0:-1]
    
    array_name = ['elec'] if '-' not in probe_info[0][-1]\
        else list(set([i[-1].split('-')[0] for i in probe_info]))
    
    array_name.sort()
    
    probegroup = ProbeGroup()
    for array in array_name:
        # create an electrode group for this shank
        probe_2d = Probe(ndim=2, si_units='um')
        probe_2d.annotate(
            name = array, 
            manufacturer="blackrock microsystem",
            escription = '96-ch Utah array x 2'
            )
        positions = []
        device_channel = []
        
        # add electrodes to the electrode table
        for ielec in probe_info:
            if array not in ielec[-1]:
                continue
        
            positions.append([float(ielec[0])*400, float(ielec[1])*400])
            device_channel.append((ord(ielec[2])-65)*32+int(ielec[3])-1)
            # test2.append(ielec[2]+ielec[3])
            
        probe_2d.set_contacts(positions=np.array(positions), 
                              shapes='circle', 
                              shape_params={'radius': 20})
        probe_2d.set_device_channel_indices(device_channel)
        probe_2d.create_auto_shape(probe_type='tip')
        probegroup.add_probe(probe_2d)
        
    write_probeinterface(os.path.join(output_path, file_name+'.json'),
                         probegroup)


#%%
parser = argparse.ArgumentParser(argument_default=None)

parser.add_argument('-mp', '--map_path', 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/SN+11386-000049.cmp')

parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()
get_probe_json(args.map_path, args.output, "Bohr_Utah_96x2")
