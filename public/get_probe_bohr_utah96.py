#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 15:42:01 2024

@author: cuilab
"""

import numpy as np
from probeinterface import Probe, ProbeGroup


def get(map_path):

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
            description = 'one 96 Utah array'
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
        
    return probegroup