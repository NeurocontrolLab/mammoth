<<<<<<< HEAD
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 22:03:01 2024

@author: cuilab
"""

import argparse
import os
import matplotlib.pyplot as plt
from get_probe_bohr_utah96 import get as get_probe
from probeinterface.plotting import plot_probe


def run(map_path, output_dir):
    probegroup = get_probe(map_path)

    fig, (ax1,ax2) = plt.subplots(2,1)
    plot_probe(probegroup.probes[0], ax=ax1)
    plot_probe(probegroup.probes[1], ax=ax2)
    ax2.title.set_text('')

    #%% save to appointed path
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    plt.savefig(os.path.join(output_dir, 'channel_map.png'), dpi=300, bbox_inches = 'tight')

#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Abel/Data_recording/20240705_centerOut_001/description', 
                    metavar='/the/path/you/want/to/save', help='output folder')
parser.add_argument('-mp', '--map_path', 
                    default='/home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Spike_sorting/SN+11386-000049.cmp')

args = parser.parse_args()

run(args.map_path, args.output)
=======
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 22:03:01 2024

@author: cuilab
"""

import argparse
import os
import matplotlib.pyplot as plt
from get_probe_bohr_utah96 import get as get_probe
from probeinterface.plotting import plot_probe


def run(map_path, output_dir):
    probegroup = get_probe(map_path)

    fig, (ax1,ax2) = plt.subplots(2,1)
    plot_probe(probegroup.probes[0], ax=ax1)
    plot_probe(probegroup.probes[1], ax=ax2)
    ax2.title.set_text('')

    #%% save to appointed path
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    plt.savefig(os.path.join(output_dir, 'channel_map.png'), dpi=300, bbox_inches = 'tight')

#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Abel/Data_recording/20240705_centerOut_001/description', 
                    metavar='/the/path/you/want/to/save', help='output folder')
parser.add_argument('-mp', '--map_path', 
                    default='/home/cuihe_lab/lichenyang/DATA_AMAX/Neucyber-NC-2023-A-01/Bohr/Spike_sorting/SN+11386-000049.cmp')

args = parser.parse_args()

run(args.map_path, args.output)
>>>>>>> 30ed201402aebf165a59b0f25615ba2e47694b03
