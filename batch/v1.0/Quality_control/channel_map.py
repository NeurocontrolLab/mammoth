import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from probeinterface import read_probeinterface,Probe, ProbeGroup
from probeinterface.plotting import plot_probe

#%% parse the input arguments
parser = argparse.ArgumentParser(argument_default=None)
parser.add_argument("-f", "--file", type=str,
                    default='/media/lenovo/Extreme Pro/Bohr data/20240201_br_utah_interception_120_semi_brain_control_4_001', 
                    metavar='/the/path/your/data/located/in', help='input folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Abel/Data_recording/20240705_centerOut_001', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()
raw_dirname = args.file

with open('/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/SN+11386-000049.cmp','r') as f:
    probe_info = f.readlines()[14::]

probe_info = [i.split() for i in probe_info][0:-1]
array_name = ['elec'] if '-' not in probe_info[0][-1]\
    else list(set([i[-1].split('-')[0] for i in probe_info]))

nchannels_per_array = len(probe_info)
electrode_counter = []
probegroup = ProbeGroup()

for array in array_name:
    # create an electrode group for this shank
    
    probe_2d = Probe(ndim=2, si_units='um')
    probe_2d.annotate(
        name = array, 
        manufacturer="blackrock microsystem",
        escription = 'one 96 Utah array'
        )
    positions = []
    device_channel = []
    # test2 = []
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

fig, (ax1,ax2) = plt.subplots(2,1)
plot_probe(probegroup.probes[0], ax=ax1)
plot_probe(probegroup.probes[1], ax=ax2)
ax2.title.set_text('')

#%% save to appointed path
if not os.path.exists(os.path.join(args.output,'description')):
    os.mkdir(os.path.join(args.output,'description'))
plt.savefig(os.path.join(args.output, 'description', 'channel_map.png'), dpi=300,bbox_inches = 'tight')