import re
import json
import pandas as pd
import argparse
import os
import joblib

# Parse the input and output folder
parser = argparse.ArgumentParser(description='extract trial info')
parser.add_argument('-i', '--input', type=str, 
                    default='/media/lenovo/新加卷/Jacob_BMI/BCI2DCenter20221025_good_ljf/001_lag/20221025_Jacob_interception_semi_80_rifit_1', 
                    metavar='/the/path/your/data/located/in', help='input folder')
parser.add_argument('-o', '--output', type=str, 
                    default='/media/lenovo/新加卷/Jacob_BMI/BCI2DCenter20221025_good_ljf/001_lag/20221025_Jacob_interception_semi_80_rifit_1/1_trial_info', 
                    metavar='/the/path/you/want/to/save', help='output folder')

args = parser.parse_args()
input_folder = args.input
output_folder = args.output

def trial_extract(single_trial):    
    
    # Extract the trial information from the single trial
    trial_d = {}
    marker = [i for i in single_trial if 'Event marker:' in i]
    frame = [i for i in single_trial if 'frame info' in i]
    trial_info = [i for i in single_trial if 'Trial info' in i]
    pos = [i for i in single_trial if 'Publishing: pos' in i]
    state = [i for i in single_trial if 'Publishing: state' in i]
    coefficient = [i for i in single_trial if 'Publishing: coefficient' in i]

    # Extract the marker information
    trial_d['marker'] = \
        [int(re.findall(r".*Event marker\: (.*)", i)[0]) for i in marker]
    
    trial_d['marker_time'] = \
        [float(re.findall(r"\[INFO\] \[(.*)\] \[", i)[0]) for i in marker]
    
    # Extract the frame information
    trial_d['frame'] = [json.loads(re.findall(r".*frame info\: (.*)", i)[0]) \
                            for i in frame]
    trial_d['frame_time'] = \
        [float(re.findall(r"\[INFO\] \[(.*)\] \[", i)[0]) for i in frame]
    
    # Extract the trial information
    trial_info_l = []
    for i in trial_info:
        trial_info_d = {}
        ext = re.findall(r".*Target position\: \[(.*)\], Trial index\: (.*)", i)[0]
        trial_info_d['Target position'] = [float(i) for i in ext[0].split()]
        trial_info_d['Trial index'] = int(ext[1])
        trial_info_l.append(trial_info_d)
    trial_d['trial_info'] = trial_info_l
    trial_d['trial_info_time'] = \
        [float(re.findall(r"\[INFO\] \[(.*)\] \[", i)[0]) for i in trial_info]

    # Extract the pos, state, and coefficient information
    pos_l = []
    for i in pos:
        pos_d = {}
        ext = re.findall(r".*Publishing\: pos\: \[(.*)\], vel\: \[(.*)\], vec_len\: (.*), decoding results\: \[(.*)\], target pos\: \[(.*)\]", i)[0]
        pos_d['pos'] = [float(i) for i in ext[0].split()]
        pos_d['vel'] = [float(i) for i in ext[1].split()]
        pos_d['vec_len'] = float(ext[2])
        pos_d['decoding results'] = [float(i) for i in ext[3].split()]
        pos_d['target pos'] = [float(i) for i in ext[4].split()]
        pos_l.append(pos_d)
        trial_d['pos'] = pos_l
    trial_d['pos_time'] = \
        [float(re.findall(r"\[INFO\] \[(.*)\] \[", i)[0]) for i in pos]

    # Extract the state information
    state_l = []
    for i in state:
        state_d = {}
        ext = re.findall(r".*Publishing\: state\: desired state\: array\('d', \[(.*)\]\), state\: array\('d', \[(.*)\]\)", i)[0]
        state_d['desired state'] = [float(i) for i in ext[0].split(', ')]
        state_d['state'] = [float(i) for i in ext[1].split(', ')]
        state_l.append(state_d)
        trial_d['state'] = state_l
    trial_d['state_time'] = \
        [float(re.findall(r"\[INFO\] \[(.*)\] \[", i)[0]) for i in state]
    
    # Extract the coefficient information
    coefficient_l = []
    for i in coefficient:
        coefficient_d = {}
        ext = re.findall(r".*Publishing\: coefficient\: decoder coefficient\: (.*), assist coefficient\: (.*), PosFlag\: (.*), GenerationSpeed\: (.*), AssistCoefficient1\: (.*), AssistCoefficient2\: (.*), moving_target_flag\: (.*), training flag\: (.*), ring flag\: (.*)", i)[0]
        coefficient_d['decoder coefficient'] = float(ext[0])
        coefficient_d['assist coefficient'] = float(ext[1])
        coefficient_d['PosFlag'] = int(ext[2])
        coefficient_d['GenerationSpeed'] = float(ext[3])
        coefficient_d['AssistCoefficient1'] = float(ext[4])
        coefficient_d['AssistCoefficient2'] = float(ext[5])
        coefficient_d['moving_target_flag'] = int(ext[6])
        coefficient_d['training flag'] = int(ext[7])
        coefficient_d['ring flag'] = int(ext[8])
        coefficient_l.append(coefficient_d)
        trial_d['coefficient'] = coefficient_l
    trial_d['coefficient_time'] = \
        [float(re.findall(r"\[INFO\] \[(.*)\] \[", i)[0]) for i in coefficient]
        
    return trial_d

# Read the behavior log file
with open(os.path.join(input_folder,'behavior.log'),'r') as f:
    bhv = f.readlines()

# Extract the trial
marker_24_pos = [i for i, x in enumerate(bhv) if 'Event marker: 24' in x]
marker_5_pos = [i for i, x in enumerate(bhv) if 'Event marker: 5' in x]

if len(marker_24_pos) != len(marker_5_pos):
    marker_24_pos = marker_24_pos[:-1]
# Check if the number of 24 and 5 markers are equal
assert len(marker_24_pos) == len(marker_5_pos), \
    'Number of 24 and 5 markers are not equal'

# Extract the trials
trial = [bhv[i:j+1] for i, j in zip(marker_24_pos, marker_5_pos)]

# Extract the trial information
trials = [trial_extract(i) for i in trial]

# Save the trial information    
with open(os.path.join(output_folder, 'trials.db'), 'wb') as f:
    joblib.dump(trials, f)
