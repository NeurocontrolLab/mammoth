#!/bin/bash

#SBATCH -J {{subject}}_formatting
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/MAMMOTH/logs/{{subject}}/formatting_job_%j.out 
# SBATCH -p q_gpu_c
# SBATCH --gres=gpu:1

#SBATCH -n 1
#SBATCH -p q_fat_2
#SBATCH -c 2

for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/{{subject}}/Data_recording/*; do

    keyword='wrong'

    if echo $dir | grep -q $keyword; then
        continue
    fi

    keyword='check'

    if echo $dir | grep -q $keyword; then
        continue
    fi

    echo $dir
    if ! find $dir/formatted_data -name "neural_data_no_sort.nwb" | read -r; then
      echo "integrate no sort neural data"
      /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_formatter/data_formatter_neural_blackrock.py -r $dir -o $dir/formatted_data -mp {{mappath}} -flag 0
    fi   

    if ! find $dir/formatted_data -name "continuous_behavior.nwb" | read -r; then
      echo "integrate continuous behavior"
      /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_formatter/data_formatter_bhv_continuous.py -r $dir -o $dir/formatted_data
    fi

done
