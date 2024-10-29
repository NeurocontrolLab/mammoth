#!/bin/bash

#SBATCH -J {{subject}}sorting
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/MAMMOTH/{{subject}}/sorting_job_%j.out 
#SBATCH -p q_gpu_c
#SBATCH --gres=gpu:1


for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/{{subject}}/Data_recording/*; do
    target_folder="sorted_data"
 
    if find $dir -type d -name $target_folder | read; then
        continue
    fi

    keyword='wrong'

    if echo $dir | grep -q $keyword; then
        continue
    fi

    ulimit -n 10240 && /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_sorter_blackrock.py -sorter 'kilosort2_5' -r $dir -mp {{mappath}} -o $dir/sorted_data -cp /AMAX/cuihe_lab/cuilab_share/sorter_container

done
