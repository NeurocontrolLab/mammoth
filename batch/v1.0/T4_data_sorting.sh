#!/bin/bash

#SBATCH -J bohrsorting
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/MAMMOTH/Bohr/sorting_job_%j.out 
#SBATCH -p q_gpu_c
#SBATCH --gres=gpu:1

# target_path="/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Data_recording"

# python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/do_list.py
# python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/sorting_script.py

for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/*; do
    target_folder="sorted_data"
 
    if find $dir -type d -name $target_folder | read; then
        continue
    fi

    keyword='wrong'

    if echo $dir | grep -q $keyword; then
        continue
    fi

    ulimit -n 10240 && /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_sorter_blackrock.py -sorter 'kilosort2_5' -r $dir -mp /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Bohr_Utah_96x2.json -o $dir/sorted_data -cp /AMAX/cuihe_lab/cuilab_share/sorter_container

done
