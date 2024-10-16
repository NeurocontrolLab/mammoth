#!/bin/bash

#SBATCH -J bohrorganizing
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/MAMMOTH/Bohr/organizing_job_%j.out 
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p q_cn_2

# python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/do_list.py
# python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/sorting_script.py

for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/*; do
    target_folder="sorted_data"
 
    if find $dir -type d -name $target_folder | read; then
        continue
    fi
    keyword="wrong"

    if echo $dir | grep -q $keyword; then
        continue
    fi

    echo $dir

    if [ -d $dir ]; then
        /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_organizer.py -r $dir
    fi
done
