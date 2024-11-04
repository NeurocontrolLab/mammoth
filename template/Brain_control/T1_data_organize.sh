#!/bin/bash

#SBATCH -J {{subject}}_organizing
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/MAMMOTH/logs/{{subject}}/organizing_job_%j.out 
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p q_fat_2

for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/{{subject}}/Brain_control/*; do
    target_folder="sorted_data"
 
    if find $dir -type d -name $target_folder | read; then
        continue
    fi
  
    keyword="wrong"
    if echo $dir | grep -q $keyword; then
        continue
    fi
    
    keyword='ndt_bmi_log'
    if echo $dir | grep -q $keyword; then
        continue
    fi

    echo $dir

    if [ -d $dir ]; then
        /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_organizer/data_organizer.py -r $dir
    fi
done
