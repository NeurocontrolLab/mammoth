#!/bin/bash

#SBATCH -J nezhaorg
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/Nezha/Code/data_organizer/job_orgnizer.out
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p q_cn

# python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/do_list.py
# python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/sorting_script.py

for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Brain_control/*; do
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
        python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/data_organizer/data_organization.py -f $dir
    fi
done
