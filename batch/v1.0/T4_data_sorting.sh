#!/bin/bash

#SBATCH -J bohrsort
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/Bohr/Code/sorting_controller/job_bc_convert.out
#SBATCH -p q_gpu_c
#SBATCH --gres=gpu:1

# target_path="/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Data_recording"

# python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/do_list.py
# python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/sorting_controller/sorting_script.py

for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Brain_control/*; do
    target_folder="sorted_data"
 
    if find $dir -type d -name $target_folder | read; then
        continue
    fi

    keyword='wrong'

    if echo $dir | grep -q $keyword; then
        continue
    fi

    ulimit -n 10240 && python /AMAX/cuihe_lab/cuilab_share/Bohr/Code/sorting_controller/sorting_script_kilo_2_5.py -dp $dir

done
