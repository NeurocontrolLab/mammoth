#!/bin/bash

#SBATCH -J bohrchecking
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/MAMMOTH/Bohr/checking_job_%j.out 
# SBATCH -p q_gpu_c
# SBATCH --gres=gpu:1

#SBATCH -n 1
#SBATCH -p q_fat_2
#SBATCH -c 4

for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Data_recording/*; do

    keyword='wrong'

    if echo $dir | grep -q $keyword; then
        continue
    fi

    keyword='check'

    if echo $dir | grep -q $keyword; then
        continue
    fi
    
    echo $dir
    if ! find $dir/description -name "channel_map.png" | read -r; then
      echo "plot channel map"
        /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_checker/data_checker_channel_map_plotting.py -mp /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Bohr_Utah_96x2.json -o $dir/description
    fi
    
    if ! find $dir/description -name "Time_consistency_check.png" | read -r; then
      echo "check time consistency"
        /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_checker/data_checker_time_consistency.py -d $dir/formatted_data -o $dir/description
    fi
    
    if ! find $dir/description -name "correlation_test_score.txt" | read -r; then
      if find $dir/formatted_data -name "neural_data_no_sort.nwb" | read -r; then
        echo "neural correlation caculating"
        # 执行 pipeline.py 脚本
        /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_checker/data_checker_neural_correlation.py -d $dir/formatted_data -o $dir/description -s $dir/description
      fi
    fi

    if ! find $dir/description -name "chn_consis_summary.png" | read -r; then
      if find $dir/formatted_data -name "neural_data.nwb" | read -r; then
        echo "channel consistency caculation"
        # 执行 pipeline.py 脚本
        /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_checker/data_checker_channel_consistency.py -r $dir -d $dir/formatted_data -o $dir/description
      fi
    fi

done
