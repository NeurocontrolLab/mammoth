#!/bin/bash

#SBATCH -J bohrqc
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/Bohr/Code/Data_recording/Quality_control/job.out
#SBATCH -n 1
#SBATCH -p q_fat_2
#SBATCH -c 4

for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Data_recording/*; do

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
        python /AMAX/cuihe_lab/cuilab_share/Bohr/Code/Quality_control/channel_map.py -f $dir -o $dir
    fi
    
    if ! find $dir/description -name "Time_consistency_check.png" | read -r; then
      echo "check time consistency"
        python /AMAX/cuihe_lab/cuilab_share/Bohr/Code/Data_recording/Quality_control/behavior_formatter.py -f $dir -o $dir
        python /AMAX/cuihe_lab/cuilab_share/Bohr/Code/Quality_control/Time_consistency_check.py -f $dir -o $dir
    fi
    
    if ! find $dir/description -name "correlation_test_score.txt" | read -r; then
      if find $dir/formatted_data -name "neural_data_no_sort.nwb" | read -r; then
        echo "neural correlation caculating"
        # 执行 pipeline.py 脚本
        python /AMAX/cuihe_lab/cuilab_share/Bohr/Code/Quality_control/neural_correlation_analysis.py -f $dir -o $dir
      fi
    fi

    if ! find $dir/description -name "chn_consis_summary.png" | read -r; then
      if find $dir/formatted_data -name "neural_data.nwb" | read -r; then
        echo "channel consistency caculation"
        # 执行 pipeline.py 脚本
        python /AMAX/cuihe_lab/cuilab_share/Bohr/Code/Quality_control/channel_consistency.py -f $dir -o $dir
      fi
    fi

done
