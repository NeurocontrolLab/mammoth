#!/bin/bash

#SBATCH -J {{subject}}_ofc
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/MAMMOTH/logs/{{subject}}/ofc_job_%j.out 
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p q_fat_2

for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/{{subject}}/Brain_control/*; do
    
    if [! -d $dir]; then
      continue
    fi  

    keyword='wrong'

    if echo $dir | grep -q $keyword; then
        continue
    fi

    keyword='check'

    if echo $dir | grep -q $keyword; then
        continue
    fi
    
    keyword='ndt_bmi_log'
    if echo $dir | grep -q $keyword; then
        continue
    fi
    
    echo $dir

    target_folder="sorted_data"
 
    if ! find $dir -type d -name $target_folder | read; then
        echo "organize"
        /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_organizer/data_organizer.py -r $dir    
    fi

    if ! find $dir/formatted_data -name "neural_data_no_sort.nwb" | read -r; then
      echo "integrate no sort neural data"
      /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_formatter/data_formatter_neural_blackrock.py -r $dir -o $dir/formatted_data -mp {{mappath}} -flag 0
    fi   

    if ! find $dir/formatted_data -name "continuous_behavior.nwb" | read -r; then
      echo "integrate continuous behavior"
      /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_formatter/data_formatter_bhv_continuous_BC.py -r $dir -o $dir/formatted_data
    fi


    if ! find $dir/description -name "channel_map.png" | read -r; then
      echo "plot channel map"
        /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_checker/data_checker_channel_map_plotting.py -mp {{mappah}} -o $dir/description
    fi
    
    if ! find $dir/description -name "Time_consistency_check.png" | read -r; then
      echo "check time consistency"
        /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_checker/data_checker_time_consistency.py -d $dir/formatted_data -o $dir/description
    fi
    
    if ! find $dir/description -name "correlation_test_score.txt" | read -r; then
      if find $dir/formatted_data -name "neural_data_no_sort.nwb" | read -r; then
        echo "neural correlation caculating"
        /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_checker/data_checker_neural_correlation.py -d $dir/formatted_data -o $dir/description -s $dir/description
      fi
    fi

    echo "O-F-C completed"

done
