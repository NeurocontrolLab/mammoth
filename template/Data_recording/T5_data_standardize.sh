#!/bin/bash

#SBATCH -J {{subject}}_structuring
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/MAMMOTH/logs/{{subject}}/structuring_job_%j.out 
# SBATCH -p q_gpu_c
# SBATCH --gres=gpu:1

#SBATCH -n 1
#SBATCH -p q_fat_2
#SBATCH -c 4


for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/{{subject}}/Data_recording/*; do
    keyword='wrong'

    if echo $dir | grep -q $keyword; then
        continue
    fi
    
    keyword='check'

    if echo $dir | grep -q $keyword; then
        continue
    fi

    target_folder="formatted_data"
 
    if ! find $dir -type d -name $target_folder | read; then
        continue
    fi
    
    target_folder="sorted_data"
 
    if ! find $dir -type d -name $target_folder | read; then
        continue
    fi
    
    echo $dir
    if ! find $dir/formatted_data -name "neural_data.nwb" | read -r; then
      echo "integrate neural data"

      s=1      
      for sdir in $dir/sorted_data/kilosort2_5_output/*; do
      
        if ! find $sdir -name "cluster_info.tsv" | read -r; then
            s=0
            break
        #    echo $sorting_files
        fi
      done
      
      if s==1; then
          /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_formatter/data_formatter_neural_blackrock.py -r $dir -o $dir/formatted_data -mp {{mappath}} -flag 1
   
      else
          echo "manual sorting not yet finished"
          continue
      fi
    fi

    
    if find $dir/formatted_data -name "neural_data.nwb" | read -r; then
      if ! find $dir/description -name "chn_consis_summary.png" | read -r; then
          echo "channel consistency caculation"
          /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_checker/data_checker_channel_consistency.py -r $dir -d $dir/formatted_data -o $dir/description
      fi
  
      if find $dir/formatted_data -name "standard_data.nwb" | read -r; then
          echo "standard_data.nwb already exists"
          continue 
      fi
  
      ulimit -n 10240 && /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_structurer/data_structurer_standardNWB.py -r $dir -mp {{mappath}} -o $dir/formatted_data
      echo "standard_data.nwb integrated"
    fi
done
