#!/bin/bash

#SBATCH -J bohrnwb
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/Bohr/Code/nwb_builder/job.%j.out 
#SBATCH -n 1
#SBATCH -p q_fat_2
#SBATCH -c 2

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
    if ! find $dir/formatted_data -name "neural_data_no_sort.nwb" | read -r; then
      echo "integrate no sort neural data"
      python /AMAX/cuihe_lab/cuilab_share/mammoth_public_2024/data_formatter_neural_blackrock.py -r $dir -o $dir/formatted_data -mp /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Bohr_Utah_96x2.json -cl ['TCR']
    fi

    if ! find $dir/formatted_data -name "neural_data.nwb" | read -r; then
      echo "integrate neural data"
      
      sorting_files=$(find $dir -type f -name "cluster_info.tsv" -exec printf "%s\n" {} +)
      if [ -e "$sorting_files" ]; then
         echo $sorting_files
      else
         continue
      fi
      
      python /AMAX/cuihe_lab/cuilab_share/mammoth_public_2024/data_formatter_neural_blackrock.py -r $dir -o $dir/formatted_data -mp /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Bohr_Utah_96x2.json -cl ['spike', 'TCR', 'LFP']
    fi

    if ! find $dir/formatted_data -name "continuous_behavior.nwb" | read -r; then
      echo "integrate continuous behavior"
      python /AMAX/cuihe_lab/cuilab_share/mammoth_public_2024/data_formatter_bhv_continuous.py -r $dir -o $dir/formatted_data
    fi

    if ! find $dir/formatted_data -name "trial_behavior.nwb" | read -r; then
      echo "integrate trial behavior"
      bhv_files=$(find $dir -type f -name "*.bhv2" -exec printf "%s\n" {} +)
      if [ -e "$bhv_files" ]; then
         echo $bhv_files
      else
         continue
      fi

      matlab -nodesktop -nosplash -r "data_path='$bhv_files';addpath('/AMAX/cuihe_lab/cuilab_share/mammoth_public_2024');data2ml;quit"
      python /AMAX/cuihe_lab/cuilab_share/mammoth_public_2024/data_formatter_trial_ML.py -r $dir -o $dir/formatted_data
    fi

done