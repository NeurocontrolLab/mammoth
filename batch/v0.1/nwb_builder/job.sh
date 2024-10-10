#!/bin/bash

#SBATCH -J bohrnwb
#SBATCH -o /AMAX/cuihe_lab/cuilab_share/Bohr/Code/nwb_builder/job.%j.out 
#SBATCH -n 1
#SBATCH -p q_fat_2
#SBATCH -c 2

for dir in /AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Bohr/Brain_control/*; do

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
      python /AMAX/cuihe_lab/cuilab_share/Bohr/Code/nwb_builder/data_integration_pipeline/neural_data/br_share/pipeline_no_sort.py -f $dir -o $dir
    fi

    if ! find $dir/formatted_data -name "neural_data.nwb" | read -r; then
      echo "integrate neural data"
      python /AMAX/cuihe_lab/cuilab_share/Bohr/Code/nwb_builder/data_integration_pipeline/neural_data/br_share/pipeline.py -f $dir -o $dir
    fi

    if ! find $dir/formatted_data -name "continuous_behavior.nwb" | read -r; then
      echo "integrate continuous behavior"
      python /AMAX/cuihe_lab/cuilab_share/Bohr/Code/nwb_builder/data_integration_pipeline/continuous_behavior/aie_share/pipeline.py -f $dir -o $dir
    fi

    if ! find $dir/formatted_data -name "trial_behavior.nwb" | read -r; then
      echo "integrate trial behavior"
      bhv_files=$(find $dir -type f -name "*.bhv2" -exec printf "%s\n" {} +)
      if [ -e "$bhv_files" ]; then
         echo $bhv_files
      else
         continue
      fi

      matlab -nodesktop -nosplash -r "data_path='$bhv_files';addpath('/AMAX/cuihe_lab/cuilab_share/Bohr/Code/nwb_builder/data_integration_pipeline/trial_behavior/monkeylogic_share');data2ml;quit"
      python /AMAX/cuihe_lab/cuilab_share/Bohr/Code/nwb_builder/data_integration_pipeline/trial_behavior/monkeylogic_share/pipeline.py -f $dir -o $dir
    fi

done
