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
      /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_formatter/data_formatter_neural_blackrock.py -r $dir -o $dir/formatted_data -mp /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Bohr/Bohr_Utah_96x2.json -flag 0
      echo "neural_data_no_sort.nwb integrated"
    fi

    if ! find $dir/formatted_data -name "neural_data.nwb" | read -r; then
      echo "integrate neural data"
      
      sorting_files=$(find $dir -type f -name "cluster_info.tsv" -exec printf "%s\n" {}+)
      if [ -n "$sorting_files" ]; then
        echo $sorting_files
        /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_formatter/data_formatter_neural_blackrock.py -r $dir -o $dir/formatted_data -mp /AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01/Abel/Abel_Utah_64x4.json -flag 1
        echo "neural_data.nwb done"
      else
        echo "not yet manually sorted"
      fi
      

    if ! find $dir/formatted_data -name "continuous_behavior.nwb" | read -r; then
      echo "integrate continuous behavior"
      /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_formatter/data_formatter_bhv_continuous.py -r $dir -o $dir/formatted_data
      echo "continuous_behavior.nwb integrated"
    fi

    if ! find $dir/formatted_data -name "trial_behavior.nwb" | read -r; then
      echo "integrate trial behavior"
      bhv_files=$(find $dir -type f -name "*.bhv2" -exec printf "%s\n" {}+)
      if [ -e "$bhv_files" ]; then
         echo $bhv_files
      else
         continue
      fi

      matlab -nodesktop -nosplash -r "data_path='$bhv_files';addpath('/AMAX/cuihe_lab/cuilab_share/mammoth_public_2024/dependencies');data2ml;quit"
      /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/MAMMOTH/mammoth_public_2024/data_formatter/data_formatter_trial_ML.py -r $dir -o $dir/formatted_data
      echo "trial_behavior.nwb integrated"
    fi
    
done