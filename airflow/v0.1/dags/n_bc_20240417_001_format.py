#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: auto
                
                
import datetime
import data_organizer as organizer

from airflow import DAG
from airflow.operators.bash_operator import BashOperator


root = format(r'/AMAX/cuihe_lab/chenyun/test')
session = format(r'/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240417_centerOut_001')


default_args = {
    "owner": "cuilab",
    "start_date": datetime.date.today().strftime('%Y-%m-%d')
}

dag = DAG(
        dag_id="n_bc_20240417_001_format", 
        description="To integrate the data as .nwb files.",
        default_args=default_args, 
        schedule="@once",
        tags=["Nezha", "Brain_control", "format"],
        is_paused_upon_creation=False)

t1 = BashOperator(task_id='integrate_neural_data', 
                  bash_command='''dir='/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240417_centerOut_001'
echo $dir
if ! find $dir -name "neural_data.nwb" | read -r; then
  echo "integrate neural data"
  /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/nwb_builder/data_integration_pipeline/neural_data/br_share/pipeline.py -f $dir -o $dir
else
  echo "neural data already integrated"
fi''',
                  dag=dag)

t2 = BashOperator(task_id='integrate_trial_behavior', 
                  bash_command='''dir='/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240417_centerOut_001'
echo $dir                                                             
if ! find $dir -name "trial_behavior.nwb" | read -r; then
  echo "integrate trial behavior"
  bhv_files=$(find $dir -type f -name "*.bhv2" -exec printf "%s
" {} +)
  if [ -e "$bhv_files" ]; then
     echo $bhv_files
     
     matlab -nodesktop -nosplash -r "data_path='$bhv_files';addpath('/AMAX/cuihe_lab/cuilab_share/Nezha/Code/nwb_builder/data_integration_pipeline/trial_behavior/monkeylogic_share/');data2ml;quit"
  /AMAX/cuihe_lab/share_rw_anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/nwb_builder/data_integration_pipeline/trial_behavior/monkeylogic_share/pipeline.py -f $dir -o $dir
  
  else
     echo "no bhv files thus no trial behavior"
  fi
else 
  echo "trial behavior already integrated"
fi
''',
                  dag=dag)


t1 >> t2
                