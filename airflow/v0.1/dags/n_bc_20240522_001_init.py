#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: airflow
                

import os                
import datetime
import data_organizer as organizer
import data_checker as checker

from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.operators.bash import BashOperator


root = format(r'/AMAX/cuihe_lab/chenyun/test')
session = format(r'/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001')


default_args = {
    "owner": "cuilab",
    "start_date": datetime.date.today().strftime('%Y-%m-%d')
}

dag = DAG(
        dag_id="n_bc_20240522_001_init", 
        description="To organize the data and verify.",
        default_args=default_args, 
        schedule="@once",
        tags=["Nezha", "Brain_control", "init"],
        is_paused_upon_creation=False)

t1 = PythonOperator(task_id='organize_the_files', 
                    python_callable=organizer.organize_file, 
                    op_kwargs={'root_dir': session},
                    dag=dag)

t2 = BashOperator(task_id='convert_rec_file', 
                  bash_command='''dir='/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001/raw_data'
if ! find $dir -type d -name "*kilosort" | read; then
  echo "convert .rec files"
  rec_files=$(find $dir -type f -name "*.rec" -exec printf "%s
" {} +)
  echo $rec_files
  ulimit -n 10240 && module load trodes/2.4.2 && trodesexport -rec $rec_files -lfp -lfplowpass 300 -dio -kilosort -spikeband -spikehighpass 300 -spikelowpass 6000 -spikes -thresh 50
else
  echo "already converted"
fi''',
                  dag=dag)

t3 = BashOperator(task_id='format_neural_data_no_sort',
                  bash_command='''dir='/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001'                                                   
echo $dir
if ! find $dir -name "neural_data_no_sort.nwb" | read -r; then
  echo "integrate neural data"
  /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/nwb_builder/data_integration_pipeline/neural_data/br_share/pipeline_no_sort.py -f $dir -o $dir
else
  echo "neural data no sort already integrated"
fi''',
                  dag=dag)

t4 = BashOperator(task_id='format_continuous_behavior',
                  bash_command='''dir='/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001'
echo $dir                               
if ! find $dir -name "continuous_behavior.nwb" | read -r; then
  echo "integrate continuous behavior"
  /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/nwb_builder/data_integration_pipeline/continuous_behavior/aie_share/pipeline.py -f $dir -o $dir
else
  echo "continuous behavior already integrated"
fi''',
                  dag=dag)

t5 = PythonOperator(task_id='gen_check_code', 
                    python_callable=checker.check_file, 
                    op_kwargs={'session': session},
                    dag=dag)

t6 = BashOperator(task_id='run_check_code',
                  bash_command='/AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/jupyter nbconvert %s --to notebook --execute --inplace --allow-errors'
                  % os.path.join(
                      os.path.join(session, 'description'),
                      'quality_control_no_sort.ipynb'))

t1 >> t2 >> t3 >> t4 >> t5 >> t6
