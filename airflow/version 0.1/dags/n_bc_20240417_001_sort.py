#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: auto
                
                
import datetime

from airflow import DAG
from airflow.operators.bash_operator import BashOperator


root = format(r'/AMAX/cuihe_lab/chenyun/test')
session = format(r'/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240417_centerOut_001')


default_args = {
    "owner": "cuilab",
    "start_date": datetime.date.today().strftime('%Y-%m-%d')
}

dag = DAG(
        dag_id="n_bc_20240417_001_sort", 
        description="To integrate the data as .nwb files.",
        default_args=default_args, 
        schedule="@once",
        tags=["Nezha", "Brain_control", "sort"])

t1 = BashOperator(task_id='kilosort', 
                  bash_command='''dir='/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240417_centerOut_001'
echo $dir                                                             
if ! find $dir -name "autokilo*.txt" | read -r; then
  echo "sort with kilosort2.5"
  module load singularity
  /AMAX/cuihe_lab/share_rw/anaconda3/envs/lcy_env/bin/python /AMAX/cuihe_lab/share_rw/airpipeline/dags/sort_kilosort_nezha_spikegadgets.py -dp $dir
else
  echo "already sorted with kilosort2.5"
fi''',
                  dag=dag)

t2 = BashOperator(task_id='mountainsort', 
                  bash_command='''dir='/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240417_centerOut_001'
echo $dir                                                             
if ! find $dir -name "automount*.txt" | read -r; then
  echo "sort with mountainsort5"
  module load singularity
  /AMAX/cuihe_lab/share_rw/anaconda3/envs/lcy_env/bin/python /AMAX/cuihe_lab/share_rw/airpipeline/dags/sort_mountsort_nezha_spikegadgets.py -dp $dir
else
  echo "already sorted with mountainsort5
fi''',
                  dag=dag)


t1 >> t2
                