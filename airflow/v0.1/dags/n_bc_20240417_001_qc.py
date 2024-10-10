#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: auto
                
import os                
import datetime
import data_qcer as qcer

from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.operators.bash import BashOperator


root = format(r'/AMAX/cuihe_lab/chenyun/test')
session = format(r'/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240417_centerOut_001')


default_args = {
    "owner": "cuilab",
    "start_date": datetime.date.today().strftime('%Y-%m-%d')
}

dag = DAG(
        dag_id="n_bc_20240417_001_qc", 
        description="To generate quality control report.",
        default_args=default_args, 
        schedule="@once",
        tags=["Nezha", "Brain_control", "qc"],
        is_paused_upon_creation=False)

t1 = PythonOperator(task_id='gen_qc_code', 
                    python_callable=qcer.check_file, 
                    op_kwargs={'session': session},
                    dag=dag)

t2 = BashOperator(task_id='run_qc_code',
                  bash_command='/AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/jupyter nbconvert %s --to notebook --execute --inplace --allow-errors'
                  % os.path.join(
                      os.path.join(session, 'description'),
                      'quality_control.ipynb'))

t1 >> t2
                