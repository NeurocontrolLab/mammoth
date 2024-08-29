#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 14:54:50 2024

@author: chenyun
"""

from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from datetime import datetime


default_args = {
    "owner": "test",
    "start_date": datetime(2024, 5, 28)
}
dag = DAG("Hello-World2", 
        description="test DAG",
        default_args=default_args, 
        schedule_interval='0 8 * * *')

t0 = BashOperator(task_id="task0", 
                  bash_command="python /AMAX/cuihe_lab/share_rw/airpipeline/dags/trytime.py ", 
                  dag=dag)

t0
