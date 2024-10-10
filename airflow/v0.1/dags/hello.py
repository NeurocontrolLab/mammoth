#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 14:11:37 2024

@author: chenyun
"""

from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from datetime import datetime


default_args = {
    "owner": "test",
    "start_date": datetime(2024, 5, 28)
}

dag = DAG("Hello-World", 
        description="test DAG",
        default_args=default_args, 
        schedule_interval='0 8 * * *')

t1 = BashOperator(task_id="hello", 
                  bash_command="echo 'Hello World, today is {{ ds }}'", 
                  dag=dag)
