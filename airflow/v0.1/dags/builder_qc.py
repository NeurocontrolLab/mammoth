#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 19:09:58 2024

@author: chenyun
"""

import os
from string import Template

code = Template('''#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: auto
                
import os                
import datetime
import data_qcer as qcer

from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.operators.bash import BashOperator


root = format(r'${root0}')
session = format(r'${session0}')


default_args = {
    "owner": "cuilab",
    "start_date": datetime.date.today().strftime('%Y-%m-%d')
}

dag = DAG(
        dag_id="${idname0}", 
        description="To generate quality control report.",
        default_args=default_args, 
        schedule="@once",
        tags=["${monkey0}", "${paradigm0}", "qc"],
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
                ''')
                
                
def gen_qc_script(root, session):
    
    # root = format(r'/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01')
    # session = format(r'/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001')
    
    monkey = session.split(root)[1].split(os.path.sep)[1]
    paradigm = session.split(root)[1].split(os.path.sep)[2]
    edate = session.split(root)[1].split(os.path.sep)[3].split('_')[0]
    eid = session.split(root)[1].split(os.path.sep)[3].split('_')[-1]
    idname = '%s_%s_%s_%s_qc' % (monkey[0].lower(), 
                               ''.join([i[0] for i in paradigm.lower().split('_')]),
                               edate,eid)
 
    gen_auto = code.substitute(root0=root, session0=session, idname0=idname,
                               monkey0=monkey, paradigm0=paradigm)
    
    # print(gen_auto)
    
    with open('/AMAX/cuihe_lab/share_rw/airpipeline/dags/%s.py' % idname, 'w') as file:
        file.write(gen_auto)
