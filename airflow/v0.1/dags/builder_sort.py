#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:47:24 2024

@author: chenyun
"""

import os
from string import Template


class BashTemplate(Template):
    delimiter = '&'


kilosort_code_dict = {}
kilosort_code_dict['nezha'] = BashTemplate('''dir='&{session0}'
echo $dir                                                             
if ! find $dir -name "autokilo*.txt" | read -r; then
  echo "sort with kilosort2.5"
  module load singularity
  /AMAX/cuihe_lab/share_rw/anaconda3/envs/lcy_env/bin/python /AMAX/cuihe_lab/share_rw/airpipeline/dags/sort_kilosort_nezha_&{system0}.py -dp $dir
else
  echo "already sorted with kilosort2.5"
fi''')


spykingcircus_code_dict = {}
spykingcircus_code_dict['nezha'] = BashTemplate('''dir='&{session0}'
echo $dir                                                             
if ! find $dir -name "autoskc*.txt" | read -r; then
  echo "sort with spykingcircus"
  module load singularity
  /AMAX/cuihe_lab/share_rw/anaconda3/envs/lcy_env/bin/python /AMAX/cuihe_lab/share_rw/airpipeline/dags/sort_spykingcircus_nezha_&{system0}.py -dp $dir
else
  echo "already sorted with spykingcircus"
fi''')

vote_template = BashTemplate('''dir='&{session0}'
echo $dir                                                             
if ! find $dir -name "autokilo_voted*.txt" | read -r; then
  echo "combine the sorting results"
  /AMAX/cuihe_lab/share_rw/anaconda3/envs/lcy_env/bin/python /AMAX/cuihe_lab/share_rw/airpipeline/dags/sort_voter.py -dp $dir
else
  echo "already voted"
fi''')


code = Template("""#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: auto
                
                
import datetime

from airflow import DAG
from airflow.operators.bash_operator import BashOperator


root = format(r'${root0}')
session = format(r'${session0}')


default_args = {
    "owner": "cuilab",
    "start_date": datetime.date.today().strftime('%Y-%m-%d')
}

dag = DAG(
        dag_id="${idname0}", 
        description="To integrate the data as .nwb files.",
        default_args=default_args, 
        schedule="@once",
        tags=["${monkey0}", "${paradigm0}", "sort"])

t1 = BashOperator(task_id='kilosort', 
                  bash_command='''${kilosort_script0}''',
                  dag=dag)

t2 = BashOperator(task_id='spykingcircus_sort', 
                  bash_command='''${spykingcircus_sort_script0}''',
                  dag=dag)

t3 = BashOperator(task_id='vote',
                  bash_command='''${vote_script0}''',
                  dag=dag)

t1 >> t2 >> t3
                """)
                
                
def gen_sort_script(root, session):
    
    #root = format(r'/AMAX/cuihe_lab/chenyun/test')
    #session = format(r'/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001')
    
    monkey = session.split(root)[1].split(os.path.sep)[1]
    paradigm = session.split(root)[1].split(os.path.sep)[2]
    edate = session.split(root)[1].split(os.path.sep)[3].split('_')[0]
    eid = session.split(root)[1].split(os.path.sep)[3].split('_')[-1]
    idname = '%s_%s_%s_%s_sort' % (monkey[0].lower(), 
                               ''.join([i[0] for i in paradigm.lower().split('_')]),
                               edate,eid)
 
    # gen_auto_bash = sh_dict[monkey.lower()].substitute(session0=session)
    # with open('/AMAX/cuihe_lab/share_rw/airpipeline/dags/%s.sh' % idname, 'w') as file:
    #     file.write(gen_auto_bash)
    
    kilosort_script = kilosort_code_dict[monkey.lower()].substitute(
        session0=session, system0='spikegadgets')
    skc_script = spykingcircus_code_dict[monkey.lower()].substitute(
        session0=session, system0='spikegadgets')
    vote_script = vote_template.substitute(session0=session)
    gen_auto_py = code.substitute(root0=root, session0=session, idname0=idname,
                                  monkey0=monkey, paradigm0=paradigm,
                                  kilosort_script0=kilosort_script,
                                  spykingcircus_sort_script0=skc_script,
                                  vote_script0=vote_script)
    
    # print(gen_auto)
    
    with open('/AMAX/cuihe_lab/share_rw/airpipeline/dags/%s.py' % idname, 'w') as file:
        file.write(gen_auto_py)
        
