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

nd_dict = {}
tb_dict = {}

nd_dict['nezha'] = BashTemplate('''dir='&{session0}'
echo $dir
if ! find $dir -name "neural_data.nwb" | read -r; then
  echo "integrate neural data"
  /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/nwb_builder/data_integration_pipeline/neural_data/br_share/pipeline.py -f $dir -o $dir
else
  echo "neural data already integrated"
fi''')


tb_dict['nezha'] = BashTemplate('''dir='&{session0}'
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
''')


code = Template("""#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: auto
                
                
import datetime
import data_organizer as organizer

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
        tags=["${monkey0}", "${paradigm0}", "format"],
        is_paused_upon_creation=False)

t1 = BashOperator(task_id='integrate_neural_data', 
                  bash_command='''${nd_script0}''',
                  dag=dag)

t2 = BashOperator(task_id='integrate_trial_behavior', 
                  bash_command='''${tb_script0}''',
                  dag=dag)


t1 >> t2
                """)
                
                
def gen_format_script(root, session):
    
    #root = format(r'/AMAX/cuihe_lab/chenyun/test')
    #session = format(r'/AMAX/cuihe_lab/chenyun/test/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001')
    
    monkey = session.split(root)[1].split(os.path.sep)[1]
    paradigm = session.split(root)[1].split(os.path.sep)[2]
    edate = session.split(root)[1].split(os.path.sep)[3].split('_')[0]
    eid = session.split(root)[1].split(os.path.sep)[3].split('_')[-1]
    idname = '%s_%s_%s_%s_format' % (monkey[0].lower(), 
                               ''.join([i[0] for i in paradigm.lower().split('_')]),
                               edate,eid)
 
    # gen_auto_bash = sh_dict[monkey.lower()].substitute(session0=session)
    # with open('/AMAX/cuihe_lab/share_rw/airpipeline/dags/%s.sh' % idname, 'w') as file:
    #     file.write(gen_auto_bash)
    
    nd_script = nd_dict[monkey.lower()].substitute(session0=session)
    tb_script = tb_dict[monkey.lower()].substitute(session0=session)
    gen_auto_py = code.substitute(root0=root, session0=session, idname0=idname,
                                  monkey0=monkey, paradigm0=paradigm,
                                  nd_script0=nd_script,
                                  tb_script0=tb_script)
    
    # print(gen_auto)
    
    with open('/AMAX/cuihe_lab/share_rw/airpipeline/dags/%s.py' % idname, 'w') as file:
        file.write(gen_auto_py)
        
