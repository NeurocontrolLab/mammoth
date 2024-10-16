#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 19:34:42 2024

@author: cuilab
"""

import os
from string import Template


#%% 
format_neural_data_no_sort_code = {}
format_continuous_behavior = {}

class BashTemplate(Template):
    delimiter = '&'

format_neural_data_no_sort_code['nezha'] = BashTemplate('''dir='&{session0}'                                                   
echo $dir
if ! find $dir -name "neural_data_no_sort.nwb" | read -r; then
  echo "integrate neural data"
  /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/nwb_builder/data_integration_pipeline/neural_data/br_share/pipeline_no_sort.py -f $dir -o $dir
else
  echo "neural data no sort already integrated"
fi''')

format_continuous_behavior['nezha'] = BashTemplate('''dir='&{session0}'
echo $dir                               
if ! find $dir -name "continuous_behavior.nwb" | read -r; then
  echo "integrate continuous behavior"
  /AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python /AMAX/cuihe_lab/cuilab_share/Nezha/Code/nwb_builder/data_integration_pipeline/continuous_behavior/aie_share/pipeline.py -f $dir -o $dir
else
  echo "continuous behavior already integrated"
fi''')

class BashTemplate2(Template):
    delimiter = '#'
convert_rec_template = BashTemplate2('''dir='#{raw0}'
if ! find $dir -type d -name "*kilosort" | read; then
  echo "convert .rec files"
  rec_files=$(find $dir -type f -name "*.rec" -exec printf "%s\n" {} +)
  echo $rec_files
  ulimit -n 10240 && module load trodes/2.4.2 && trodesexport -rec $rec_files -lfp -lfplowpass 300 -dio -kilosort -spikeband -spikehighpass 300 -spikelowpass 6000 -spikes -thresh 50
else
  echo "already converted"
fi''')

code = Template("""#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: airflow
                

import os                
import datetime
import data_organizer as organizer
import data_checker as checker

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
        description="To organize the data and verify.",
        default_args=default_args, 
        schedule="@once",
        tags=["${monkey0}", "${paradigm0}", "init"],
        is_paused_upon_creation=False)

t1 = PythonOperator(task_id='organize_the_files', 
                    python_callable=organizer.organize_file, 
                    op_kwargs={'root_dir': session},
                    dag=dag)

t2 = BashOperator(task_id='convert_rec_file', 
                  bash_command='''${convert_rec_code0}''',
                  dag=dag)

t3 = BashOperator(task_id='format_neural_data_no_sort',
                  bash_command='''${no_sort_format0}''',
                  dag=dag)

t4 = BashOperator(task_id='format_continuous_behavior',
                  bash_command='''${con_bhv_format0}''',
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
""")
                
                
def gen_init_script(root, session):
    
    # root = format(r'/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01')
    # session = format(r'/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Brain_control/20240522_centerOutKalmanNet_threshold70_001')
    
    monkey = session.split(root)[1].split(os.path.sep)[1]
    paradigm = session.split(root)[1].split(os.path.sep)[2]
    edate = session.split(root)[1].split(os.path.sep)[3].split('_')[0]
    eid = session.split(root)[1].split(os.path.sep)[3].split('_')[-1]
    idname = '%s_%s_%s_%s_init' % (monkey[0].lower(), 
                               ''.join([i[0] for i in paradigm.lower().split('_')]),
                               edate,eid)
    no_sort_format = format_neural_data_no_sort_code[monkey.lower()].substitute(
        session0=session)
 
    con_bhv_format = format_continuous_behavior[monkey.lower()].substitute(
        session0=session)
    convert_rec_code = convert_rec_template.substitute(
        raw0=os.path.join(session, 'raw_data'))
    gen_auto = code.substitute(root0=root, session0=session, idname0=idname,
                         monkey0=monkey, paradigm0=paradigm,
                         convert_rec_code0=convert_rec_code,
                         no_sort_format0=no_sort_format,
                         con_bhv_format0=con_bhv_format)
    
    # print(gen_auto)
    
    with open('/AMAX/cuihe_lab/share_rw/airpipeline/dags/%s.py' % idname, 'w') as file:
        file.write(gen_auto)