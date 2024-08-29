#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 3 12:49:11 2024

@author: yunchen
"""

import pandas as pd
import os
import time
import builder_init
import builder_format
import builder_qc
import builder_sort
from data_organizer import target_folder_name_list 

from airflow import DAG
from airflow.operators.python import PythonOperator
from datetime import datetime


#%% setup path
# root0 = format(r'/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01')
# metapath0 = format(r'/AMAX/cuihe_lab/share_rw/metadata')
root0 = format(r'/AMAX/cuihe_lab/chenyun/test')
metapath0 = format(r'/AMAX/cuihe_lab/chenyun/test')


#%% define functions
def get_all_sessions(root_path, meta_path):
    data_menu = pd.DataFrame(columns=['name', 'type', 'date', 'path'])
    for dir_path, dir_names, files in os.walk(root_path, topdown=True):
        for dir_name_i in dir_names:
            dir_path_i = os.path.join(dir_path, dir_name_i)
            relative_path = dir_path_i.split(root_path)[1].split(os.sep)
            relative_depth = len(relative_path)-1
            if relative_depth == 3:          
                data_menu = data_menu._append({'name': relative_path[1],
                                               'type': relative_path[2],
                                               'date': relative_path[3],
                                               'path': dir_path_i},
                                               ignore_index=True)
            else:
                break                                 
    
    data_menu.to_csv(os.path.join(meta_path, '%s_all_sessions_list.csv'
                     % time.strftime("%Y%m%d", time.localtime(time.time()))))
    print('Already got all session folders.')


def get_newest_file(path, name):
    lists = os.listdir(path)
    lists = [i for i in lists 
             if (name in i) and (i[-4:]=='.csv')]
    lists.sort()
    return os.path.join(path, lists[-1])
    

def get_and_scan_fit_sessions(meta_path):
    f_session_df = pd.read_csv(get_newest_file(meta_path,'all_sessions'))
    
    f_session_df = f_session_df.drop(
        f_session_df[f_session_df['type']=='backup'].index)
    f_session_df = f_session_df.drop(
        f_session_df[f_session_df['type']=='Monkey_training'].index)
    f_session_df = f_session_df.drop(
        [i for i in f_session_df.index 
         if 'wrong' in f_session_df.loc[i, 'date'].lower()])
    
    f_session_df = f_session_df.reset_index(drop=False)
    
    for i in f_session_df.index:
        spath = f_session_df.loc[i, 'path']
        
        dir_list = next(os.walk(spath))[1]
        
        # organized
        f_session_df.loc[i, 'organized'] = int(
            len([dd for dd in dir_list if dd not in target_folder_name_list])==0)
       
        # checked
        f_session_df.loc[i, 'checked'] = (
             0 if not os.path.exists(os.path.join(spath, 'description')) 
             else int('quality_control_no_sort.ipynb' in 
                      os.listdir(os.path.join(spath, 'description'))))
        
        # sorted
        if not os.path.exists(os.path.join(spath, 'sorted_data')):
            f_session_df.loc[i, 'sorted'] = '0'
            continue
        
        kspath = os.path.join(spath, 'sorted_data/kilosort2_5_output') 
        marktxt = next(os.walk(kspath))[2]
       
        if 'verified' in marktxt:
            f_session_df.loc[i, 'sorted'] = marktxt[0].split('.')[0]
            continue
        if 'wrong' in marktxt:
            f_session_df.loc[i, 'sorted'] = marktxt[0].split('.')[0]
            continue
        
        wrong_flag = len([i for i in next(os.walk(kspath))[1] 
                          if 'wrong' in i])
        if wrong_flag!=0:
            if len(marktxt)>0:
                os.rename(os.path.join(kspath, marktxt), 
                          os.path.join(kspath, 'wrong_shanks_%d.txt' % wrong_flag))
            else:
                with open(
                        os.path.join(kspath, 
                                     'wrong_shanks_%d.txt' % wrong_flag), 
                        'w') as file:
                    file.write("")  
            continue
       
        shank_to_verify = []
        for shank in sorted(next(os.walk(kspath))[1]):
            if 'cluster_info.tsv' not in os.listdir(
                    os.path.join(kspath, shank+'/sorter_output')):
                shank_to_verify.append(shank[-1])
            else:
                ks_info = pd.read_csv(os.path.join(
                    kspath, shank+'/sorter_output/cluster_info.tsv'),
                    sep='\t')
                groupnum = len([j for j in ks_info['group'] 
                                if type(j)==str])
                labelnum = len(ks_info['KSLabel'])
                grouprate = groupnum/labelnum
                if grouprate < 0.5:
                    shank_to_verify.append(shank[-1])
            
        if len(shank_to_verify)>0:
            add = '_toverify%s' % (''.join(shank_to_verify))
        else:
            add = '_verified'
        
        if len(marktxt)>0:
            if 'verif' in marktxt[0]:
                dev = 'toverif' if '_toverif' in marktxt[0] else 'verif' 
                newtxtname = ''.join([marktxt[0].split(dev)[0], add, '.txt'])
                os.rename(os.path.join(kspath, marktxt[0]), 
                          os.path.join(kspath, newtxtname))
            
            else:
                newtxtname = ''.join([marktxt[0].split('.')[0], add, '.txt'])
                os.rename(os.path.join(kspath, marktxt[0]), 
                          os.path.join(kspath, newtxtname))
        else:
            with open(os.path.join(kspath, 'autokilo%s.txt'%add), 'w') as file:
                file.write("")
        
        marktxtnew = next(os.walk(kspath))[2]
        f_session_df.loc[i, 'sorted'] = marktxtnew[0].split('.')[0]
    
                    
        # formatted
        nwbcount = ([f for f in os.listdir(os.path.join(spath, 'formatted_data'))
              if f[-4:]=='.nwb'] 
                    if os.path.exists(os.path.join(spath, 'formatted_data'))
                    else [])
        f_session_df.loc[i, 'formatted'] = ('0' if len(nwbcount)==0 else ','.join(nwbcount))
        
        
        # quality control
        f_session_df.loc[i, 'quality_control'] = (
             0 if not os.path.exists(os.path.join(spath, 'description')) 
             else int('quality_control.ipynb' in 
                      os.listdir(os.path.join(spath, 'description'))))
            
    f_session_df.to_csv(os.path.join(meta_path, '%s_scan_fit_sessions.csv'
                     % time.strftime("%Y%m%d", time.localtime(time.time()))))
    print('Already scan fit sessions.')    


def gen_init_script(root_path, meta_path):
    f_session_df = pd.read_csv(get_newest_file(meta_path, 'scan_fit'))
    for i in f_session_df[(f_session_df['organized']==0)
                          |(f_session_df['checked']==0)].index:
        print('init: %s' % f_session_df.loc[i, 'path'])
        builder_init.gen_init_script(root_path, f_session_df.loc[i, 'path'])


def gen_sort_script(root_path, meta_path):
    f_session_df = pd.read_csv(get_newest_file(meta_path, 'scan_fit'))
    for i in f_session_df[(f_session_df['organized']==1)
                          &(f_session_df['checked']==1)].index:
        if 'verified' not in f_session_df.loc[i, 'sorted']:
            if 'voted' not in f_session_df.loc[i, 'sorted']:
                print('sort: %s' % f_session_df.loc[i, 'path'])
                builder_sort.gen_sort_script(root_path, f_session_df.loc[i, 'path'])


def gen_format_script(root_path, meta_path):
    f_session_df = pd.read_csv(get_newest_file(meta_path, 'scan_fit'))
    for i in f_session_df[
            (f_session_df['organized']==1)
            &(f_session_df['checked']==1)].index:
        if 'verified' in f_session_df.loc[i, 'sorted']:
            print('format: %s' % f_session_df.loc[i, 'path'])
            builder_format.gen_format_script(root_path, f_session_df.loc[i, 'path'])


def gen_qc_script(root_path, meta_path):
    f_session_df = pd.read_csv(get_newest_file(meta_path, 'scan_fit'))
    for i in f_session_df[(f_session_df['checked']==1)
                          &(f_session_df['quality_control']==0)].index:
        if 'neural_data.nwb' in f_session_df.loc[i, 'formatted']:
            print('quality_control: %s' % f_session_df.loc[i, 'path'])
            builder_qc.gen_qc_script(root_path, f_session_df.loc[i, 'path'])

#%% define DAG
default_args = {
    "owner": "cuilab",
    "start_date": datetime(2024, 5, 29)
}


dag = DAG(
        dag_id="Dataset_scanning", 
        description="To scan the database",
        default_args=default_args, 
        # schedule="daily",
        schedule_interval='0 8 * * *')


t1 = PythonOperator(task_id='get_all_sessions', 
                    python_callable=get_all_sessions, 
                    op_kwargs={'root_path': root0,
                             'meta_path': metapath0},
                    dag=dag)
                    
t2 = PythonOperator(task_id='get_and_scan_fit_sessions', 
                    python_callable=get_and_scan_fit_sessions, 
                    op_kwargs={'root_path': root0,
                             'meta_path': metapath0},
                    dag=dag)

t_init = PythonOperator(task_id='batch_init_sessions', 
                    python_callable=gen_init_script, 
                    op_kwargs={'root_path': root0,
                               'meta_path': metapath0},
                    dag=dag)

t_sort = PythonOperator(task_id='batch_sort_sessions',
                    python_callable=gen_sort_script,
                    op_kwargs={'root_path': root0,
                               'meta_path': metapath0},
                    dag=dag)

t_format = PythonOperator(task_id='batch_format_sessions',
                    python_callable=gen_format_script,
                    op_kwargs={'root_path': root0,
                               'meta_path': metapath0},
                    dag=dag)

t_qc = PythonOperator(task_id='batch_qc_sessions',
                    python_callable=gen_qc_script,
                    op_kwargs={'root_path': root0,
                               'meta_path': metapath0},
                    dag=dag)


t1 >> t2 >> t_init >> t_sort >> t_format >> t_qc
