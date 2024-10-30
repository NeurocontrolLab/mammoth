import streamlit as st
import pandas as pd
import os
import subprocess

#%% style and title of the page
st.markdown(':mammoth: MAMMOTH @CuiLab')

st.title("Overview")

#%%
with st.sidebar:
    database = st.selectbox('Choose dataset', ["Neucyber-NC-2024-A-01"], 0, key='s_database')
    if not database:
        st.error("Please select one database.")

#%% dataset setup
# database_dir = '/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01'
database_dir = os.path.join('/AMAX/cuihe_lab/share_rw', database)
code_dir = '/AMAX/cuihe_lab/cuilab_share/MAMMOTH'
curr_overview = [f for f in os.listdir(database_dir) if ('.csv' in f) and ('dataset_overview' in f)][0]   
overview_df = pd.read_csv(os.path.join(database_dir, curr_overview))
overview_df = overview_df.iloc[:, 1:]


#%% overview part
if 'date' not in st.session_state:
    st.session_state['date'] = curr_overview.split('.')[0][-8:]

st.write("#### Last update: %s" % st.session_state['date'])

if st.button("Scan now!", key="scan"):
    scan_code = os.path.join(code_dir, 'scan_dataset.py')
    result = subprocess.run(["/AMAX/cuihe_lab/share_rw/anaconda3/envs/smartneo_env/bin/python", 
                    "%s/scan_dataset.py" % code_dir,
                    "-r", database_dir, 
                    "-o", database_dir], 
                    capture_output=True, text=True)
    
    curr_overview = [f for f in os.listdir(database_dir) if ('.csv' in f) and ('dataset_overview' in f)][0]   
    overview_df = pd.read_csv(os.path.join(database_dir, curr_overview))
    overview_df = overview_df.iloc[:, 1:]

    st.session_state['date'] = curr_overview.split('.')[0][-8:]
    
    st.write(result.stdout)
    if result.stderr:
        st.write(result.stderr)


data = overview_df.copy()

subjects = st.multiselect("Choose subjects", list(pd.unique(data["subject"]))+["All"], ["Abel"])
if not subjects:
    st.error("Please select at least one subject.")
elif (len(subjects)==1) and (subjects[0]=="All"):
    pass
else:
    data = data[data["subject"].isin(subjects)]

types = st.multiselect("Choose types", list(pd.unique(data["type"]))+["All"], ["Data_recording"])
if not types:
    st.error("Please choose at least one type.")
elif (len(types)==1) and (types[0]=="All"):
    pass
else:
    data = data[data["type"].isin(types)]

filters = st.multiselect("Choose filters", ["none", "wrong", "check"], ["none"])
if not filters:
    pass
elif (len(filters)==1) and (types[0]=="none"):
    pass
else:
    for filter in filters:
        data = data[~data["session"].str.contains(filter)]

st.write("", data.sort_index())


