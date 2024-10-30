import streamlit as st
import pandas as pd
import os


#%% dataset setup
database_dir = '/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01'
code_dir = '/AMAX/cuihe_lab/cuilab_share/MAMMOTH'
curr_overview = [f for f in os.listdir(database_dir) if ('.csv' in f) and ('dataset_overview' in f)][0]   
overview_df = pd.read_csv(os.path.join(database_dir, curr_overview))
overview_df = overview_df.iloc[:, 1:]


#%% style and title of the page
st.markdown('MAMMOTH @CuiLab')

st.title("Neucyber-NC-2024-A-01")


#%% overview part
st.header("Overview")

data = overview_df.copy()

subjects = st.multiselect("Choose subjects", list(pd.unique(data["subject"])), ["Abel"])
if not subjects:
    st.error("Please select at least one subject.")
else:
    data = data[data["subject"].isin(subjects)]

types = st.multiselect("Choose type", list(pd.unique(data["type"])), ["Data_recording"])
if not types:
    st.error("Please choose at least one type.")
else:
    data = data[data["type"].isin(types)]

st.write("Overview", data.sort_index())
