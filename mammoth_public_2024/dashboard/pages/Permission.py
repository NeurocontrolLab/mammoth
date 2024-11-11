import streamlit as st
import pandas as pd
import os
import subprocess

#%% style and title of the page
st.markdown(':mammoth: MAMMOTH @CuiLab')

st.title("Permission")

#%%
with st.sidebar:
    database = st.selectbox('Choose dataset', ["Neucyber-NC-2024-A-01"], 0, key='s_database')
    if not database:
        st.error("Please select one database.")

#%% dataset setup
# database_dir = '/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01'
database_dir = os.path.join('/AMAX/cuihe_lab/share_rw', database)


#%% permission part
target_dir = st.text_input('Please enter the directory:', database_dir)

user = ''
users = st.selectbox('Choose users', ['All', 'Someone'], 0, key='select_users')
if users == 'Someone':
    user = st.text_input('Please enter the userID:')
elif users == 'All':
    user = users

if st.button("Permit now!", key="permit"):
    if user == 'All':
        result = subprocess.run(["find", database_dir, "-user", os.getlogin(), "-exec", "chmod", "770", "{}", "+"], 
                        capture_output=True, text=True)
        st.write(result.stdout)

    else:
        result = subprocess.run(["find", database_dir, "-user", os.getlogin(), "-exec", "setfacl", "-m", "u:%s:rwx" % user, "{}", "+"],
                        capture_output=True, text=True)
        st.write(result.stdout)

    if result.stderr:
        st.write(result.stderr)
    





