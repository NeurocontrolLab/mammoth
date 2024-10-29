import streamlit as st
import pandas as pd
import subprocess
import os
from io import StringIO

st.markdown('MAMMOTH @CuiLab')

st.title("Neucyber-NC-2024-A-01")

database_dir = '/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01'
code_dir = '/AMAX/cuihe_lab/cuilab_share/MAMMOTH'
curr_overview = [f for f in os.listdir(database_dir) if ('.csv' in f) and ('dataset_overview' in f)][0]   
overview_df = pd.read_csv(os.path.join(database_dir, curr_overview))
overview_df = overview_df.iloc[:, 1:]

#%% Run
st.header("Run the pipeline")

snote = ['Choose subjects', '', '', '', '']
tnote = ['Choose type', '', '', '', '']
button_names = ['Organize', 'Format', 'Check', 'Sort', 'Standardize']
buttons = []
mappath_dict = {}
for subject in list(pd.unique(overview_df["subject"])):
    subject_dir = os.path.join(database_dir, subject)
    mapfile = [f for f in os.listdir(subject_dir) if ('.json' in f) and (subject.lower() in f.lower())]
    if len(mapfile)==1:
        mappath_dict[subject] = os.path.join(subject_dir, mapfile[0])

col1, col2, col3 = st.columns(3)    
i = 0
with col1: 
    curr_subject = st.selectbox(snote[i], list(pd.unique(overview_df["subject"])), 0, key='subject_ms%d' % i)
    if not curr_subject:
        st.error("Please select at least one subject.")

with col2:
    curr_type = st.selectbox(tnote[i], list(pd.unique(overview_df["type"])), 0, key='type_ms%d' % i)
    if not curr_type:
        st.error("Please choose at least one type.")

with col3:
    st.empty()
    curr_button = st.button(button_names[i], key='btn%d' % i)
    button_names.append(curr_button)

if curr_button:
    st.write('Organize: %s, %s' % (curr_subject, curr_type))

    folder_name = f"scripts/{curr_subject}_{curr_type}"
    if not os.path.exists(os.path.join(code_dir, folder_name)):
        os.makedirs(os.path.join(code_dir, folder_name))

    script_file_path = os.path.join(code_dir, folder_name, '%s.sh' % button_names[i])

    if os.path.exists(script_file_path):
        st.write("Find the existed script, running...")
    else:
        st.write("Generating a new script...")
        
        template_dir = os.path.join(code_dir, 'template')
        template_file = [f for f in os.listdir(template_dir) if ('.' in f) and (str(i+1) in f)][0]
        with open(os.path.join(template_dir, template_file), "r") as template:
            script_content = template.read()
            script_content = script_content.replace("{{subject}}", curr_subject)
            script_content = script_content.replace("{{mappath}}", mappath_dict[curr_subject])

        with open(script_file_path, "w") as new:
            new.write(script_content)
            os.chmod(script_file_path, 0o755)

    result = subprocess.run([script_file_path], capture_output=True, text=True)

    job_id = subprocess.run(["job"], capture_output=True, text=True)
    job_info = StringIO(job_id.stdout)
    job_df = pd.read_csv(job_info, delim_whitespace=True)
    job_df.set_index("JOBID")
    st.table(job_df)

    st.write("Output: ")
    st.write(result.stdout)

    if result.stderr:
        st.write("Error: ")
        st.write(result.stderr)



# if buttons[1]:
#     st.write('Organize for Abel and Data!')


# st.button("Just go till Check!")

#%%
uploaded_file = st.file_uploader("Upload your shell script", type="sh")

if uploaded_file is not None:
    with open("uploaded_script.sh", "wb") as f:

        f.write(uploaded_file.getbuffer())

        if st.button("Run"):
            result = subprocess.run(["/bin/bash", "uploaded_script.sh"],
                                    capture_output=True,
                                    text=True)
            st.text("Output:")
            st.text(result.stdout)
            st.text("Errors:")
            st.text(result.stderr)


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



