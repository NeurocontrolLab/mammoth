import streamlit as st
import pandas as pd
import subprocess
import os
from io import StringIO
import datetime

#%% dataset setup
database_dir = '/AMAX/cuihe_lab/share_rw/Neucyber-NC-2024-A-01'
code_dir = '/AMAX/cuihe_lab/cuilab_share/MAMMOTH'
curr_overview = [f for f in os.listdir(database_dir) if ('.csv' in f) and ('dataset_overview' in f)][0]   
overview_df = pd.read_csv(os.path.join(database_dir, curr_overview))
overview_df = overview_df.iloc[:, 1:]

mappath_dict = {}
for subject in list(pd.unique(overview_df["subject"])):
    subject_dir = os.path.join(database_dir, subject)
    mapfile = [f for f in os.listdir(subject_dir) if ('.json' in f) and (subject.lower() in f.lower())]
    if len(mapfile)==1:
        mappath_dict[subject] = os.path.join(subject_dir, mapfile[0])

log_path = os.path.join(code_dir, 'log.txt')

#%% style and title of the page
st.markdown('MAMMOTH @CuiLab')

st.title("Neucyber-NC-2024-A-01")


#%% pipeline section
st.header("Pipeline")

snote = ['Choose subjects', '', '', '', '', '']
tnote = ['Choose type', '', '', '', '', '']
step_names = ['Organize', 'Format', 'Check', 'O-F-C', 'Sort', 'Standardize']

tabs = st.tabs(step_names)

global log
if os.path.exists(log_path):
    log = pd.read_csv(log_path)
else:
    log = pd.DataFrame(columns=['JOBID', 'script_path', 'datetime'])

for i, sn in enumerate(step_names):
    with tabs[i]:
        col1, col2, col3 = st.columns(3)    
    
        with col1: 
            curr_subject = st.selectbox('Choose subjects', list(pd.unique(overview_df["subject"])), 0, key='subject_sb%d' % i)
            if not curr_subject:
                st.error("Please select at least one subject.")

        with col2:
            curr_type = st.selectbox('Choose type', list(pd.unique(overview_df[overview_df['subject']==curr_subject]['type'])), 0, key='type_sb%d' % i)
            if not curr_type:
                st.error("Please choose at least one type.")

        with col3:
            st.text("Click to start '%s'." % sn)
            curr_button = st.button("Run", key='btn%d' % i)
        
        if curr_button:
            st.write('%s the data for: %s, %s' % (sn, curr_subject, curr_type))

            folder_name = f"scripts/{curr_subject}_{curr_type}"
            if not os.path.exists(os.path.join(code_dir, folder_name)):
                os.makedirs(os.path.join(code_dir, folder_name))

            script_file_path = os.path.join(code_dir, folder_name, '%s.sh' % sn)

            if os.path.exists(script_file_path):
                st.write("Find the existed script, running...")
            else:
                st.write("Generating a new script...")
                
                template_dir = os.path.join(code_dir, 'template')
                template_file = [f for f in os.listdir(template_dir) if ('.' in f) and (sn.lower() in f.lower())][0]
                with open(os.path.join(template_dir, template_file), "r") as template:
                    script_content = template.read()
                    script_content = script_content.replace("{{subject}}", curr_subject)
                    script_content = script_content.replace("{{mappath}}", mappath_dict[curr_subject])

                with open(script_file_path, "w") as new:
                    new.write(script_content)
                    os.chmod(script_file_path, 0o755)

            result = subprocess.run(["sbatch", script_file_path], capture_output=True, text=True)
            curr_job = result.stdout.split("job ")[-1]
            log.loc[len(log)] = [curr_job, script_file_path, datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')]
            log.to_csv(log_path)

            
if st.button("See current jobs", key='show_job_status'):
    job_id = subprocess.run(["job"], capture_output=True, text=True)
    job_info = StringIO(job_id.stdout)
    job_df = pd.read_csv(job_info, delim_whitespace=True)
    job_df.set_index("JOBID", inplace=True)
    job_df.index.name = "JOBID"
    st.table(job_df)


previous_jobs = ["job %s for %s with %s at %s" % (i["JOBID"], i["script_path"].split(os.sep)[-2], i["script_path"].split(os.sep)[-1], i["datetime"])
                 for _, i in log.iterrows()]
previous_job = st.selectbox('See previous jobs', previous_jobs, 0, key='select_previous_job')
if previous_job:
    with st.expander("See job output ↓"):
        sdir = os.path.join(code_dir, previous_job.split("for ")[1].split("_")[0])
        sjob = previous_job.split(" ")[1]
        job_out_file = [f for f in os.listdir(sdir) if ('.out' in f) and (sjob in f)][0]

        with open(os.path.join(sdir, job_out_file), "r") as f:
            output_content = f.read()
            st.write(output_content)          


# if buttons[1]:
#     st.write('Organize for Abel and Data!')


# st.button("Just go till Check!")
