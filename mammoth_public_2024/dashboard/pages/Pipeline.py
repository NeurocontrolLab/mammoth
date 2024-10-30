import streamlit as st
import pandas as pd
import subprocess
import os
from io import StringIO
import datetime

#%% style and title of the page
st.markdown(':mammoth: MAMMOTH @CuiLab')

st.title("Pipeline")


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

mappath_dict = {}
for subject in list(pd.unique(overview_df["subject"])):
    subject_dir = os.path.join(database_dir, subject)
    mapfile = [f for f in os.listdir(subject_dir) if ('.json' in f) and (subject.lower() in f.lower())]
    if len(mapfile)==1:
        mappath_dict[subject] = os.path.join(subject_dir, mapfile[0])

log_path = os.path.join(code_dir, 'log.txt')


#%% pipeline section
st.header("Start new jobs")

step_names = ['Organize', 'Format', 'Check', 'O-F-C', 'Sort', 'Standardize']

tabs = st.tabs(["#%s" % i for i in step_names])

global log
if os.path.exists(log_path):
    log = pd.read_csv(log_path)
else:
    log = pd.DataFrame(columns=['JOBID', 'script_path', 'datetime'])

log.drop(list(log.filter(regex='Unnamed')), axis=1, inplace=True)

for i, sn in enumerate(step_names):
    with tabs[i]:
        col1, col2, col3 = st.columns(3)    
    
        with col1: 
            st.write('Choose subject')
            curr_subject = st.selectbox('Choose subject', list(pd.unique(overview_df["subject"])), 0, key='subject_sb%d' % i, label_visibility="collapsed")
            if not curr_subject:
                st.error("Please select at least one subject.")

        with col2:
            st.write('Choose type')
            curr_type = st.selectbox('Choose type', list(pd.unique(overview_df[overview_df['subject']==curr_subject]['type'])), 0, key='type_sb%d' % i, label_visibility="collapsed")
            if not curr_type:
                st.error("Please choose at least one type.")

        with col3:
            st.write("Click to start **'%s'**." % sn)
            curr_button = st.button("   :point_right: Run :point_left:   ", key='btn%d' % i)
        
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


#%%
st.header("See current jobs")            
if st.button("Show me!", key='show_job_status'):
    job_id = subprocess.run(["job"], capture_output=True, text=True)
    job_info = StringIO(job_id.stdout)
    job_df = pd.read_csv(job_info, delim_whitespace=True)
    job_df.set_index("JOBID", inplace=True)
    job_df.index.name = "JOBID"
    st.table(job_df)


#%%
st.header("Check previous jobs")
previous_jobs = ["job %s for %s with %s at %s" % (i["JOBID"], i["script_path"].split(os.sep)[-2], i["script_path"].split(os.sep)[-1], i["datetime"])
                 for _, i in log.iterrows()]
previous_job = st.selectbox('Select from latest 20 jobs', [" "]+list(reversed(previous_jobs[:20])), 0, key='select_previous_job')
if previous_job and previous_job!=" ":
    with st.expander("Check job output ↓"):
        sdir = os.path.join(code_dir, "logs", previous_job.split("for ")[1].split("_")[0])
        sjob = previous_job.split(" ")[1]
        job_out_file = [f for f in os.listdir(sdir) if ('.out' in f) and (sjob in f)][0]

        try:
            with open(os.path.join(sdir, job_out_file), "r") as f:
                output_content = f.read()     
        except IOError:
            st.error("No such file")
        else:
            st.write(output_content)     
