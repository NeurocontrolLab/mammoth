#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 10:23:14 2024

@author: cuilab
"""

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
import csv
import copy
from datetime import datetime


#%%
"""
Created on Wed Aug 28 14:10:11 2024

@author: cuilab
"""
#%%
import pandas as pd
import os
import time
import pprint
import json
import re
import numpy as np
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
from tkcalendar import DateEntry
import csv
import copy
from datetime import datetime


#%%
# overview all files in the root directory
def get_all_files(root_dir):
    data_menu = pd.DataFrame(columns=['path', 'name', 'create_time', 'last_edit_time'])
    for dir_path, dir_names, files in os.walk(root_dir):
        if len(files) > 0:
            for file in files:
                data_menu = data_menu._append({'path': dir_path,
                                               'name': file,
                                               'create_time': time.ctime(os.path.getctime(os.path.join(dir_path, file))),
                                               'last_edit_time': time.ctime(os.path.getmtime(os.path.join(dir_path, file)))},
                                              ignore_index=True)
    return data_menu


# get the hierarchy of directories relative to the root directory
def get_hierarchy(curr_dir, root_dir=None):
    if root_dir is None:
        root_dir = curr_dir
    assert os.path.exists(curr_dir), "This path does not exist."
    v_dm = get_all_files(curr_dir)
    for i in range(len(v_dm)):
        relative_path = root_dir.split(os.sep)[-1] + v_dm.loc[i, 'path'].split(root_dir, 1)[1]
        for hierarchy, folder in enumerate(relative_path.split(os.sep)):
            v_dm.loc[i, 'h%d' % (hierarchy+1)] = folder

    return v_dm


# check folders' names
def check_hierarchy(root_path, layer, formats=None, values=None):
    # formats = r'20[0-9][0-9][0-1][0-9][0-3][0-9]_.+_[0][1-9]'
    v_dm = get_hierarchy(root_path)
    v_dm_check = v_dm[layer].unique()

    assert formats or values, "Need formats or values."

    f_warning_list = []
    if formats:
        for e in v_dm_check:
            if not re.search(formats, e):
                paths = [i.split(e)[0] for i in v_dm.loc[v_dm[layer] == e, 'path']]
                for path in np.unique(paths):
                    f_warning_list.append(path+e)

    if len(f_warning_list) > 0:
        print('Wrong formats.')
        pprint.pprint(f_warning_list)

    v_warning_list = []
    if values:
        for e in v_dm_check:
            if e not in values:
                paths = [i.split(e)[0] for i in v_dm.loc[v_dm[layer] == e, 'path']]
                for path in np.unique(paths):
                    v_warning_list.append(path+e)
    if len(v_warning_list) > 0:
        print('Wrong folder names.')
        pprint.pprint(v_warning_list)

    if (len(f_warning_list) == 0) & (len(v_warning_list) == 0):
        print("Checked! It's OK now!")

    return f_warning_list, v_warning_list


# provide folders with correct names
def view_good_dir(root_path, checklist, layer='h4'):
    # checklist is a list of dict {'layer':, 'format':, 'values':}
    v_dm = get_hierarchy(root_path)

    bad_dirs = []
    for checkpoint in checklist:
        (dir1, dir2) = check_hierarchy(root_path, **checkpoint)
        bad_dirs += dir1 + dir2

    all_dirs = v_dm['path'].unique()
    good_dirs = [i for i in all_dirs if len([j for j in bad_dirs if j in i]) == 0]

    num = 5 - int(layer[-1])
    good_dirs = np.unique([i.rsplit(os.sep, num)[0] for i in good_dirs]).tolist()

    pprint.pprint(good_dirs)

    return good_dirs


# check the tree of given root directory
def view_tree(root_dir, show_files_or_not=0):
    tree_list = []
    v_dm = get_hierarchy(root_dir)
    hs = [i for i in v_dm.columns if ('h' in i) & (str.isdigit(i.split('h')[-1]))]
    # for key in v_dm[hs[0]].unique():
    #     curr = v_dm.loc[v_dm[hs[0]] == key]
    #     print(' ' + key)
    #     for key in curr[hs[1]].unique():
    #         curr = curr.loc[curr[hs[1]] == key]
    #         print('  ' + key)
    #         for key in curr[hs[2]].unique():
    #             curr = curr.loc[curr[hs[2]] == key]
    #             print('   ' + key)
    #             for j in range(len(curr)):
    #                 curr = curr.reset_index(drop=True)
    #                 print('     ' + curr.loc[j, 'name'])

    def print_nested_values(df, layers, level=0):
        if level >= len(layers):
            if show_files_or_not:
                for j in range(len(df)):
                    df = df.reset_index(drop=True)
                    tree_list.append('|' + '-' * level * 4 + df.loc[j, 'name'])
                return
            else:
                return

        h = hs[level]
        for key in [i for i in df[h].unique() if i != 'nan']:
            subset = df[df[h] == key].copy()
            tree_list.append('|' + '-' * level * 4 + key)
            print_nested_values(subset, hs, level + 1)

    print_nested_values(v_dm, hs)

    return tree_list


# check the updates of given directory
def view_update(root_dir, nums=10):
    v_dm = get_hierarchy(root_dir)
    v_dm = v_dm.sort_values(by='last_edit_time', ascending=False)
    v_dm = v_dm.reset_index(drop=True)
    print(v_dm.loc[:nums, [i for i in v_dm.columns if i not in ['path', 'create_time']]])


# check metadata entry
def view_entry(path, print_or_not=1):
    if 'description' in path:
        with open(path + r'\\' + 'metadata.json', encoding='utf-8') as f:
            curr_entry = json.load(f)
    else:
        with open(path + r'\\description\\' + 'metadata.json', encoding='utf-8') as f:
            curr_entry = json.load(f)

    if print_or_not:
        pprint.pprint(curr_entry)

    return curr_entry


# metadata table
def view_entry_table(path):
    v_dm = get_all_files(path)
    v_dm_with_md = v_dm.loc[v_dm['name'] == 'metadata.json', 'path']

    type_list = ['recording', 'training', 'bci']
    entry_type = type_list[[i in path for i in type_list] is True]

    mt_dm = eval('hd.entry_' + entry_type)

    for k, v in v_dm_with_md.items():
        mt_entry = view_entry(format(v.split(r'\description')[0]), print_or_not=0)
        mt_dm = mt_dm._append(mt_entry, ignore_index=True)

    return mt_dm


# query on demand
def query_entry_table(path, print_or_not=1, condition=None, show_columns='all'):

    # condition format: [column, sign, value]
    # condition example: ['name', '==', 'Qianqian']
    # condition example: [['name', '==', 'Qianqian'], ['date', '==', '20231212']]
    mt_dm = view_entry_table(path)

    v_mt_dm = mt_dm.copy()
    if condition is not None:
        assert isinstance(condition, list), \
            "condition should be a list as [column, sign, value], or a list of such lists"
        if isinstance(condition[0], str):
            if condition[-1] in ['{}', 'None']:
                v_mt_dm = eval("mt_dm[mt_dm['%s'] %s %s]" % tuple(condition))
            else:
                v_mt_dm = eval("mt_dm[mt_dm['%s'] %s '%s']" % tuple(condition))
        elif isinstance(condition[0], list):
            v_mt_dm = mt_dm.copy()
            for c in condition:
                assert isinstance(c, list), \
                    "condition should be a list as [column, sign, value], or a list of such lists"
                if c[-1] in ['{}', 'None']:
                    v_mt_dm = eval("v_mt_dm[v_mt_dm['%s'] %s %s]" % tuple(c))
                    # print("v_mt_dm[v_mt_dm['%s'] %s %s]" % tuple(c))
                else:
                    v_mt_dm = eval("v_mt_dm[v_mt_dm['%s'] %s '%s']" % tuple(c))
                    # print("v_mt_dm[v_mt_dm['%s'] %s '%s']" % tuple(c))
                # pprint.pprint(v_mt_dm)

    if show_columns != 'all':
        v_mt_dm = v_mt_dm[show_columns]

    if print_or_not:
        pprint.pprint(v_mt_dm)

    return v_mt_dm


#%%
# main window
root = tk.Tk()
root.title("Metadata")

# set row height
s = ttk.Style()

s.configure("Treeview", rowheight=25)

var = tk.StringVar()

# get directory
def choose_directory():
    global curr_dir, hierarchy_print
    curr_dir = filedialog.askdirectory()
    print(curr_dir)
    if not curr_dir:
        return
    else:
        var.set(curr_dir)
    hierarchy_print = view_tree(curr_dir)
    insert_tree(hierarchy_print)
    

def insert_tree(tree_list):
    for i in range(len(tree_list)):
        tree_listbox.insert(i+1, tree_list[i])


# add element
columns = ["Tree", "View"]
treeview = ttk.Treeview(root, columns=columns, show="headings", height=30)

treeview.column("Tree", width=200)
treeview.column("View", width=800)

treeview.heading("Tree", text="")
treeview.heading("View", text="")


def show_metadata():
    selected_indices = tree_listbox.curselection()
    for id in selected_indices:
        bottom_dir = hierarchy_print[id]
        # print(bottom_dir)
        dir_path = bottom_dir.split('-')[-1]
        for i in np.arange(bottom_dir.count('-')-4, -1, -4):
            cids = np.array([j for j, k in enumerate(hierarchy_print) if k.count('-')==i])
            
            cid1 = cids[np.argwhere(cids<id).flatten()]
            
            cid2 = cid1[np.argmin(abs(cid1-id))]
            # print(cid2)
            dir_path = os.path.join(hierarchy_print[cid2].split('-')[-1], dir_path)

        print('full path = {}'.format(dir_path.split('|')[-1]))



top_frame = ttk.Frame(root)
top_frame.pack(side=tk.TOP, fill=tk.X, pady=10)

load_button = ttk.Button(top_frame, text="Choose directory", command=choose_directory)
load_button.pack(side=tk.LEFT, padx=10)

dir_label = tk.Label(top_frame, textvariable=var, font=('Arial', 8), width=70, height=3)
dir_label.pack(side=tk.LEFT, padx=10)

tree_listbox = tk.Listbox(top_frame, width=90, height=30, selectmode=tk.MULTIPLE)
tree_listbox.pack(side=tk.LEFT, padx=10, pady=30)

view_button = ttk.Button(top_frame, text="View metadata", command=show_metadata)
view_button.pack(side=tk.LEFT, padx=10, pady=40)

root.mainloop()
# %%
