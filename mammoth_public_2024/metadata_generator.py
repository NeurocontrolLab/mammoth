<<<<<<< HEAD
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 10:23:14 2024

@author: cuilab
"""

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
from tkcalendar import DateEntry
import csv
import copy
from datetime import datetime

metadata = {}
initial_metadata = None

# main window
root = tk.Tk()
root.title("Metadata")

# add element
columns = ["Key", "Value"]
treeview = ttk.Treeview(root, columns=columns, show="headings", height=30)

treeview.column("Key", width=200)
treeview.column("Value", width=800)

treeview.heading("Key", text="Key")
treeview.heading("Value", text="Value")

# set row height
s = ttk.Style()
s.configure("Treeview", rowheight=25)

# read csv
def load_metadata_from_csv():
    global metadata, initial_metadata
    file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
    if not file_path:
        return

    metadata = {}
    for item in treeview.get_children():
        treeview.delete(item)

    with open(file_path, mode='r', newline='') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            if len(row) != 2:
                continue
            key, value = row
            update_metadata_from_key_value(key, value)

    initial_metadata = copy.deepcopy(metadata)
    insert_metadata()

# update dictionary
def update_metadata_from_key_value(key, value):
    keys = key.split('/')
    ref = metadata
    for i, k in enumerate(keys[:-1]):
        try:
            k = int(k)
            if isinstance(ref, list):
                while len(ref) <= k:
                    ref.append({})
                ref = ref[k]
            else:
                raise ValueError(f"Expected list but got {type(ref).__name__}")
        except ValueError:
            if isinstance(ref, dict):
                if k not in ref:
                    if i < len(keys) - 1 and keys[i + 1].isdigit():
                        ref[k] = []
                    else:
                        ref[k] = {}
                ref = ref[k]
            else:
                raise ValueError(f"Expected dict but got {type(ref).__name__}")

    last_key = keys[-1]
    try:
        last_key = int(last_key)
        if isinstance(ref, list):
            while len(ref) <= last_key:
                ref.append(None)
            ref[last_key] = value
        else:
            raise ValueError(f"Expected list for index but got {type(ref).__name__}")
    except ValueError:
        ref[last_key] = value


# insert metadata
def insert_metadata():
    for item in treeview.get_children():
        treeview.delete(item)

    for main_key, main_value in metadata.items():
        insert_data(main_key, main_value)

# interation
def insert_data(parent_key, data):
    if isinstance(data, dict):
        for key, value in data.items():
            full_key = f"{parent_key}/{key}"
            insert_data(full_key, value)
    elif isinstance(data, list):
        for idx, item in enumerate(data):
            full_key = f"{parent_key}/{idx}"
            if isinstance(item, dict):
                insert_data(full_key, item)
            else:
                treeview.insert("", "end", values=(full_key, item))
    else:
        treeview.insert("", "end", values=(parent_key, data))

# edit
def set_metadata_value(event):
    selected_item = treeview.selection()[0]
    item_values = treeview.item(selected_item, "values")
    full_key = item_values[0]
    current_value = item_values[1]

    bbox = treeview.bbox(selected_item, column="#2")
    if not bbox:
        return

    # datetime editor
    def saveedit(new_value):
        treeview.set(selected_item, column="#2", value=new_value)
        update_metadata_from_key_value(full_key, new_value)
        entryedit.destroy()

    if "session_start_time" in full_key:
        date_win = tk.Toplevel(root)
        date_win.geometry(f"+{bbox[0]}+{bbox[1]}")

        # date option
        entryedit = DateEntry(date_win, date_pattern='yyyy-mm-dd', background='darkblue', foreground='white', borderwidth=2)
        entryedit.grid(row=0, column=0)
        entryedit.set_date(datetime.strptime(current_value, '%Y-%m-%d %H:%M:%S').date())

        # time option
        time_frame = ttk.Frame(date_win)
        time_frame.grid(row=0, column=1, padx=10)
        hour_var = tk.StringVar(value=current_value.split()[1].split(':')[0])
        minute_var = tk.StringVar(value=current_value.split()[1].split(':')[1])
        second_var = tk.StringVar(value=current_value.split()[1].split(':')[2])
        ttk.Label(time_frame, text="Hour:").grid(row=0, column=0)
        hour_entry = ttk.Entry(time_frame, width=3, textvariable=hour_var)
        hour_entry.grid(row=0, column=1)
        ttk.Label(time_frame, text="Minute:").grid(row=1, column=0)
        minute_entry = ttk.Entry(time_frame, width=3, textvariable=minute_var)
        minute_entry.grid(row=1, column=1)
        ttk.Label(time_frame, text="Second:").grid(row=2, column=0)
        second_entry = ttk.Entry(time_frame, width=3, textvariable=second_var)
        second_entry.grid(row=2, column=1)

        def on_date_time_select():
            selected_date = entryedit.get_date()
            selected_hour = int(hour_var.get())
            selected_minute = int(minute_var.get())
            selected_second = int(second_var.get())
            selected_datetime = datetime(selected_date.year, selected_date.month, selected_date.day, selected_hour, selected_minute, selected_second)
            formatted_datetime = selected_datetime.strftime('%Y-%m-%d %H:%M:%S')
            saveedit(formatted_datetime)
            date_win.destroy()

        confirm_btn = ttk.Button(date_win, text="OK", command=on_date_time_select)
        confirm_btn.grid(row=3, column=0, columnspan=2, pady=10)

    else:
        # txt edition
        entryedit = ttk.Entry(root)
        entryedit.insert(0, current_value)
        entryedit.place(x=bbox[0], y=bbox[1], width=bbox[2], height=bbox[3])
        entryedit.focus_set()
        entryedit.bind("<Return>", lambda e: saveedit(entryedit.get()))
        entryedit.bind("<FocusOut>", lambda e: saveedit(entryedit.get()))

treeview.bind("<Double-1>", set_metadata_value)

# reset metadata
def reset_metadata():
    if messagebox.askokcancel("Reset Data", " Do you want to save data"):
        global metadata
        metadata = copy.deepcopy(initial_metadata)
        insert_metadata()

# clear data
def clear_data():
    if messagebox.askokcancel("Clear Data", " Do you want to clear data"):
        global metadata
        metadata = {}
        for item in treeview.get_children():
            treeview.delete(item)

# save metadata as CSV 
def save_metadata_to_csv():
    if messagebox.askokcancel("Save Data", " Do you want to save data"):
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if not file_path:
            return
        try:
            with open(file_path, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(["Key", "Value"])
                for main_key, main_value in metadata.items():
                    for key, value in flatten_metadata(main_key, main_value):
                        writer.writerow([key, value])
            messagebox.showinfo("Save Data", "Data saved successfully")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save csv:{e}")

def flatten_metadata(parent_key, data):
    flat_data = []
    if isinstance(data, dict):
        for key, value in data.items():
            full_key = f"{parent_key}/{key}"
            flat_data.extend(flatten_metadata(full_key, value))
    elif isinstance(data, list):
        for idx, item in enumerate(data):
            if isinstance(item, dict):
                for key, value in item.items():
                    full_key = f"{parent_key}/{idx}/{key}"
                    flat_data.extend(flatten_metadata(full_key,value))
            else:
                full_key = f"{parent_key}/{idx}"
                flat_data.extend(flatten_metadata(full_key, item))
    else:
        flat_data.append((parent_key, data))
    return flat_data

# add buttons
bottom_frame = ttk.Frame(root)
bottom_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=10)

load_button = ttk.Button(bottom_frame, text="Load CSV", command=load_metadata_from_csv)
load_button.pack(side=tk.LEFT, padx=10)

reset_button = ttk.Button(bottom_frame, text="Reset", command=reset_metadata)
reset_button.pack(side=tk.LEFT, padx=10)

clear_button = ttk.Button(bottom_frame, text="Clear Data", command=clear_data)
clear_button.pack(side=tk.LEFT, padx=10)

save_button = ttk.Button(bottom_frame, text="Save to CSV", command=save_metadata_to_csv)
save_button.pack(side=tk.LEFT, padx=10)

# final setting
treeview.pack(fill="both", expand=True)
root.mainloop()
=======
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 10:23:14 2024

@author: cuilab
"""

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
from tkcalendar import DateEntry
import csv
import copy
from datetime import datetime

metadata = {}
initial_metadata = None

# main window
root = tk.Tk()
root.title("Metadata")

# add element
columns = ["Key", "Value"]
treeview = ttk.Treeview(root, columns=columns, show="headings", height=30)

treeview.column("Key", width=200)
treeview.column("Value", width=800)

treeview.heading("Key", text="Key")
treeview.heading("Value", text="Value")

# set row height
s = ttk.Style()
s.configure("Treeview", rowheight=25)

# read csv
def load_metadata_from_csv():
    global metadata, initial_metadata
    file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
    if not file_path:
        return

    metadata = {}
    for item in treeview.get_children():
        treeview.delete(item)

    with open(file_path, mode='r', newline='') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            if len(row) != 2:
                continue
            key, value = row
            update_metadata_from_key_value(key, value)

    initial_metadata = copy.deepcopy(metadata)
    insert_metadata()

# update dictionary
def update_metadata_from_key_value(key, value):
    keys = key.split('/')
    ref = metadata
    for i, k in enumerate(keys[:-1]):
        try:
            k = int(k)
            if isinstance(ref, list):
                while len(ref) <= k:
                    ref.append({})
                ref = ref[k]
            else:
                raise ValueError(f"Expected list but got {type(ref).__name__}")
        except ValueError:
            if isinstance(ref, dict):
                if k not in ref:
                    if i < len(keys) - 1 and keys[i + 1].isdigit():
                        ref[k] = []
                    else:
                        ref[k] = {}
                ref = ref[k]
            else:
                raise ValueError(f"Expected dict but got {type(ref).__name__}")

    last_key = keys[-1]
    try:
        last_key = int(last_key)
        if isinstance(ref, list):
            while len(ref) <= last_key:
                ref.append(None)
            ref[last_key] = value
        else:
            raise ValueError(f"Expected list for index but got {type(ref).__name__}")
    except ValueError:
        ref[last_key] = value


# insert metadata
def insert_metadata():
    for item in treeview.get_children():
        treeview.delete(item)

    for main_key, main_value in metadata.items():
        insert_data(main_key, main_value)

# interation
def insert_data(parent_key, data):
    if isinstance(data, dict):
        for key, value in data.items():
            full_key = f"{parent_key}/{key}"
            insert_data(full_key, value)
    elif isinstance(data, list):
        for idx, item in enumerate(data):
            full_key = f"{parent_key}/{idx}"
            if isinstance(item, dict):
                insert_data(full_key, item)
            else:
                treeview.insert("", "end", values=(full_key, item))
    else:
        treeview.insert("", "end", values=(parent_key, data))

# edit
def set_metadata_value(event):
    selected_item = treeview.selection()[0]
    item_values = treeview.item(selected_item, "values")
    full_key = item_values[0]
    current_value = item_values[1]

    bbox = treeview.bbox(selected_item, column="#2")
    if not bbox:
        return

    # datetime editor
    def saveedit(new_value):
        treeview.set(selected_item, column="#2", value=new_value)
        update_metadata_from_key_value(full_key, new_value)
        entryedit.destroy()

    if "session_start_time" in full_key:
        date_win = tk.Toplevel(root)
        date_win.geometry(f"+{bbox[0]}+{bbox[1]}")

        # date option
        entryedit = DateEntry(date_win, date_pattern='yyyy-mm-dd', background='darkblue', foreground='white', borderwidth=2)
        entryedit.grid(row=0, column=0)
        entryedit.set_date(datetime.strptime(current_value, '%Y-%m-%d %H:%M:%S').date())

        # time option
        time_frame = ttk.Frame(date_win)
        time_frame.grid(row=0, column=1, padx=10)
        hour_var = tk.StringVar(value=current_value.split()[1].split(':')[0])
        minute_var = tk.StringVar(value=current_value.split()[1].split(':')[1])
        second_var = tk.StringVar(value=current_value.split()[1].split(':')[2])
        ttk.Label(time_frame, text="Hour:").grid(row=0, column=0)
        hour_entry = ttk.Entry(time_frame, width=3, textvariable=hour_var)
        hour_entry.grid(row=0, column=1)
        ttk.Label(time_frame, text="Minute:").grid(row=1, column=0)
        minute_entry = ttk.Entry(time_frame, width=3, textvariable=minute_var)
        minute_entry.grid(row=1, column=1)
        ttk.Label(time_frame, text="Second:").grid(row=2, column=0)
        second_entry = ttk.Entry(time_frame, width=3, textvariable=second_var)
        second_entry.grid(row=2, column=1)

        def on_date_time_select():
            selected_date = entryedit.get_date()
            selected_hour = int(hour_var.get())
            selected_minute = int(minute_var.get())
            selected_second = int(second_var.get())
            selected_datetime = datetime(selected_date.year, selected_date.month, selected_date.day, selected_hour, selected_minute, selected_second)
            formatted_datetime = selected_datetime.strftime('%Y-%m-%d %H:%M:%S')
            saveedit(formatted_datetime)
            date_win.destroy()

        confirm_btn = ttk.Button(date_win, text="OK", command=on_date_time_select)
        confirm_btn.grid(row=3, column=0, columnspan=2, pady=10)

    else:
        # txt edition
        entryedit = ttk.Entry(root)
        entryedit.insert(0, current_value)
        entryedit.place(x=bbox[0], y=bbox[1], width=bbox[2], height=bbox[3])
        entryedit.focus_set()
        entryedit.bind("<Return>", lambda e: saveedit(entryedit.get()))
        entryedit.bind("<FocusOut>", lambda e: saveedit(entryedit.get()))

treeview.bind("<Double-1>", set_metadata_value)

# reset metadata
def reset_metadata():
    if messagebox.askokcancel("Reset Data", " Do you want to save data"):
        global metadata
        metadata = copy.deepcopy(initial_metadata)
        insert_metadata()

# clear data
def clear_data():
    if messagebox.askokcancel("Clear Data", " Do you want to clear data"):
        global metadata
        metadata = {}
        for item in treeview.get_children():
            treeview.delete(item)

# save metadata as CSV 
def save_metadata_to_csv():
    if messagebox.askokcancel("Save Data", " Do you want to save data"):
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if not file_path:
            return
        try:
            with open(file_path, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(["Key", "Value"])
                for main_key, main_value in metadata.items():
                    for key, value in flatten_metadata(main_key, main_value):
                        writer.writerow([key, value])
            messagebox.showinfo("Save Data", "Data saved successfully")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save csv:{e}")

def flatten_metadata(parent_key, data):
    flat_data = []
    if isinstance(data, dict):
        for key, value in data.items():
            full_key = f"{parent_key}/{key}"
            flat_data.extend(flatten_metadata(full_key, value))
    elif isinstance(data, list):
        for idx, item in enumerate(data):
            if isinstance(item, dict):
                for key, value in item.items():
                    full_key = f"{parent_key}/{idx}/{key}"
                    flat_data.extend(flatten_metadata(full_key,value))
            else:
                full_key = f"{parent_key}/{idx}"
                flat_data.extend(flatten_metadata(full_key, item))
    else:
        flat_data.append((parent_key, data))
    return flat_data

# add buttons
bottom_frame = ttk.Frame(root)
bottom_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=10)

load_button = ttk.Button(bottom_frame, text="Load CSV", command=load_metadata_from_csv)
load_button.pack(side=tk.LEFT, padx=10)

reset_button = ttk.Button(bottom_frame, text="Reset", command=reset_metadata)
reset_button.pack(side=tk.LEFT, padx=10)

clear_button = ttk.Button(bottom_frame, text="Clear Data", command=clear_data)
clear_button.pack(side=tk.LEFT, padx=10)

save_button = ttk.Button(bottom_frame, text="Save to CSV", command=save_metadata_to_csv)
save_button.pack(side=tk.LEFT, padx=10)

# final setting
treeview.pack(fill="both", expand=True)
root.mainloop()
>>>>>>> 30ed201402aebf165a59b0f25615ba2e47694b03
