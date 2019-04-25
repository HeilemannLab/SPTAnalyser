# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 13:23:16 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Jupyter Notebook Widgets for Merge HDF5 (trc->hdf5 and hmm->archive).
"""

import tkinter as tk 
from tkinter.filedialog import askopenfilename
from ipywidgets import widgets
from IPython.display import clear_output


class WidgetMergeHdf5():
    def __init__(self):
        self.load_hdf5_button = self.create_load_hdf5_button()
        self.load_hdf5_box = self.create_load_hdf5_box()
        self.hdf5_file_name = ""
        self.load_archive_button = self.create_load_archive_button()
        self.load_archive_box = self.create_load_archive_box()
        self.archive_file_name = ""
        self.merge_button = self.create_merge_button()
        self.merged_name_box = self.create_merged_name_box()
        
    def create_load_hdf5_button(self):
        """
        Button to load the file.
        """
        button = widgets.Button(
                description='browse',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='browse for file')
                #icon='check')
        return button
    
    def create_load_hdf5_box(self, val = "", desc = "Insert path"):
        """
        Box for inserting the file path of the fitting hdf5 file for the HMM analysis (analysis was done with .trc file).
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def open_hdf5_file(self, b):  # b = ???
        """
        Give the file button opening powers.
        """
        root = tk.Tk()  # window class
        root.withdraw()  # close the window 
        root.update()  # close the window
        root.name = askopenfilename(title="Import hdf5 file", filetypes=(("text files", "*.h5"),("all files", "*.*")))
        self.hdf5_file_name = root.name
        root.update()
        root.destroy()
        self.load_hdf5_box.value = self.hdf5_file_name
        
    def change_hdf5_box(self, change):
        self.hdf5_file_name = self.load_hdf5_box.value  
    
    def create_load_archive_button(self):
        """
        Button to load the file.
        """
        button = widgets.Button(
                description='browse',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='browse for file')
                #icon='check')
        return button
    
    def create_load_archive_box(self, val = "", desc = "Insert path"):
        """
        Box for inserting the file path of the archive file created for each cell in the HMM analysis.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_merged_name_box(self, val = "merged", desc = "Name"):
        """
        Box for inserting the file path of the archive file created for each cell in the HMM analysis.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def open_archive_file(self, b):
        """
        Give the file button opening powers.
        """
        root = tk.Tk()  # window class
        root.withdraw()  # close the window 
        root.update()  # close the window
        root.name = askopenfilename(title="Import archive file", filetypes=(("text files", "*.h5"),("all files", "*.*")))
        self.archive_file_name = root.name
        root.update()
        root.destroy()
        self.load_archive_box.value = self.archive_file_name
        
    def change_archive_box(self, change):
        self.archive_file_name = self.load_archive_box.value  
    
    def create_merge_button(self):
        """
        Button to merge hdf5 and archive file.
        """
        button = widgets.Button(
                description='merge & save',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='merge & save files')
                #icon='check')
        return button
    

    def create_clear_output(self):
        clear_output()
