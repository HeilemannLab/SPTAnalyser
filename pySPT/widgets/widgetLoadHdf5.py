# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 09:14:39 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Jupyter Notebook widget handling for loadHdf5 class.
"""

import tkinter as tk 
import os
import os.path
import tkinter.filedialog as fd
from ipywidgets import widgets
from IPython.display import display
from IPython.display import clear_output

class WidgetLoadHdf5():
    def __init__(self):
        self.data_set_dir = ""
        self.file_names = []
        self.suffix = ".h5"
        self.dir_button = self.create_dir_button()
        self.dir_box = self.create_dir_box()
        self.dir_name = ""
        self.init_button = self.create_init_button()
        
        
    def search_sub_folders(self, dirName):
        if (dirName):
            self.data_set_dir = dirName
            for root, dirs, files in os.walk(self.data_set_dir):
                self.extend_list(root, files)
                
    def extend_list(self, root, files):
        for name in files:
            #if name.endswith(self.suffix) and "background" not in name:
            if name.endswith(self.suffix):
                self.file_names.append(os.path.join(root, name))
                
    def create_dir_button(self):
        """
        Button to load a directory as search platform.
        """
        button = widgets.Button(
                description='browse',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='browse for directory')
                #icon='check')
        return button    
    
    def open_dir(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.name = fd.askdirectory(initialdir=os.getcwd(),title='Please select a directory') 
        root.update()
        root.destroy()
        self.dir_name = root.name
        self.dir_box.value=self.dir_name
        self.got_dir = True
        
    def create_dir_box(self, val = "directory to be searched in", desc = "directory"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the directory with description, alternative for dir loading button.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def change_dir_box(self, change):
        self.dir_name = self.dir_box.value  
        self.got_dir = True
        
    def create_init_button(self):
        """
        Initialize objects
        """
        button = widgets.Button(
                description='initialize',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='initialize objects')
                #icon='check')
        return button    
