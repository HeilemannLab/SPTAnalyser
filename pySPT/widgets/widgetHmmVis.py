# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 17:28:20 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import tkinter as tk 
import os
import os.path
import tkinter.filedialog as fd
from ipywidgets import widgets
from IPython.display import display
from IPython.display import clear_output
import datetime


class WidgetHmmVis():
    def __init__(self):
        # loading
        self.load_dir_box = self.create_load_dir_box()
        self.load_dir_button = self.create_load_dir_button()
        self.load_dir_name = ""
        self.file_names = []
        self.suffix = ".h5"
        # plotting
        self.plot_button = self.create_plot_button()
        # saving
        self.save_plots_checkbox = self.create_save_plots_checkbox()
        self.save_dir_box = self.create_save_dir_box()
        self.save_dir_button = self.create_save_dir_button()
        self.save_dir_name = ""
        self.save_folder_name_box = self.create_save_folder_name_box() 
        self.save_button = self.create_save_button()
        
    def create_load_dir_box(self, val = "path", desc = "Complete path"):  # val = in box, desc = infront of box
        """
        Box for inserting the path with description, alternative for file loading button.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_load_dir_button(self):
        """
        Button to load the file.
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
        self.load_dir_name = root.name
        self.load_dir_box.value = self.load_dir_name
        
    def change_load_dir_box(self, change):
        self.load_dir_name = self.load_dir_box.value  
    
    # plotting
    
    def create_plot_button(self):
        """
        Button to plot.
        """
        button = widgets.Button(
                description='plot',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='click to plot')
                #icon='check')
        return button
    
    def search_sub_folders(self, dir_name):
        if dir_name:
            for root, dirs, files in os.walk(dir_name):
                self.extend_list(root, files)
 
    def extend_list(self, root, files):
        for name in files:
            if name.endswith(self.suffix):
                self.file_names.append(os.path.join(root, name))
   
    # saving
    
    def create_save_plots_checkbox(self):
        """
        If true, all plots will be saved in a folder.
        """
        checkbox = widgets.Checkbox(value=True,
                         description='Save plots',
                         disabled=False)
        return checkbox
    
    
    def create_save_dir_box(self, val = "directory", desc = "Insert directory"):  # val = in box, desc = infront of box
        """
        Box for inserting the path with description, alternative for file loading button.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_save_dir_button(self):
        """
        Button to load the file.
        """
        button = widgets.Button(
                description='browse',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='browse for directory')
                #icon='check')
        return button
    
    def save_open_dir(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.name = fd.askdirectory(initialdir=os.getcwd(),title='Please select a directory') 
        root.update()
        root.destroy()
        self.save_dir_name = root.name
        self.save_dir_box.value = self.save_dir_name
        
    def change_save_dir_box(self, change):
        self.save_dir_name = self.save_dir_box.value  
    
    def create_save_folder_name_box(self, val = "name", desc = "Folder name"):  # val = in box, desc = infront of box
        """
        Box for inserting the folder name in which the statistics file and the plots are saved.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text

    def create_save_button(self):
        """
        Button to load the file.
        """
        button = widgets.Button(
                description='save',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='click to save')
                #icon='check')
        return button
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
# =============================================================================
#     def __init__(self):
#         self.file_names = []  # list of file names for cell files
#         self.suffix = ".h5"
#         self.dir_button = self.create_dir_button()
#         self.dir_box = self.create_dir_box()
#         self.dir_name = ""  # input for directory
#         self.open_dir = ""
#         self.test_button = self.create_test_button()
#         self.plot_button = self.create_plot_button()
#         self.save_button = self.create_save_button()
#         self.save_dir_box = self.create_save_dir_box()
#         self.save_folder_name = ""
#         self.save_dir_name = ""
#         self.save_dir_browse_button = self.create_save_dir_browse_button()
#         self.save_folder_name_box = self.create_save_folder_name_box()
#         self.save_plots_checkbox = self.create_save_plots_checkbox()
# =============================================================================
        
        
        
        
        
        
        
        
        
        
# =============================================================================
#     def create_plot_button(self):
#         button = widgets.Button(
#                 description='plot',
#                 disabled=False,
#                 button_style='', # 'success', 'info', 'warning', 'danger' or ''
#                 tooltip='browse for directory')
#                 #icon='check')
#         return button  
#             
#     def create_test_button(self):
#         button = widgets.Button(
#                 description='browse',
#                 disabled=False,
#                 button_style='', # 'success', 'info', 'warning', 'danger' or ''
#                 tooltip='browse for directory')
#                 #icon='check')
#         return button  
#     
#     def create_dir_box(self, val = "", desc = "directory"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
#         """
#         Box for inserting the directory with description, alternative for dir loading button.
#         """
#         style = {'description_width': 'initial'}  # display too long desc
#         text = widgets.Text(value=str(val), placeholder='directory for filtered data', description=str(desc), disabled=False, style = style)
#         return text
# 
#     def search_sub_folders(self, dir_name):
#         if dir_name:
#             for root, dirs, files in os.walk(dir_name):
#                 self.extend_list(root, files)
# 
#     def extend_list(self, root, files):
#         for name in files:
#             if name.endswith(self.suffix):
#                 self.file_names.append(os.path.join(root, name))
#                
#     def create_dir_button(self):
#         """
#         Button to load a directory as search platform.
#         """
#         button = widgets.Button(
#                 description='browse',
#                 disabled=False,
#                 button_style='', # 'success', 'info', 'warning', 'danger' or ''
#                 tooltip='browse for directory')
#                 #icon='check')
#         return button   
#     
#     def open_dir(self, b):
#         print("xxx")
#         root = tk.Tk()
#         root.withdraw()
#         root.update()
#         root.name = fd.askdirectory(initialdir=os.getcwd(),title='Please select a directory') 
#         root.update()
#         root.destroy()
#         self.open_dir = root.name
#         self.dir_box.value=self.open_dir
#         
#     def change_dir_box(self, change):
#         self.open_dir = self.dir_box.value  
#         
#     
#     # saving  
#     
#     def create_save_button(self):
#         button = widgets.Button(
#                 description='save',
#                 disabled=False,
#                 button_style='', # 'success', 'info', 'warning', 'danger' or ''
#                 tooltip='save results')
#                 #icon='check')
#         return button 
#     
#     
#     def create_save_dir_browse_button(self):
#         button = widgets.Button(
#                 description='save',
#                 disabled=False,
#                 button_style='', # 'success', 'info', 'warning', 'danger' or ''
#                 tooltip='save results')
#                 #icon='check')
#         return button
#         
#     def create_save_dir_box(self, val = "", desc = "directory"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
#         """
#         Box for inserting the directory with description, alternative for dir loading button.
#         """
#         style = {'description_width': 'initial'}  # display too long desc
#         text = widgets.Text(value=str(val), placeholder='directory for filtered data', description=str(desc), disabled=False, style = style)
#         return text
#     
#     def change_save_dir_box(self, change):
#         self.dir_save = self.save_dir_box.value  
#         
#     def create_save_folder_name_box(self, val = "folder name", desc = "directory"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
#         """
#         Box for inserting the folder name for saving.
#         """
#         style = {'description_width': 'initial'}  # display too long desc
#         text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
#         return text
#     
# 
#     def create_save_plots_checkbox(self):
#         checkbox = widgets.Checkbox(value=True,
#                          description='Save plots',
#                          disabled=False)
#         return checkbox
#     
#     
#     def save_open_dir(self, b):
#         root = tk.Tk()
#         root.withdraw()
#         root.update()
#         root.name = fd.askdirectory(initialdir=os.getcwd(),title='Please select a directory') 
#         root.update()
#         root.destroy()
#         self.save_dir_name = root.name
#         self.save_dir_box.value=self.dir_save
# 
# =============================================================================
