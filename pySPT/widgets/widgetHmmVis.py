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
        self.add_graphviz_bin_box = self.create_add_graphviz_bin_box()
        self.add_graphviz_tmp_box = self.create_add_graphviz_tmp_box()
        self.load_dir_box = self.create_load_dir_box()
        self.load_dir_button = self.create_load_dir_button()
        self.load_dir_name = ""
        self.file_names = []
        self.suffix = ".h5"
        # plotting
        self.plot_button = self.create_plot_button()
        self.cell_options = []
        self.trajectory_options = []
        self.drop_down_cells = self.create_drop_down_cells()
        self.drop_down_trajectories = self.create_drop_down_trajectories()
        self.plot_trajectory_button = self.create_plot_trajectory_button()
        # saving
        self.save_plots_checkbox = self.create_save_plots_checkbox()
        self.save_dir_box = self.create_save_dir_box()
        self.save_dir_button = self.create_save_dir_button()
        self.save_dir_name = ""
        self.save_folder_name_box = self.create_save_folder_name_box() 
        self.save_button = self.create_save_button()
        
    def create_add_graphviz_bin_box(self, val = "C:\\Program Files (x86)\\Graphviz2.38\\bin", desc = "Graphviz bin path"):  # val = in box, desc = infront of box 
        """
        Box for inserting the path to the bin folder of the graphviz installation (in programs).
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_add_graphviz_tmp_box(self, val = "C:\\Users\\pcoffice37\\Documents\\graphviz_tmp", desc = "Graphviz tmp path"):  # val = in box, desc = infront of box 
        """
        Box for inserting the path to the bin folder of the graphviz installation (in programs).
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
        
    def create_load_dir_box(self, val = "directory", desc = "Directory"):  # val = in box, desc = infront of box
        """
        Box for inserting the directory from which all .h5 files will be loaded.
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
                
    def create_drop_down_cells(self):
        drop_down_cells = widgets.Dropdown(
                options=self.cell_options,
                description='Number:',
                disabled=False)
        return drop_down_cells
    
    def get_cell_names(self, cells):
        cell_names = []
        for cell in cells:
            cell_names.append(cell.hmm_cell_name)
        self.drop_down_cells.options = sorted(cell_names)
        return sorted(cell_names)
    
    def create_drop_down_trajectories(self):
        drop_down_trajectories = widgets.Dropdown(
                options=self.trajectory_options,
                description='Number:',
                disabled=False)
        return drop_down_trajectories
    
    def create_plot_trajectory_button(self):
        """
        Button to plot single trajectories
        """
        button = widgets.Button(
                description='plot',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='plot chosen trajectory')
                #icon='check')
        return button
   
    # saving
    
    def create_save_plots_checkbox(self):
        """
        If true, all plots will be saved in a folder.
        """
        checkbox = widgets.Checkbox(value=True,
                         description='Save plots',
                         disabled=False)
        return checkbox
    
    
    def create_save_dir_box(self, val = "", desc = "Insert directory"):  # val = in box, desc = infront of box
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
        
    def calc_date(self):
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        day = str(now.day)
        if len(month) == 1:
            month = str(0) + month
        if len(day) == 1:
            day = str(0) + day
        date = str(year + month + day)
        return date
    
    def create_save_folder_name_box(self, desc = "Folder name"):
        """
        Box for inserting the raw base name for statistics h5 file.
        """
        current_date = self.calc_date()
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(current_date + "_hmm_results"), placeholder='name of folder', description=str(desc), disabled=False, style = style)
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
    
    def create_clear_output(self):
        clear_output()
    
    