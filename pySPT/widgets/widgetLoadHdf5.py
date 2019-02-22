# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 09:14:39 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Handling widgets of trackStatistics JNB.
"""

import tkinter as tk 
import os
import os.path
import tkinter.filedialog as fd
from ipywidgets import widgets
from IPython.display import display
from IPython.display import clear_output
import math

class WidgetLoadHdf5():
    def __init__(self):
        self.data_set_dir = ""
        self.file_names = []
        self.suffix = ".h5"
        self.dir_button = self.create_dir_button()
        self.dir_box = self.create_dir_box()
        self.dir_name = ""
        self.init_button = self.create_init_button()
        self.chosen_cell = ""
        self.cell_options = []
        self.trajectory_options = []
        self.drop_down_cells = self.create_drop_down_cells()
        self.drop_down_trajectories = self.create_drop_down_trajectories()
        self.plot_button = self.create_plot_button()
        self.filter_button = self.create_filter_button()
        self.min_length_box = self.create_min_length_box()
        self.max_length_box = self.create_max_length_box()
        self.min_D_box = self.create_min_D_box()
        self.max_D_box = self.create_max_D_box()
        self.immob_type_check_box = self.create_immob_type_check_box()
        self.confined_type_check_box = self.create_confined_type_check_box()
        self.free_type_check_box = self.create_free_type_check_box()
        self.analyse_successful_check_box = self.create_analyse_successful_check_box()
        self.analyse_not_successful_check_box = self.create_analyse_not_successful_check_box()
        self.plot_diffusions_button = self.create_plot_diffusions_button()
        # Plot diffusion histogram
        self.bin_size_box = self.create_bin_size_box()
        
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

    def create_drop_down_cells(self):
        drop_down_cells = widgets.Dropdown(
                options=self.cell_options,
                description='Number:',
                disabled=False)
        return drop_down_cells
    
    def create_drop_down_trajectories(self):
        drop_down_trajectories = widgets.Dropdown(
                options=self.trajectory_options,
                description='Number:',
                disabled=False)
        return drop_down_trajectories
    
    def get_trajectory_numbers(self, cell, cell_trajectories):
        trajectory_numbers = []
        for trajectory in cell_trajectories[cell]:
            trajectory_numbers.append(trajectory)
        self.drop_down_trajectories.options = trajectory_numbers
        return trajectory_numbers
    
    def get_cell_names(self, cells):
        cell_names = []
        for cell in cells:
            cell_names.append(cell.name)
        self.drop_down_cells.options = cell_names
        return cell_names
    
    def create_plot_button(self):
        button = widgets.Button(
                description="plot",
                disabled=False,
                button_style="",
                tooltip = "plot chosen trajectory")
        return button
    
    def create_clear_output(self):
        clear_output()
    
    def create_filter_button(self):
        button = widgets.Button(
                description="apply filter",
                disabled=False,
                button_style="",
                tooltip = "apply filter")
        return button
    
    def create_min_length_box(self, val = "min length" , desc = "Trajectory"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the minimum length of a trajectory.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_max_length_box(self, val = "max length", desc = "Trajectory"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the max length of a trajectory.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_min_D_box(self, val = "min value" , desc = "Diffusion coefficient"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the minimum D value.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_max_D_box(self, val = "max value", desc = "Diffusion coefficient"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the max D value.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_immob_type_check_box(self):
        """
        True -> check box is already selected; False -> check box is not selected.
        """
        checkbox = widgets.Checkbox(value=True,
                         description='Immobile',
                         disabled=False)
        return checkbox
    
    def create_confined_type_check_box(self):
        checkbox = widgets.Checkbox(value=True,
                         description='Confined',
                         disabled=False)
        return checkbox
    
    def create_free_type_check_box(self):
        checkbox = widgets.Checkbox(value=True,
                         description='Free',
                         disabled=False)
        return checkbox
    
    def create_analyse_successful_check_box(self):
        checkbox = widgets.Checkbox(value=True,
                         description='Analyse successful',
                         disabled=False)
        return checkbox

    def create_analyse_not_successful_check_box(self):
        checkbox = widgets.Checkbox(value=True,
                         description='Analyse not successful',
                         disabled=False)
        return checkbox
    
    def create_plot_diffusions_button(self):
        button = widgets.Button(
                description="plot",
                disabled=False,
                button_style="",
                tooltip = "plot diffusion coefficients")
        return button
    
    # Plot diffusion histogram
    
    def create_bin_size_box(self, val = "0.1" , desc = "bin size"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the bin size for log10(D) histogram.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='size for log10(D) histogram', description=str(desc), disabled=False, style = style)
        return text
    
    