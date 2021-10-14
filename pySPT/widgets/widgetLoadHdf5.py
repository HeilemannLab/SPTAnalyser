"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Widget handling for trackStatistics.ipynb.
"""

import tkinter as tk 
import os
import os.path
import tkinter.filedialog as fd
from ipywidgets import widgets
from IPython.display import clear_output
import datetime


class WidgetLoadHdf5():
    def __init__(self, bin_size, x_range, y_range):
        # Load cell.h5 files
        self.file_names = []  # list of file names for cell files
        self.suffix = ".h5"
        self.dir_button = self.create_dir_button()
        self.dir_box = self.create_dir_box()
        self.dir_name = ""  # input for directory
        self.init_cells_button = self.create_init_cells_button()
        # Load bg.h5 files
        self.init_background_button = self.create_init_background_button()
        self.dir_name_bg = ""  # directory of bg file search
        self.dir_button_bg = self.create_dir_button_bg()
        self.dir_box_bg = self.create_dir_box_bg()
        self.file_names_bg = []  # list of file names for bg files
        self.chosen_cell = ""
        self.cell_options = []  # list of filtered cells for drop down menu
        self.trajectory_options = []  # list of filtered trajectories for drop down menu
        self.drop_down_cells = self.create_drop_down_cells()
        self.drop_down_trajectories = self.create_drop_down_trajectories()
        self.plot_button = self.create_plot_button()
        self.filter_button = self.create_filter_button()
        self.min_length_box = self.create_min_length_box()
        self.min_length = 0  
        self.max_length_box = self.create_max_length_box()
        self.min_D_box = self.create_min_D_box()
        self.max_D_box = self.create_max_D_box()
        self.immob_type_check_box = self.create_immob_type_check_box()
        self.confined_type_check_box = self.create_confined_type_check_box()
        self.free_type_check_box = self.create_free_type_check_box()
        self.analyse_not_successful_check_box = self.create_analyse_not_successful_check_box()
        # Calculate the dynamic localization error
        self.calc_sigma_dyn_button = self.create_calc_sigma_dyn_button()
        # Plot diffusion histogram
        self.bin_size_box = self.create_bin_size_box(val=bin_size)
        self.MSD_delta_t_n = self.create_MSD_delta_t_n(val=x_range)
        self.MSD_y_lim = self.create_MSD_y_lim(val=y_range)
        # Save statistics
        self.save_dir_button = self.create_save_dir_button()
        self.save_dir_box = self.create_save_dir_box()
        self.save_name_box = self.create_save_raw_base_name_box()
        self.dir_save = ""
        self.filtered_dataset_checkbox = self.create_filtered_dataset_checkbox()
        self.hmm_trc_checkbox = self.create_hmm_trc_checkbox()
        self.Dplot_checkbox = self.create_Dplot_checkbox()
        self.save_button = self.create_save_button()
        self.save_folder_name_box = self.create_save_folder_name_box()
    
    def search_sub_folders(self, dir_name, is_cell=True):
        if dir_name:
            for root, dirs, files in os.walk(dir_name):
                self.extend_list(root, files, is_cell)

    def extend_list(self, root, files, is_cell=True):
        for name in files:
            if name.endswith(self.suffix):
                if is_cell:
                    self.file_names.append(os.path.join(root, name))
                else:
                    self.file_names_bg.append(os.path.join(root, name))
                
    def create_dir_button(self):
        button = widgets.Button(
                description="browse",
                disabled=False,
                button_style="",
                tooltip="browse for directory")
        return button    
    
    def open_dir(self, b):
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        root.lift()
        root.update()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title="Please select a directory")
        root.update()
        root.destroy()
        self.dir_name = root.name
        self.dir_box.value=self.dir_name
        
    def create_dir_box(self, val="", desc="Directory"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="directory to be searched in", description=str(desc), disabled=False, style=style)
        return text
    
    def change_dir_box(self, change):
        self.dir_name = self.dir_box.value
        
    def create_init_cells_button(self):
        button = widgets.Button(
                description="initialize",
                disabled=False,
                button_style="",
                tooltip="initialize objects")
        return button    
    
    # Load bg.h5 files
    
    def create_dir_button_bg(self):
        button = widgets.Button(
                description="browse",
                disabled=False,
                button_style="",
                tooltip="browse for directory")
        return button    
    
    def open_dir_bg(self, b):
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        root.lift()
        root.update()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title='Please select a directory')
        root.update()
        root.destroy()
        self.dir_name_bg = root.name
        self.dir_box_bg.value = self.dir_name_bg
        
    def create_dir_box_bg(self, val = "", desc = "Directory"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="directory to be searched in", description=str(desc), disabled=False, style=style)
        return text
    
    def change_dir_box_bg(self, change):
        self.dir_name_bg = self.dir_box_bg.value
    
    def create_init_background_button(self):
        button = widgets.Button(
                description="initialize",
                disabled=False,
                button_style="",
                tooltip="initialize objects")
        return button  
    
    def create_drop_down_cells(self):
        drop_down_cells = widgets.Dropdown(
                options=self.cell_options,
                description="Cell:",
                disabled=False)
        return drop_down_cells
    
    def create_drop_down_trajectories(self):
        drop_down_trajectories = widgets.Dropdown(
                options=self.trajectory_options,
                description="Trajectory:",
                disabled=False)
        return drop_down_trajectories
    
    def get_trajectory_numbers(self, cell, cell_trajectories):
        trajectory_numbers = []
        for trajectory in cell_trajectories[cell]:
            trajectory_numbers.append(trajectory)
        self.drop_down_trajectories.options = trajectory_numbers
        return trajectory_numbers
    
    def get_cell_names(self, cells, filtered_cell_trajectories):
        """
        :param cells: list of cell objects from coverslip.
        :param cell_trajectories_filtered: list with trajectories (for each cell 1 list).
        """
        cell_names = []
        for cell in cells:
            if filtered_cell_trajectories[cells.index(cell)]:  # if cell has trajectory entries 
                cell_names.append(cell.name)
        self.drop_down_cells.options = sorted(cell_names)  # alphabetically sorted
        return sorted(cell_names)
    
    def create_plot_button(self):
        button = widgets.Button(
                description="plot",
                disabled=False,
                button_style="",
                tooltip="plot chosen trajectory")
        return button
    
    def create_clear_output(self):
        clear_output()
    
    def create_filter_button(self):
        button = widgets.Button(
                description="filter & plot",
                disabled=False,
                button_style="",
                tooltip="apply filter and display global information")
        return button
    
    def create_min_length_box(self, val="", desc="Trajectory"):
        style = {"description_width": "initial"}
        text_min = widgets.Text(value=str(val), placeholder="min length", description=str(desc), disabled=False, style=style)
        return text_min
    
    def create_max_length_box(self, val="", desc="Trajectory"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="max length", description=str(desc), disabled=False, style=style)
        return text
    
    def create_min_D_box(self, val="", desc="Diffusion coefficient"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="min value", description=str(desc), disabled=False, style=style)
        return text
    
    def create_max_D_box(self, val="", desc="Diffusion coefficient"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="max value", description=str(desc), disabled=False, style=style)
        return text
    
    def create_immob_type_check_box(self):
        checkbox = widgets.Checkbox(value=True,
                         description="Immobile",
                         disabled=False)
        return checkbox
    
    def create_confined_type_check_box(self):
        checkbox = widgets.Checkbox(value=True,
                         description="Confined",
                         disabled=False)
        return checkbox
    
    def create_free_type_check_box(self):
        checkbox = widgets.Checkbox(value=True,
                         description="Free",
                         disabled=False)
        return checkbox

    def create_analyse_not_successful_check_box(self):
        checkbox = widgets.Checkbox(value=True,
                         description="Type determination not successful",
                         disabled=False)
        return checkbox
    
    def create_calc_sigma_dyn_button(self):
        """
        Button to calculate the dynamic localization error per cell.
        """
        button = widgets.Button(
                description="calc",
                disabled=False,
                button_style="",
                tooltip="dynamic localization error")
        return button  
    
    # Plot diffusion histogram
    
    def create_bin_size_box(self, val="0.1", desc="bin size"):
        """
        Box for inserting the bin size for log10(D) histogram.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="size for log10(D) histogram", description=str(desc), disabled=False, style=style)
        return text

    def create_MSD_delta_t_n(self, val="None", desc="x range"):
        """
        Box for inserting the x range in seconds for MSD plot.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="x range in seconds", description=str(desc), disabled=False, style=style)
        return text

    def create_MSD_y_lim(self, val="None", desc="y range"):
        """
        Box for inserting the y range for MSD plot.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="y limit of MSD plot", description=str(desc), disabled=False, style=style)
        return text

    # Save h5 statistics

    def create_save_dir_button(self):
        button = widgets.Button(
                description="browse",
                disabled=False,
                button_style="",
                tooltip="browse for directory")
        return button  
    
    def save_open_dir(self, b):
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        root.lift()
        root.update()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title="Please select a directory")
        root.update()
        root.destroy()
        self.dir_save = root.name
        self.save_dir_box.value=self.dir_save
        
    def create_save_dir_box(self, val="", desc="directory"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="directory for filtered data", description=str(desc), disabled=False, style=style)
        return text
    
    def change_save_dir_box(self, change):
        self.dir_save = self.save_dir_box.value   
    
    def create_save_raw_base_name_box(self, val="statistics", desc="file name"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="name for statistics .h5 file", description=str(desc), disabled=False, style=style)
        return text
    
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
    
    def create_save_folder_name_box(self, desc="folder name"):
        current_date = self.calc_date()
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(current_date + "_filtered"), placeholder="name of folder", description=str(desc), disabled=False, style=style)
        return text
        
    def create_filtered_dataset_checkbox(self):
        checkbox = widgets.Checkbox(value=True,
                         description="Save filtered dataset",
                         disabled=False)
        return checkbox
    
    def create_Dplot_checkbox(self):
        checkbox = widgets.Checkbox(value=True,
                         description="Save global plots",
                         disabled=False)
        return checkbox
    
    def create_hmm_trc_checkbox(self):
        checkbox = widgets.Checkbox(value=False,
                         description="Save filtered .trc files",
                         disabled=False)
        return checkbox
    
    def create_save_button(self):
        button = widgets.Button(
                description="save",
                disabled=False,
                button_style="",
                tooltip="save statistics")
        return button
