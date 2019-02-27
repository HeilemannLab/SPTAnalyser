# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 16:21:52 2019

@author: pcoffice37

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

class WidgetTrackAnalysis():
    def __init__(self):
        self.data_set_dir = ""
        self.file_names = []
        self.suffix = ".trc"
        self.dir_button = self.create_dir_button()
        self.dir_box = self.create_dir_box()
        self.dir_name = ""
        self.roi_name = ""
        self.roi_box = self.create_roi_box()
        self.roi_button = self.create_roi_button()
        self.got_dir = False
        self.got_roi = False
        self.camera_pixel_size_box = self.create_camera_pixel_size_box()
        self.camera_pixel_amount_box = self.create_camera_pixel_amount_box()
        self.camera_integration_time_box = self.create_camera_integration_time_box()
        self.min_track_length_box = self.create_min_track_length_box()
        self.dof_box = self.create_dof_box()
        self.D_min_box = self.create_D_min_box()
        self.run_button = self.create_run_button()
        self.chosen_cell = ""
        self.cell_options = []
        self.trajectory_options = []
        self.drop_down_cells = self.create_drop_down_cells()
        self.drop_down_trajectories = self.create_drop_down_trajectories()
        self.plot_button = self.create_plot_button()
        self.save_button = self.create_save_button()
        # Plot diffusion histogram
        self.bin_size_box = self.create_bin_size_box()
        self.plot_diff_button = self.create_plot_diff_button()
   
    def searchSubFolders(self, dirName):
        if (dirName):
            self.data_set_dir = dirName
            for root, dirs, files in os.walk(self.data_set_dir):
                self.extendList(root, files)
                
    def extendList(self, root, files):
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

    def create_roi_box(self, val = "path of roi", desc = "roi"):  # val = "F:\\Marburg\\single_colour_tracking\\resting\\roi.log"
        """
        Box for inserting the roi file, alternative for roi loading button.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text

    def change_roi_box(self, change):
        self.roi_name = self.roi_box.value
        self.got_roi = True
    
    def create_roi_button(self):
        """
        Button to load a roi file for cell size.
        """
        button = widgets.Button(
                description='browse',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='browse for roi')
                #icon='check')
        return button        

    def open_roi(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.name = fd.askopenfilename(title="Import roi.log file", filetypes=(("log files", "*.log"),("all files", "*.*")))
        root.update()
        root.destroy()
        self.roi_name = root.name
        self.roi_box.value=self.roi_name
        self.got_roi = True     
        
    def create_camera_pixel_size_box(self, val = "158", desc = "Pixel size [nm]"):
        """
        Box for inserting the pixel size in nm of the camera.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
        
    def create_camera_pixel_amount_box(self, val = "65536", desc = "Amount of pixel on the camera"):
        """
        Box for inserting the amount of pixel on the camera
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_min_track_length_box(self, val = "20", desc = "Min track length"):
        """
        Box for inserting the minimal track length for tau threshold calculation.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_camera_integration_time_box(self, val = "0.02", desc = "Camera integration time [s]"):
        """
        Box for inserting the camera integration time for tau threshold calculation.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_D_min_box(self, val = "0.0065", desc = "Minimal detectable D  [\u03BCm\u00b2/s]"):
        """
        Box for inserting the camera integration time for tau threshold calculation.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_dof_box(self, val = "4", desc = "degree of freedom of D"):
        """
        Box for inserting the camera integration time for tau threshold calculation.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_run_button(self):
        """
        Button to load a roi file for cell size.
        """
        button = widgets.Button(
                description='run',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='run the analysis')
                #icon='check')
        return button             
        
    def printFileNames(self):
        for name in self.fileNames:
            print("Found .trc file: %s" %(name))
            
    def float_progress(self):
        float_progress = widgets.FloatProgress(
                value=7.5,
                min=0,
                max=10.0,
                step=0.1,
                description='Loading:',
                bar_style='info',
                orientation='horizontal')
        return float_progress
    
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
        self.drop_down_cells.options = sorted(cell_names)
        return sorted(cell_names)
    
    def create_plot_button(self):
        button = widgets.Button(
                description="plot",
                disabled=False,
                button_style="",
                tooltip = "plot chosen trajectory")
        return button

    def create_clear_output(self):
        clear_output()
        
    def warning_trc_file(self):
        print("No *.trc files were loaded.")
        
    def create_save_button(self):
        button = widgets.Button(
                description="save",
                disabled=False,
                button_style="",
                tooltip = "save entire analysis")
        return button
    
    # Plot diffusion histogram
    
    def create_bin_size_box(self, val = "0.1" , desc = "bin size"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the bin size for log10(D) histogram.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='size for log10(D) histogram', description=str(desc), disabled=False, style = style)
        return text
    
    def create_plot_diff_button(self):
        button = widgets.Button(
                description="plot",
                disabled=False,
                button_style="",
                tooltip = "plot diffusion histogram")
        return button


        
            
def main():
    communicator = Mol2Judi()
    communicator.searchSubFolders("E:/Receptor-signaling/met/rtk/resting")
    communicator.printFileNames()
    
if __name__ == '__main__':
    main()
    