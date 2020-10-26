# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 16:21:52 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import tkinter as tk 
import os
import os.path
import tkinter.filedialog as fd
from ipywidgets import widgets
from IPython.display import clear_output

class WidgetTrackAnalysis():
    def __init__(self, pixel_size, area_camera, camera_dt, n_points_D, fit_area_MSD, dof_D, min_D, min_length,
                 hmm_min_length, hmm_float, bin_size, x_range, y_range):
        if int(n_points_D) >= int(min_length):
            min_length = str(int(n_points_D)+1)
        if int(n_points_D) >= int(hmm_min_length):
            hmm_min_length = str(int(n_points_D)+1)
        self.data_set_dir = ""
        self.ignore_words_box = self.create_ignore_words_box()
        self.masked_words = []
        self.software_button = self.create_software_button()
        self.file_names = []
        self.suffix = ""
        self.dir_button = self.create_dir_button()
        self.dir_box = self.create_dir_box()
        self.dir_name = ""
        self.roi_name = ""
        self.roi_box = self.create_roi_box()
        self.roi_button = self.create_roi_button()
        self.got_dir = False
        self.got_roi = False
        self.camera_pixel_size_box = self.create_camera_pixel_size_box(val=pixel_size)
        self.camera_pixel_amount_box = self.create_camera_pixel_amount_box(val=area_camera)
        self.camera_integration_time_box = self.create_camera_integration_time_box(val=camera_dt)
        self.min_track_length_box = self.create_min_track_length_box(val=min_length)
        self.rossier_fit_area_box = self.create_rossier_fit_area_box(val=fit_area_MSD)
        self.dof_box = self.create_dof_box(val=dof_D)
        self.D_min_box = self.create_D_min_box(val=min_D)
        self.points_D_fit_box = self.create_points_D_fit_box(val=n_points_D)
        self.hmm_check_box = self.create_hmm_check_box()
        self.hmm_trc_float_precision_box = self.create_hmm_trc_float_precision_box(val=hmm_float)
        self.microscope_check_box = self.create_microscope_check_box()
        self.min_track_length_hmm_box = self.create_min_track_length_hmm_box(val=hmm_min_length)
        self.run_button = self.create_run_button()
        self.chosen_cell = ""
        self.cell_options = []
        self.trajectory_options = []
        self.drop_down_cells = self.create_drop_down_cells()
        self.drop_down_trajectories = self.create_drop_down_trajectories()
        self.plot_button = self.create_plot_button()
        self.save_button = self.create_save_button()
        self.trajectory_id_button = self.create_trajectory_id_button()
        # Plot diffusion histogram
        self.bin_size_box = self.create_bin_size_box(val=bin_size)
        self.MSD_delta_t_n = self.create_MSD_delta_t_n(val=x_range)
        self.MSD_y_lim = self.create_MSD_y_lim(val=y_range)
        self.plot_diff_button = self.create_plot_diff_button()
        
        
    def create_rossier_fit_area_box(self, val = "0.6", desc = "Fit area MSD"):
        """
        Box for inserting amount of MSD values to be fitted by rossier [0..1].
        """        
        style = {"description_width": "initial"}
        text = widgets.Text(value = str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_hmm_trc_float_precision_box(self, val = "10", desc = "Trc float precision"):  
        """
        Box for inserting the min track length for hmm.
        """        
        style = {"description_width": "initial"}
        text = widgets.Text(value = str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_hmm_check_box(self):
        """
        True -> Save hmm.trc file in pySPT/hmm folder.
        """
        checkbox = widgets.Checkbox(value=True,
                         description='Save .trc file',
                         disabled=False)
        return checkbox
    
    def create_microscope_check_box(self):
        """
        True -> Save microscope file (containing dt, pixel size, sigma dyn) in pySPT/hmm folder.
        """
        checkbox = widgets.Checkbox(value=True,
                         description='Save .microscope file',
                         disabled=False)
        return checkbox
    
    def create_min_track_length_hmm_box(self, val = "20", desc = "Min track length"):  # val = "F:\\Marburg\\single_colour_tracking\\resting\\roi.log"
        """
        Box for inserting the min track length for hmm.
        """        
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text #value=str(self.calc_min_track_length_hmm()
        
    def create_software_button(self):
        """
        Radiobutton to choose from PALMTracer, rapidSTORM+swift or ThunderSTORM+swift.
        PALMTracer has .trc files, rapidSTORM .txt and ThunderSTORM .csv files that will be loaded.
        """
        button = widgets.RadioButtons(
                options = ["ThunderSTORM", "rapidSTORM", "PALMTracer"],
                disabled = False)
        return button
        
    def create_ignore_words_box(self, val = "", desc = "Mask words"):
        """
        The string in the box contains words that lead to not loading a file if one of the words is contained.
        Commaseparate mask words.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='type something', description=str(desc), disabled=False, style = style)
        return text       
        
    def create_trajectory_id_button(self):
        """
        Radiobutton to choose between seg id or track id,
        id is used for diffusion type analysis.
        """
        button = widgets.RadioButtons(
                options = ["seg id", "track id"],
                disabled = False)
        return button
   
    def searchSubFolders(self, dirName):
        if (dirName):
            self.data_set_dir = dirName
            for root, dirs, files in os.walk(self.data_set_dir):
                self.determine_suffix()
                self.extendList(root, files)

    def determine_suffix(self):
        """
        Depending on the chosen software, the file ending of the target files differs.
        """
        if self.software_button.value == "ThunderSTORM":
            self.suffix = "tracked.csv"
        elif self.software_button.value == "rapidSTORM":
            self.suffix = ".tracked.txt"
        elif self.software_button.value == "PALMTracer":
            self.suffix = ".trc"
                
    def extendList(self, root, files):
        # create mask word list
        ignore_words = self.ignore_words_box.value.replace(" ", "")
        self.masked_words = ignore_words.split(",")          
        if self.masked_words == [""]:
            self.masked_words = [] 
        for name in files:
            #if name.endswith(self.suffix) and "background" not in name:
            if name.endswith(self.suffix):
                if self.masked_words:
                    if not any(x in name for x in self.masked_words):  #does not work with empty string
                        self.file_names.append(os.path.join(root, name))
                else:
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
        
    def create_dir_box(self, val = "", desc = "Directory"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the directory with description, alternative for dir loading button.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='directory to be searched in', description=str(desc), disabled=False, style = style)
        return text
    
    def change_dir_box(self, change):
        self.dir_name = self.dir_box.value  
        self.got_dir = True

    def create_roi_box(self, val = "", desc = "roi"):  # val = "F:\\Marburg\\single_colour_tracking\\resting\\roi.log"
        """
        Box for inserting the roi file, alternative for roi loading button.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='path of roi', description=str(desc), disabled=False, style = style)
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
    
    def create_points_D_fit_box(self, val = "4", desc = "Number of points fitted for D"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
        
    def create_D_min_box(self, val = "0.0038", desc = "Minimal detectable D  [\u03BCm\u00b2/s]"):
        """
        Box for inserting the camera integration time for tau threshold calculation.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_dof_box(self, val = "4", desc = "Degree of freedom of D"):
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

    def create_MSD_delta_t_n(self, val = "None" , desc = "x range"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the bin size for log10(D) histogram.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='number of MSD values shown', description=str(desc), disabled=False, style = style)
        return text

    def create_MSD_y_lim(self, val = "None" , desc = "y range"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the bin size for log10(D) histogram.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='y limit of MSD plot', description=str(desc), disabled=False, style = style)
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