# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 10:36:58 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Handle widgets for the localization uncertainty JNB.
"""

import tkinter as tk 
import os
from tkinter.filedialog import askopenfilename
import tkinter.filedialog as fd
from ipywidgets import widgets
from IPython.display import display
from IPython.display import clear_output
#from ..preAnalysis import expDisplacement


class WidgetPrecision():
    def __init__(self):
        self.software_button = self.create_software_button()
        # precision per folder
        self.got_dir = False
        self.dir_name = ""
        self.dir_button = self.create_dir_button()
        self.dir_box = self.create_dir_box()
        self.run_button_folder = self.create_run_button_folder()
        self.save_button_folder = self.create_save_button_folder()
        self.save_fig_checkbox = self.create_save_fig_checkbox()
        self.got_dir_save = False
        self.box_foldername = self.create_box_foldername()
        self.dir_button_save = self.create_dir_button()
        self.dir_box_save = self.create_dir_box(placeholder="Directory to save")
        self.save_figures_checkbox_folder = self.create_save_figures_checkbox_folder()
        self.check_microscope_folder = self.create_check_microscope_folder()
        # precision per file
        self.file_name = ""
        self.got_file_name = False
        self.file_text_box = self.create_file_box()
        self.file_button = self.create_file_button()
        self.run_button = self.create_run_button()
        self.camera_pixel_size_box = self.create_camera_pixel_size_box()
        self.camera_integration_time_box = self.create_camera_integration_time_box()
        self.check_microscope = self.create_check_microscope()
        self.save_button = self.create_save_button()
        self.save_figures_checkbox = self.create_save_figures_checkbox()

    # precison per folder

    def create_software_button(self):
        """
        Radiobutton to choose between rapidSTORM and thunderSTORM.
        """
        button = widgets.RadioButtons(
                options = ["ThunderSTORM", "rapidSTORM"],
                disabled = False)
        return button

    def create_dir_button(self):
        """
        Button to load a directory as search platform.
        """
        button = widgets.Button(
            description='browse',
            disabled=False,
            button_style='',  # 'success', 'info', 'warning', 'danger' or ''
            tooltip='browse for directory')
        # icon='check')
        return button

    def create_dir_box(self, val="",
                       desc="Directory", placeholder='directory to be searched in'):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the directory with description, alternative for dir loading button.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder=placeholder, description=str(desc),
                            disabled=False, style=style)
        return text

    def open_dir(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title='Please select a directory')
        root.update()
        root.destroy()
        self.dir_name = root.name
        self.dir_box.value = self.dir_name
        self.got_dir = True

    def change_dir_box(self, change):
        self.dir_name = self.dir_box.value
        self.got_dir = True

    def create_file_button(self):
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

    def create_run_button_folder(self):
        """
        Button for running the analysis, has an on click event.
        """
        button = widgets.Button(
                description='run',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='run the analysis')
                #icon='check')
        return button

    def create_save_button_folder(self):
        """
        Button to save the results, has an on click event.
        """
        button = widgets.Button(
            description='save',
            disabled=False,
            button_style='',  # 'success', 'info', 'warning', 'danger' or ''
            tooltip='save the results')
        # icon='check')
        return button

    def create_save_fig_checkbox(self):
        check_box = widgets.Checkbox(
                value=True,
                description='Save plot',
                disabled=False)
        return check_box

    def open_dir_save(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title='Please select a directory')
        root.update()
        root.destroy()
        self.dir_name_save = root.name
        self.dir_box_save.value = self.dir_name_save
        self.got_dir_save = True

    def change_dir_box_save(self, change):
        self.dir_name_save = self.dir_box_save.value
        self.got_dir_save = True

    def create_save_box(self, val="", desc="Complete path"):  # val = in box, desc = infront of box
        """
        Box for inserting the path with description, alternative for file loading button.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='insert path', description=str(desc), disabled=False,
                            style=style)
        return text

    def create_box_foldername(self, val="precision",
                       desc="Foldername"):  # val = in box, desc = infront of box; val = "C:\\Users\\pcoffice37\\Documents\\testing_file_search"
        """
        Box for inserting the directory with description, alternative for dir loading button.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='name of folder', description=str(desc),
                            disabled=False, style=style)
        return text

    def create_save_figures_checkbox_folder(self):
        check_box = widgets.Checkbox(
                value=True,
                description='Save plots',
                disabled=False)
        return check_box

    def create_check_microscope_folder(self):
        """
        For HMM-analysis a microscope file with camera px size, integration time and localization
        precision [nm] is needed. If True, file will be saved.
        """
        check_box = widgets.Checkbox(
                value=True,
                description='Save microscope file for HMM analysis?',
                disabled=False)
        return check_box

    # precision per file
        
    def open_file(self, b):  # b = ???
        """
        Give the file button opening powers.
        """
        root = tk.Tk()  # window class
        root.withdraw()  # close the window 
        root.update()  # close the window
        if self.software_button.value == "ThunderSTORM":
            root.name = askopenfilename(title="Import .csv file", filetypes=(("csv files", "*.csv"),("all files", "*.*")))
        else:
            root.name = askopenfilename(title="Import .txt file", filetypes=(("text files", "*.txt"),("all files", "*.*")))
        self.file_name = root.name
        root.update()
        root.destroy()
        self.file_text_box.value=self.file_name
        if os.path.isfile(self.file_name):
            self.got_file_name = True

    def create_file_box(self, val = "", desc = "Complete path"):  # val = in box, desc = infront of box
        """
        Box for inserting the path with description, alternative for file loading button.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='insert path', description=str(desc), disabled=False, style = style)
        return text
    
    def change_file_box(self, change):
        self.file_name = self.file_text_box.value  
        if os.path.isfile(self.file_name):
            self.got_file_name = True
    
    def create_run_button(self):
        """
        Button for running the analysis, has an on click event.
        """
        button = widgets.Button(
                description='run',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='run the analysis')
                #icon='check')
        return button
    
        self.camera_pixel_size_box = self.create_camera_pixel_size_box()
        self.camera_pixel_amount_box = self.create_camera_pixel_amount_box()
        self.camera_integration_time_box = self.create_camera_integration_time_box()
        
    def create_camera_integration_time_box(self, val = "0.02", desc = "Camera integration time [s]"):
        """
        Box for inserting the camera integration time for tau threshold calculation.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_camera_pixel_size_box(self, val = "158", desc = "Pixel size [nm]"):
        """
        Box for inserting the pixel size in nm of the camera.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_check_microscope(self):
        """
        For HMM-analysis a microscope file with camera px size, integration time and localization
        precision [nm] is needed. If True, file will be saved.
        """
        check_box = widgets.Checkbox(
                value=True,
                description='Save microscope file for HMM analysis?',
                disabled=False)
        return check_box

    def create_save_button(self):
        """
        Button to save the results, has an on click event.
        """
        button = widgets.Button(
                description='save',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='save the results')
                #icon='check')
        return button

    def warning_wrong_file_path(self):
        print("The file path is empty or does not exist.")
        
    def warning_wrong_file(self):
        print("A file with false columns was loaded.")
        
    def create_clear_output(self):
        clear_output()
        
    def create_save_figures_checkbox(self):
        check_box = widgets.Checkbox(
                value=True,
                description='Save plots',
                disabled=False)
        return check_box
        
        
def main():
    widget_precision = WidgetPrecision()
    widget_precision.open_file()
    print(widget_precision.file_name, widget_precision.folder_name, widget_precision.base_name, sep="\n")
    
    
if __name__ == "__main__":
    main()
    