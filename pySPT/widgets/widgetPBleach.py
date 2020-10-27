"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Widget handling for pBleach.ipynb.
"""

import tkinter as tk 
import os
from tkinter.filedialog import askopenfilename
from ipywidgets import widgets
from IPython.display import clear_output


class WidgetPBleach():
    def __init__(self, k, camera_dt, n_points):
        self.software_button = self.create_software_button()
        self.file_name = ""
        self.got_file_name = False
        self.file_text_box = self.create_file_box()
        self.file_button = self.create_file_button()
        self.run_button = self.create_run_button()
        self.save_button = self.create_save_button()
        self.init_k = self.create_init_k_box(val=k)
        self.integration_time = self.create_init_integration_time_box(val=camera_dt)
        self.ignore_points = self.create_ignore_points_box(val=n_points)
        self.save_plot_checkbox = self.create_save_plot_checkbox()
        
    def create_software_button(self):
        button = widgets.RadioButtons(
                options=["ThunderSTORM", "rapidSTORM"],
                disabled=False)
        return button
        
    def create_file_button(self):
        button = widgets.Button(
                description="browse",
                disabled=False,
                button_style="",
                tooltip="browse for file")
        return button

    def open_file(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        if self.software_button.value == "ThunderSTORM":
            root.name = askopenfilename(title="Import .tracked.csv file", filetypes=(("csv files", "*.csv"), ("all files", "*.*")))
        elif self.software_button.value == "rapidSTORM":
            root.name = askopenfilename(title="Import .tracked.txt file", filetypes=(("text files", "*.txt"), ("all files", "*.*")))
        self.file_name = root.name
        root.update()
        root.destroy()
        self.file_text_box.value = self.file_name
        if os.path.isfile(self.file_name):
            self.got_file_name = True
        
    def create_file_box(self, val="", desc="Complete path"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Insert path', description=str(desc), disabled=False, style=style)
        return text
    
    def create_ignore_points_box(self, val="0", desc="Number of points"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Type something", description=str(desc), disabled=False, style=style)
        return text
    
    def change_file_box(self, change):
        self.file_name = self.file_text_box.value  
        if os.path.isfile(self.file_name):
            self.got_file_name = True
        
    def create_init_k_box(self, val="0.5", desc="Initial k [s\u207B\u00B9]"):
        """
        Box for inserting the initial exponential decay rate k in 1/s.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Type something", description=str(desc), disabled=False, style=style)
        return text
    
    def create_init_integration_time_box(self, val="0.02", desc="Camera integration time [s]"):
        """
        Box for inserting the camera integration time in s.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="insert path", description=str(desc), disabled=False, style=style)
        return text    
        
    def create_run_button(self):
        button = widgets.Button(
                description="run",
                disabled=False,
                button_style="",
                tooltip="run the analysis")
        return button        

    def create_save_button(self):
        button = widgets.Button(
                description="save",
                disabled=False,
                button_style="",
                tooltip="save the results")
        return button    
    
    def create_save_plot_checkbox(self):
        check_box = widgets.Checkbox(
            value=True,
            description="Save plot",
            disabled=False)
        return check_box
        
    def create_clear_output(self):
        clear_output()
        
    def warning_wrong_file_path(self):
        print("The file path is empty or does not exist.")
        
    def warning_wrong_file(self):
        print("A file with false columns was loaded.")
