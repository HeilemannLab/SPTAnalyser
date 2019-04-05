# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 14:23:51 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""
        

import tkinter as tk 
from tkinter.filedialog import askopenfilename
from ipywidgets import widgets
from IPython.display import clear_output
#from ..preAnalysis import expDisplacement


class WidgetTrcFormat():
    def __init__(self):
        self.software_button = self.create_software_button()
        self.file_name = ""
        self.got_file_name = False
        self.file_text_box = self.create_file_box()
        self.file_button = self.create_file_button()
        self.pixel_size_box = self.create_pixel_size_box()
        self.min_track_length_box = self.create_min_track_length_box()
        self.run_button = self.create_run_button()
        self.save_button = self.create_save_button()
        
    def create_min_track_length_box(self, val = "2", desc = "Min track length"):
        """
        Box for inserting the minimal track length for tau threshold calculation.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
        
    def create_software_button(self):
        """
        Radiobutton to choose between rapidSTORM and thunderSTORM.
        """
        button = widgets.RadioButtons(
                options = ["thunderSTORM", "rapidSTORM"],
                disabled = False)
        return button
    
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
        
    def open_file(self, b):  # b = ???
        """
        Give the file button opening powers.
        """
        root = tk.Tk()  # window class
        root.withdraw()  # close the window 
        root.update()  # close the window
        if self.software_button.value == "thunderSTORM":
            root.name = askopenfilename(title="Import tracked.seg file", filetypes=(("csv files", "*.csv"),("all files", "*.*")))
        elif self.software_button.value == "rapidSTORM":
            root.name = askopenfilename(title="Import tracked.seg file", filetypes=(("text files", "*.txt"),("all files", "*.*")))
        self.file_name = root.name
        root.update()
        root.destroy()
        self.file_text_box.value=self.file_name
        self.got_file_name = True

    def create_file_box(self, val = "path", desc = "Complete path"):  # val = in box, desc = infront of box
        """
        Box for inserting the path with description, alternative for file loading button.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def change_file_box(self, change):
        self.file_name = self.file_text_box.value  
        self.got_file_name = True
        
    def create_pixel_size_box(self, val = "158", desc = "pixel size in nm"): 
        """
        Box for inserting the pixel size in nm.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
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
        print("This file path does not exist.")
        
    def warning_wrong_file(self):
        print("A file with false columns was loaded.")

    def create_clear_output(self):
        clear_output()
        
        
def main():
    widget_trc_format = WidgetTrcFormat()
    
    
if __name__ == "__main__":
    main()
