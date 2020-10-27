"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Widget handling for MergeHdf5.ipynb.
"""

import tkinter as tk 
from tkinter.filedialog import askopenfilename
from ipywidgets import widgets
from IPython.display import clear_output


class WidgetMergeHdf5():
    def __init__(self, load_path, archive_path, save_path):
        self.load_hdf5_button = self.create_load_hdf5_button()
        self.load_hdf5_box = self.create_load_hdf5_box(val=load_path)
        self.hdf5_file_name = ""
        self.load_archive_button = self.create_load_archive_button()
        self.load_archive_box = self.create_load_archive_box(val=archive_path)
        self.archive_file_name = ""
        self.save_path_button = self.create_save_path_button()
        self.save_path_box = self.create_save_path_box(val=save_path)
        self.save_file_name = ""
        self.merge_button = self.create_merge_button()
        
    def create_load_hdf5_button(self):
        button = widgets.Button(
                description="browse",
                disabled=False,
                button_style="",
                tooltip="browse for file")
        return button
    
    def create_load_hdf5_box(self, val="", desc="Insert path"):
        """
        Box for inserting the path to the *.txt file with fitting h5 file paths from SPTAnalyser.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Type something", description=str(desc), disabled=False, style=style)
        return text
    
    def open_hdf5_file(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.name = askopenfilename(title="Import hdf5 path info", filetypes=(("text files", "*.txt"), ("all files", "*.*")))
        self.hdf5_file_name = root.name
        root.update()
        root.destroy()
        self.load_hdf5_box.value = self.hdf5_file_name
        
    def change_hdf5_box(self, change):
        self.hdf5_file_name = self.load_hdf5_box.value  
    
    def create_load_archive_button(self):
        button = widgets.Button(
                description="browse",
                disabled=False,
                button_style="",
                tooltip="browse for file")
        return button
    
    def create_load_archive_box(self, val="", desc="Insert path"):
        """
        Box for inserting the path to the *.txt file with archive.h5 file paths from ermine.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Type something", description=str(desc), disabled=False, style=style)
        return text

    def create_save_path_button(self):
        button = widgets.Button(
                description="browse",
                disabled=False,
                button_style="",
                tooltip="browse for file")
        return button
    
    def create_save_path_box(self, val="", desc="Insert path"):
        """
        Box for inserting the path to the *.txt file with save paths for merged *.h5 files.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Type something", description=str(desc), disabled=False, style=style)
        return text
    
    def change_save_box(self, change):
        self.save_file_name = self.save_path_box.value  
    
    def open_save_paths(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.name = askopenfilename(title="Import save path info", filetypes=(("text files", "*.txt"), ("all files", "*.*")))
        self.save_file_name = root.name
        root.update()
        root.destroy()
        self.save_path_box.value = self.save_file_name
    
    def open_archive_file(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.name = askopenfilename(title="Import archive path info", filetypes=(("text files", "*.txt"), ("all files", "*.*")))
        self.archive_file_name = root.name
        root.update()
        root.destroy()
        self.load_archive_box.value = self.archive_file_name
        
    def change_archive_box(self, change):
        self.archive_file_name = self.load_archive_box.value  
    
    def create_merge_button(self):
        button = widgets.Button(
                description="merge & save",
                disabled=False,
                button_style="",
                tooltip="merge & save files")
        return button

    def create_clear_output(self):
        clear_output()
