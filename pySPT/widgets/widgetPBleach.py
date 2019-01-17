# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 17:14:19 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import os
import tkinter as tk 
from tkinter.filedialog import askopenfilename
from ipywidgets import widgets
from IPython.display import display
from IPython.display import clear_output

class WidgetPBleach():
    def __init__(self):
        self.file_name = ""
        self.folder_name = ""
        self.base_name = ""
        self.file_name_text = self.create_text_str()
        self.file_dialog_button = self.create_file_dialog()
        self.run_button = self.create_run_button()
        self.save_button = self.create_save_button()
        self.clear_output = self.create_clear_output()
        
        
    def create_text_str(self, val = "path", desc = "Complete path"):  # val = in box, desc = infront of box
        """
        Box for inserting the path with description.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_file_dialog(self):
        """
        Loading button.
        """
        button = widgets.Button(
                description='browse',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='browse for file')
                #icon='check')
        return button
    
    def open_file(self, b):  # b = ???
        root = tk.Tk()  # window class
        root.withdraw()  # close the window 
        root.update()  # close the window
        root.name = askopenfilename(title="Import tracked.seg file", filetypes=(("text files", "*.txt"),("all files", "*.*")))
        self.file_name = root.name
        self.folder_name = os.path.dirname(self.file_name)  # returns the head of path /foo/bar/item the directory name of path = item -> return /foo/bar
        self.base_name = os.path.basename(self.file_name)[:-4]  # returns the tail of path -> item, [:-4] delete the ending .txt 
        root.update()
        root.destroy()
        self.file_name_text.value=self.file_name  # inserts path in box for inserting the path
        
    def create_run_button(self):
        button = widgets.Button(
                description='run',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='run the analysis')
                #icon='check')
        return button        
    
    def create_clear_output(self):
        clear_output()
    
    def create_save_button(self):
        pass
    

        
        
def main():
    pass
        
if __name__ == "__main__":
    main()
    
    