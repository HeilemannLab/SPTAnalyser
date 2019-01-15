# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 15:40:08 2019

@author: Johanna Rahm, Research Group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt am Main.
"""

import os
import tkinter as tk 
import datetime
from tkinter.filedialog import askopenfilename
from ipywidgets import widgets
from IPython.display import display

class WidgetExpDisp():
    def __init__(self):
        self.file_name = ""
        self.folder_name = ""
        self.base_name = ""
        self.file_name_text = self.create_text_str()
        self.file_dialog_button = self.create_file_dialog()
        
    def create_text_str(self, val = "path", desc = "Complete path"):  # val = in box, desc = infront of box
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def create_file_dialog(self):
        button = widgets.Button(
                description='browse',
                disabled=False,
                button_style='', # 'success', 'info', 'warning', 'danger' or ''
                tooltip='browse for file')
                #icon='check')
        return button
    
    def print_file_name(self):
        #print(self.file_name_text.value)
        display(self.file_name_text.value)  # testing purposes, rather display than print
        
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
        self.file_name_text.value=self.file_name
        

def main():
    widget_exp_disp = WidgetExpDisp()
    widget_exp_disp.open_file()
    print(widget_exp_disp.file_name, widget_exp_disp.folder_name, widget_exp_disp.base_name, sep="\n")
    
    
if __name__ == "__main__":
    main()
    
    