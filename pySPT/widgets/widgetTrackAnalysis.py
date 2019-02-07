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
        
    def searchSubFolders(self, dirName):
        if (dirName):
            self.data_set_dir = dirName
            for root, dirs, files in os.walk(self.data_set_dir):
                self.extendList(root, files)
                
    def extendList(self, root, files):
        for name in files:
            if name.endswith(self.suffix) and "background" not in name:
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
        
    def create_dir_box(self, val = "directory to be searched in", desc = "directory"):  # val = in box, desc = infront of box
        """
        Box for inserting the directory with description, alternative for dir loading button.
        """
        style = {'description_width': 'initial'}  # display too long desc
        text = widgets.Text(value=str(val), placeholder='Type something', description=str(desc), disabled=False, style = style)
        return text
    
    def change_dir_box(self, change):
        self.dir_name = self.dir_box.value  
        self.got_dir = True

    def create_roi_box(self, val = "path of roi", desc = "roi"):
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
        self.roi_box.value=self.dir_name
        self.got_roi = True      
      
    def printFileNames(self):
        for name in self.fileNames:
            print("Found .trc file: %s" %(name))
            
def main():
    communicator = Mol2Judi()
    communicator.searchSubFolders("E:/Receptor-signaling/met/rtk/resting")
    communicator.printFileNames()
    
if __name__ == '__main__':
    main()
    