# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 16:52:13 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import os

class WidgetDirStructure():
    def __init__(self):
        self.file_name = ""
        self.folder_name = ""
        self.base_name = ""
        self.pre_analysis_dir = ""
        
    def name_handling(self, file_name):
        self.file_name = file_name
        self.folder_name = os.path.dirname(self.file_name)  # returns the head of path /foo/bar/item the directory name of path = item -> return /foo/bar
        self.base_name = os.path.basename(self.file_name)#[:-4] # returns the tail of path -> item, [:-4] delete the ending .seg    
        #print(self.file_name, self.folder_name, self.base_name, sep="\n")
        
    def create_folder(self):
        self.pySPT_dir = self.folder_name + "\pySPT_" + self.base_name[:-20]
        if not os.path.exists(self.pySPT_dir):
            os.makedirs(self.pySPT_dir)
        self.pre_analysis_dir = self.pySPT_dir + "\preAnalysis"
        if not os.path.exists(self.pre_analysis_dir):
            os.makedirs(self.pre_analysis_dir)
            
        
def main():
    widget_dir_structure = WidgetDirStructure()
    file_name = "F:/Marburg/single_colour_tracking/resting/160404_CS5_Cell1/cell_1_MMStack_Pos0.ome.tif.tracked.seg"
    widget_dir_structure.name_handling(file_name)
    widget_dir_structure.create_folder()
        
if __name__ == "__main__":
    main()
    
    