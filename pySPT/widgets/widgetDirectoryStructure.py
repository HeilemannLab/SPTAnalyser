"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Create a folder and subfolders in same directory as analysed file to structure output files.
SPTAnalyser + filename -> preAnalysis, analysis, hmm.
"""

import os


class WidgetDirStructure():
    def __init__(self):
        self.file_name = ""  # full name (path)
        self.folder_name = ""  # path to folder
        self.base_name = ""  # tail of path
        self.raw_base_name = ""  # tail of path without .endings
        self.sub_folder = ""  # name of desired subfolder (\\preAnalysis, \\analysis ...)
        self.sub_folder_dir = ""  # directory of subfolder

    def name_handling(self, file_name):
        self.file_name = file_name
        # returns the head of path /foo/bar/item the directory name of path = item -> return /foo/bar
        self.folder_name = os.path.dirname(self.file_name)
        # returns the tail of path -> item
        self.base_name = os.path.basename(self.file_name)
        
    def create_raw_base_name(self):
        """
        Delete all .xxx.yyy.zzz in a base_name.
        Example: cell_1_MMStack_Pos0.ome.tif.tracked.seg.txt -> cell_1_MMStack_Pos0
        """
        self.raw_base_name = ""
        for i in self.base_name:
            if i == ".":
                break
            else:
                self.raw_base_name += i  
            
    def create_folder(self):
        """
        Create folder SPTAnalyser + raw base name in same directory as target file.
        Create sub folder with chosen name in SPTAnalyser folder.
        """
        self.pySPT_dir = self.folder_name + "\SPTAnalyser_" + self.raw_base_name
        if not os.path.exists(self.pySPT_dir):
            os.makedirs(self.pySPT_dir)
        self.sub_folder_dir = self.pySPT_dir + self.sub_folder
        if not os.path.exists(self.sub_folder_dir):
            os.makedirs(self.sub_folder_dir)
