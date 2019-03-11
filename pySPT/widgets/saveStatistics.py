# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 14:19:23 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

TrackStatistics: Save statistics as hdf5
"""

import h5py
import tkinter as tk
import os
import numpy as np

class SaveStatistics():
    def __init__(self):
        self.h5_file = []  # h5 file
        self.grp00 = []  # groups for structure
        self.grp01 = []
        self.grp02 = []
        self.grp03 = []
        self.grp04 = []
        self.trc_file_hdf5 = ""  # path of file with .h5 ending
        self.D_length = 0  # col length from log hist D
    
    def create_h5(self, path):
        self.create_h5_name(path)
        self.create_h5_file()
        self.groups()
        
    def create_h5_file(self):
        self.h5_file = h5py.File(self.trc_file_hdf5, "w")  # w- or x = Create file, fail if exists

    def create_h5_name(self, path):     
        # splitext -> tupel with path split from .* ending. It splits at the last dot in name.
        self.trc_file_hdf5 = os.path.splitext(path)[0] + ".h5"    
        
    def groups(self):
        self.grp00 = self.h5_file.create_group("cellInfo")
        self.grp01 = self.h5_file.create_group("cellCounts")
        self.grp02 = self.h5_file.create_group("backgroundInfo")
        self.grp03 = self.h5_file.create_group("backgroundCounts")
        self.grp04 = self.h5_file.create_group("relativeFrequencies")
        
    def cells(self, data):
        my_datatype = np.dtype([("cell name", h5py.special_dtype(vlen=str)),
                                ("cell size [\u03BCm\u00b2]", float)])
        dset = self.grp00.create_dataset("cells", (np.shape(data)[0],), dtype = my_datatype)
        data_array = np.array(data, dtype=my_datatype)
        dset[...] = data_array
        
        
    def backgrounds(self, data):
        my_datatype = np.dtype([("background name", h5py.special_dtype(vlen=str)),
                                ("background size [\u03BCm\u00b2]", float)])
        dset = self.grp02.create_dataset("backgrounds", (np.shape(data)[0],), dtype = my_datatype)
        data_array = np.array(data, dtype=my_datatype)
        dset[...] = data_array
        self.h5_file.close()
        


