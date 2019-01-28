# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 14:14:32 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Convert a tracked rapidSTORM file into .trc format (PALMTracer)
"""


import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import datetime


class TrcFormat():
    def __init__(self):
        self.file_name = ""
        self.loaded_file = []
        self.trc_file = []
        self.trc_file_sorted = []
        
    def load_localization_file(self):
        """
        Array with columns:
        col0 = sed_id (seg_id = track_id if segmentation is False)
        col1 = frame
        col2 = x
        col3 = y
        col4 = intensity
        """
        seg_id_index = list(self.column_order.keys())[list(self.column_order.values()).index('"seg_id"')]
        self.loaded_file = np.loadtxt(self.file_name, usecols = (10, 4, 0, 2, 5)) 


            
    def create_trc_file(self):
        self.trc_file = np.zeros([np.size(self.loaded_file[:,0]),6])
        
        print(self.loaded_file[:,0])
        #print(self.trc_file, len(self.trc_file[0]))
        seg_id = self.loaded_file[:,0]
        frame = self.loaded_file[:,1]
        position_x = self.loaded_file[:,2]
        position_y = self.loaded_file[:,3]
        intensity = self.loaded_file[:,4]
        self.trc_file[:,0] = seg_id
        self.trc_file[:,1] = frame
        self.trc_file[:,2] = position_x
        self.trc_file[:,3] = position_y
        self.trc_file[:,5] = intensity
        print(self.trc_file)
        
    def sort_trc_file(self):
        dtype = [("seg_id", int), ("frame", int), ("position_x", float), ("position_y", float), ("0", int), ("intensity", float)]
        values = list(map(tuple,self.trc_file))  # convert np.ndarrays to tuples
        structured_array = np.array(values, dtype=dtype)  # create structured array
        self.trc_file_sorted = np.sort(structured_array, order=["seg_id", "frame"])  # sort by dtype name
        
    def save_trc_file(self):
        out_file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\pySPT_cell_1_MMStack_Pos0\\preAnalysis\\sorted.txt"
        np.savetxt(out_file_name, 
                   X=self.trc_file_sorted,
                   fmt = ("%i","%i", "%.3f", "%.3f", "%i", "%.3f"))

        
def main():
    trc_format= TrcFormat()
    trc_format.file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.tif.tracked.txt"
    trc_format.load_localization_file()
    trc_format.create_trc_file()
    trc_format.sort_trc_file()
    trc_format.save_trc_file()

if __name__ == "__main__":
    main()
    