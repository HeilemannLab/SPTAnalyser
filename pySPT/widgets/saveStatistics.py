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
        self.grp04 = self.h5_file.create_group("filterInfo")
        self.grp05 = self.h5_file.create_group("diffusionHistogram")
        self.grp06 = self.h5_file.create_group("statistics")

    def groups_bg(self):
        self.grp02 = self.h5_file.create_group("backgroundInfo")
        self.grp03 = self.h5_file.create_group("backgroundCounts")
        
    def cells(self, data):
        my_datatype = np.dtype([("cell name", h5py.special_dtype(vlen=str)),
                                ("cell size [\u03BCm\u00b2]", float)])
        dset = self.grp00.create_dataset("cells", (np.shape(data)[0],), dtype = my_datatype)
        data_array = np.array(data, dtype=my_datatype)
        dset[...] = data_array
        self.h5_file.close()
        
    def backgrounds(self, data):
        my_datatype = np.dtype([("background name", h5py.special_dtype(vlen=str)),
                                ("background size [\u03BCm\u00b2]", float)])
        dset = self.grp02.create_dataset("backgrounds", (np.shape(data)[0],), dtype = my_datatype)
        data_array = np.array(data, dtype=my_datatype)
        dset[...] = data_array
        
        
    def cell_counts(self, cell_name, number, diffusion_coeff, counts):
        dset = self.grp01.create_dataset(cell_name, (number, ), dtype = np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                         ("counts / area [1/\u03BCm\u00b2]", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion_coeff
        dset["counts / area [1/\u03BCm\u00b2]"] = counts
        
    def bg_counts(self, bg_name, number, diffusion_coeff, counts):
        dset = self.grp03.create_dataset(bg_name, (number, ), dtype = np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                         ("counts / area [1/\u03BCm\u00b2]", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion_coeff
        dset["counts / area [1/\u03BCm\u00b2]"] = counts
        
    def filter_info(self, filter_settings,filter_thresholds_values):
        dset = self.grp04.create_dataset("filters", (1,1), dtype = np.dtype([("filter immobile", int),
                                                         ("filter confined", int),
                                                         ("filter free", int),
                                                         ("type determination not successful", int),
                                                         ("min trajectory length [frames]", int),
                                                        ("max trajectory length [frames]", int),
                                                        ("min diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                         ("max diffusion coefficient [\u03BCm\u00b2/s]", float)]))
        dset["filter immobile"] = filter_settings[0]
        dset["filter confined"] = filter_settings[1]
        dset["filter free"] = filter_settings[2]
        dset["type determination not successful"] = filter_settings[3]
        dset["min trajectory length [frames]"] = filter_thresholds_values[0]
        dset["max trajectory length [frames]"] = filter_thresholds_values[1]
        dset["min diffusion coefficient [\u03BCm\u00b2/s]"] = filter_thresholds_values[2]
        dset["max diffusion coefficient [\u03BCm\u00b2/s]"] = filter_thresholds_values[3]
        
    def statistics(self, immobile, confined, free, notype, trajectories_included, trajectories_excluded, mean_Ds, mean_dDs, mean_lengths, mean_dlengths):
        dset = self.grp06.create_dataset("statistics", (1,1), dtype = np.dtype([("immobile [%]", float),
                                                         ("confined [%]", float),
                                                         ("free [%]", float),
                                                         ("no type [%]", float),
                                                         ("trajectories included", int),
                                                         ("trajectories excluded", int),
                                                         ("mean D immobile [\u03BCm\u00b2/s]", float),
                                                         ("mean D confined [\u03BCm\u00b2/s]", float),
                                                         ("mean D free [\u03BCm\u00b2/s]", float),
                                                         ("mean D no type [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D immobile [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D confined [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D free [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D no type [\u03BCm\u00b2/s]", float),
                                                         ("mean length immobile [frames]", float),
                                                         ("mean length confined [frames]", float),
                                                         ("mean length free [frames]", float),
                                                         ("mean length no type [frames]", float),
                                                         ("\u0394 mean length immobile [frames]", float),
                                                         ("\u0394 mean length confined [frames]", float),
                                                         ("\u0394 mean length free [frames]", float),
                                                         ("\u0394 mean length no type [frames]", float)]))

        dset["immobile [%]"] = immobile
        dset["confined [%]"] = confined
        dset["free [%]"] = free
        dset["no type [%]"] = notype
        dset["trajectories included"] = trajectories_included
        dset["trajectories excluded"] = trajectories_excluded
        dset["mean D immobile [\u03BCm\u00b2/s]"] = mean_Ds[0]
        dset["mean D confined [\u03BCm\u00b2/s]"] = mean_Ds[1]
        dset["mean D free [\u03BCm\u00b2/s]"] = mean_Ds[2]
        dset["mean D no type [\u03BCm\u00b2/s]"] = mean_Ds[3]
        dset["\u0394 mean D immobile [\u03BCm\u00b2/s]"] = mean_dDs[0]
        dset["\u0394 mean D confined [\u03BCm\u00b2/s]"] = mean_dDs[1]
        dset["\u0394 mean D free [\u03BCm\u00b2/s]"] = mean_dDs[2]
        dset["\u0394 mean D no type [\u03BCm\u00b2/s]"] = mean_dDs[3]
        dset["mean length immobile [frames]"] = mean_lengths[0]
        dset["mean length confined [frames]"] = mean_lengths[1]
        dset["mean length free [frames]"] = mean_lengths[2]
        dset["mean length no type [frames]"] = mean_lengths[3]
        dset["\u0394 mean length immobile [frames]"] = mean_dlengths[0]
        dset["\u0394 mean length confined [frames]"] = mean_dlengths[1]
        dset["\u0394 mean length free [frames]"] = mean_dlengths[2]
        dset["\u0394 mean length no type [frames]"] = mean_dlengths[3]
        
    def diffusion_plot(self, number, diffusion, mean_cell, dmean_cell):
        dset = self.grp05.create_dataset("histogram values", (number,), dtype = np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                 ("mean frequency cells", float),
                                                 ("\u0394 mean frequency cells", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion    
        dset["mean frequency cells"] = mean_cell
        dset["\u0394 mean frequency cells"] = dmean_cell
        
    def diffusion_plot_bg(self, number, diffusion, mean_cell, dmean_cell, mean_bg, dmean_bg, mean_cell_corr, dmean_cell_corr):
        dset = self.grp05.create_dataset("histogram values", (number,), dtype = np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                              ("mean frequency cells", float),
                                                              ("\u0394 mean frequency cells", float),
                                                              ("mean frequency background", float),
                                                              ("\u0394 mean frequency background", float),
                                                              ("mean frequency corrected", float),
                                                              ("\u0394 mean frequency corrected", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion
        dset["mean frequency cells"] = mean_cell
        dset["\u0394 mean frequency cells"] = dmean_cell
        dset["mean frequency background"] = mean_bg
        dset["\u0394 mean frequency background"] = dmean_bg
        dset["mean frequency corrected"] = mean_cell_corr
        dset["\u0394 mean frequency corrected"] = dmean_cell_corr
    
    def diffusion_plot_normalized(self, number, diffusion, mean_cell_percent, dmean_cell_percent):
        """
        Normalized values if bg is inserted."
        """
        dset = self.grp05.create_dataset("histogram values normalized", (number,), dtype = np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                 ("mean frequency cells [%]", float),
                                                 ("\u0394 mean frequency cells [%]", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion    
        dset["mean frequency cells [%]"] = mean_cell_percent
        dset["\u0394 mean frequency cells [%]"] = dmean_cell_percent
        
    def diffusion_plot_bg_normalized(self, number, diffusion, mean_cell_percent, dmean_cell_percent, mean_cell_corr_percent, dmean_cell_corr_percent):
        """
        Diffusion plot without bg correction, normalized.
        """
        dset = self.grp05.create_dataset("histogram values normalized", (number,), dtype = np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                 ("mean frequency cells [%]", float),
                                                 ("\u0394 mean frequency cells [%]", float),
                                                 ("mean frequency corrected [%]", float),
                                                 ("\u0394 mean frequency corrected [%]", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion
        dset["mean frequency cells [%]"] = mean_cell_percent
        dset["\u0394 mean frequency cells [%]"] = dmean_cell_percent
        dset["mean frequency corrected [%]"] = mean_cell_corr_percent
        dset["\u0394 mean frequency corrected [%]"] = dmean_cell_corr_percent
    
    def diffusion_bin_size(self, bin_size):
        dset = self.grp05.create_dataset("bin size", (1,1), dtype = np.dtype([("bin size", float)]))
        dset["bin size"] = bin_size
        
        


