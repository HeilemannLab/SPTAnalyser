# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:47:24 2019

@author: Johanna Rahm

Save hmmVis statistics & plot information .h5.
"""

import h5py
import tkinter as tk
import os
import numpy as np

class SaveHmmVis():
    def __init__(self, directory, folder_name):
        self.save_path = directory + "\\" + folder_name + "\\hmm_vis.h5"
        self.h5_file = self.create_h5_file()
        
    def create_h5_file(self):
        return h5py.File(self.save_path, "w")  # w- or x = Create file, fail if exists
    
    def groups(self):
        self.grp00 = self.h5_file.create_group("cellInfo")  # cell names, aic value, cell size, localizations, density
        self.grp01 = self.h5_file.create_group("statistics")  # mean aic, mean %, mean tp, mean D
        self.grp02 = self.h5_file.create_group("states")  # single Ds and populations
        self.grp03 = self.h5_file.create_group("transitionProbabilities")  # single tps
        self.grp04 = self.h5_file.create_group("stateTransitionDiagram")  # infos about the state transition diagram
        
    def cell_names(self, data):
        my_datatype = np.dtype([("cell name", h5py.special_dtype(vlen=str))])
        dset = self.grp00.create_dataset("cells", (np.shape(data)[0],), dtype = my_datatype)
        data_array = np.array(data, dtype=my_datatype)
        dset[...] = data_array
        
    def cell_info(self, data):
        my_datatype = np.dtype([("cell name", h5py.special_dtype(vlen=str)),
                                ("cell size [\u03BCm\u00b2]", float),
                                ("localizations", int),
                                ("density [locs/cell size]", float),
                                ("aic value", float),
                                ("bic value", float),
                                ("log likelihood", float)])
        dset = self.grp00.create_dataset("cell infos", (np.shape(data)[0],), dtype = my_datatype)
        data_array = np.array(data, dtype=my_datatype)
        dset[...] = data_array        
        
# =============================================================================
#     def mean_aic_value(self, mean_aic):
#         dset = self.grp01.create_dataset("mean AIC value", (1,1), dtype = np.dtype([("mean AIC value", float)]))
#         dset["mean AIC value"] = mean_aic
#         
# =============================================================================
    
        
    def mean_states(self, population, dpopulation, D_coeff, dD_coeff):
        dset = self.grp01.create_dataset("mean states", (np.shape(population)), dtype = np.dtype([("mean population", float),
                                                             ("\u0394 mean population", float),
                                                             ("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                             ("\u0394 diffusion coefficient [\u03BCm\u00b2/s]", float)]))
        dset["mean population"] = population
        dset["\u0394 mean population"] = dpopulation
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = D_coeff
        dset["\u0394 diffusion coefficient [\u03BCm\u00b2/s]"] = dD_coeff
        
    def mean_tps(self, tps):
        dset = self.grp01.create_dataset("mean transition probabilities", (np.shape(tps)))
        dset[...] = tps
        
    def single_states(self, D, cell_name, population):
        dset = self.grp02.create_dataset(cell_name, np.shape(D), dtype = np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                             ("population", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = D
        dset["population"] = population
        
    def single_tps(self, tp, cell_name):
        dset = self.grp03.create_dataset(cell_name, np.shape(tp))
        dset[...] = tp
        
    def edge_sizes(self, edge_sizes):
        dset = self.grp04.create_dataset("edge sizes", (np.shape(edge_sizes)))
        dset[...] = edge_sizes
    
    def node_sizes(self, node_sizes):
        dset = self.grp04.create_dataset("node sizes", (np.shape(node_sizes)))
        dset[...] = node_sizes
        
    def dmean_tps(self, dtps):
        dset = self.grp01.create_dataset("\u0394 mean transition probabilities", (np.shape(dtps)))
        dset[...] = dtps
        self.h5_file.close()
        

        
        
        
