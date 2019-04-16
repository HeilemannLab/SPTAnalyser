# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 10:01:02 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import numpy as np
from . import cell
import time
import os
from tqdm import tqdm_notebook as tqdm

# =============================================================================
# from pySPT.analysis import cell
# from pySPT.analysis import trajectoryStatistics
# from pySPT.analysis import coverSlip
# =============================================================================

class CoverSlip():
    def __init__(self):
        self.cells = []  # cell objects
        self.cell_trajectories = []  # lists of trajectory object, 1 list = 1 cell
        self.roi_file = []  # full path to roi file
        self.cell_files = []  # list with full paths to cell files
        self.backgrounds = []  # backgroud objects
        self.background_trajectories = []  # lists of trajectory objects, 1 list = 1 background
        self.background_files = []  # list with full paths of bg files
        # initialization parameters for cell (only from trc -> h5)
        self.pixel_size = 0.0   # [nm], hand down to cell
        self.pixel_amount = 0.0  # amount of pixels of detector (eg. 256*256), hand down to cell
        self.tau_threshold_min_length = 0.0  
        self.dt = 0.0   # hand down to cell -> trajectory
        self.tau_threshold = 0.0  # hand down to cell -> trajectory
        self.dof = 0.0  # hand down to cell -> trajectory
        self.D_min = 0.0  # hand down to cell -> trajectory

    def calc_tau_threshold(self):
        self.tau_threshold = float(self.tau_threshold_min_length)*float(self.dt)*0.6*0.5

    def create_cells(self):
        """
        Initialize cell analysis, also load roi file for cell size. Cell & trajectory
        objects from e.g. CoverSlip Class.
        """
        start = time.time()
        self.calc_tau_threshold()
        if self.roi_file:  # if no roi file path was inserted, no file can be loaded  
            roi_file = np.genfromtxt(self.roi_file, dtype=None, delimiter=",", skip_header=3, encoding=None)
        for file_name in tqdm(self.cell_files):
            one_cell = cell.Cell()  # create cell object for each file
            base_name = os.path.basename(file_name)
            raw_base_name = ""
            for i in base_name:
                if i == ".":
                    break
                else:
                    raw_base_name += i   
            one_cell.name = raw_base_name
            if self.roi_file:
                for file in roi_file:
                    if one_cell.name in file[0]:
                        one_cell.size = file[1]            
            trc_file = np.loadtxt(file_name, usecols = (0, 1, 2, 3, 4, 5))  # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = place holder (int), col 5 = intensity
            #self.trc_file = trc_file ????????????
            one_cell.trc_file = trc_file
            one_cell.tau_threshold_min_length = self.tau_threshold_min_length
            one_cell.tau_threshold = float(self.tau_threshold)
            one_cell.dt = float(self.dt)
            one_cell.pixel_amount = float(self.pixel_amount)
            one_cell.pixel_size = float(self.pixel_size)
            one_cell.dof = float(self.dof)
            one_cell.D_min = float(self.D_min)
            one_cell.create_trajectories()
            one_cell.cell_size()
            one_cell.analyse_trajectories()
            self.cell_trajectories.append(one_cell.analysed_trajectories)
            self.cells.append(one_cell)
        print("Analysis took {} s".format(time.time()-start))
        
    def plot_trajectory(self, cell_name, trajectory):
        """
        Plot a trajectory.
        :param cell_name: name of cell, index will be created to cell name.
        :param trajectory: number of trajectory -> index-1.
        """
        for cell in self.cells:
            if cell.name == cell_name:
                cell_index = self.cells.index(cell)
        self.cell_trajectories[cell_index][int(trajectory)-1].plot_particle()

        

def main():
    pass


if __name__ == "__main__":
    main()
    