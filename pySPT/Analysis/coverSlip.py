# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 10:01:02 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import numpy as np
from . import trajectory
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
        #self.background = []
        #self.background_name = [] 
        #self.background_log = []
        self.cells = []  # cell objects
        self.cell_trajectories = []  # lists of trajectory object, 1 list = 1 cell
        self.roi_file = []  # full path to roi file
        self.cell_files = []  # list with full paths to cell files
        #self.trc_file= []   # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = place holder (int), col 5 = intensity ??????
        
# =============================================================================
#     def log_transform(self):
#         """
#         Create list with log10(D) for each trajectory in each background of list.
#         """
#         for background in range(0, len(self.background)):
#             for trajectory in range(0, len(self.background[background])):
#                 self.background_log.append(np.log10(self.background[background][trajectory].D))
#         print(self.background_log)
# =============================================================================
        
# =============================================================================
#     def count_background(self):
#         """
#         Create histogram with bins = mjd_n and frequencies as np.ndarray.
#         Multiply the bins with camera integration time -> [s].
#         """
#         
#         #for background in self.background:
#             
#         max_bin = self.background_log.max()  # max mjd_n value
#         bin_size = int(max_bin)  # divides the bin range in sizes -> desired bin = max_bin/bin_size
#         hist = np.histogram(self.background_log,
#                             range = (-5, max_bin),
#                             bins = bin_size,
#                             density = True)
#         self.background_histogram = np.zeros([np.size(hist[0]),5])
#         self.background_histogram[:,0] = hist[1][:-1] # col0 = frames
#         self.background_histogram[:,1] = hist[0][:]  # col1 = frequencies
#         #self.normalized_mjd_ns()  # normalize the histogram by the sum
# =============================================================================


    def create_cells(self):
        """
        Initialize cell analysis, also load roi file for cell size. Cell & trajectory
        objects from e.g. CoverSlip Class.
        """
        start = time.time()
        self.cell_trajectories = []  # contains trajectories for each cell in separate lists
        self.cells = []  # contains cell objects
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
            for file in roi_file:
                if raw_base_name in file[0]:
                    one_cell.size = file[1]            
            trc_file = np.loadtxt(file_name, usecols = (0, 1, 2, 3, 4, 5))  # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = place holder (int), col 5 = intensity
            #self.trc_file = trc_file ????????????
            one_cell.trc_file = trc_file
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
    