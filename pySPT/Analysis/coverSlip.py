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

# =============================================================================
# from pySPT.analysis import cell
# from pySPT.analysis import trajectoryStatistics
# from pySPT.analysis import coverSlip
# =============================================================================

class CoverSlip():
    def __init__(self):
        self.background = []
        self.background_name = [] 
        self.background_log = []
        self.cells = []
        self.cell_trajectories = []
        
    def log_transform(self):
        """
        Create list with log10(D) for each trajectory in each background of list.
        """
        for background in range(0, len(self.background)):
            for trajectory in range(0, len(self.background[background])):
                self.background_log.append(np.log10(self.background[background][trajectory].D))
        print(self.background_log)
        
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

    def create_cells(self, file_names):
        # create list with file names for cells
        for file_name in file_names:
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
            trc_file = np.loadtxt(file_name, usecols = (0, 1, 2, 3, 5)) # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = intensity
            one_cell.trc_file = trc_file
            one_cell.create_trajectories()
            one_cell.analyse_trajectories()
            print("size", one_cell.size)
            cell_trajectories.append(one_cell.analysed_trajectories)
            cells.append(one_cell)
        for i in cells:
            print(i.name, i.size)
            i.plot_trajectory(4)
        

def main():
    pass


if __name__ == "__main__":
    main()
    