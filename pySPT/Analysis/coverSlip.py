# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 10:01:02 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import numpy as np
from . import cell
import time
import os
from tqdm import tqdm_notebook as tqdm
from . import trcFormat

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
        self.dt = 0.0   # hand down to cell -> trajectory
        self.tau_threshold = 0.0  # hand down to cell -> trajectory
        self.dof = 0.0  # hand down to cell -> trajectory
        self.D_min = 0.0  # hand down to cell -> trajectory
        self.points_fit_D = 4  # hand down to cell -> trajectory
        self.seg_id = True  # hand down to cell
        self.software = ""
        self.min_track_length_type = 0
        self.min_track_length_hmm = 0
        self.column_orders = []
        
    def calc_tau_threshold(self):
        self.tau_threshold = float(self.min_track_length_type)*float(self.dt)*0.6*0.5
        
    def calc_min_track_lengths(self):
        """
        The min track length has to be points_fit_D+1
        """
        if int(self.min_track_length_hmm) <= int(self.points_fit_D):
            self.min_track_length_hmm = int(self.points_fit_D) + 1
        if int(self.min_track_length_type) <= int(self.points_fit_D):
            self.min_track_length_type = int(self.points_fit_D) + 1
    
    def seg_id_boolean(self):
        """
        The return of the radio button is a string, convert it to boolean.
        """
        if self.seg_id == "seg id":
            self.seg_id = True
        elif self.seg_id == "track id":
            self.seg_id = False
            
    def check_PT_trajectory(self):
        """
        If PALMTracer is used as software, seg_id will automatically be set to False, because no segments exist.
        """
        if self.software == "PALMTracer":
            self.seg_id = False

    def create_cells(self):
        """
        Initialize cell analysis, also load roi file for cell size. Cell & trajectory
        objects from e.g. CoverSlip Class.
        """        
        start = time.time()
        self.seg_id_boolean()
        self.check_PT_trajectory()
        self.calc_min_track_lengths()
        self.calc_tau_threshold()
        #print("min track hmm, min track type, tau threshold", self.min_track_length_hmm, self.min_track_length_type, self.tau_threshold)
        if self.roi_file:  # if no roi file path was inserted, no file can be loaded  
            roi_file = np.genfromtxt(self.roi_file, dtype=None, delimiter=",", skip_header=3, encoding=None)
        for file_name in tqdm(self.cell_files):
            cell_idx = self.cell_files.index(file_name)
            one_cell = cell.Cell()  # create cell object for each file
            base_name = os.path.basename(file_name)
            raw_base_name = ""
            for i in base_name:
                if i == ".":
                    break
                else:
                    raw_base_name += i   
            one_cell.name = raw_base_name
            if self.roi_file:  # if the raw base name of the cell equals the raw roi_file [0] entry, the cell gets its size
                for file in roi_file:  # raw_cell_name = cell01, roi_file [0] = cell01.tracked.csv, cell01.tracked, cell01 -> all would be ok
                    raw_file = ""  
                    for i in file[0]:
                        if i == ".":
                            break
                        else:
                            raw_file += i
                    raw_file = raw_file.replace('"', "")
                    if one_cell.name == raw_file:
                        print("roi check", file)
                        one_cell.size = file[1]     
            # in PT the column order is set and not necessary.
            if self.software == "PALMTracer":
                trc_format = trcFormat.TrcFormat(self.software, file_name, self.pixel_size, self.min_track_length_type,
                                                 self.min_track_length_hmm, self.seg_id)
            else:
                trc_format = trcFormat.TrcFormat(self.software, file_name, self.pixel_size, self.min_track_length_type,
                                 self.min_track_length_hmm, self.seg_id, column_order=self.column_orders[cell_idx])
            trc_format.run()
        
# =============================================================================
#             # testing purpose
#             print("cs pixel size", self.pixel_size)   # [nm], hand down to cell
#             print("px amount", self.pixel_amount)  # amount of pixels of detector (eg. 256*256), hand down to cell 
#             print("dt", self.dt)   # hand down to cell -> trajectory
#             print("tau threshold", self.tau_threshold)  # hand down to cell -> trajectory
#             print("dof", self.dof )  # hand down to cell -> trajectory
#             print("dmin", self.D_min)  # hand down to cell -> trajectory
#             print("points to fit", self.points_fit_D)  # hand down to cell -> trajectory
#             print("seg id", self.seg_id)  # hand down to cell
#             print("software", self.software) 
#             print("min track length type", self.min_track_length_type)
#             print("min track length hmm", self.min_track_length_hmm)
#             print("column orders", self.column_orders)
# =============================================================================
            
            trc_file_type = trc_format.trc_file_type_filtered
            trc_file_hmm = trc_format.trc_file_hmm_filtered
            if trc_file_type and trc_file_hmm:
                one_cell.trc_file_type = trc_file_type
                one_cell.trc_file_hmm = trc_file_hmm
                one_cell.seg_id = self.seg_id
                one_cell.min_track_length_type = self.min_track_length_type
                one_cell.min_track_length_hmm = self.min_track_length_hmm
                one_cell.tau_threshold = float(self.tau_threshold)
                one_cell.dt = float(self.dt)
                one_cell.pixel_amount = float(self.pixel_amount)
                one_cell.pixel_size = float(self.pixel_size)
                one_cell.dof = float(self.dof)
                one_cell.D_min = float(self.D_min)
                one_cell.points_fit_D = int(self.points_fit_D)
                one_cell.run_analysis()
                self.cell_trajectories.append(one_cell.analysed_trajectories)
                self.cells.append(one_cell)
            else:
                print("All trajecoties are shorter as the minimum trajectory length inserted, please select a smaller minimum threshold.")
            print("Analysis took {} s".format(time.time()-start))
            
    def plot_trajectory(self, cell_name, trajectory_idx):
        """
        Plot a trajectory.
        :param cell_name: name of cell, index will be created to cell name.
        :param trajectory: number of trajectory -> index-1.
        """
        for cell in self.cells:
            if cell.name == cell_name:
                cell_index = self.cells.index(cell)
                for trajectory in cell.analysed_trajectories:
                    if trajectory_idx == trajectory.trajectory_number:
                        target_trajectory = cell.analysed_trajectories.index(trajectory)
        self.cell_trajectories[cell_index][target_trajectory].plot_particle()

        

def main():
    pass


if __name__ == "__main__":
    main()
    