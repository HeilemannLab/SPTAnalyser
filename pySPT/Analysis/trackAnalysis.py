# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 09:17:05 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Helping class for same called Jupyter Notebook, create cells and initialize analysis,
helping functions for .h5 saving (gathering informations from analyse classes).
"""

import os
import numpy as np
from . import trajectory
from . import cell

class TrackAnalysis():
    def __init__(self):
        #self.cell_trajectories = []
        #self.cells = []
        #self.file_names = []
        self.diffusion_info = []
        self.number_of_trajectories = 0
        self.rossier_info = []
        self.diff_plot = []
        self.rossier_plot = []
        
# =============================================================================
#     def create_cell_trajectories(self):
#         """
#         Initialize cell analysis, also load roi file for cell size. Cell & trajectory
#         objects from e.g. CoverSlip Class.
#         """
#         self.cell_trajectories = []  # contains trajectories for each cell in separate lists
#         self.cells = []  # contains cell objects
#         for file_name in self.file_names:
#             one_cell = cell.Cell()  # create cell object for each file
#             base_name = os.path.basename(file_name)
#             raw_base_name = ""
#             for i in base_name:
#                 if i == ".":
#                     break
#                 else:
#                     raw_base_name += i   
#             one_cell.name = raw_base_name
#             for file in roi_file:
#                 if raw_base_name in file[0]:
#                     one_cell.size = file[1]            
#             trc_file = np.loadtxt(file_name, usecols = (0, 1, 2, 3, 5)) # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = intensity
#             one_cell.trc_file = trc_file
#             one_cell.create_trajectories()
#             one_cell.analyse_trajectories()
#             print("size", one_cell.size)
#             self.cell_trajectories.append(one_cell.analysed_trajectories)
#             self.cells.append(one_cell)
#         for i in self.cells:
#             print(i.name, i.size)
#             i.plot_trajectory(4)
# =============================================================================
            
    def save_diff(self, trajectories):
        """
        diffusion info: col0 = id, col1= D, col2 = dD, col3 = MSD(0), col4 = chi2, col5 = length
        """
        self.diffusion_info =  np.zeros([np.size(trajectories), 6])
        self.number_of_trajectories = np.shape(self.diffusion_info)[0]
        i = 0
        for trajectory in trajectories:
            self.diffusion_info[i,0] = trajectory.trajectory_number 
            self.diffusion_info[i,1] = trajectory.D
            self.diffusion_info[i,2] = trajectory.dD
            self.diffusion_info[i,3] = trajectory.MSD_0
            self.diffusion_info[i,4] = trajectory.chi_D
            self.diffusion_info[i,5] = trajectory.length_trajectory
            i += 1
    
    def save_plots(self, trajectory):
        trajectory_number = trajectory.trajectory_number
        dt_D = trajectory.MSD_D[:,0]
        MSD_D = trajectory.MSD_D[:,1]
        fit_D = trajectory.MSD_D[:,2]
        residues_D = trajectory.MSD_D[:,3]
        dt_r = trajectory.MSD_fit[:,0]
        MSD_r = trajectory.MSD_fit[:,1]
        fit_r = trajectory.MSD_fit[:,2]
        residues_r = trajectory.MSD_fit[:,3]
        return trajectory_number, dt_D, MSD_D, fit_D, residues_D, dt_r, MSD_r, fit_r, residues_r     

    def save_rossier(self, trajectories):
        """
        rossier info: col0 = id, col1 = type immobile, col2 = type confined, col3 = type free
        col4 = analyse successful. col5 = tau, col6 = dtau, col7 = r, col8 = dr
        col9 = diffusion confined, col10 = d diffusion confined
        """
        self.rossier_info = np.zeros([np.size(trajectories), 11])
        i = 0
        for trajectory in trajectories:
            self.rossier_info[i,0] = trajectory.trajectory_number 
            # if trajectory is immobile the analysis was never made, all analysis output is 0 by default
            if trajectory.immobility:
                self.rossier_info[i,1] = trajectory.immobility
                self.rossier_info[i,2] = False
                self.rossier_info[i,3] = False
                self.rossier_info[i,4] = False  # no analyse was made -> False
                self.rossier_info[i,5] = 0
                self.rossier_info[i,6] = 0
                self.rossier_info[i,7] = 0
                self.rossier_info[i,8] = 0
                self.rossier_info[i,9] = 0
                self.rossier_info[i,10] = 0
            # if trajectory is confined it is not immobile and not free
            elif trajectory.confined:
                self.rossier_info[i,1] = trajectory.immobility
                self.rossier_info[i,2] = trajectory.confined
                self.rossier_info[i,3] = not trajectory.confined
            # if trajectory is free it is not confined and not immobile
            elif not trajectory.confined:
                self.rossier_info[i,1] = trajectory.immobility
                self.rossier_info[i,2] = trajectory.confined
                self.rossier_info[i,3] = not trajectory.confined
            # analysis is made for trajectories not immobile -> if analysis was successful, output gets values
            if trajectory.analyse_successful and not trajectory.immobility:
                self.rossier_info[i,4] = trajectory.analyse_successful
                self.rossier_info[i,5] = trajectory.tau
                self.rossier_info[i,6] = trajectory.dtau
                self.rossier_info[i,7] = trajectory.r
                self.rossier_info[i,8] = trajectory.dr
                self.rossier_info[i,9] = trajectory.D_conf 
                self.rossier_info[i,10] = trajectory.chi_MSD_fit
            # if analysis was not successful -> output is 0 by default
            elif not trajectory.analyse_successful and not trajectory.immobility:
                self.rossier_info[i,4] = trajectory.analyse_successful
                self.rossier_info[i,5] = 0
                self.rossier_info[i,6] = 0
                self.rossier_info[i,7] = 0
                self.rossier_info[i,8] = 0
                self.rossier_info[i,9] = 0
                self.rossier_info[i,10] = 0
            i += 1
 
    
def main():
    pass
    

if __name__ == "__main__":
    main()
    