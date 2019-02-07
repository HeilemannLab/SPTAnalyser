# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 15:32:28 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import time
import os
import numpy as np
from pySPT.analysis import cell
from pySPT.analysis import trajectoryStatistics
from pySPT.analysis import coverSlip

def main():
    #load roi
    file_name_roi = "F:\\Marburg\\single_colour_tracking\\resting\\roi.log"
    roi_file = np.genfromtxt(file_name_roi, dtype=None, delimiter=",", skip_header=3, encoding=None)
    print("Roi file", roi_file)
    #print(roi_file[1][0])
    #load cell files
    file_names = []
    file_name01 = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.MIA\\tracking\\cell_1_MMStack_Pos0.ome_MIA.trc"
    file_name02 = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell2\\cell_2_MMStack_Pos0.ome.MIA\\tracking\\cell_2_MMStack_Pos0.ome_MIA.trc"
    #file_name03 = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell3\\cell_3_MMStack_Pos0.ome.MIA\\tracking\\cell_3_MMStack_Pos0.ome_MIA.trc"
    #file_name10 = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell10\\cell_10_MMStack_Pos0.ome.MIA\\tracking\\cell_10_MMStack_Pos0.ome_MIA.trc"
    file_names.append(file_name01)
    file_names.append(file_name02)
    #file_names.append(file_name03)
    #file_names.append(file_name10)
    
# =============================================================================
#     cell_trajectories_name = []
#     raw_base_name = ["cell_1_MMStack_Pos0", "cell_2_MMStack_Pos0"]
#     cell_trajectories_name.append(raw_base_name)
#     for name in raw_base_name:
#         print("name")
#         for file in roi_file:
#             if name in file[0]:
#                 cell_size = file[1]
# =============================================================================

                

            
    
# =============================================================================
#     # load background files
#     file_nameb01 = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Background1\\background_1_MMStack_Pos0.ome.MIA\\tracking\\background_1_MMStack_Pos0.ome_MIA.trc"
#     file_nameb02 = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Background2\\background_2_MMStack_Pos0.ome.MIA\\tracking\\background_2_MMStack_Pos0.ome_MIA.trc"
#     file_nameb03 = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Background3\\background_3_MMStack_Pos0.ome.MIA\\tracking\\background_3_MMStack_Pos0.ome_MIA.trc"
#     file_names_bg = []
#     file_names_bg.append(file_nameb01)
#     file_names_bg.append(file_nameb02)
#     file_names_bg.append(file_nameb03)
#     
# =============================================================================

    # create list with file names for cells
    cell_trajectories = []  # contains trajectories for each cell in separate lists
    cells = []  # contains cell objects
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
    
    
    
# =============================================================================
#     # create list with file names for background
#     background_trajectories = []
#     background_trajectories_name = []
#     for file_name in file_names_bg:
#         one_cell = cell.Cell()
#         trc_file = np.loadtxt(file_name, usecols = (0, 1, 2, 3, 5)) # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = intensity
#         one_cell.trc_file = trc_file
#         one_cell.create_trajectories()
#         one_cell.analyse_trajectories()
#         background_trajectories.append(one_cell.analysed_trajectories)
#         base_name = os.path.basename(file_name)
#         raw_base_name = ""
#         for i in base_name:
#             if i == ".":
#                 break
#             else:
#                 raw_base_name += i
#         background_trajectories_name.append(raw_base_name)
# =============================================================================
# =============================================================================
#     cover_slip = coverSlip.CoverSlip()
#     cover_slip.background = background_trajectories
#     cover_slip.background_name = background_trajectories_name
#     cover_slip.log_transform()    
# =============================================================================
    
# =============================================================================
#     # initialize trajectory statics
#     trajectory_stats = trajectoryStatistics.TrajectoryStatistics()
#     trajectory_stats.cell_trajectories = cell_trajectories
#     trajectory_stats.cell_trajectories_name = cell_trajectories_name
#     trajectory_stats.get_index()
#     trajectory_stats.create_init_filter_lst()
#     trajectory_stats.filter_length(20)
#     trajectory_stats.filter_type(["immobile", "free", "confined"])  # get trajectories with diffusion types [x, y, z]
#     trajectory_stats.filter_D()
#     trajectory_stats.plot_trajectory(1,1)
# 
# =============================================================================

    
    #trajectory_stats.background_frequencies()

        
# =============================================================================
    #checking certain MSD values and D ...
#     print(trajectory_stats.cell_trajectories[0][0].localizations)
#     print(trajectory_stats.cell_trajectories[0][0].D)
#     print(trajectory_stats.cell_trajectories[0][0].MSDs)
# =============================================================================
                    
# =============================================================================
    #print a certain trajectory 
#     trajectory_stats.plot_trajectory(1,218)  # cell, trajectory
# =============================================================================
    
# =============================================================================
    #comparing ML fit results with py fit results
#     print(trajectory_stats.trajectory_list[0][247].analyse_successful)
#     print("Diffusion coeff:", trajectory_stats.trajectory_list[0][247].D)
#     print("D_conf", trajectory_stats.trajectory_list[0][247].D_conf)
#     print("r", trajectory_stats.trajectory_list[0][247].r)
#     print("tau", trajectory_stats.trajectory_list[0][247].tau)
#     print("ML_D_conf", trajectory_stats.trajectory_list[0][247].D_conf_ML)
#     print("ML_r", trajectory_stats.trajectory_list[0][247].r_ML)
#     print("ML_tau", trajectory_stats.trajectory_list[0][247].tau_ML)
# =============================================================================
    
    
if __name__ == "__main__":
    main()
    