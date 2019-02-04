# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 15:32:28 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import time
import numpy as np
from pySPT.analysis import cell
from pySPT.analysis import trajectoryStatistics

def main():
    one_cell = cell.Cell()
    start = time.time()
    one_cell.load_file()
    end = time.time()
    print("time load file, create objects", end-start)
    #cell.get_trajectories()
    start2=time.time()
    #one_cell.run_processes()
    one_cell.analyse_trajectories()  # analyse trajectories without multiprocessing
    end2 = time.time()
    print("Object analysis took {} seconds".format(end2-start2))
    #cell.plot_trajectorie(1)
    trajectory_stats = trajectoryStatistics.TrajectoryStatistics()
    trajectory_stats.trajectory_list.append(one_cell.analysed_trajectories)
    #print(trajectory_stats.trajectory_list)
    trajectory_stats.plot_trajectory(1,218)  # cell, trajectory

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
    