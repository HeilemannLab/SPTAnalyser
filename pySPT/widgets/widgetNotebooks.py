# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 14:02:19 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Large juncs of code that would scare the user of the JNB :) 
"""

#from .analysis import trajectory
#from pySPT.analysis import trajectory
#from pySPT.analysis import cell
from . import loadHdf5  # one dot if in same directory
from ..analysis import cell  # two dots for switching the folder
from ..analysis import trajectory
import time

def init_filter_notebook(cover_slip, widget_load_hdf5, load_hdf5):
    """
    :param: Objects created in JNB.
    """
    start = time.time()
    widget_load_hdf5.search_sub_folders(widget_load_hdf5.dir_name)
    load_hdf5.file_names = widget_load_hdf5.file_names
    load_hdf5.run_load_hdf5()
    for cell_index in range(0, load_hdf5.cell_numbers):
        one_cell = cell.Cell()
        one_cell.trc_file = load_hdf5.trc_files[cell_index]
        one_cell.pixel_size = load_hdf5.pixel_sizes[cell_index]
        one_cell.pixel_amount = load_hdf5.pixel_amounts[cell_index]
        one_cell.size = load_hdf5.sizes[cell_index]
        one_cell.name = load_hdf5.names[cell_index]
        for trajectory_index in range(0, load_hdf5.trajectory_numbers[cell_index]):
            one_trajectory = trajectory.Trajectory(load_hdf5.locs[cell_index][trajectory_index])
            one_trajectory.trajectory_number = trajectory_index+1
            one_trajectory.MSDs = load_hdf5.cells_trajectories_MSDs[cell_index][trajectory_index]
            one_trajectory.times = load_hdf5.cells_trajectories_times[cell_index][trajectory_index]
            one_trajectory.MSD_fit = load_hdf5.cells_trajectories_MSD_fit[cell_index][trajectory_index]
            one_trajectory.MSD_D = load_hdf5.cells_trajectories_MSD_D[cell_index][trajectory_index]
            one_trajectory.dt = load_hdf5.dts[cell_index]
            one_trajectory.dof = load_hdf5.dofs[cell_index]
            one_trajectory.D_min = load_hdf5.D_mins[cell_index]
            one_trajectory.length_trajectory = load_hdf5.cells_lengths_trajectories[cell_index][trajectory_index][0]
            one_trajectory.length_MSD = load_hdf5.cells_lengths_MSDs[cell_index][trajectory_index][0]
            one_trajectory.D = load_hdf5.cells_trajectories_D[cell_index][trajectory_index][0]
            one_trajectory.dD = load_hdf5.cells_trajectories_dD[cell_index][trajectory_index][0]
            one_trajectory.chi_D = load_hdf5.cells_trajectories_chi2_D[cell_index][trajectory_index][0]
            one_trajectory.chi_MSD_fit = load_hdf5.cells_trajectories_chi2_rossier[cell_index][trajectory_index][0]
            one_trajectory.MSD_0 = load_hdf5.cells_trajectories_MSD0[cell_index][trajectory_index][0]
            one_trajectory.fit_area = load_hdf5.fit_areas[cell_index]
            one_trajectory.tau = load_hdf5.cells_trajectories_tau[cell_index][trajectory_index][0]
            one_trajectory.dtau = load_hdf5.cells_trajectories_dtau[cell_index][trajectory_index][0]
            one_trajectory.D_conf = load_hdf5.cells_trajectories_Dconf[cell_index][trajectory_index][0]
            one_trajectory.r = load_hdf5.cells_trajectories_r[cell_index][trajectory_index][0]
            one_trajectory.dr = load_hdf5.cells_trajectories_dr[cell_index][trajectory_index][0]
            one_trajectory.tau_threshold = load_hdf5.tau_thresholds[cell_index]
            one_trajectory.immobility = bool(load_hdf5.cells_trajectories_type[cell_index][trajectory_index][0])
            one_trajectory.confined = bool(load_hdf5.cells_trajectories_type[cell_index][trajectory_index][1])
            one_trajectory.analyse_successful = bool(load_hdf5.cells_trajectories_analyse_successful[cell_index][trajectory_index][0])
            one_cell.analysed_trajectories.append(one_trajectory)
        cover_slip.cells.append(one_cell)
        cover_slip.cell_trajectories.append(one_cell.analysed_trajectories)
        cover_slip.cell_files = load_hdf5.file_names
    print("Initialization took {} s".format(time.time()-start))