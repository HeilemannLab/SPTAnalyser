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
from . import hdf5
from . import saveFiltered
from ..analysis import cell  # two dots for switching the folder
from ..analysis import trajectory
from ..analysis import trackAnalysis
import time
import numpy as np 
import os
from ipywidgets import HBox, VBox
from tqdm import tqdm_notebook as tqdm


def init_filter_notebook(cover_slip, widget_load_hdf5, load_hdf5, is_cell=True):
    """
    JNB: TrackStatistics
    :param: Objects created in JNB.
    :param is_cell: if true -> cell objects distributed to cell attributes of coverslip, 
    if false -> bg object distributed to bg attributes of coverslip.
    """
    start = time.time()
    load_hdf5.clear()
    if is_cell:
        widget_load_hdf5.search_sub_folders(widget_load_hdf5.dir_name, is_cell)
        load_hdf5.file_names = widget_load_hdf5.file_names
    else:
        widget_load_hdf5.search_sub_folders(widget_load_hdf5.dir_name_bg, is_cell)
        load_hdf5.file_names = widget_load_hdf5.file_names_bg
    load_hdf5.run_load_hdf5()  # initialize the loadHdf5 class -> fill all needed attributes with values.
    for cell_index in range(load_hdf5.cell_numbers):  # distribute the attributes to the objects
        one_cell = cell.Cell()
        one_cell.trc_file = load_hdf5.trc_files[cell_index]
        one_cell.pixel_size = load_hdf5.pixel_sizes[cell_index]
        one_cell.pixel_amount = load_hdf5.pixel_amounts[cell_index]
        one_cell.size = load_hdf5.cell_sizes[cell_index]
        one_cell.name = load_hdf5.names[cell_index]
        one_cell.tau_threshold = load_hdf5.tau_thresholds[cell_index]
        one_cell.dt = load_hdf5.dts[cell_index]
        one_cell.dof = load_hdf5.dofs[cell_index]
        one_cell.D_min = load_hdf5.D_mins[cell_index]
        one_cell.tau_threshold_min_length = load_hdf5.tau_min_trajectory_lengths[cell_index]
        for trajectory_index in range(0, load_hdf5.trajectory_numbers[cell_index]):
            one_trajectory = trajectory.Trajectory(load_hdf5.locs[cell_index][trajectory_index], one_cell.tau_threshold, one_cell.dt, one_cell.dof, one_cell.D_min)
            one_trajectory.trajectory_number = trajectory_index+1
            one_trajectory.MSDs = load_hdf5.cells_trajectories_MSDs[cell_index][trajectory_index]
            one_trajectory.times = load_hdf5.cells_trajectories_times[cell_index][trajectory_index]
            one_trajectory.MSD_fit = load_hdf5.cells_trajectories_MSD_fit[cell_index][trajectory_index]
            one_trajectory.MSD_D = load_hdf5.cells_trajectories_MSD_D[cell_index][trajectory_index]
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
            one_trajectory.immobility = bool(load_hdf5.cells_trajectories_type[cell_index][trajectory_index][0])
            one_trajectory.confined = bool(load_hdf5.cells_trajectories_type[cell_index][trajectory_index][1])
            one_trajectory.analyse_successful = bool(load_hdf5.cells_trajectories_analyse_successful[cell_index][trajectory_index][0])
            one_cell.analysed_trajectories.append(one_trajectory)
        # distinguish between cells & bg 
        if is_cell:
            cover_slip.cells.append(one_cell)
            cover_slip.cell_trajectories.append(one_cell.analysed_trajectories)
            cover_slip.cell_files = load_hdf5.file_names
        else:
            cover_slip.backgrounds.append(one_cell)
            cover_slip.background_trajectories.append(one_cell.analysed_trajectories)
            cover_slip.background_files = load_hdf5.file_names
    print("Initialization took {} s".format(time.time()-start))
    
    
def init_track_stats_widget_arrangement(widget11, widget21, widget31, widget41, widget51, widget12, widget22, widget32, widget42):
    """
    JNB Track Statistics.
    HBox are line arrangements, col/row
    return VBox with its inserted lines.
    """
    first_line = HBox([widget11, widget12])
    second_line = HBox([widget21, widget22])
    third_line = HBox([widget31, widget32])
    fourth_line = HBox([widget41, widget42])
    fifth_line = HBox([widget51])
    return VBox([first_line, second_line, third_line, fourth_line, fifth_line])
    

def init_save_track_analysis(cover_slip, cell_index, track_analysis):
    """
    JNB: track Analysis, saving.
    :param: Objects created in JNB.
    """       
    h5 = hdf5.Hdf5()
    
    h5.create_h5(cover_slip.cell_files[cell_index])

    cell = cover_slip.cells[cell_index]
    one_trajectory = cover_slip.cell_trajectories[cell_index][0]  # get trajectory attributes, that are the same for every trajectory
    h5.data_settings(cell.dt, cell.pixel_size, cell.pixel_amount, cell.size, cell.tau_threshold, cover_slip.tau_threshold_min_length, one_trajectory.fit_area, cell.dof, cell.D_min)
    
    h5.statistics(track_analysis.cell_type_count[cell_index][0], track_analysis.cell_type_count[cell_index][1], track_analysis.cell_type_count[cell_index][2], track_analysis.total_trajectories_cell[cell_index])


    h5.trc(np.shape(cell.trc_file), cell.trc_file[:,0], cell.trc_file[:,1], cell.trc_file[:,2], cell.trc_file[:,3], cell.trc_file[:,4], cell.trc_file[:,5])
    
    for trajectory in cover_slip.cell_trajectories[cell_index]:
        plot = track_analysis.save_plots(trajectory)
        h5.data_diffusion_plots(plot[0], plot[1], plot[2], plot[3], plot[4])
    for trajectory in cover_slip.cell_trajectories[cell_index]:
        plot = track_analysis.save_plots(trajectory)
        h5.data_rossier_plots(plot[0], plot[5], plot[6], plot[7], plot[8])
    for trajectory in cover_slip.cell_trajectories[cell_index]:
        h5.msd(trajectory.trajectory_number, trajectory.times, trajectory.MSDs)
    
    track_analysis.save_diff(cover_slip.cell_trajectories[cell_index])
    diff_info = track_analysis.diffusion_info
    h5.data_diffusion_info(track_analysis.number_of_trajectories, diff_info[:,0], diff_info[:,1], diff_info[:,2], diff_info[:,3], diff_info[:,4], diff_info[:,5])

    track_analysis.save_rossier(cover_slip.cell_trajectories[cell_index])
    rossier_info = track_analysis.rossier_info
    h5.data_rossier_info(track_analysis.number_of_trajectories, rossier_info[:,0],  rossier_info[:,1],  rossier_info[:,2],  rossier_info[:,3],  rossier_info[:,4],  rossier_info[:,5],  rossier_info[:,6],  rossier_info[:,7],  rossier_info[:,8],  rossier_info[:,9], rossier_info[:,10])
    
    
def init_save_filtered_analysis(cover_slip, cell_index, track_stats, directory, folder_name):
    """
    JNB Track Statistics, create filtered h5 files for cells that are not excluded by filters.
    """
    h5_filtered = saveFiltered.SaveFiltered()
    track_analysis = trackAnalysis.TrackAnalysis()
    cell_name = track_stats.cells[cell_index].name
    h5_filtered.create_h5(directory + "\\" + folder_name + "\\" + cell_name)
    cell = track_stats.cells[cell_index]
    one_trajectory = track_stats.cell_trajectories_filtered[cell_index][0]  # get trajectory attributes, that are the same for every trajectory
    h5_filtered.data_settings(cell.dt, cell.pixel_size, cell.pixel_amount, cell.size, cell.tau_threshold, cell.tau_threshold_min_length, one_trajectory.fit_area, cell.dof, cell.D_min)
    
    h5_filtered.statistics(track_stats.cell_type_count[cell_index][0], track_stats.cell_type_count[cell_index][1], track_stats.cell_type_count[cell_index][2], track_stats.total_trajectories_filtered_cell[cell_index], (track_stats.total_trajectories_cell[cell_index]-track_stats.total_trajectories_filtered_cell[cell_index]))
       
    h5_filtered.trc(np.shape(track_stats.filtered_trc_files[cell_index]), track_stats.filtered_trc_files[cell_index][:,0], track_stats.filtered_trc_files[cell_index][:,1], track_stats.filtered_trc_files[cell_index][:,2], track_stats.filtered_trc_files[cell_index][:,3], track_stats.filtered_trc_files[cell_index][:,4], track_stats.filtered_trc_files[cell_index][:,5])
    
    for trajectory in track_stats.cell_trajectories_filtered[cell_index]:
        plot = track_analysis.save_plots(trajectory)
        h5_filtered.data_diffusion_plots(plot[0], plot[1], plot[2], plot[3], plot[4])
    for trajectory in track_stats.cell_trajectories_filtered[cell_index]:
        plot = track_analysis.save_plots(trajectory)
        h5_filtered.data_rossier_plots(plot[0], plot[5], plot[6], plot[7], plot[8])
    for trajectory in track_stats.cell_trajectories_filtered[cell_index]:
        h5_filtered.msd(trajectory.trajectory_number, trajectory.times, trajectory.MSDs)
        
    h5_filtered.filter_info(track_stats.filter_settings, track_stats.filter_thresholds_values)
    
    track_analysis.save_diff(track_stats.cell_trajectories_filtered[cell_index])
    diff_info = track_analysis.diffusion_info
    h5_filtered.data_diffusion_info(track_analysis.number_of_trajectories, diff_info[:,0], diff_info[:,1], diff_info[:,2], diff_info[:,3], diff_info[:,4], diff_info[:,5])
    
    track_analysis.save_rossier(track_stats.cell_trajectories_filtered[cell_index])
    rossier_info = track_analysis.rossier_info
    h5_filtered.data_rossier_info(track_analysis.number_of_trajectories, rossier_info[:,0],  rossier_info[:,1],  rossier_info[:,2],  rossier_info[:,3],  rossier_info[:,4],  rossier_info[:,5],  rossier_info[:,6],  rossier_info[:,7],  rossier_info[:,8],  rossier_info[:,9], rossier_info[:,10])
    
    
def init_save_track_stats(h5_stats, track_stats, directory, folder_name, name):
    """
    JNB Track Statistics, create h5 file for all statistics.
    :param path: head path for statistics h5 file.
    :param name: raw base name for statistics h5 file.
    """
    if not os.path.exists(directory + "\\" + folder_name):
        os.makedirs(directory + "\\" + folder_name)
    
    h5_stats.create_h5(directory + "\\" + folder_name + "\\" + name)
    
    # if background files were loaded    
    if track_stats.background_trajectories:
        h5_stats.groups_bg()
        
        h5_stats.diffusion_plot_bg_normalized(len(track_stats.hist_diffusion), track_stats.hist_diffusion, track_stats.mean_frequencies_percent, track_stats.mean_error_percent, track_stats.corrected_frequencies_percent, track_stats.corrected_frequencies_error_percent)

        h5_stats.diffusion_plot_bg(len(track_stats.hist_diffusion), track_stats.hist_diffusion, track_stats.mean_frequencies,
                                   track_stats.mean_error, track_stats.mean_frequencies_bg, track_stats.mean_error_bg,
                                   track_stats.corrected_frequencies, track_stats.corrected_frequencies_error)        
        for bg in track_stats.backgrounds:
            bg_name = bg.name
            bg_index = track_stats.backgrounds.index(bg)
            h5_stats.bg_counts(bg_name, len(track_stats.hist_diffusion), track_stats.hist_diffusion, track_stats.diffusion_frequencies_bg[:,bg_index])
        background_info = []
        for i in range(len(track_stats.bg_sizes)):
            one_bg = (track_stats.backgrounds[i].name, track_stats.bg_sizes[i])
            background_info.append(one_bg)
        h5_stats.backgrounds(background_info)
    
    h5_stats.filter_info(track_stats.filter_settings, track_stats.filter_thresholds_values)
    
    h5_stats.statistics(track_stats.type_percentage()[0], track_stats.type_percentage()[1],
                         track_stats.type_percentage()[2], track_stats.total_trajectories_filtered, (track_stats.total_trajectories - track_stats.total_trajectories_filtered))
    
    h5_stats.diffusion_bin_size(track_stats.bin_size)
    
    # cell files are always loaded
    for cell in track_stats.cells:
        cell_name = cell.name
        cell_index = track_stats.cells.index(cell)
        h5_stats.cell_counts(cell_name, len(track_stats.hist_diffusion), track_stats.hist_diffusion, track_stats.diffusion_frequencies[:,cell_index])
    cell_info = []  # name, size
    for i in range(len(track_stats.cell_sizes)):
        one_cell = (track_stats.cells[i].name, track_stats.cell_sizes[i])
        cell_info.append(one_cell)
    
    # if no bg file was loaded only mean cell frequencies are determined
    if not track_stats.background_trajectories:
        h5_stats.diffusion_plot_normalized(len(track_stats.hist_diffusion), track_stats.hist_diffusion, track_stats.mean_frequencies_percent, track_stats.mean_error_percent)
        
        h5_stats.diffusion_plot(len(track_stats.hist_diffusion), track_stats.hist_diffusion, track_stats.mean_frequencies,
                                   track_stats.mean_error)
        
    h5_stats.cells(cell_info)

        