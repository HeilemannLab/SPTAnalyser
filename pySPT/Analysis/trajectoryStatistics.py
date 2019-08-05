# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 11:45:55 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""


import numpy as np
import copy
import math
import matplotlib.pyplot as plt
import time


class TrajectoryStatistics():
    def __init__(self):
        self.cells = []  # contains lists of cell objects
        self.cell_trajectories = [] # [[],[]] contains list of cells, cells contain trajectories
        self.cell_trajectories_filtered = []  # filtered for thresholds & booleans
        self.cell_trajectories_filtered_thresholds = []  # filtered for thresholds only
        self.cell_trajectories_filtered_index = [] 
        self.backgrounds = []  # background objects
        self.background_trajectories = [] # [[],[]] contains list of cells, cells contain trajectories
        self.background_trajectories_filtered = []  # filtered for thresholds & booleans
        self.background_trajectories_filtered_thresholds = []  # filtered for thresholds only
        self.background_trajectories_filtered_index = []  # deep copy of original cell trajectories index
        self.filtered_trc_files = []  # contains filtered trc type files of cells
        self.trc_files_hmm = []
        self.filter_settings = []  # list with boolean values for immob, confined, free, analyse success, not success
        self.filter_thresholds_values = []  # values for minD, maxD, minlength, maxlength
        self.min_D = math.inf
        self.max_D = - math.inf
        self.total_trajectories = 0  # amount of trajectories in data set
        self.total_trajectories_cell = []  # list with amount of total trajectories per cell 
        self.total_trajectories_filtered = 0  # amount of trajectories in data set after filter
        self.total_trajectories_filtered_cell = []  # list with amount of total trajectories per cell after filter
        self.cell_type_count = []  # tupels with (immobile, confined, free %) for each cell
        self.bin_size = 0.1  # bin size for freq vs D plot
        self.cell_sizes = []  # filled by jnb -> cell.size for cell in cover_slip.cells
        self.bg_sizes = []  # filled by jnb -> background.size for bg in cover_slip.bg
        self.hist_log_Ds = []  # histograms (logD vs freq) from all cells as np arrays in this list
        self.hist_log_Ds_bg = []
        self.diffusion_frequencies = []  # only freq (divided by cell size) of cells
        self.diffusion_frequencies_bg = []  # only freq (divided my bg size) of bg
        self.hist_diffusion = []  # diffusions from histogram calculation, transformed back -> 10^-(log10(D))
        self.mean_frequencies = []  # mean frequencies
        self.mean_frequencies_percent = []  # mean frequencies, size corrected, sum = 1 in percent
        self.mean_error = []  # standard error of mean value
        self.mean_error_percent = []  # standard error of mean value
        self.mean_frequencies_bg = [] # mean frequencies
        self.mean_error_bg = [] # standard error of mean value
        self.normalization_factor = 0.0  # 100/sum of all mean frequencies
        self.normalization_factor_corrected = 0.0  # 100/sum of all mean frequencies
        self.corrected_frequencies_error = []
        self.corrected_frequencies_error_percent = []
        self.tau_threshold_min_length = float(math.inf)  # from all cells, the min rossier length will be set as the default value
        self.corrected_frequencies = []  # mean cell frequencies - mean bg frequencies
        self.corrected_frequencies_percent = []  # mean cell frequencies - mean bg frequencies   
        self.sigma_dyns = []  # dynamic localization error, based on filtered trajectories         
        self.diff_fig = []  # log diffusion plot
        self.fig_name = []
        self.D_mean_types = []  # containes averaged D for immobile, confined and free diffusion
        self.dD_mean_types = []
        self.length_mean_types = []
        self.dlength_mean_types = []
        self.D_cell_types = []  
        self.dD_cell_types = []
        self.length_cell_types = []
        self.dlength_cell_types = []

    def calc_mean_statistics(self):
        immobile_tracks = []
        confined_tracks = []
        free_tracks = []
        notype_tracks = []
        for cell in self.cell_trajectories_filtered:
            immobile_tracks_cell = []
            confined_tracks_cell = []
            free_tracks_cell = []
            notype_tracks_cell = []
            for i in cell:
                if i.immobility and not i.analyse_successful and not i.confined:
                    immobile_tracks.append(i)
                    immobile_tracks_cell.append(i)
                elif not i.immobility and i.confined and i.analyse_successful:
                    confined_tracks.append(i)
                    confined_tracks_cell.append(i)
                elif not i.immobility and not i.confined and i.analyse_successful:
                    free_tracks.append(i)
                    free_tracks_cell.append(i)
                elif not i.immobility and not i.analyse_successful:
                    notype_tracks.append(i)
            immob_trajectories_cell = [track.D for track in immobile_tracks_cell]
            immobile_D_cell = np.mean(immob_trajectories_cell)
            immobile_dD_cell = np.std(immob_trajectories_cell, ddof=1)/math.sqrt(len(immob_trajectories_cell))
            conf_trajectories_cell = [track.D for track in confined_tracks_cell]
            confined_D_cell = np.mean(conf_trajectories_cell)
            confined_dD_cell = np.std(conf_trajectories_cell, ddof=1)/math.sqrt(len(conf_trajectories_cell))
            free_trajectories_cell = [track.D for track in free_tracks_cell]
            free_D_cell = np.mean(free_trajectories_cell)
            free_dD_cell = np.std(free_trajectories_cell, ddof=1)/math.sqrt(len(free_trajectories_cell))
            notype_trajectories_cell = [track.D for track in notype_tracks_cell]
            notype_D_cell = np.mean(notype_trajectories_cell)
            notype_dD_cell = np.std(notype_trajectories_cell, ddof=1)/math.sqrt(len(notype_trajectories_cell))
            D_type_cell = (immobile_D_cell, confined_D_cell, free_D_cell, notype_D_cell)
            dD_type_cell = (immobile_dD_cell, confined_dD_cell, free_dD_cell, notype_dD_cell)
            self.D_cell_types.append(D_type_cell)
            self.dD_cell_types.append(dD_type_cell)
            immob_trajectories_l_cell = [track.length_trajectory for track in immobile_tracks_cell]
            immobile_length_cell = np.mean(immob_trajectories_l_cell)
            immobile_dlength_cell = np.std(immob_trajectories_l_cell, ddof=1)/math.sqrt(len(immob_trajectories_l_cell))
            conf_trajectories_l_cell = [track.length_trajectory for track in confined_tracks_cell]
            confined_length_cell = np.mean(conf_trajectories_l_cell)
            confined_dlength_cell = np.std(conf_trajectories_l_cell, ddof=1)/math.sqrt(len(conf_trajectories_l_cell))
            free_trajectories_l_cell = [track.length_trajectory for track in free_tracks_cell]
            free_length_cell = np.mean(free_trajectories_l_cell)
            free_dlength_cell = np.std(free_trajectories_l_cell, ddof=1)/math.sqrt(len(free_trajectories_l_cell))
            notype_trajectories_l_cell = [track.length_trajectory for track in notype_tracks_cell]
            notype_length_cell = np.mean(notype_trajectories_l_cell)
            notype_dlength_cell = np.std(notype_trajectories_l_cell, ddof=1)/math.sqrt(len(notype_trajectories_l_cell))
            length_type_cell = (immobile_length_cell, confined_length_cell, free_length_cell, notype_length_cell)
            dlength_type_cell = (immobile_dlength_cell, confined_dlength_cell, free_dlength_cell, notype_dlength_cell)            
            self.length_cell_types.append(length_type_cell)
            self.dlength_cell_types.append(dlength_type_cell)
        immob_D = [track.D for track in immobile_tracks]
        immobile_D_mean = np.mean(immob_D)
        immobile_dD_mean = np.std(immob_D, ddof=1)/math.sqrt(len(immob_D))
        conf_D = [track.D for track in confined_tracks]
        confined_D_mean = np.mean(conf_D)
        confined_dD_mean = np.std(conf_D, ddof=1)/math.sqrt(len(conf_D))
        free_D = [track.D for track in free_tracks]
        free_D_mean = np.mean(free_D)
        free_dD_mean = np.std(free_D, ddof=1)/math.sqrt(len(free_D))
        notype_D = [track.D for track in notype_tracks]
        notype_D_mean = np.mean(notype_D)
        notype_dD_mean = np.std(notype_D, ddof=1)/math.sqrt(len(notype_D))
        self.D_mean_types.append(immobile_D_mean)
        self.D_mean_types.append(confined_D_mean)
        self.D_mean_types.append(free_D_mean)
        self.D_mean_types.append(notype_D_mean)
        self.dD_mean_types.append(immobile_dD_mean)
        self.dD_mean_types.append(confined_dD_mean)
        self.dD_mean_types.append(free_dD_mean)
        self.dD_mean_types.append(notype_dD_mean)
        immob_length = [track.length_trajectory for track in immobile_tracks]
        immobile_length_mean = np.mean(immob_length)
        immobile_dlength_mean = np.std(immob_length, ddof=1)/math.sqrt(len(immob_length))
        conf_length = [track.length_trajectory for track in confined_tracks]
        confined_length_mean = np.mean(conf_length)
        confined_dlength_mean = np.std(conf_length, ddof=1)/math.sqrt(len(conf_length))
        free_length = [track.length_trajectory for track in free_tracks]
        free_length_mean = np.mean(free_length)
        free_dlength_mean = np.std(free_length, ddof=1)/math.sqrt(len(free_length))
        notype_length = [track.length_trajectory for track in notype_tracks]
        notype_length_mean = np.mean(notype_length)
        notype_dlength_mean = np.std(notype_length, ddof=1)/math.sqrt(len(notype_length))
        self.length_mean_types.append(immobile_length_mean)
        self.length_mean_types.append(confined_length_mean)
        self.length_mean_types.append(free_length_mean)
        self.length_mean_types.append(notype_length_mean)
        self.dlength_mean_types.append(immobile_dlength_mean)
        self.dlength_mean_types.append(confined_dlength_mean)
        self.dlength_mean_types.append(free_dlength_mean)
        self.dlength_mean_types.append(notype_dlength_mean)
      

    def calc_min_rossier_length(self):
        #self.tau_threshold_min_length = float(math.inf)
        for cell in self.cells:
            print("before cell", cell.tau_threshold_min_length, type(cell.tau_threshold_min_length))
            print("before self", self.tau_threshold_min_length, type(self.tau_threshold_min_length))
            if cell.tau_threshold_min_length < self.tau_threshold_min_length:
                self.tau_threshold_min_length = cell.tau_threshold_min_length
                print(type(self.tau_threshold_min_length))
        self.tau_threshold_min_length = str(self.tau_threshold_min_length)
        print(self.tau_threshold_min_length, type(self.tau_threshold_min_length))
        
    def create_filtered_framework(self):
        """
        JNB needs a list in the shape of the resulting cell trajecories filtered list.
        """
        for cell in self.cell_trajectories:
            self.cell_trajectories_filtered.append([])
        
    def default_statistics(self):
        self.cell_trajectories_filtered = []  # deep copy of original cell trajectories
        self.cell_trajectories_filtered_thresholds = []
        self.cell_trajectories_filtered_index = []
        self.background_trajectories_filtered = []
        self.background_trajectories_filtered_thresholds = []
        self.background_trajectories_filtered_index = []
        self.total_trajectories = 0  # amount of trajectories in data set
        self.total_trajectories_cell = []  # list with amount of total trajectories per cell 
        self.cell_type_count = []
        self.sigma_dyns = []
        #######self.filtered_trc_files = []
        self.trc_files_hmm = []
        self.D_mean_types = []  
        self.dD_mean_types = []
        self.length_mean_types = []
        self.dlength_mean_types = []
        self.D_cell_types = []  
        self.dD_cell_types = []
        self.length_cell_types = []
        self.dlength_cell_types = []
        
    def run_statistics(self, min_length, max_length, min_D, max_D, filter_immob, filter_confined,
                       filter_free, filter_analyse_not_successful):
        """
        Initialize filtered lists, filter functions and percentage function. Filters are applied to the cell & background dataset.
        :param min/max: Min/max thresholds for diffusion coefficient and trajectory length. Only trajectories in the closed interval are included.
        :param filter: if checked in JNB -> True. Filter for types and for type determination not successful. 
        """
        start = time.time()
        self.default_statistics()
        self.filter_settings = [filter_immob, filter_confined, filter_free, filter_analyse_not_successful]
        self.calc_amount_trajectories()
        self.get_trc_files_hmm()
        try:
            min_length = int(min_length)
            print("min trajectory length:", min_length)
        except ValueError:
            min_length = self.get_min_length()
            print("min trajectory length:", min_length)
        try:
            max_length = int(max_length)
            print("max trajectory length:", max_length)
        except ValueError:
            max_length = self.get_max_length()
            print("max trajectory length:", max_length)
        try:
            min_D = float(min_D)
            print("min diffusion coefficient: {} [\u03BCm\u00b2/s]".format(min_D))
        except ValueError:
            min_D = self.get_min_D()
            print("min diffusion coefficient: {} [\u03BCm\u00b2/s]".format(min_D))
        try:
            max_D = float(max_D)
            print("max diffusion coefficient: {} [\u03BCm\u00b2/s]".format(max_D))
        except ValueError:
            max_D = self.get_max_D()
            print("max diffusion coefficient: {} [\u03BCm\u00b2/s]".format(max_D))
        self.filter_thresholds_values = [min_length, max_length, min_D, max_D]
        self.filter_thresholds(min_length, max_length, min_D, max_D)
        self.filter_type(filter_immob, filter_confined, filter_free, filter_analyse_not_successful)
        self.sort_filtered_trajectories()
        self.create_index_lst()
        self.calc_sigma_dyns()
        self.calc_mean_statistics()
        print("Initialization took {} s".format(time.time()-start))
        if filter_immob:
            print("Filter for immobile.")
        if filter_confined:
            print("Filter for confined.")
        if filter_free:
            print("Filter for free.")
        if not filter_analyse_not_successful:
            print("Filter for type determination successful only.")
        elif filter_analyse_not_successful:
            print("Include type determination not successful.")
        type_percentage = self.type_percentage()
        print("%.2f %% are immobile, mean diffusion coefficient = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(type_percentage[0], self.D_mean_types[0], self.length_mean_types[0]))
        print("%.2f %% are confined, mean diffusion coefficient = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(type_percentage[1], self.D_mean_types[1], self.length_mean_types[1]))
        print("%.2f %% are free, mean diffusion coefficient = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(type_percentage[2], self.D_mean_types[2], self.length_mean_types[2]))      
        print("%.2f %% could not be analysed, mean diffusion coefficient = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(type_percentage[3], self.D_mean_types[3], self.length_mean_types[3]))
        if self.total_trajectories_filtered == 0:
            print("The selection excludes all data.")
        print("Trajectories included:", self.total_trajectories_filtered)
        print("Trajectories excluded:", self.total_trajectories - self.total_trajectories_filtered)
        
# =============================================================================
        # filter testing
        #print("cell types %", self.cell_type_count)
#         print("cell count", self.total_trajectories_filtered_cell)
#         print("index", self.cell_trajectories_filtered_index)
#         print("cell count", len(self.cell_trajectories_filtered_index[0]), len(self.cell_trajectories_filtered_index[1]))
#         print("background index", self.background_trajectories_filtered_index)
#         print("bg count", len(self.background_trajectories_filtered_index[0]), len( self.background_trajectories_filtered_index[1]))
# =============================================================================
        
    def plot_trajectory(self, cell, number):
        cell = int(cell) - 1
        number = int(number) -1
        self.cell_trajectories[cell][number].plot_particle()

    def get_max_D(self):
        """
        Max diffusion coefficient of trajectory.
        :return: int of max diffusion coefficient.
        """
        max_D = - math.inf
        for cell in range(0, len(self.cell_trajectories)):
            for trajectory in self.cell_trajectories[cell]:
                if max_D < trajectory.D:
                    max_D = trajectory.D
        return float(max_D)
    
    def get_min_D(self):
        """
        Min diffusion coefficient of trajectory. 
        :return: int of min diffusion coefficient.
        """
        min_D = math.inf
        for cell in range(0, len(self.cell_trajectories)):
            for trajectory in self.cell_trajectories[cell]:
                if min_D > trajectory.D:
                    min_D = trajectory.D
        return float(min_D)

    def get_max_length(self):
        """
        Max length of trajectory.
        :return: int of max length.
        """
        max_length = - math.inf
        for cell in range(0, len(self.cell_trajectories)):
            for trajectory in self.cell_trajectories[cell]:
                if max_length < trajectory.length_trajectory:
                    max_length = trajectory.length_trajectory
        return int(max_length)

    def get_min_length(self):
        """
        Max length of trajectory.
        :return: int of max length.
        """
        min_length = math.inf
        for cell in range(0, len(self.cell_trajectories)):
            for trajectory in self.cell_trajectories[cell]:
                if min_length > trajectory.length_trajectory:
                    min_length = trajectory.length_trajectory
        return int(min_length)
            
    def filter_thresholds(self, min_length, max_length, min_diff, max_diff):
        for cell_index in range(len(self.cell_trajectories)):
            filtered_cell = [trajectory for trajectory in self.cell_trajectories[cell_index] if trajectory.length_trajectory >= min_length 
                        and trajectory.length_trajectory <= max_length
                        and trajectory.D >= min_diff and trajectory.D <= max_diff]
            self.cell_trajectories_filtered_thresholds.append(filtered_cell)
        if self.background_trajectories:
            for bg_index in range(len(self.background_trajectories)):
                filtered_bg = [trajectory for trajectory in self.background_trajectories[bg_index] if trajectory.length_trajectory >= min_length 
                            and trajectory.length_trajectory <= max_length
                            and trajectory.D >= min_diff and trajectory.D <= max_diff]
                self.background_trajectories_filtered_thresholds.append(filtered_bg)
        
    def filter_type(self, filter_immob, filter_confined, filter_free, filter_type_not_successful):
        for cell_index in range(len(self.cell_trajectories)):
            filtered_cell = []
            if filter_immob:
                filtered_cell_immob = [trajectory for trajectory in self.cell_trajectories_filtered_thresholds[cell_index] if trajectory.immobility
                                   and not trajectory.confined and not trajectory.analyse_successful]
                filtered_cell.extend(filtered_cell_immob)
            if filter_confined:
                filtered_cell_confined = [trajectory for trajectory in self.cell_trajectories_filtered_thresholds[cell_index] if trajectory.confined
                                   and not trajectory.immobility and trajectory.analyse_successful]
                filtered_cell.extend(filtered_cell_confined)
            if filter_free:
                filtered_cell_free = [trajectory for trajectory in self.cell_trajectories_filtered_thresholds[cell_index] if not trajectory.confined
                                   and not trajectory.immobility and trajectory.analyse_successful]
                filtered_cell.extend(filtered_cell_free)
            if filter_type_not_successful:
                filtered_cell_type_not_successful = [trajectory for trajectory in self.cell_trajectories_filtered_thresholds[cell_index] if not trajectory.analyse_successful
                                                     and not trajectory.immobility]
                filtered_cell.extend(filtered_cell_type_not_successful)
            self.cell_trajectories_filtered.append(filtered_cell)
        if self.background_trajectories:
            for bg_index in range(len(self.background_trajectories)):
                filtered_bg = []
                if filter_immob:
                    filtered_cell_immob = [trajectory for trajectory in self.background_trajectories_filtered_thresholds[bg_index] if trajectory.immobility
                                       and not trajectory.confined and not trajectory.analyse_successful]
                    filtered_bg.extend(filtered_cell_immob)
                if filter_confined:
                    filtered_cell_confined = [trajectory for trajectory in self.background_trajectories_filtered_thresholds[bg_index] if trajectory.confined
                                       and not trajectory.immobility and trajectory.analyse_successful]
                    filtered_bg.extend(filtered_cell_confined)
                if filter_free:
                    filtered_cell_free = [trajectory for trajectory in self.background_trajectories_filtered_thresholds[bg_index] if not trajectory.confined
                                       and not trajectory.immobility and trajectory.analyse_successful]
                    filtered_bg.extend(filtered_cell_free)
                if filter_type_not_successful:
                    filtered_cell_type_not_successful = [trajectory for trajectory in self.background_trajectories_filtered_thresholds[bg_index] if not trajectory.analyse_successful
                                                         and not trajectory.immobility]
                    filtered_bg.extend(filtered_cell_type_not_successful)
                self.background_trajectories_filtered.append(filtered_bg)

    def sort_filtered_trajectories(self):
        for cell_index in range(len(self.cell_trajectories)):
            self.cell_trajectories_filtered[cell_index].sort(key=lambda trajectory: trajectory.trajectory_number)
        for bg_index in range(len(self.background_trajectories)):
            self.background_trajectories_filtered[bg_index].sort(key=lambda trajectory: trajectory.trajectory_number)
    
    def filter_cell_trc(self):
        """
        For each cell a filtered trc file is created. Original trc file is copied.
        If first column of trc file (index) is in the filtered cell index list, the column will be in the filtered trc file.
        """
        self.filtered_trc_files = []
        for cell_index in range(len(self.cell_trajectories)):
            if self.cells[cell_index].seg_id:
                trajectory_seg_id = 6
            else:
                trajectory_seg_id = 0
            trc_file = copy.deepcopy(self.cells[cell_index].converted_trc_file_type)
            trc_idx = np.isin(trc_file[:,trajectory_seg_id], self.cell_trajectories_filtered_index[cell_index])
            trc_file=trc_file[trc_idx,:]
            self.filtered_trc_files.append(trc_file)     
            
    def get_trc_files_hmm(self):
        for cell_index in range(len(self.cell_trajectories)):
            self.trc_files_hmm.append(self.cells[cell_index].filtered_trc_file_hmm)
        
    def create_index_lst(self):
        """
        Create a list which contains trajectory object numbers for each cell -> [[],[]].
        If bg -> create a list for bg as well.
        """
        for cell_index in range(len(self.cell_trajectories)):
            trajectory_index = []
            for trajectory in self.cell_trajectories_filtered[cell_index]:
                trajectory_index.append(trajectory.trajectory_number)
            self.cell_trajectories_filtered_index.append(trajectory_index)
        if self.background_trajectories:
            for bg_index in range(len(self.background_trajectories)):
                trajectory_index_bg = []
                for trajectory in self.background_trajectories_filtered[bg_index]:
                    trajectory_index_bg.append(trajectory.trajectory_number)
                self.background_trajectories_filtered_index.append(trajectory_index_bg)
                
    def calc_amount_trajectories(self):
        for cell_index in range(0, len(self.cell_trajectories)):
            self.total_trajectories_cell.append(len(self.cell_trajectories[cell_index]))
        self.total_trajectories = np.sum(self.total_trajectories_cell)
        
    def type_percentage(self):
        """
        (immob/confined -> true/false=immob, false/true=conf, false/false=free)
        Calculate percentage of immobile free and confined based on total number of trajectories in all cells.
        If no trajectory exists (total_trajectories = 0) percentages will be set to zero, no calculation will be made.
        """
        data_selected = True
        self.total_trajectories_filtered = 0
        self.total_trajectories_filtered_cell = []
        self.cell_type_count = []
        for cell_index in range(0, len(self.cell_trajectories_filtered)):
            self.total_trajectories_filtered_cell.append(len(self.cell_trajectories_filtered[cell_index]))
        self.total_trajectories_filtered = np.sum(self.total_trajectories_filtered_cell)
        if self.total_trajectories_filtered == 0:
            data_selected = False
        if data_selected:
            # overview over the entire dataset
            count_immobile = 0
            count_confined = 0
            count_free = 0
            count_not_successful = 0
            for cell in self.cell_trajectories_filtered:
                # overview for one cell
                count_immobile_cell = 0
                count_confined_cell = 0
                count_free_cell = 0
                count_not_successful_cell = 0
                for trajectory in cell:
                    if trajectory.immobility and not trajectory.confined and not trajectory.analyse_successful:
                        count_immobile_cell += 1
                        count_immobile += 1
                    if trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                        count_confined_cell += 1
                        count_confined += 1
                    # has to be not confined AND not immobile (otherwise it will count the immobile particles as well)
                    if not trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                        count_free_cell += 1
                        count_free +=1
                    if not trajectory.analyse_successful and not trajectory.immobility:
                        count_not_successful_cell += 1
                        count_not_successful += 1
                cell_index = self.cell_trajectories_filtered.index(cell)
                if self.total_trajectories_filtered_cell[cell_index]:
                    ratio_immobile_cell = count_immobile_cell/self.total_trajectories_filtered_cell[cell_index]*100
                    ratio_confined_cell = count_confined_cell/self.total_trajectories_filtered_cell[cell_index]*100
                    ratio_free_cell = count_free_cell/self.total_trajectories_filtered_cell[cell_index]*100
                    ratio_not_successful_cell = count_not_successful_cell/self.total_trajectories_filtered_cell[cell_index]*100
                    cell_types_percent = (ratio_immobile_cell, ratio_confined_cell, ratio_free_cell, ratio_not_successful_cell)
                else:
                    cell_types_percent = (0.0,0.0,0.0,0.0)  # if all trajectories from a cell are filtered, the type % are 0
                self.cell_type_count.append(cell_types_percent)
            ratio_immobile = count_immobile/self.total_trajectories_filtered*100
            ratio_confined = count_confined/self.total_trajectories_filtered*100
            ratio_free = count_free/self.total_trajectories_filtered*100
            ratio_not_successful = count_not_successful/self.total_trajectories_filtered*100
        else:  # if all trajectories are filtered the type % are 0 
            ratio_immobile = 0
            ratio_confined = 0
            ratio_free = 0
            ratio_not_successful = 0
        return ratio_immobile, ratio_confined, ratio_free, ratio_not_successful
    
    # plot diffusion vs frequencies.
    
    def run_diffusion_histogram(self, desired_bin_size, plot=True):
        # bin size can only be something that can be converted to float (integer or float, comma separated)
        try:
            float(desired_bin_size)
        except:
            print("Insert a dot separated float or integer as bin size (e.g. 0.1).")
        # bin size can not be 0
        else:
            if float(desired_bin_size) != 0.0:
                self.clear_attributes()
                self.bin_size = desired_bin_size
                self.determine_max_min_diffusion()
                self.diffusions_log(float(desired_bin_size))
                self.calc_nonlogarithmic_diffusions()
                self.determine_mean_frequency()
                self.calc_mean_error()
# =============================================================================
#                 print("normalization factor", self.normalization_factor)
#                 print("mean cell freq", self.mean_frequencies)
#                 print("mean cell freq %", self.mean_frequencies_percent)
#                 print("mean freq error", self.mean_error_percent)
#                 print("mean freq error%", self.mean_error_percent)
# =============================================================================
                if self.background_trajectories:
                    self.diffusions_log_bg(float(desired_bin_size))
                    self.determine_mean_frequency(is_cell=False)
                    self.calc_mean_error(is_cell=False)
                    self.calc_bg_corrected_freq()
# =============================================================================
#                     print("D", self.hist_diffusion)
#                     print("normalization factor corr", self.normalization_factor_corrected)
#                     print("corr freq", self.corrected_frequencies)
#                     print("corr freq %", self.corrected_frequencies_percent)
#                     print("dcorr freq", self.corrected_frequencies_error)
#                     print("dcorr freq%", self.corrected_frequencies_error_percent)
#                     print("mean bg freq", self.mean_frequencies_bg)
#                     print("cell", self.hist_log_Ds)
#                     print("BG", self.hist_log_Ds_bg)
#                     print("diffusion freq", self.diffusion_frequencies)
# =============================================================================
                if plot:
                    self.plot_bar_log_bins()
                    if self.background_trajectories:
                        self.plot_bar_log_bins_bg_corrected()
            else:
                print("Bin size can not be zero.")
        
    def clear_attributes(self):
        """
        If one restarts filtering, these attributes are empy (otherwise they would append).
        """
        self.hist_log_Ds = []  # histograms (logD vs freq) from all cells as np arrays in this list
        self.hist_diffusion = []  # diffusions from histogram calculation, transformed back -> 10^-(log10(D))
        self.mean_frequencies = []  # mean frequencies, size corrected
        self.mean_error = []
        self.hist_log_Ds_bg = []
        self.mean_frequencies_bg = []
        self.mean_error_bg = []
        
    def determine_max_min_diffusion(self):
        """
        Create np array with log10(D) and D. -> min and max values can be determined over that. 
        The min value has to be positive, because logarithm of value <= 0 are not defined.
        """
        for cell in self.cell_trajectories:
            for trajectory in cell:
                if trajectory.D < self.min_D and trajectory.D > 0:
                    self.min_D = trajectory.D
                if trajectory.D > self.max_D:
                    self.max_D = trajectory.D
                    
    def diffusions_log_bg(self, desired_bin_size):
        """
        For each cell initialize histogram with cell size and target array.
        """
        #desired_bin_size = 0.05
        for bg_index in range(0, len(self.background_trajectories_filtered)):
            log_Ds = np.zeros(len(self.background_trajectories_filtered[bg_index]))
            bg_size = self.bg_sizes[bg_index]
            for trajectory_index in range(0, len(self.background_trajectories_filtered[bg_index])):
                log_Ds[trajectory_index] =  np.log10(self.background_trajectories_filtered[bg_index][trajectory_index].D)
            self.calc_diffusion_frequencies(log_Ds, desired_bin_size, bg_size, is_cell=False)
                    
    def diffusions_log(self, desired_bin_size):
        """
        For each cell initialize histogram with cell size and target array.
        """
        #desired_bin_size = 0.05
        for cell_index in range(0, len(self.cell_trajectories_filtered)):
            log_Ds = np.zeros(len(self.cell_trajectories_filtered[cell_index]))
            cell_size = self.cell_sizes[cell_index]
            for trajectory_index in range(0, len(self.cell_trajectories_filtered[cell_index])):
                if self.cell_trajectories_filtered[cell_index][trajectory_index].D > 0:
                    log_Ds[trajectory_index] =  np.log10(self.cell_trajectories_filtered[cell_index][trajectory_index].D)
            self.calc_diffusion_frequencies(log_Ds, desired_bin_size, cell_size) 
        
    def calc_diffusion_frequencies(self, log_diff, desired_bin, size, is_cell=True):
        """
        :param log_diff: np array with log10(D) of one cell.
        :param size: cell size.
        :param desired_bin: bin size.
        """
        # min & max determined by diffusions_log_complete function
        min_bin = np.ceil(-6/desired_bin)*desired_bin
        max_bin = np.ceil(2/desired_bin)*desired_bin 
        bin_size = int(np.ceil((max_bin - min_bin)/desired_bin))
        #print(max_bin, min_bin, bin_size)
        hist = np.histogram(log_diff,
                            range = (min_bin, max_bin),
                            bins = bin_size)
        log_diffusion_hist = np.zeros([np.size(hist[0]),2])
        log_diffusion_hist[:,0] = hist[1][:-1]  # log(D)
        log_diffusion_hist[:,1] = hist[0][:]  # freq
        log_diffusion_hist[:,1] = log_diffusion_hist[:,1]  / size
        if is_cell:
            self.hist_log_Ds.append(log_diffusion_hist)
        else:
            self.hist_log_Ds_bg.append(log_diffusion_hist)
        
    def calc_nonlogarithmic_diffusions(self):
        """
        Calculate the nonlogarithmic diffusion coefficients from log10(D) from histogram.
        """
        self.hist_diffusion = 10**self.hist_log_Ds[0][:,0]
        
    def determine_mean_frequency(self, is_cell=True):
        """
        Mean frequency will be calculated based on all frequencies of cells.
        Normalize an array (sum of elements = 1) and represent is in percent (*100).
        """
        if is_cell:
            self.diffusion_frequencies = self.create_np_array(np.shape(self.hist_log_Ds)[1], len(self.cell_trajectories_filtered))
            for i in range (0, len(self.cell_trajectories_filtered)):
                self.diffusion_frequencies[:,i] = self.hist_log_Ds[i][:,1]
            self.mean_frequencies = self.calc_mean_frequencies(self.diffusion_frequencies)
            self.normalization_factor = 100/np.sum(self.mean_frequencies)
            self.mean_frequencies_percent = self.mean_frequencies * self.normalization_factor
        else:
            self.diffusion_frequencies_bg = self.create_np_array(np.shape(self.hist_log_Ds_bg)[1], len(self.background_trajectories_filtered))
            for i in range(0, len(self.background_trajectories_filtered)):
                self.diffusion_frequencies_bg[:,i] = self.hist_log_Ds_bg[i][:,1]
            self.mean_frequencies_bg = self.calc_mean_frequencies(self.diffusion_frequencies_bg)
            #self.normalization_factor_bg = 100/np.sum(self.mean_frequencies_bg)
            #self.mean_frequencies_bg *= self.normalization_factor_bg
        
    def calc_mean_error(self, is_cell=True):
        """
        Standard deviation (N-1) divided by square root of number of elements.
        Normalize an array (sum of elements = 1) and represent is in percent (*100).
        """
        if is_cell:
            self.mean_error =  np.std(self.diffusion_frequencies, axis=1, ddof=1)/(np.shape(self.diffusion_frequencies)[1])**(1/2) 
            self.mean_error_percent = self.mean_error * self.normalization_factor 
        else:
            self.mean_error_bg = np.std(self.diffusion_frequencies_bg, axis=1, ddof=1)/(np.shape(self.diffusion_frequencies_bg)[1])**(1/2)
            
    def calc_bg_corrected_freq(self):
        """
        Substract background frequencies of cell frequencies. if < 0 -> set to 0.
        Calc error.
        """
        self.corrected_frequencies = np.zeros(np.shape(self.mean_frequencies)[0],)
        self.corrected_frequencies_error = self.corrected_frequencies
        self.corrected_frequencies = np.subtract(self.mean_frequencies, self.mean_frequencies_bg)
        for i in range(np.shape(self.corrected_frequencies)[0]):
            if self.corrected_frequencies[i] < 0:
                self.corrected_frequencies[i] = 0
            self.corrected_frequencies_error[i] = ((self.mean_error[i])**2+(self.mean_error_bg[i])**2)**(1/2)
        self.normalization_factor_corrected = 100/np.sum(self.corrected_frequencies)
        self.corrected_frequencies_error_percent = self.corrected_frequencies_error * self.normalization_factor_corrected
        self.corrected_frequencies_percent = self.corrected_frequencies * self.normalization_factor_corrected
        
    def plot_bar_log_bins(self):
        self.diff_fig.append(plt.figure())
        self.fig_name.append("diffusion_hist")
        plt.subplot(111, xscale="log")
        (_, caps, _) = plt.errorbar(self.hist_diffusion, self.mean_frequencies_percent, yerr=self.mean_error_percent, capsize=4, label="relative frequency")  # capsize length of cap
        for cap in caps:
            cap.set_markeredgewidth(1)  # markeredgewidth thickness of cap (vertically)
        plt.xlim(self.min_D, self.max_D)
        plt.legend()
        plt.title("Distribution of diffusion coefficients")
        plt.ylabel("normalized relative occurence [%]")
        plt.xlabel("D [\u03BCm\u00b2/s]")
        plt.show() 
    
    def save_diff_fig(self, directory, folder_name):
        for fig, name in zip(self.diff_fig, self.fig_name):
            fig.savefig(directory + "\\"+ folder_name +  "\\" + name + ".pdf", format="pdf", transparent=True)

    def plot_bar_log_bins_bg_corrected(self):
        self.diff_fig.append(plt.figure())
        self.fig_name.append("diffusion_hist_background")
        plt.subplot(111, xscale="log")
        (_, caps, _) = plt.errorbar(self.hist_diffusion, self.corrected_frequencies_percent, yerr=self.corrected_frequencies_error_percent, capsize=4, label="relative frequency")  # capsize length of cap
        for cap in caps:
            cap.set_markeredgewidth(1)  # markeredgewidth thickness of cap (vertically)
        plt.xlim(self.min_D, self.max_D)
        plt.legend()
        plt.title("Distribution of diffusion coefficients (background corrected)")
        plt.ylabel("Normalized relative occurence [%]")
        plt.xlabel("D [\u03BCm\u00b2/s]")
        plt.show() 
        
    def calc_mean_frequencies(self, np_array):
        """
        Determine mean value over the horizonal entries of an np.array.
        -> [2,2][6,4] = [4,3].
        :param np_array: Np array to build mean over.
        """
        #mean_frequencies = np.zeros(np.size(self.hist_log_Ds[0]))
        mean_frequencies = np_array.mean(axis=1)
        return mean_frequencies
            
    def create_np_array(self, length, columns=1):
        """
        :param length: number of np arrays.
        :param columns: amount of entries in a numpy array.
        :return: return a numpy array.
        """
        np_array = np.zeros((length,columns))
        return np_array
    
    def calc_sigma_dyns(self):
        """
        Based on filtered trajectories, a new dynamic localization error is calculated.
        Parameters: Mean D, mean MSD_0 per cell & loc, dt.
        """
        for cell_idx in range(len(self.cell_trajectories_filtered)):
            if len(self.cell_trajectories_filtered[cell_idx]) == 0:
                self.sigma_dyns.append("no trajectories included")
            else:
                dof = self.cell_trajectories_filtered[cell_idx][0].dof
                dt = self.cell_trajectories_filtered[cell_idx][0].dt
                mean_D = np.mean([track.D for track in self.cell_trajectories_filtered[cell_idx]])
                mean_MSD_0 = np.mean([track.MSD_0 for track in self.cell_trajectories_filtered[cell_idx]])
                sigma_dyn = math.sqrt((mean_MSD_0+(4/3)*mean_D*dt)/dof)
                self.sigma_dyns.append(sigma_dyn)      
 
    
def main():
    pass
        
    
if __name__ == "__main__":
    main()
        