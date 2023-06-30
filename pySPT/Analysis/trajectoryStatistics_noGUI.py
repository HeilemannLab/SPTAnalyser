"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Load analysed *.h5 files from trackAnalysis, filter them, calculate dynamic localization errors,
globally visualize them, store filtered dataset.
"""

import numpy as np
import pandas as pd
import copy
import math
import matplotlib.pyplot as plt
import seaborn as sns
import time


class TrajectoryStatistics():
    def __init__(self):
        self.cells = []  # contains lists of cell objects
        self.cell_trajectories = [] # [[],[]] contains list of cells, cells contain trajectories
        self.cell_trajectories_filtered = []  # filtered for thresholds & booleans
        self.trajectories_immob_cells_filtered = []  # filtered for thresholds & booleans & split into types
        self.trajectories_conf_cells_filtered = []
        self.trajectories_free_cells_filtered = []
        self.trajectories_notype_cells_filtered = []
        self.trajectories_immob_notype_cells_filtered = []
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
        self.bin_size = 0.1  # bin size for freq vs D plot
        self.cell_sizes = []  # filled by jnb -> cell.size for cell in cover_slip.cells
        self.bg_sizes = []  # filled by jnb -> background.size for bg in cover_slip.bg
        self.hist_log_Ds = []  # histograms (logD vs freq) from all cells as np arrays in this list
        self.hist_log_Ds_bg = []
        self.diffusion_frequencies = []  # only freq (divided by cell size) of cells
        self.diffusion_frequencies_norm = []  # only freq (divided by cell size) of cells and norm to 100%
        self.diffusion_frequencies_bg = []  # only freq (divided my bg size) of bg
        self.hist_diffusion = []  # diffusions from histogram calculation, transformed back -> 10^-(log10(D))
        self.hist_diffusion_immob = []  # mean histogram frequencies of diffusions type immobile & error
        self.hist_diffusion_conf = []
        self.hist_diffusion_free = []
        self.hist_diffusion_notype = []
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
        self.diff_fig_types = []  # log diffusion plot per type
        self.diff_fig_types_merge = []
        self.fig_name = []
        self.percentages_cell_types = []
        self.type_percentages_mean = []  # mean values of type % immob, conf, free, notype, immob+notype
        self.type_percentages_error = []  # error of mean of type % immob, conf, free, notype, immob+notype
        self.D_cell_types = []
        self.dD_cell_types = []
        self.D_mean_types = []  # containes averaged D for immobile, confined and free diffusion
        self.dD_mean_types = []
        self.length_cell_types = []
        self.dlength_cell_types = []
        self.length_mean_types = []
        self.dlength_mean_types = []
        self.delta_t_types = []  # MSD plot x axis per type
        self.mean_MSDs_immob = []  # mean MSD values of all immobile trajectories
        self.mean_error_MSDs_immob = []  # # mean error of MSD values of all immobile trajectories
        self.mean_MSDs_conf = []
        self.mean_error_MSDs_conf = []
        self.mean_MSDs_free = []
        self.mean_error_MSDs_free = []
        self.mean_MSDs_notype = []
        self.mean_error_MSDs_notype = []
        self.mean_MSDs_immob_notype = []
        self.mean_error_MSDs_immob_notype = []
        self.MSD_fig_types = []
        self.MSD_fig_types_merge = []
        self.cell_sizes_fig = []
        self.n_segments_fig = []

    def create_filtered_framework(self):
        """
        JNB needs a list in the shape of the resulting cell trajecories filtered list.
        """
        for cell in self.cell_trajectories:
            self.cell_trajectories_filtered.append([])
        
    def default_statistics(self):
        self.cell_trajectories_filtered = []  # deep copy of original cell trajectories
        self.trajectories_immob_cells_filtered = []
        self.trajectories_conf_cells_filtered = []
        self.trajectories_free_cells_filtered = []
        self.trajectories_notype_cells_filtered = []
        self.trajectories_immob_notype_cells_filtered = []
        self.cell_trajectories_filtered_thresholds = []
        self.cell_trajectories_filtered_index = []
        self.background_trajectories_filtered = []
        self.background_trajectories_filtered_thresholds = []
        self.background_trajectories_filtered_index = []
        self.total_trajectories = 0  # amount of trajectories in data set
        self.total_trajectories_cell = []  # list with amount of total trajectories per cell
        self.sigma_dyns = []
        self.trc_files_hmm = []
        self.percentages_cell_types = []
        self.type_percentages_mean = []
        self.type_percentages_error = []
        self.D_mean_types = []  
        self.dD_mean_types = []
        self.D_cell_types = []  
        self.dD_cell_types = []
        self.length_mean_types = []
        self.dlength_mean_types = []
        self.length_cell_types = []
        self.dlength_cell_types = []
        self.type_percentages_cell_fig = []
        self.diff_cell_fig = []
        self.lengths_cell_fig = []
        
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
        self.D_average()
        self.length_average()
        self.type_percentage()
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
        if self.total_trajectories_filtered == 0:
            print("The selection excludes all data.")
        else:
            print("Trajectories included:", self.total_trajectories_filtered)
            print("Trajectories excluded:", self.total_trajectories - self.total_trajectories_filtered)
            print("%.2f %% are immobile, mean diffusion coefficient = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(self.type_percentages_mean[0], self.D_mean_types[0], self.length_mean_types[0]))
            print("%.2f %% are confined, mean diffusion coefficient = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(self.type_percentages_mean[1], self.D_mean_types[1], self.length_mean_types[1]))
            print("%.2f %% are free, mean diffusion coefficient = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(self.type_percentages_mean[2], self.D_mean_types[2], self.length_mean_types[2]))
            print("%.2f %% could not be analysed, mean diffusion coefficient = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(self.type_percentages_mean[3], self.D_mean_types[3], self.length_mean_types[3]))
            print("%.2f %% are immobile+notype, mean diffusion coefficient = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(self.type_percentages_mean[4], self.D_mean_types[4], self.length_mean_types[4]))

    def plot_trajectory(self, cell, number):
        cell = int(cell) - 1
        number = int(number) - 1
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
                self.trajectories_immob_cells_filtered.append(filtered_cell_immob)
            if not filter_immob:
                self.trajectories_immob_cells_filtered.append([])
            if filter_confined:
                filtered_cell_confined = [trajectory for trajectory in self.cell_trajectories_filtered_thresholds[cell_index] if trajectory.confined
                                   and not trajectory.immobility and trajectory.analyse_successful]
                filtered_cell.extend(filtered_cell_confined)
                self.trajectories_conf_cells_filtered.append(filtered_cell_confined)
            if not filter_confined:
                self.trajectories_conf_cells_filtered.append([])
            if filter_free:
                filtered_cell_free = [trajectory for trajectory in self.cell_trajectories_filtered_thresholds[cell_index] if not trajectory.confined
                                   and not trajectory.immobility and trajectory.analyse_successful]
                filtered_cell.extend(filtered_cell_free)
                self.trajectories_free_cells_filtered.append(filtered_cell_free)
            if not filter_free:
                self.trajectories_free_cells_filtered.append([])
            if filter_type_not_successful:
                filtered_cell_type_not_successful = [trajectory for trajectory in self.cell_trajectories_filtered_thresholds[cell_index] if not trajectory.analyse_successful
                                                     and not trajectory.immobility]
                filtered_cell.extend(filtered_cell_type_not_successful)
                self.trajectories_notype_cells_filtered.append(filtered_cell_type_not_successful)
            if not filter_type_not_successful:
                self.trajectories_notype_cells_filtered.append([])
            self.cell_trajectories_filtered.append(filtered_cell)

        for cell_idx in range(len(self.cells)):
            merge_sub = []
            merge_sub.extend(self.trajectories_immob_cells_filtered[cell_idx])
            merge_sub.extend(self.trajectories_notype_cells_filtered[cell_idx])
            self.trajectories_immob_notype_cells_filtered.append(merge_sub)

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
            trc_idx = np.isin(trc_file[:, trajectory_seg_id], self.cell_trajectories_filtered_index[cell_index])
            trc_file = trc_file[trc_idx, :]
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

    # calc statistics

    def D_average(self):
        """
        Calculate average diffusion coefficient and mean error per type for 3 and 4 type combination, average over cells
        """
        immobile_tracks, confined_tracks, free_tracks, notype_tracks = [], [], [], []
        for cell in self.cell_trajectories_filtered:
            immobile_tracks_cell, confined_tracks_cell, free_tracks_cell, notype_tracks_cell = [], [], [], []
            for i in cell:
                if i.immobility and not i.analyse_successful and not i.confined:
                    immobile_tracks.append(i)
                    immobile_tracks_cell.append(i)
                if not i.immobility and i.confined and i.analyse_successful:
                    confined_tracks.append(i)
                    confined_tracks_cell.append(i)
                if not i.immobility and not i.confined and i.analyse_successful:
                    free_tracks.append(i)
                    free_tracks_cell.append(i)
                if not i.immobility and not i.analyse_successful:
                    notype_tracks.append(i)
                    notype_tracks_cell.append(i)
            immob_trajectories_cell = [track.D for track in immobile_tracks_cell]
            immobile_D_cell = np.mean(immob_trajectories_cell)
            immobile_dD_cell = np.std(immob_trajectories_cell, ddof=1) / math.sqrt(len(immob_trajectories_cell))
            conf_trajectories_cell = [track.D for track in confined_tracks_cell]
            confined_D_cell = np.mean(conf_trajectories_cell)
            confined_dD_cell = np.std(conf_trajectories_cell, ddof=1) / math.sqrt(len(conf_trajectories_cell))
            free_trajectories_cell = [track.D for track in free_tracks_cell]
            free_D_cell = np.mean(free_trajectories_cell)
            free_dD_cell = np.std(free_trajectories_cell, ddof=1) / math.sqrt(len(free_trajectories_cell))
            notype_trajectories_cell = [track.D for track in notype_tracks_cell]
            notype_D_cell = np.mean(notype_trajectories_cell)
            notype_dD_cell = np.std(notype_trajectories_cell, ddof=1) / math.sqrt(len(notype_trajectories_cell))
            immob_notype_trajectories_cell = immob_trajectories_cell + notype_trajectories_cell
            immob_notype_D_cell = np.mean(immob_notype_trajectories_cell)
            immob_notype_dD_cell = np.std(immob_notype_trajectories_cell, ddof=1) / math.sqrt(len(immob_notype_trajectories_cell))
            D_type_cell = (immobile_D_cell, confined_D_cell, free_D_cell, notype_D_cell, immob_notype_D_cell)
            dD_type_cell = (immobile_dD_cell, confined_dD_cell, free_dD_cell, notype_dD_cell, immob_notype_dD_cell)
            self.D_cell_types.append(D_type_cell)
            self.dD_cell_types.append(dD_type_cell)
        self.D_mean_types = np.nanmean(self.D_cell_types, axis=0)
        self.dD_mean_types = np.nanmean(self.dD_cell_types, axis=0)

    def length_average(self):
        """
        Calculate average track length and mean error per type for 3 and 4 type combination, average over cells
        """
        immobile_tracks, confined_tracks, free_tracks, notype_tracks = [], [], [], []
        for cell in self.cell_trajectories_filtered:
            immobile_tracks_cell, confined_tracks_cell, free_tracks_cell, notype_tracks_cell = [], [], [], []
            for i in cell:
                if i.immobility and not i.analyse_successful and not i.confined:
                    immobile_tracks.append(i)
                    immobile_tracks_cell.append(i)
                if not i.immobility and i.confined and i.analyse_successful:
                    confined_tracks.append(i)
                    confined_tracks_cell.append(i)
                if not i.immobility and not i.confined and i.analyse_successful:
                    free_tracks.append(i)
                    free_tracks_cell.append(i)
                if not i.immobility and not i.analyse_successful:
                    notype_tracks.append(i)
                    notype_tracks_cell.append(i)
            immob_trajectories_l_cell = [track.length_trajectory for track in immobile_tracks_cell]
            immobile_length_cell = np.mean(immob_trajectories_l_cell)
            immobile_dlength_cell = np.std(immob_trajectories_l_cell, ddof=1)/math.sqrt(len(immob_trajectories_l_cell))
            conf_trajectories_l_cell = [track.length_trajectory for track in confined_tracks_cell]
            confined_length_cell = np.mean(conf_trajectories_l_cell)
            confined_dlength_cell = np.std(conf_trajectories_l_cell, ddof=1) / math.sqrt(len(conf_trajectories_l_cell))
            free_trajectories_l_cell = [track.length_trajectory for track in free_tracks_cell]
            free_length_cell = np.mean(free_trajectories_l_cell)
            free_dlength_cell = np.std(free_trajectories_l_cell, ddof=1) / math.sqrt(len(free_trajectories_l_cell))
            notype_trajectories_l_cell = [track.length_trajectory for track in notype_tracks_cell]
            notype_length_cell = np.mean(notype_trajectories_l_cell)
            notype_dlength_cell = np.std(notype_trajectories_l_cell, ddof=1) / math.sqrt(len(notype_trajectories_l_cell))
            immob_notype_trajectories_l_cell = immob_trajectories_l_cell + notype_trajectories_l_cell
            immob_notype_length_cell = np.mean(immob_notype_trajectories_l_cell)
            immob_notype_dlength_cell = np.std(immob_notype_trajectories_l_cell, ddof=1)/math.sqrt(len(immob_notype_trajectories_l_cell))
            length_type_cell = (immobile_length_cell, confined_length_cell, free_length_cell, notype_length_cell, immob_notype_length_cell)
            dlength_type_cell = (immobile_dlength_cell, confined_dlength_cell, free_dlength_cell, notype_dlength_cell, immob_notype_dlength_cell)
            self.length_cell_types.append(length_type_cell)
            self.dlength_cell_types.append(dlength_type_cell)
        self.length_mean_types = np.nanmean(self.length_cell_types, axis=0)
        self.dlength_mean_types = np.nanmean(self.dlength_cell_types, axis=0)

    def type_percentage(self):
        """
        Calculate average type in % and mean error for 3 and 4 type combination, average over cells
        """
        data_selected = True
        self.total_trajectories_filtered = 0
        self.total_trajectories_filtered_cell = []
        for cell_index in range(0, len(self.cell_trajectories_filtered)):
            self.total_trajectories_filtered_cell.append(len(self.cell_trajectories_filtered[cell_index]))
        self.total_trajectories_filtered = np.sum(self.total_trajectories_filtered_cell)
        if self.total_trajectories_filtered == 0:
            data_selected = False
        if data_selected:
            for cell_index, cell in enumerate(self.cell_trajectories_filtered):
                # overview for one cell
                count_immobile_cell, count_confined_cell, count_free_cell, count_not_successful_cell = 0, 0, 0, 0
                for trajectory in cell:
                    if trajectory.immobility and not trajectory.confined and not trajectory.analyse_successful:
                        count_immobile_cell += 1
                    if trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                        count_confined_cell += 1
                    # has to be not confined AND not immobile (otherwise it will count the immobile particles as well)
                    if not trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                        count_free_cell += 1
                    if not trajectory.analyse_successful and not trajectory.immobility:
                        count_not_successful_cell += 1
                count_immob_notype = count_immobile_cell + count_not_successful_cell
                if self.total_trajectories_filtered_cell[cell_index]:
                    ratio_immobile_cell = count_immobile_cell/self.total_trajectories_filtered_cell[cell_index]*100
                    ratio_confined_cell = count_confined_cell/self.total_trajectories_filtered_cell[cell_index]*100
                    ratio_free_cell = count_free_cell/self.total_trajectories_filtered_cell[cell_index]*100
                    ratio_not_successful_cell = count_not_successful_cell/self.total_trajectories_filtered_cell[cell_index]*100
                    ratio_immob_notype = count_immob_notype/self.total_trajectories_filtered_cell[cell_index]*100
                    type_percentages_cell = (ratio_immobile_cell, ratio_confined_cell, ratio_free_cell, ratio_not_successful_cell, ratio_immob_notype)
                else:
                    type_percentages_cell = (0.0, 0.0, 0.0, 0.0, 0.0)  # if all trajectories from a cell are filtered, the type % are 0
                self.percentages_cell_types.append(type_percentages_cell)
            self.type_percentages_mean = np.nanmean(self.percentages_cell_types, axis=0)
            self.type_percentages_error = np.nanstd(self.percentages_cell_types, axis=0, ddof=1)/math.sqrt(len(self.percentages_cell_types))

    # plot diffusion vs frequencies.

    def plot_diffusion_hist_types(self, mean_log_hist_immob, error_log_hist_immob, mean_log_hist_conf,
                                  error_log_hist_conf, mean_log_hist_free, error_log_hist_free,
                                  mean_log_hist_fit_fail, error_log_hist_fit_fail,
                                  mean_log_hist_immob_fil_fail, error_log_hist_immob_fit_fail, merge):
        x = plt.figure()
        plt.subplot(111, xscale="log")
        if merge:
            (_, caps, _) = plt.errorbar(self.hist_diffusion, mean_log_hist_immob_fil_fail, yerr=error_log_hist_immob_fit_fail,
                                        capsize=4, label="relative frequency immobile+notype", ecolor="#4169e1", color="#4169e1")
            for cap in caps:
                cap.set_markeredgewidth(1)
        else:

            (_, caps, _) = plt.errorbar(self.hist_diffusion, mean_log_hist_immob, yerr=error_log_hist_immob, capsize=4,
                                        label="relative frequency immobile", ecolor="#4169e1", color="#4169e1")
            for cap in caps:
                cap.set_markeredgewidth(1)
            (_, caps, _) = plt.errorbar(self.hist_diffusion, mean_log_hist_fit_fail, yerr=error_log_hist_fit_fail,
                                        capsize=4,
                                        label="relative frequency no type", ecolor="#8b008b",
                                        color="#8b008b")
            for cap in caps:
                cap.set_markeredgewidth(1)

        (_, caps, _) = plt.errorbar(self.hist_diffusion, mean_log_hist_conf, yerr=error_log_hist_conf, capsize=4,
                                    label="relative frequency confined", ecolor="#228b22", color="#228b22")
        for cap in caps:
            cap.set_markeredgewidth(1)

        (_, caps, _) = plt.errorbar(self.hist_diffusion, mean_log_hist_free, yerr=error_log_hist_free, capsize=4,
                                    label="relative frequency free", ecolor="#ff8c00", color="#ff8c00")
        for cap in caps:
            cap.set_markeredgewidth(1)

        for c, i in enumerate(self.mean_frequencies_percent):
            if i != 0 and c == 0:
                xlim_min = self.hist_diffusion[c]
                break
            elif i != 0:
                xlim_min = self.hist_diffusion[c-1]
                break
        for c, i in enumerate(np.flip(self.mean_frequencies_percent)):
            if i != 0 and c == 0:
                xlim_max = self.hist_diffusion[len(self.hist_diffusion) - c]
                break
            elif i != 0:
                xlim_max = self.hist_diffusion[len(self.hist_diffusion) - c + 1]
                break

        plt.xlim(xlim_min, xlim_max)
        plt.legend()
        plt.title("Distribution of diffusion coefficients per type")
        plt.ylabel("Normalized relative occurence [%]")
        plt.xlabel("D [\u03BCm\u00b2/s]")
        return x

    def diffusion_hist_types(self, desired_bin_size, merge):
        log_Ds_immob = [[] for i in self.cell_sizes]
        for c, cell in enumerate(self.trajectories_immob_cells_filtered):
            log_Ds_immob_cell = [np.log10(trajectory.D) for trajectory in cell if trajectory.D > 0]
            log_Ds_immob[c] = log_Ds_immob_cell

        log_Ds_conf = [[] for i in self.cell_sizes]
        for c, cell in enumerate(self.trajectories_conf_cells_filtered):
            log_Ds_conf_cell = [np.log10(trajectory.D) for trajectory in cell if trajectory.D > 0]
            log_Ds_conf[c] = log_Ds_conf_cell

        log_Ds_free = [[] for i in self.cell_sizes]
        for c, cell in enumerate(self.trajectories_free_cells_filtered):
            log_Ds_free_cell = [np.log10(trajectory.D) for trajectory in cell if trajectory.D > 0]
            log_Ds_free[c] = log_Ds_free_cell

        log_Ds_notype = [[] for i in self.cell_sizes]
        for c, cell in enumerate(self.trajectories_notype_cells_filtered):
            log_Ds_notype_cell = [np.log10(trajectory.D) for trajectory in cell if trajectory.D > 0]
            log_Ds_notype[c] = log_Ds_notype_cell

        log_Ds_immob_notype = [[] for i in self.cell_sizes]
        for c, cell in enumerate(self.trajectories_immob_notype_cells_filtered):
            log_Ds_immob_notype_cell = [np.log10(trajectory.D) for trajectory in cell if trajectory.D > 0]
            log_Ds_immob_notype[c] = log_Ds_immob_notype_cell

        hist_immob, hist_conf, hist_free, hist_notype, hist_immob_notype = [], [], [], [], []
        for i, cell_size in enumerate(self.cell_sizes):
            hist_immob_cell = self.calc_diffusion_frequencies(log_Ds_immob[i], desired_bin_size, cell_size)
            hist_conf_cell = self.calc_diffusion_frequencies(log_Ds_conf[i], desired_bin_size, cell_size)
            hist_free_cell = self.calc_diffusion_frequencies(log_Ds_free[i], desired_bin_size, cell_size)
            hist_notype_cell = self.calc_diffusion_frequencies(log_Ds_notype[i], desired_bin_size, cell_size)
            hist_immob_notype_cell = self.calc_diffusion_frequencies(log_Ds_immob_notype[i], desired_bin_size, cell_size)
            hist_immob.append(hist_immob_cell[:,1])  # only frequency counts
            hist_conf.append(hist_conf_cell[:, 1])
            hist_free.append(hist_free_cell[:, 1])
            hist_notype.append(hist_notype_cell[:, 1])
            hist_immob_notype.append(hist_immob_notype_cell[:, 1])

        mean_log_hist_immob = np.nanmean(hist_immob, axis=0) * self.normalization_factor
        error_log_hist_immob = np.nanstd(hist_immob, axis=0, ddof=1) / math.sqrt(len(hist_immob)) * self.normalization_factor
        mean_log_hist_conf = np.nanmean(hist_conf, axis=0) * self.normalization_factor
        error_log_hist_conf = np.nanstd(hist_conf, axis=0, ddof=1) / math.sqrt(len(hist_conf)) * self.normalization_factor
        mean_log_hist_free = np.nanmean(hist_free, axis=0) * self.normalization_factor
        error_log_hist_free = np.nanstd(hist_free, axis=0, ddof=1) / math.sqrt(len(hist_free)) * self.normalization_factor
        mean_log_hist_notype = np.nanmean(hist_notype, axis=0) * self.normalization_factor
        error_log_hist_notype = np.nanstd(hist_notype, axis=0, ddof=1) / math.sqrt(len(hist_notype)) * self.normalization_factor
        mean_log_hist_immob_notype = np.nanmean(hist_immob_notype, axis=0) * self.normalization_factor
        error_log_hist_immob_notype = np.nanstd(hist_immob_notype, axis=0, ddof=1) / math.sqrt(len(hist_immob_notype)) * self.normalization_factor

        self.hist_diffusion_immob = [mean_log_hist_immob, error_log_hist_immob]
        self.hist_diffusion_conf = [mean_log_hist_conf, error_log_hist_conf]
        self.hist_diffusion_free = [mean_log_hist_free, error_log_hist_free]
        self.hist_diffusion_notype = [mean_log_hist_notype, error_log_hist_notype]
        self.hist_diffusion_immob_notype = [mean_log_hist_immob_notype, error_log_hist_immob_notype]
        self.diff_fig_types = self.plot_diffusion_hist_types(mean_log_hist_immob, error_log_hist_immob, mean_log_hist_conf,
                                       error_log_hist_conf, mean_log_hist_free, error_log_hist_free,
                                       mean_log_hist_notype, error_log_hist_notype, mean_log_hist_immob_notype,
                                       error_log_hist_immob_notype, merge=False)
        self.diff_fig_types_merge = self.plot_diffusion_hist_types(mean_log_hist_immob, error_log_hist_immob, mean_log_hist_conf,
                                       error_log_hist_conf, mean_log_hist_free, error_log_hist_free,
                                       mean_log_hist_notype, error_log_hist_notype, mean_log_hist_immob_notype,
                                       error_log_hist_immob_notype, merge)

    def run_diffusion_histogram(self, desired_bin_size, x_lim="None", y_lim="None", y_min=0, merge=True):
        """
        :param merge: If True, immobile and notype are handled as one type in a plot and separately in a second plot.
        """
        try:
            float(desired_bin_size)
        except:
            print("Insert a dot separated float or integer as bin size (e.g. 0.1).")
        else:
            if float(desired_bin_size) != 0.0:
                self.clear_attributes()
                self.bin_size = desired_bin_size
                self.determine_max_min_diffusion()
                self.diffusions_log(float(desired_bin_size))
                self.calc_nonlogarithmic_diffusions()
                self.determine_mean_frequency()
                self.calc_mean_error()
                if self.background_trajectories:
                    self.diffusions_log_bg(float(desired_bin_size))
                    self.determine_mean_frequency(is_cell=False)
                    self.calc_mean_error(is_cell=False)
                    self.calc_bg_corrected_freq()
                self.plot_bar_log_bins()
                self.diffusion_hist_types(float(desired_bin_size), merge)
                if self.background_trajectories:
                    self.plot_bar_log_bins_bg_corrected()
                # plot MSD type
                x_lim = None if x_lim == "None" else float(x_lim)
                y_lim = None if y_lim == "None" else float(y_lim)
                self.MSD_types(x_lim, y_lim, y_min, merge)
                # plot cell distributions
                self.plot_sizes_ntracks()
                self.plot_type_distributions()
                self.plot_D_distributions()
                self.plot_length_distributions()
            else:
                print("Bin size can not be zero.")

    # plot MSD per type

    def plot_MSD_types(self, MSD_delta_t_n, y_lim, y_min, merge=False):
        """
        Plot the mean MSD curve per diffusion type. If merge, handle immobile and notype as one class.
        """
        camera_time = self.cell_trajectories[0][0].dt  # camera integration time in s
        delta_t_immob = [camera_time * i for i in range(1, len(self.mean_MSDs_immob) + 1)]
        delta_t_conf = [camera_time * i for i in range(1, len(self.mean_MSDs_conf) + 1)]
        delta_t_free = [camera_time * i for i in range(1, len(self.mean_MSDs_free)+1)]
        delta_t_fit_fail = [camera_time * i for i in range(1, len(self.mean_MSDs_notype) + 1)]
        delta_t_immob_fit_fail = [camera_time * i for i in range(1, len(self.mean_MSDs_immob_notype)+1)]
        self.delta_t_types.append(delta_t_immob); self.delta_t_types.append(delta_t_conf)
        self.delta_t_types.append(delta_t_free); self.delta_t_types.append(delta_t_fit_fail)
        self.delta_t_types.append(delta_t_immob_fit_fail)

        x = plt.figure()
        if merge:
            (_, caps, _) = plt.errorbar(delta_t_immob_fit_fail, self.mean_MSDs_immob_notype, yerr=self.mean_error_MSDs_immob_notype,
                                        capsize=4, label="MSD immob+notype", ecolor="#4169e1", color="#4169e1")  # capsize length of cap
            for cap in caps:
                cap.set_markeredgewidth(1)  # markeredgewidth thickness of cap (vertically)
        else:
            (_, caps, _) = plt.errorbar(delta_t_immob, self.mean_MSDs_immob, yerr=self.mean_error_MSDs_immob, capsize=4,
                                        label="MSD immobile", ecolor="#4169e1", color="#4169e1")  # capsize length of cap
            for cap in caps:
                cap.set_markeredgewidth(1)  # markeredgewidth thickness of cap (vertically)
            (_, caps, _) = plt.errorbar(delta_t_fit_fail, self.mean_MSDs_notype, yerr=self.mean_error_MSDs_notype, capsize=4,
                                        label="MSD no type", ecolor="#8b008b", color="#8b008b")  # capsize length of cap
            for cap in caps:
                cap.set_markeredgewidth(1)  # markeredgewidth thickness of cap (vertically)

        (_, caps, _) = plt.errorbar(delta_t_conf, self.mean_MSDs_conf, yerr=self.mean_error_MSDs_conf, capsize=4,
                                    label="MSD confined", ecolor="#228b22", color="#228b22")  # capsize length of cap
        for cap in caps:
            cap.set_markeredgewidth(1)  # markeredgewidth thickness of cap (vertically)

        (_, caps, _) = plt.errorbar(delta_t_free, self.mean_MSDs_free, yerr=self.mean_error_MSDs_free, capsize=4,
                                    label="MSD free", ecolor="#ff8c00", color="#ff8c00")  # capsize length of cap
        for cap in caps:
            cap.set_markeredgewidth(1)  # markeredgewidth thickness of cap (vertically)

        x_lim_max = None if MSD_delta_t_n is None else MSD_delta_t_n
        plt.xlim(0, x_lim_max)
        plt.ylim(y_min, y_lim)
        plt.legend()
        plt.title("Mean MSD values per type")
        plt.ylabel("Mean MSD [\u03BCm\u00b2]")
        plt.xlabel("Time step [s]")
        return x

    def plot_violine(self, lst, title, x_axis, y_axis, color="cornflowerblue"):
        df = pd.DataFrame(lst)
        fig = plt.figure(figsize=(3, 5))
        ax = sns.violinplot(data=df, color=color, scale="width", cut=0, inner="quartile")
        ax = sns.stripplot(data=df, color="0.25")
        ax.set_ylabel(y_axis)
        ax.set_xlabel(x_axis)
        ax.set_title(title)
        ax.set(xticklabels=[])
        return fig

    def plot_violine_subplots(self, lsts, title, x_axiss, y_axis, colors):
        fig, axs = plt.subplots(1, len(lsts), sharey=True)
        for c, (lst, x_axis, color) in enumerate(zip(lsts, x_axiss, colors)):
            df = pd.DataFrame(lst)
            sns.violinplot(data=df, color=color, scale="width", cut=0, inner="quartile", ax=axs[c])
            sns.stripplot(data=df, color=color, edgecolor=(19/255,27/255,30/255), linewidth=1, ax=axs[c])
            axs[c].set_xlabel(x_axis)
            axs[c].set(xticklabels=[])
        fig.text(0.02, 0.5, y_axis, va="center", rotation="vertical")
        plt.suptitle(title)
        return fig

    def plot_D_distributions(self):
        # average D per cell, global and per type
        all_Ds = []
        for cell in self.cell_trajectories_filtered:
             all_Ds.append(np.mean([trajectory.D for trajectory in cell]))
        # get mean D values per cell and type
        immob_notype_Ds = []
        for cell in self.trajectories_immob_notype_cells_filtered:
            immob_notype_Ds.append(np.mean([trajectory.D for trajectory in cell]))
        conf_Ds = []
        for cell in self.trajectories_conf_cells_filtered:
            conf_Ds.append(np.mean([trajectory.D for trajectory in cell]))
        free_Ds = []
        for cell in self.trajectories_free_cells_filtered:
            free_Ds.append(np.mean([trajectory.D for trajectory in cell]))
        immob_Ds = []
        for cell in self.trajectories_immob_cells_filtered:
             immob_Ds.append(np.mean([trajectory.D for trajectory in cell]))
        notype_Ds = []
        for cell in self.trajectories_notype_cells_filtered:
            notype_Ds.append(np.mean([trajectory.D for trajectory in cell]))
        lsts = [all_Ds, immob_notype_Ds, conf_Ds, free_Ds, immob_Ds, notype_Ds]
        title = "mean diffusion coefficients per cell"
        x_axiss = ["all", "immob+notype", "conf", "free", "immobile", "notype"]
        y_axis = "diffusion coefficient [\u00B5m\u00B2/s]"
        colors = ["cornflowerblue", "#4169e1", "#228b22", "#ff8c00", "#4169e1", "#8b008b"]
        self.diff_cell_fig = self.plot_violine_subplots(lsts, title, x_axiss, y_axis, colors)

    def plot_sizes_ntracks(self):
        # cell sizes
        target = "segments" if self.cells[0].seg_id else "trajectories"  # SEG OR TRACK
        self.cell_sizes_fig = self.plot_violine(self.cell_sizes, "cell sizes", "cells", "cell size [\u00B5m\u00B2]")
        self.n_segments_fig = self.plot_violine([len(tracks) for tracks in self.cell_trajectories_filtered], "number of " + target + " per cell", "cells", "number of " + target)

    def plot_length_distributions(self):
        # average lengths per cell, global and per type
        target = "segment" if self.cells[0].seg_id else "trajectory"  # SEG OR TRACK
        all_lengths = []
        for cell in self.cell_trajectories_filtered:
             all_lengths.append(np.mean([trajectory.length_trajectory for trajectory in cell]))
        # get mean D values per cell and type
        immob_notype_lengths = []
        for cell in self.trajectories_immob_notype_cells_filtered:
            immob_notype_lengths.append(np.mean([trajectory.length_trajectory for trajectory in cell]))
        conf_lengths = []
        for cell in self.trajectories_conf_cells_filtered:
            conf_lengths.append(np.mean([trajectory.length_trajectory for trajectory in cell]))
        free_lengths = []
        for cell in self.trajectories_free_cells_filtered:
            free_lengths.append(np.mean([trajectory.length_trajectory for trajectory in cell]))
        immob_lengths = []
        for cell in self.trajectories_immob_cells_filtered:
             immob_lengths.append(np.mean([trajectory.length_trajectory for trajectory in cell]))
        notype_lengths = []
        for cell in self.trajectories_notype_cells_filtered:
            notype_lengths.append(np.mean([trajectory.length_trajectory for trajectory in cell]))

        lsts = [all_lengths, immob_notype_lengths, conf_lengths, free_lengths, immob_lengths, notype_lengths]
        title = "mean " + target + " lengths per cell"
        x_axiss = ["all", "immob+notype", "conf", "free", "immobile", "notype"]
        y_axis = "length [frames]"
        colors = ["cornflowerblue", "#4169e1", "#228b22", "#ff8c00", "#4169e1", "#8b008b"]
        self.lengths_cell_fig = self.plot_violine_subplots(lsts, title, x_axiss, y_axis, colors)

    def plot_type_distributions(self):
        immob_notype_percentages = [cell[4] for cell in self.percentages_cell_types]
        conf_percentages = [cell[1] for cell in self.percentages_cell_types]
        free_percentages = [cell[2] for cell in self.percentages_cell_types]
        immob_percentages = [cell[0] for cell in self.percentages_cell_types]
        notype_percentages = [cell[3] for cell in self.percentages_cell_types]

        lsts = [immob_notype_percentages, conf_percentages, free_percentages, immob_percentages, notype_percentages]
        title = "mean type percentages per cell"
        x_axiss = ["immob+notype", "conf", "free", "immobile", "notype"]
        y_axis = "type [\u0025]"
        colors = ["#4169e1", "#228b22", "#ff8c00", "#4169e1", "#8b008b"]
        self.type_percentages_cell_fig = self.plot_violine_subplots(lsts, title, x_axiss, y_axis, colors)

    def calc_mean_error_different_lengths(self, arrays):
        """
        Calculate the mean and mean error of a list of arrays with varying array lengths.
        [[1,2], [1,3,4]] -> [1,2.5,4]
        :param arrays: List of arrays
        :return: List of mean and mean error
        """
        # determine the maximum array length
        max_length = 0
        for i in arrays:
            max_length = len(i) if len(i) > max_length else max_length

        arrays_sorted = []
        for i in range(max_length):
            step = []
            for j in arrays:
                try:
                    step.append(j[i])
                except IndexError:
                    pass
            arrays_sorted.append(step)
        mean = [np.nanmean(i) for i in arrays_sorted]
        mean_error = [np.nanstd(i, ddof=1) / math.sqrt(len(i)) for i in arrays_sorted]
        return mean, mean_error

    def MSD_types(self, MSD_delta_t_n, y_lim, y_min, merge):
        average_MSDs_immob = []
        average_MSDs_immob_error = []
        for cell in self.trajectories_immob_cells_filtered:
            MSDs_cell = [trajectory.MSDs for trajectory in cell]
            mean_MSDs_cell, error_MSDs_cell = self.calc_mean_error_different_lengths(MSDs_cell)
            average_MSDs_immob.append(mean_MSDs_cell)
            average_MSDs_immob_error.append(error_MSDs_cell)
        self.mean_MSDs_immob, _ = self.calc_mean_error_different_lengths(average_MSDs_immob)
        self.mean_error_MSDs_immob, _ = self.calc_mean_error_different_lengths(average_MSDs_immob_error)

        average_MSDs_conf = []
        average_MSDs_conf_error = []
        for cell in self.trajectories_conf_cells_filtered:
            MSDs_cell = [trajectory.MSDs for trajectory in cell]
            mean_MSDs_cell, error_MSDs_cell = self.calc_mean_error_different_lengths(MSDs_cell)
            average_MSDs_conf.append(mean_MSDs_cell)
            average_MSDs_conf_error.append(error_MSDs_cell)
        self.mean_MSDs_conf, _ = self.calc_mean_error_different_lengths(average_MSDs_conf)
        self.mean_error_MSDs_conf, _ = self.calc_mean_error_different_lengths(average_MSDs_conf_error)

        average_MSDs_free = []
        average_MSDs_free_error = []
        for cell in self.trajectories_free_cells_filtered:
            MSDs_cell = [trajectory.MSDs for trajectory in cell]
            mean_MSDs_cell, error_MSDs_cell = self.calc_mean_error_different_lengths(MSDs_cell)
            average_MSDs_free.append(mean_MSDs_cell)
            average_MSDs_free_error.append(error_MSDs_cell)
        self.mean_MSDs_free, _ = self.calc_mean_error_different_lengths(average_MSDs_free)
        self.mean_error_MSDs_free, _ = self.calc_mean_error_different_lengths(average_MSDs_free_error)

        average_MSDs_notype = []
        average_MSDs_notype_error = []
        for cell in self.trajectories_notype_cells_filtered:
            MSDs_cell = [trajectory.MSDs for trajectory in cell]
            mean_MSDs_cell, error_MSDs_cell = self.calc_mean_error_different_lengths(MSDs_cell)
            average_MSDs_notype.append(mean_MSDs_cell)
            average_MSDs_notype_error.append(error_MSDs_cell)
        self.mean_MSDs_notype, _ = self.calc_mean_error_different_lengths(average_MSDs_notype)
        self.mean_error_MSDs_notype, _ = self.calc_mean_error_different_lengths(average_MSDs_notype_error)

        average_MSDs_immob_notype = []
        average_MSDs_immob_notype_error = []
        for cell_immob, cell_notype in zip(self.trajectories_immob_cells_filtered, self.trajectories_notype_cells_filtered):
            MSDs_cell_immob = [trajectory.MSDs for trajectory in cell_immob]
            MSDs_cell_notype = [trajectory.MSDs for trajectory in cell_notype]
            MSDs_cell = MSDs_cell_immob + MSDs_cell_notype
            mean_MSDs_cell, error_MSDs_cell = self.calc_mean_error_different_lengths(MSDs_cell)
            average_MSDs_immob_notype.append(mean_MSDs_cell)
            average_MSDs_immob_notype_error.append(error_MSDs_cell)
        self.mean_MSDs_immob_notype, _ = self.calc_mean_error_different_lengths(average_MSDs_immob_notype)
        self.mean_error_MSDs_immob_notype, _ = self.calc_mean_error_different_lengths(average_MSDs_immob_notype_error)

        if merge:
            self.MSD_fig_types = self.plot_MSD_types(MSD_delta_t_n, y_lim, y_min, merge=False)
            self.MSD_fig_types_merge = self.plot_MSD_types(MSD_delta_t_n, y_lim, y_min, merge=merge)
        else:
            self.plot_MSD_types(MSD_delta_t_n, y_lim, y_min)

    def clear_attributes(self):
        """
        If one restarts filtering, these attributes are empty (otherwise they would append).
        """
        self.hist_log_Ds = []  # histograms (logD vs freq) from all cells as np arrays in this list
        self.hist_diffusion = []  # diffusions from histogram calculation, transformed back -> 10^-(log10(D))
        self.hist_diffusion_immob = []  # first entry frequencies, 2nd entry errors of diffusion hist per type
        self.hist_diffusion_conf = []
        self.hist_diffusion_free = []
        self.hist_diffusion_notype = []
        self.mean_frequencies = []  # mean frequencies, size corrected
        self.mean_error = []
        self.hist_log_Ds_bg = []
        self.mean_frequencies_bg = []
        self.mean_error_bg = []
        self.mean_MSDs_immob = []
        self.mean_error_MSDs_immob = []
        self.mean_MSDs_conf = []
        self.mean_error_MSDs_conf = []
        self.mean_MSDs_free = []
        self.mean_error_MSDs_free = []
        self.mean_MSDs_notype = []
        self.mean_error_MSDs_notype = []
        self.mean_MSDs_immob_notype = []
        self.mean_error_MSDs_immob_notype = []
        self.delta_t_types = []

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
        for bg_index in range(0, len(self.background_trajectories_filtered)):
            log_Ds = np.zeros(len(self.background_trajectories_filtered[bg_index]))
            bg_size = self.bg_sizes[bg_index]
            for trajectory_index in range(0, len(self.background_trajectories_filtered[bg_index])):
                log_Ds[trajectory_index] = np.log10(self.background_trajectories_filtered[bg_index][trajectory_index].D)
            self.calc_diffusion_frequencies(log_Ds, desired_bin_size, bg_size, is_cell=False)
                    
    def diffusions_log(self, desired_bin_size):
        """
        For each cell initialize histogram with cell size and target array.
        """
        for cell_index in range(0, len(self.cell_trajectories_filtered)):
            log_Ds = np.zeros(len(self.cell_trajectories_filtered[cell_index]))
            cell_size = self.cell_sizes[cell_index]
            for trajectory_index in range(0, len(self.cell_trajectories_filtered[cell_index])):
                if self.cell_trajectories_filtered[cell_index][trajectory_index].D > 0:
                    log_Ds[trajectory_index] = np.log10(self.cell_trajectories_filtered[cell_index][trajectory_index].D)
            log_Ds = np.delete(log_Ds, np.where(log_Ds == 0))
            self.calc_diffusion_frequencies(log_Ds, desired_bin_size, cell_size) 
        
    def calc_diffusion_frequencies(self, log_diff, desired_bin, size, is_cell=True):
        """
        :param log_diff: np array with log10(D) of one cell.
        :param size: cell size.
        :param desired_bin: bin size.
        :param is_cell: If True, treat as cell, else as background.
        """
        # min & max determined by diffusions_log_complete function
        min_bin = np.ceil(-6/desired_bin)*desired_bin
        max_bin = np.ceil(2/desired_bin)*desired_bin 
        bin_size = int(np.ceil((max_bin - min_bin)/desired_bin))
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
        return log_diffusion_hist
        
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
            self.diffusion_frequencies = self.create_np_array(np.shape(self.hist_log_Ds)[1],
                                                              len(self.cell_trajectories_filtered))
            for i in range (0, len(self.cell_trajectories_filtered)):
                self.diffusion_frequencies[:,i] = self.hist_log_Ds[i][:,1]
            self.mean_frequencies = self.calc_mean_frequencies(self.diffusion_frequencies)
            self.normalization_factor = 100/np.sum(self.mean_frequencies)
            self.diffusion_frequencies_norm = self.diffusion_frequencies * self.normalization_factor
            self.mean_frequencies_percent = self.mean_frequencies * self.normalization_factor
        else:
            self.diffusion_frequencies_bg = self.create_np_array(np.shape(self.hist_log_Ds_bg)[1],
                                                                 len(self.background_trajectories_filtered))
            for i in range(0, len(self.background_trajectories_filtered)):
                self.diffusion_frequencies_bg[:,i] = self.hist_log_Ds_bg[i][:,1]
            self.mean_frequencies_bg = self.calc_mean_frequencies(self.diffusion_frequencies_bg)
        
    def calc_mean_error(self, is_cell=True):
        """
        Standard deviation (N-1) divided by square root of number of elements.
        Normalize an array (sum of elements = 1) and represent is in percent (*100).
        """
        if is_cell:
            self.mean_error =  np.std(self.diffusion_frequencies, axis=1,
                                      ddof=1)/(np.shape(self.diffusion_frequencies)[1])**(1/2)
            self.mean_error_percent = self.mean_error * self.normalization_factor 
        else:
            self.mean_error_bg = np.std(self.diffusion_frequencies_bg, axis=1,
                                        ddof=1)/(np.shape(self.diffusion_frequencies_bg)[1])**(1/2)
            
    def calc_bg_corrected_freq(self):
        """
        Subtract background frequencies of cell frequencies. if < 0 -> set to 0.
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
        (_, caps, _) = plt.errorbar(self.hist_diffusion, self.mean_frequencies_percent, yerr=self.mean_error_percent,
                                    capsize=4, label="relative frequency")
        for cap in caps:
            cap.set_markeredgewidth(1)

        for c, i in enumerate(self.mean_frequencies_percent):
            if i != 0 and c == 0:
                xlim_min = self.hist_diffusion[c]
                break
            elif i != 0:
                xlim_min = self.hist_diffusion[c-1]
                break
        for c, i in enumerate(np.flip(self.mean_frequencies_percent)):
            if i != 0 and c == 0:
                xlim_max = self.hist_diffusion[len(self.hist_diffusion) - c]
                break
            elif i != 0:
                xlim_max = self.hist_diffusion[len(self.hist_diffusion) - c + 1]
                break

        plt.xlim(xlim_min, xlim_max)
        plt.legend()
        plt.title("Distribution of diffusion coefficients")
        plt.ylabel("normalized relative occurence [%]")
        plt.xlabel("D [\u03BCm\u00b2/s]")
    
    def save_diff_fig(self, directory, folder_name):
        for fig, name in zip(self.diff_fig, self.fig_name):  # diffusion hist / bg hist
            fig.savefig(directory + "\\"+ folder_name +  "\\" + name + ".svg")
        # diffusion hist type
        self.diff_fig_types.savefig(directory + "\\" + folder_name + "\\" + "diffusion_hist_type" + ".svg")
        self.diff_fig_types_merge.savefig(directory + "\\" + folder_name + "\\" + "diffusion_hist_type_merge" + ".svg")
        # MSD type
        self.MSD_fig_types.savefig(directory + "\\" + folder_name + "\\" + "MSD_type" + ".svg")
        self.MSD_fig_types_merge.savefig(directory + "\\" + folder_name + "\\" + "MSD_type_merge" + ".svg")
        self.cell_sizes_fig.savefig(directory + "\\" + folder_name + "\\" + "cell_sizes" + ".svg", bbox_inches="tight")
        target = "segments" if self.cells[0].seg_id else "trajectories"  # SEG OR TRACK
        self.n_segments_fig.savefig(directory + "\\" + folder_name + "\\" + "number_" + target + ".svg", bbox_inches="tight")
        self.type_percentages_cell_fig.savefig(directory + "\\" + folder_name + "\\" + "type_percentages_cell" + ".svg", bbox_inches="tight")
        self.diff_cell_fig.savefig(directory + "\\" + folder_name + "\\" + "diffusion_cell" + ".svg", bbox_inches="tight")
        self.lengths_cell_fig.savefig(directory + "\\" + folder_name + "\\" + target + "_lengths_cell" + ".svg", bbox_inches="tight")
        #self.lengths_cell_fig.savefig(directory + "\\" + folder_name + "\\" + target + "_lengths_cell" + ".pdf", format="pdf", transparent=True, bbox_inches="tight")

    def plot_bar_log_bins_bg_corrected(self):
        self.diff_fig.append(plt.figure())
        self.fig_name.append("diffusion_hist_background")
        plt.subplot(111, xscale="log")
        (_, caps, _) = plt.errorbar(self.hist_diffusion, self.corrected_frequencies_percent,
                                    yerr=self.corrected_frequencies_error_percent,
                                    capsize=4, label="relative frequency")
        for cap in caps:
            cap.set_markeredgewidth(1)
        plt.xlim(self.min_D, self.max_D)
        plt.legend()
        plt.title("Distribution of diffusion coefficients (background corrected)")
        plt.ylabel("Normalized relative occurence [%]")
        plt.xlabel("D [\u03BCm\u00b2/s]")

    def calc_mean_frequencies(self, np_array):
        """
        Determine mean value over the horizontal entries of an np.array.
        -> [2,2][6,4] = [4,3].
        :param np_array: Np array to build mean over.
        """
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
                if (mean_MSD_0+(4/3)*mean_D*dt)/dof > 0:
                    sigma_dyn = math.sqrt((mean_MSD_0+(4/3)*mean_D*dt)/dof)
                else:
                    sigma_dyn = "negative argument in square root, error not calculable for this target"
                self.sigma_dyns.append(sigma_dyn)      
