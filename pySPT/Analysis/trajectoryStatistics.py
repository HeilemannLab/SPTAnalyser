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
        self.cell_trajectories_filtered = []  # deep copy of original cell trajectories
        self.cell_trajectories_filtered_new = []
        self.cell_trajectories_filtered_new_new = []
        self.cell_trajectories_index = []
        self.cell_trajectories_index_new = []
        self.cell_trajectories_filtered_index = []  # deep copy of original cell trajectories index
        self.backgrounds = []  # background objects
        self.background_trajectories = [] # [[],[]] contains list of cells, cells contain trajectories
        self.background_trajectories_filtered = []  # deep copy of original cell trajectories
        self.background_trajectories_filtered_new = []
        self.background_trajectories_filtered_new_new = []
        self.background_trajectories_index = []
        self.background_trajectories_index_new = []
        self.background_trajectories_filtered_index = []  # deep copy of original cell trajectories index
        self.filter_settings = []  # list with boolean values for immob, confined, free, analyse success, not success
        self.filtered_trc_files = []  # contains filtered trc files of cells
        self.min_D = math.inf
        self.max_D = - math.inf
        self.bin_size = 0.1  # bin size for freq vs D plot
        self.total_trajectories = 0  # amount of trajectories in data set
        self.total_trajectories_cell = []  # list with amount of total trajectories per cell 
        self.total_trajectories_filtered = 0  # amount of trajectories in data set after filter
        self.total_trajectories_filtered_cell = []  # list with amount of total trajectories per cell after filter
        self.cell_type_count = []  # tupels with (immobile, confined, free %) for each cell
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
                                        # for the filter setting trajectory min length
                                        
        self.filtered_trajectories = []

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
        
    def default_statistics(self):
        self.cell_trajectories_filtered = []  # deep copy of original cell trajectories
        self.cell_trajectories_filtered_new = []
        self.cell_trajectories_filtered_new_new = []
        self.cell_trajectories_index = []
        self.cell_trajectories_index_new = []
        self.cell_trajectories_filtered_index = []  # deep copy of original cell trajectories index
        self.background_trajectories_filtered_new = []
        self.background_trajectories_filtered_new_new = []
        self.background_trajectories_index = []
        self.background_trajectories_index_new = []
        self.background_trajectories_filtered = []  # deep copy of original cell trajectories
        self.background_trajectories_filtered_index = []  # deep copy of original cell trajectories index
        self.total_trajectories = 0  # amount of trajectories in data set
        self.total_trajectories_cell = []  # list with amount of total trajectories per cell 
        
    def run_statistics(self, min_length, max_length, min_D, max_D, filter_immob, filter_confined,
                       filter_free, filter_analyse_successful, filter_analyse_not_successful):
        """
        Initialize filtered lists, filter functions and percentage function.
        """
        start = time.time()
        self.default_statistics()
        self.get_index()
        self.create_init_filter_lst()
        self.filter_settings = [filter_immob, filter_confined, filter_free, filter_analyse_successful, filter_analyse_not_successful]
        self.calc_amount_trajectories()
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
        end1 = time.time()
        if self.background_trajectories:
            self.filter_thresholds(min_length, max_length, min_D, max_D, bg=True)
        else:
            self.filter_thresholds(min_length, max_length, min_D, max_D)
        self.filter_type_new(filter_immob, filter_confined, filter_free,
                         filter_analyse_successful, filter_analyse_not_successful)
        self.create_index_lst()
        print("New initialization took {} s".format(time.time()-start))
        start2 = time.time()
        self.filter_length(min_length, max_length)
        self.filter_D(min_D, max_D)
        self.filter_type(filter_immob, filter_confined, filter_free,
                         filter_analyse_successful, filter_analyse_not_successful)
        self.filter_cell_trc()
        print("Initialization took {} s".format(end1-start+time.time()-start2))

        if filter_analyse_successful and not filter_analyse_not_successful:
            print("Filter for type determination successful only.")
        elif filter_analyse_not_successful and not filter_analyse_successful:
            print("Filter for type determination not successful only.")
        print("%.1f %% are immobile" %(self.type_percentage()[0]))
        print("%.1f %% are confined" %(self.type_percentage()[1]))
        print("%.1f %% are free" %(self.type_percentage()[2]))        
        if self.total_trajectories_filtered == 0:
            print("The selection excludes all data.")
        print("Trajectories included:", self.total_trajectories_filtered)
        print("Trajectories excluded:", self.total_trajectories - self.total_trajectories_filtered)
        print("cell types %", self.cell_type_count)
        print("cell count", self.total_trajectories_filtered_cell)
        print("index", self.cell_trajectories_filtered_index)
        print("background index", self.background_trajectories_filtered_index)
        
    def create_init_filter_lst(self):
        """
        Create copy of initial cell trajectories & index list.
        """
        self.cell_trajectories_filtered = copy.deepcopy(self.cell_trajectories)
        self.cell_trajectories_filtered_index = copy.deepcopy(self.cell_trajectories_index)
        if self.background_trajectories:
            self.background_trajectories_filtered = copy.deepcopy(self.background_trajectories)
            self.background_trajectories_filtered_index = copy.deepcopy(self.background_trajectories_index)

    def get_index(self):
        """
        Each cell is a list, create i lists for the cells. In the lists append the numbers of elements = trajectory numbers.
        Starting with 1.
        """
        for cell in range(0, len(self.cell_trajectories)):
            i = []
            for trajectory in range(0, len(self.cell_trajectories[cell])):
                i.append(trajectory+1)  # trajectory numbering starts with 1
            self.cell_trajectories_index.append(i)
        if self.background_trajectories:
            for background in range(0, len(self.background_trajectories)):
                i = []
                for trajectory in range(0, len(self.background_trajectories[background])):
                    i.append(trajectory+1)  # trajectory numbering starts with 1
                self.background_trajectories_index.append(i)
        
    def plot_trajectory(self, cell, number):
        cell = int(cell) - 1
        number = int(number) -1
        self.cell_trajectories[cell][number].plot_particle()

    def get_max_length(self):
        """
        Max length of trajectory.
        :return: int of max length.
        """
        max_length = - math.inf
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if max_length < trajectory.length_trajectory:
                    max_length = trajectory.length_trajectory
        return int(max_length)

    def get_min_length(self):
        """
        Max length of trajectory.
        :return: int of max length.
        """
        min_length = math.inf
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if min_length > trajectory.length_trajectory:
                    min_length = trajectory.length_trajectory
        return int(min_length)
        
    def filter_length(self, min_length, max_length): # max_length=default
        """
        Filter by length of trajectory (PALMTracer min_length = 20). Trajectories in the range between
        min_length and max_length will be filtered -> min_length <= trajectory <= max_length.
        :param min_length: minimal length of trajectories to be considered.
        :param max_length: by default None -> None will be evaluated as the max length trajectory
        """
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if trajectory.length_trajectory < int(min_length) or trajectory.length_trajectory > int(max_length):
                    self.crop_lst(cell, trajectory)   
        if self.background_trajectories_filtered:
            for background in range(0, len(self.background_trajectories_filtered)):
                for trajectory in self.background_trajectories_filtered[background]:
                    if trajectory.length_trajectory < int(min_length) or trajectory.length_trajectory > int(max_length):
                        self.crop_lst(background, trajectory, is_cell=False) 

    def crop_lst(self, cell, trajectory, is_cell=True):
        """
        Delete trajectory and its index.
        """
        if not is_cell:
            trajectory_index = self.background_trajectories_filtered[cell].index(trajectory)
            self.background_trajectories_filtered[cell] = self.background_trajectories_filtered[cell][0:trajectory_index] + \
            self.background_trajectories_filtered[cell][trajectory_index+1:]
            self.background_trajectories_filtered_index[cell] = self.background_trajectories_filtered_index[cell][0:trajectory_index] + \
            self.background_trajectories_filtered_index[cell][trajectory_index+1:] 
        else:
            trajectory_index = self.cell_trajectories_filtered[cell].index(trajectory)
            self.cell_trajectories_filtered[cell] = self.cell_trajectories_filtered[cell][0:trajectory_index] + \
            self.cell_trajectories_filtered[cell][trajectory_index+1:]
            self.cell_trajectories_filtered_index[cell] = self.cell_trajectories_filtered_index[cell][0:trajectory_index] + \
            self.cell_trajectories_filtered_index[cell][trajectory_index+1:] 
            
    def filter_thresholds(self, min_length, max_length, min_diff, max_diff, bg=False):
        for cell_index in range(len(self.cell_trajectories)):
            filtered_cell = [trajectory for trajectory in self.cell_trajectories[cell_index] if trajectory.length_trajectory >= min_length 
                        and trajectory.length_trajectory <= max_length
                        and trajectory.D >= min_diff and trajectory.D <= max_diff]
            self.cell_trajectories_filtered_new.append(filtered_cell)
        if bg:
            for bg_index in range(len(self.background_trajectories)):
                filtered_bg = [trajectory for trajectory in self.background_trajectories[bg_index] if trajectory.length_trajectory >= min_length 
                            and trajectory.length_trajectory <= max_length
                            and trajectory.D >= min_diff and trajectory.D <= max_diff]
                self.background_trajectories_filtered_new.append(filtered_bg)
        
    def filter_type_new(self, filter_immob, filter_confined, filter_free,
                    filter_type_successful, filter_type_not_successful):
        """
        :param filter_immob: if checked in JNB -> True -> filter data includes immob trajectories, else False -> will be neglected.
        """
        # if filter was not selected, value is False
        self.cell_trajectories_filtered_new_new = []
        for cell_index in range(len(self.cell_trajectories)):
            filtered_cell = []
            filtered_cell_type = []
            if filter_immob:
                filtered_cell_immob = [trajectory for trajectory in self.cell_trajectories_filtered_new[cell_index] if trajectory.immobility
                                   and not trajectory.confined and not trajectory.analyse_successful]
                filtered_cell_type.extend(filtered_cell_immob)
            if filter_confined:
                filtered_cell_confined = [trajectory for trajectory in self.cell_trajectories_filtered_new[cell_index] if trajectory.confined
                                   and not trajectory.immobility and trajectory.analyse_successful]
                filtered_cell_type.extend(filtered_cell_confined)
            if filter_free:
                filtered_cell_free = [trajectory for trajectory in self.cell_trajectories_filtered_new[cell_index] if not trajectory.confined
                                   and not trajectory.immobility and trajectory.analyse_successful]
                filtered_cell_type.extend(filtered_cell_free)
            if filter_type_successful:
                filtered_cell_type_successful = [trajectory for trajectory in filtered_cell_type if trajectory.analyse_successful
                                                 or trajectory.immobility]
                filtered_cell.extend(filtered_cell_type_successful)
            if filter_type_not_successful:
                filtered_cell_type_not_successful = [trajectory for trajectory in filtered_cell_type if not trajectory.analyse_successful
                                                     and not trajectory.immobility]
                filtered_cell.extend(filtered_cell_type_not_successful)
            self.cell_trajectories_filtered_new_new.append(filtered_cell)
        #print("cell count new", len(self.cell_trajectories_filtered_new_new[0]), len(self.cell_trajectories_filtered_new_new[1]))
        if self.background_trajectories:
            self.background_trajectories_filtered_new_new = []
            for bg_index in range(len(self.background_trajectories)):
                filtered_bg = []
                filtered_bg_type = []
                if filter_immob:
                    filtered_cell_immob = [trajectory for trajectory in self.background_trajectories_filtered_new[bg_index] if trajectory.immobility
                                       and not trajectory.confined and not trajectory.analyse_successful]
                    filtered_bg_type.extend(filtered_cell_immob)
                if filter_confined:
                    filtered_cell_confined = [trajectory for trajectory in self.background_trajectories_filtered_new[bg_index] if trajectory.confined
                                       and not trajectory.immobility and trajectory.analyse_successful]
                    filtered_bg_type.extend(filtered_cell_confined)
                if filter_free:
                    filtered_cell_free = [trajectory for trajectory in self.background_trajectories_filtered_new[bg_index] if not trajectory.confined
                                       and not trajectory.immobility and trajectory.analyse_successful]
                    filtered_bg_type.extend(filtered_cell_free)
                if filter_type_successful:
                    filtered_cell_type_successful = [trajectory for trajectory in filtered_bg_type if trajectory.analyse_successful
                                                     or trajectory.immobility]
                    filtered_bg.extend(filtered_cell_type_successful)
                if filter_type_not_successful:
                    filtered_cell_type_not_successful = [trajectory for trajectory in filtered_bg_type if not trajectory.analyse_successful
                                                         and not trajectory.immobility]
                    filtered_bg.extend(filtered_cell_type_not_successful)
                self.background_trajectories_filtered_new_new.append(filtered_bg)
        #print("bg count new", len(self.background_trajectories_filtered_new_new[0]), len(self.background_trajectories_filtered_new_new[1]))
            
# =============================================================================
#         for cell_index in range(len(self.cell_trajectories)):
#             for trajectory in range(len(self.cell_trajectories_filtered_new_new[cell_index])):
#                 one_trajectory = self.cell_trajectories_filtered_new_new[cell_index][trajectory]
#                 print(cell_index, one_trajectory.trajectory_number, one_trajectory.immobility, one_trajectory.confined, not one_trajectory.confined,
#                       one_trajectory.analyse_successful, not one_trajectory.analyse_successful)
# =============================================================================
        
    def create_index_lst(self):
        for cell_index in range(len(self.cell_trajectories)):
            trajectory_index = []
            for trajectory in self.cell_trajectories_filtered_new_new[cell_index]:
                trajectory_index.append(trajectory.trajectory_number)
            trajectory_index_sorted = trajectory_index.sort()
# =============================================================================
#             print("index", trajectory_index, type(trajectory_index))
#             print("sorted", trajectory_index_sorted)
# =============================================================================
            self.cell_trajectories_index_new.append(trajectory_index)
# =============================================================================
#         print("new cell index", self.cell_trajectories_index_new)
# =============================================================================
        if self.background_trajectories:
            for bg_index in range(len(self.background_trajectories)):
                trajectory_index_bg = []
                for trajectory in self.background_trajectories_filtered_new_new[bg_index]:
                    trajectory_index_bg.append(trajectory.trajectory_number)
                #trajectory_index_bg_sorted = trajectory_index_bg.sort()
                self.background_trajectories_index_new.append(trajectory_index_bg)
# =============================================================================
#             print("new bg index", self.background_trajectories_index_new)
# =============================================================================
        
    def filter_type(self, filter_immob, filter_confined, filter_free,
                    filter_type_successful, filter_type_not_successful):
        """
        :param filter_immob: if checked in JNB -> True -> filter data includes immob trajectories, else False -> will be neglected.
        """
        # if filter was not selected, value is False
        if not filter_immob:
            for cell in range(0, len(self.cell_trajectories_filtered)):
                for trajectory in self.cell_trajectories_filtered[cell]:
                    # meaning that one wants to get rid of all immobile particles -> crop them out of the lists.
                    if trajectory.immobility and not trajectory.confined and not trajectory.analyse_successful:
                        self.crop_lst(cell, trajectory)
        if not filter_confined:
            for cell in range(0, len(self.cell_trajectories_filtered)):
                for trajectory in self.cell_trajectories_filtered[cell]:
                    if trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                        self.crop_lst(cell, trajectory)
        if not filter_free:
            for cell in range(0, len(self.cell_trajectories_filtered)):
                for trajectory in self.cell_trajectories_filtered[cell]:
                    if not trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                        self.crop_lst(cell, trajectory)
        if not filter_type_successful:  # only without successful type det -> immob False & analyse success False
            for cell in range(0, len(self.cell_trajectories_filtered)):
                for trajectory in self.cell_trajectories_filtered[cell]:
                    if trajectory.analyse_successful or trajectory.immobility:  # -> crop immob true or analyse successful true
                        self.crop_lst(cell, trajectory)
        if not filter_type_not_successful:  # only with successful type det -> immob = True or analyse successful
            for cell in range(0, len(self.cell_trajectories_filtered)):
                for trajectory in self.cell_trajectories_filtered[cell]:
                    if not trajectory.analyse_successful and not trajectory.immobility:  # -> crop analyse success false; immobility false
                        self.crop_lst(cell, trajectory)
        
        if self.background_trajectories_filtered:
            if not filter_immob:
                for background in range(0, len(self.background_trajectories_filtered)):
                    for trajectory in self.background_trajectories_filtered[background]:
                        # meaning that one wants to get rid of all immobile particles -> crop them out of the lists.
                        if trajectory.immobility and not trajectory.confined and not trajectory.analyse_successful:
                            self.crop_lst(background, trajectory, is_cell=False)
            if not filter_confined:
                for background in range(0, len(self.background_trajectories_filtered)):
                    for trajectory in self.background_trajectories_filtered[background]:
                        if trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                            self.crop_lst(background, trajectory, is_cell=False)
            if not filter_free:
                for background in range(0, len(self.background_trajectories_filtered)):
                    for trajectory in self.background_trajectories_filtered[background]:
                        if not trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                            self.crop_lst(background, trajectory, is_cell=False)
            if not filter_type_successful:
                for background in range(0, len(self.background_trajectories_filtered)):
                    for trajectory in self.background_trajectories_filtered[background]:
                        if trajectory.analyse_successful or trajectory.immobility:
                            self.crop_lst(background, trajectory, is_cell=False)
            if not filter_type_not_successful:
                for background in range(0, len(self.background_trajectories_filtered)):
                    for trajectory in self.background_trajectories_filtered[background]:
                        if not trajectory.analyse_successful and not trajectory.immobility:
                            self.crop_lst(background, trajectory, is_cell=False)
    
    def get_max_D(self):
        max_D = - math.inf
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if max_D < trajectory.D:
                    max_D = trajectory.D
        return float(max_D)
    
    def get_min_D(self):
        min_D = math.inf
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if min_D > trajectory.D:
                    min_D = trajectory.D
        return float(min_D)
        
    def filter_D(self, min_D, max_D):
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if trajectory.D < min_D or trajectory.D > max_D:
                    self.crop_lst(cell, trajectory)
        if self.background_trajectories_filtered:
            for background in range(0, len(self.background_trajectories_filtered)):
                for trajectory in self.background_trajectories_filtered[background]:
                    if trajectory.D < min_D or trajectory.D > max_D:
                        self.crop_lst(background, trajectory, is_cell=False)
        
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
            count_immobile = 0
            count_confined = 0
            count_free = 0
            for cell in self.cell_trajectories_filtered:
                count_immobile_cell = 0
                count_confined_cell = 0
                count_free_cell = 0
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
                cell_index = self.cell_trajectories_filtered.index(cell)
                if self.total_trajectories_filtered_cell[cell_index]:
                    ratio_immobile_cell = count_immobile_cell/self.total_trajectories_filtered_cell[cell_index]*100
                    ratio_confined_cell = count_confined_cell/self.total_trajectories_filtered_cell[cell_index]*100
                    ratio_free_cell = count_free_cell/self.total_trajectories_filtered_cell[cell_index]*100
                    cell_types_percent = (ratio_immobile_cell, ratio_confined_cell, ratio_free_cell)
                else:
                    cell_types_percent = (0.0,0.0,0.0)  # if all trajectories from a cell are filtered, the type % are 0
                self.cell_type_count.append(cell_types_percent)
            ratio_immobile = count_immobile/self.total_trajectories_filtered*100
            ratio_confined = count_confined/self.total_trajectories_filtered*100
            ratio_free = count_free/self.total_trajectories_filtered*100
        else:  # if all trajectories are filtered the type % are 0 
            ratio_immobile = 0
            ratio_confined = 0
            ratio_free = 0
        return ratio_immobile, ratio_confined, ratio_free
    
    def filter_cell_trc(self):
        self.filtered_trc_files = []
        for cell_index in range(len(self.cell_trajectories)):
            trc_file = copy.deepcopy(self.cells[cell_index].trc_file)
            bool_arr = np.array([self.filter_cell_trc_function(row, cell_index) for row in trc_file])
            self.filtered_trc_files.append(trc_file[bool_arr])

    def filter_cell_trc_function(self, row, cell_index):
        """
        If first trc index (numbering) in filtered trajectory index -> true.
        """
        return row[0] in self.cell_trajectories_filtered_index[cell_index]
    
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
        """
        for cell in self.cell_trajectories:
            for trajectory in cell:
                if trajectory.D < self.min_D:
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
                log_Ds[trajectory_index] =  np.log10(self.cell_trajectories_filtered[cell_index][trajectory_index].D)
            self.calc_diffusion_frequencies(log_Ds, desired_bin_size, cell_size) 
        
    def calc_diffusion_frequencies(self, log_diff, desired_bin, size, is_cell=True):
        """
        :param log_diff: np array with log10(D) of one cell.
        :param size: cell size.
        :param desired_bin: bin size.
        """
        # min & max determined by diffusions_log_complete function
        min_bin = np.ceil(np.log10(self.min_D)/desired_bin)*desired_bin
        max_bin = np.ceil(np.log10(self.max_D)/desired_bin)*desired_bin 
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
            print("freq, no %", self.mean_frequencies)
            print("SUMME", np.sum(self.mean_frequencies))
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
        
    def plot_bar_log_bins_bg_corrected(self):
        plt.subplot(111, xscale="log")
        (_, caps, _) = plt.errorbar(self.hist_diffusion, self.corrected_frequencies_percent, yerr=self.corrected_frequencies_error_percent, capsize=4, label="relative frequency")  # capsize length of cap
        for cap in caps:
            cap.set_markeredgewidth(1)  # markeredgewidth thickness of cap (vertically)
        plt.xlim(self.min_D, self.max_D)
        plt.legend()
        plt.title("Distribution of diffusion coefficients (background corrected)")
        plt.ylabel("normalized relative occurence [%]")
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
        



        
 
    
def main():
    pass
        
    
if __name__ == "__main__":
    main()
        