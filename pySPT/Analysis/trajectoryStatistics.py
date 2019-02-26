# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 11:45:55 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

from . import trajectory
from . import cell
import numpy as np
import pylab as pl
import copy
import math
import matplotlib.pyplot as plt

class TrajectoryStatistics():
    def __init__(self):
        self.cell_trajectories = [] # [[],[]] contains list of cells, cells contain trajectories
        self.cell_trajectories_filtered = []  # deep copy of original cell trajectories
        self.cell_trajectories_index = []
        self.cell_trajectories_filtered_index = []  # deep copy of original cell trajectories index
        self.background_trajectories = [] # [[],[]] contains list of cells, cells contain trajectories
        self.background_trajectories_filtered = []  # deep copy of original cell trajectories
        self.background_trajectories_index = []
        self.background_trajectories_filtered_index = []  # deep copy of original cell trajectories index
        self.min_D = math.inf
        self.max_D = - math.inf
        self.total_trajectories = 0  # amount of trajectories in data set
        self.cell_sizes = []
        self.hist_log_Ds = []  # histograms (logD vs freq) from all cells as np arrays in this list
        self.diffusion_frequencies = []  # only freq (divided by cell size) of cells
        self.hist_diffusion = []  # diffusions from histogram calculation, transformed back -> 10^-(log10(D))
        self.mean_frequencies = []  # mean frequencies, size corrected
        self.mean_error = []  # standard error of mean value
        self.normalization_factor = 0.0  # 100/sum of all mean frequencies
                   
    def run_statistics(self, min_length, max_length, min_D, max_D, filter_immob, filter_confined,
                       filter_free, filter_analyse_successful, filter_analyse_not_successful):
        """
        Initialize filtered lists, filter functions and percentage function.
        """
        self.get_index()
        self.create_init_filter_lst()
        try:
            max_length = int(max_length)
            print("max trajectory length:", max_length)
        except ValueError:
            max_length = self.get_max_length()
            print("max trajectory length:", max_length)
        try:
            min_length = int(min_length)
            print("min trajectory length:", min_length)
        except ValueError:
            min_length = self.get_min_length()
            print("min trajectory length:", min_length)
        try:
            max_D = float(max_D)
            print("max diffusion coefficient: {} [\u03BCm\u00b2/s]".format(max_D))
        except ValueError:
            max_D = self.get_max_D()
            print("max diffusion coefficient: {} [\u03BCm\u00b2/s]".format(max_D))
        try:
            min_D = float(min_D)
            print("min diffusion coefficient: {} [\u03BCm\u00b2/s]".format(min_D))
        except ValueError:
            min_D = self.get_min_D()
            print("min diffusion coefficient: {} [\u03BCm\u00b2/s]".format(min_D))
        self.filter_length(min_length, max_length)
        self.filter_D(min_D, max_D)
        self.filter_type(filter_immob, filter_confined, filter_free,
                         filter_analyse_successful, filter_analyse_not_successful)
        print("%.1f %% are immobile" %(self.type_percentage()[0]))
        print("%.1f %% are confined" %(self.type_percentage()[1]))
        print("%.1f %% are free" %(self.type_percentage()[2]))     
        if self.type_percentage()[0] + self.type_percentage()[1] + self.type_percentage()[2] == 0:
            print("The selection excludes all data.")
        print("Total trajectories:", self.total_trajectories)
        
    def run_statistics_no_filter(self, bin_size):
        """
        If stared from JNB trackAnalysis, no filter are applied, slight deviation in function calls therefore.
        """
        self.get_index()
        self.create_init_filter_lst()
        self.type_percentage_pre()
        self.run_plot_diffusion_histogram(bin_size)
        
    def create_init_filter_lst(self):
        """
        Create copy of initial cell trajectories & index list.
        """
        self.cell_trajectories_filtered = copy.deepcopy(self.cell_trajectories)
        self.cell_trajectories_filtered_index = copy.deepcopy(self.cell_trajectories_index)

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

    def crop_lst(self, cell, trajectory):
        """
        Delete trajectory and its index.
        """
        trajectory_index = self.cell_trajectories_filtered[cell].index(trajectory)
        self.cell_trajectories_filtered[cell] = self.cell_trajectories_filtered[cell][0:trajectory_index] + \
        self.cell_trajectories_filtered[cell][trajectory_index+1:]
        self.cell_trajectories_filtered_index[cell] = self.cell_trajectories_filtered_index[cell][0:trajectory_index] + \
        self.cell_trajectories_filtered_index[cell][trajectory_index+1:] 
        
    def filter_type(self, filter_immob, filter_confined, filter_free,
                    filter_analyse_successful, filter_analyse_not_successful):
        """
        :param filter_immob: if checked in JNB -> True -> filter data includes immob trajectories, else False -> will be neglected.
        """
        # if filter was not selected, value is False
        if not filter_immob:
            for cell in range(0, len(self.cell_trajectories_filtered)):
                for trajectory in self.cell_trajectories_filtered[cell]:
                    # meaning that one wants to get rid of all immobile particles -> crop them out of the lists.
                    if trajectory.immobility and not trajectory.confined:
                        self.crop_lst(cell, trajectory)
        if not filter_confined:
            for cell in range(0, len(self.cell_trajectories_filtered)):
                for trajectory in self.cell_trajectories_filtered[cell]:
                    if trajectory.confined and not trajectory.immobility:
                        self.crop_lst(cell, trajectory)
        if not filter_free:
            for cell in range(0, len(self.cell_trajectories_filtered)):
                for trajectory in self.cell_trajectories_filtered[cell]:
                    if not trajectory.confined and not trajectory.immobility:
                        self.crop_lst(cell, trajectory)
        if not filter_analyse_successful:
            for cell in range(0, len(self.cell_trajectories_filtered)):
                for trajectory in self.cell_trajectories_filtered[cell]:
                    if trajectory.analyse_successful:
                        self.crop_lst(cell, trajectory)
        if not filter_analyse_not_successful:
            for cell in range(0, len(self.cell_trajectories_filtered)):
                for trajectory in self.cell_trajectories_filtered[cell]:
                    if not trajectory.analyse_successful:
                        self.crop_lst(cell, trajectory)
    
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
        
    def type_percentage(self):
        """
        (immob/confined -> true/false=immob, false/true=conf, false/false=free)
        Calculate percentage of immobile free and confined based on total number of trajectories in all cells.
        If no trajectory exists (total_trajectories = 0) percentages will be set to zero, no calculation will be made.
        """
        data_selected = True
        self.total_trajectories = 0
        for cell_index in range(0, len(self.cell_trajectories_filtered)):
            self.total_trajectories += len(self.cell_trajectories_filtered[cell_index])
        if self.total_trajectories == 0:
            data_selected = False
        if data_selected:
            count_immobile = 0
            count_confined = 0
            count_free = 0
            for cell in self.cell_trajectories_filtered:
                for trajectory in cell:
                    if trajectory.immobility and not trajectory.confined:
                        count_immobile += 1
                    if trajectory.confined and not trajectory.immobility:
                        count_confined += 1
                    # has to be not confined AND not immobile (otherwise it will count the immobile particles as well)
                    if not trajectory.confined and not trajectory.immobility:
                        count_free +=1
            ratio_immobile = count_immobile/self.total_trajectories*100
            ratio_confined = count_confined/self.total_trajectories*100
            ratio_free = count_free/self.total_trajectories*100
        else:
            ratio_immobile = 0
            ratio_confined = 0
            ratio_free = 0
        return ratio_immobile, ratio_confined, ratio_free
    
    def type_percentage_pre(self):
        """
        Calculation before saving as hdf5 (immob/confined -> true/true=immob, false/true=conf, false/false=free)
        Calculate percentage of immobile free and confined based on total number of trajectories in all cells.
        If no trajectory exists (total_trajectories = 0) percentages will be set to zero, no calculation will be made.
        """
        data_selected = True
        self.total_trajectories = 0
        for cell_index in range(0, len(self.cell_trajectories_filtered)):
            self.total_trajectories += len(self.cell_trajectories_filtered[cell_index])
        if self.total_trajectories == 0:
            data_selected = False
        if data_selected:
            count_immobile = 0
            count_confined = 0
            count_free = 0
            for cell in self.cell_trajectories_filtered:
                for trajectory in cell:
                    if trajectory.immobility and trajectory.confined:
                        count_immobile += 1
                    if trajectory.confined and not trajectory.immobility:
                        count_confined += 1
                    # has to be not confined AND not immobile (otherwise it will count the immobile particles as well)
                    if not trajectory.confined and not trajectory.immobility:
                        count_free +=1
            ratio_immobile = count_immobile/self.total_trajectories*100
            ratio_confined = count_confined/self.total_trajectories*100
            ratio_free = count_free/self.total_trajectories*100
        else:
            ratio_immobile = 0
            ratio_confined = 0
            ratio_free = 0
        print("%.1f %% are immobile" %(ratio_immobile))
        print("%.1f %% are confined" %(ratio_confined))
        print("%.1f %% are free" %(ratio_free)) 
        print("Total trajectories:", self.total_trajectories)
    
    # plot diffusion vs frequencies.
    
    def run_plot_diffusion_histogram(self, desired_bin_size):
        self.clear_attributes()
        self.determine_max_min_diffusion()
        self.diffusions_log(float(desired_bin_size))
        self.calc_nonlogarithmic_diffusions()
        self.determine_mean_frequency()
        self.calc_mean_error()
        self.plot_bar_log_bins()
        
    def clear_attributes(self):
        """
        If one restarts filtering, these attributes are empy (otherwise they would append).
        """
        self.hist_log_Ds = []  # histograms (logD vs freq) from all cells as np arrays in this list
        self.hist_diffusion = []  # diffusions from histogram calculation, transformed back -> 10^-(log10(D))
        self.mean_frequencies = []  # mean frequencies, size corrected
        self.mean_error = []
        
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
        
    def calc_diffusion_frequencies(self, log_diff, desired_bin, size):
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
        log_diffusion_hist[:,1] = log_diffusion_hist[:,1]  /size
        self.hist_log_Ds.append(log_diffusion_hist)
        
    def calc_nonlogarithmic_diffusions(self):
        """
        Calculate the nonlogarithmic diffusion coefficients from log10(D) from histogram.
        """
        self.hist_diffusion = 10**self.hist_log_Ds[0][:,0]
        
    def determine_mean_frequency(self):
        """
        Mean frequency will be calculated based on all frequencies of cells.
        Normalize an array (sum of elements = 1) and represent is in percent (*100).
        """
        self.diffusion_frequencies = self.create_np_array(np.shape(self.hist_log_Ds)[1], len(self.cell_trajectories_filtered))
        for i in range (0, len(self.cell_trajectories_filtered)):
            self.diffusion_frequencies[:,i] = self.hist_log_Ds[i][:,1]
        self.mean_frequencies = self.calc_mean_frequencies(self.diffusion_frequencies)
        
        self.normalization_factor = 100/np.sum(self.mean_frequencies)
        self.mean_frequencies = self.mean_frequencies * self.normalization_factor
        
    def calc_mean_error(self):
        """
        Standard deviation (N-1) divided by square root of number of elements.
        Normalize an array (sum of elements = 1) and represent is in percent (*100).
        """
        self.mean_error =  np.std(self.diffusion_frequencies, axis=1, ddof=1)/(np.shape(self.diffusion_frequencies)[1])**(1/2) 
        self.mean_error = self.mean_error * self.normalization_factor 
        
    def plot_bar_log_bins(self):
        plt.subplot(111, xscale="log")
        (_, caps, _) = plt.errorbar(self.hist_diffusion, self.mean_frequencies, yerr=self.mean_error, capsize=4, label="relative frequency")  # capsize length of cap
        for cap in caps:
            cap.set_markeredgewidth(1)  # markeredgewidth thickness of cap (vertically)
        plt.xlim(self.min_D, self.max_D)
        plt.legend()
        plt.title("Distribution of diffusion coefficients")
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
        
    def normalize_hist(self, normalized_col):
        """
        Normalize a column and return it.
        """
        normalized_col = normalized_col / np.sum(normalized_col)
        return normalized_col
       
        
 
    
def main():
    pass
        
    
if __name__ == "__main__":
    main()
        