# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 09:17:05 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

For JNB "trackAnalysis": ...
"""


import numpy as np
import copy
import math
import matplotlib.pyplot as plt


class TrackAnalysis():
    def __init__(self):
        self.cell_trajectories = [] # [[],[]] contains list of cells, cells contain trajectories
        self.cell_trajectories_filtered = []  # deep copy of original cell trajectories
        self.cell_trajectories_index = []
        self.cell_trajectories_filtered_index = []  # deep copy of original cell trajectories index
        self.min_D = math.inf
        self.max_D = - math.inf
        self.total_trajectories = 0  # amount of trajectories in data set
        self.total_trajectories_cell = []  # amount of trajectories per cell
        self.cell_type_count = []  # tupel with percentage of types (immob, confined, free %) per cell
        self.cell_sizes = []
        self.hist_log_Ds = []  # histograms (logD vs freq) from all cells as np arrays in this list
        self.diffusion_frequencies = []  # only freq (divided by cell size) of cells
        self.hist_diffusion = []  # diffusions from histogram calculation, transformed back -> 10^-(log10(D))
        self.mean_frequencies = []  # mean frequencies, size corrected
        self.mean_error = []  # standard error of mean value
        self.normalization_factor = 0.0  # 100/sum of all mean frequencies
        # Save
        self.diffusion_info = []
        self.number_of_trajectories = 0
        self.rossier_info = []
        self.diff_plot = []
        self.rossier_plot = []
        self.diff_fig = []  # log diffusion plot

    def run_statistics_no_filter(self):
        """
        If stared from JNB trackAnalysis, no filter are applied, slight deviation in function calls therefore.
        """
        self.get_index()
        self.create_init_filter_lst()
        self.type_percentage_pre()
        
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
    
    def type_percentage_pre(self):
        """
        Calculation before saving as hdf5 (immob/confined -> true/true=immob, false/true=conf, false/false=free)
        Calculate percentage of immobile free and confined based on total number of trajectories in all cells.
        If no trajectory exists (total_trajectories = 0) percentages will be set to zero, no calculation will be made.
        """
        data_selected = True
        self.total_trajectories = 0
        self.total_trajectories_cell = []
        self.cell_type_count = []
        for cell_index in range(0, len(self.cell_trajectories_filtered)):
            self.total_trajectories_cell.append(len(self.cell_trajectories_filtered[cell_index]))
        self.total_trajectories = np.sum(self.total_trajectories_cell)
        if self.total_trajectories == 0:
            data_selected = False
        if data_selected:
            count_immobile = 0
            count_confined = 0
            count_free = 0
            count_not_successful = 0
            for cell in self.cell_trajectories_filtered:
                count_immobile_cell = 0
                count_confined_cell = 0
                count_free_cell = 0
                count_not_successful_cell = 0
                for trajectory in cell:
                    if trajectory.immobility and trajectory.confined and trajectory.analyse_successful:
                        count_immobile_cell += 1
                        count_immobile += 1
                    if trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                        count_confined_cell += 1
                        count_confined += 1
                    # has to be not confined AND not immobile (otherwise it will count the immobile particles as well)
                    if not trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                        count_free_cell += 1
                        count_free +=1
                    if not trajectory.analyse_successful:
                        count_not_successful_cell += 1
                        count_not_successful += 1
                cell_index = self.cell_trajectories_filtered.index(cell)
                ratio_immobile_cell = count_immobile_cell/self.total_trajectories_cell[cell_index]*100
                ratio_confined_cell = count_confined_cell/self.total_trajectories_cell[cell_index]*100
                ratio_free_cell = count_free_cell/self.total_trajectories_cell[cell_index]*100
                ratio_not_successful_cell = count_not_successful_cell/self.total_trajectories_cell[cell_index]*100
                cell_types_percent = (ratio_immobile_cell, ratio_confined_cell, ratio_free_cell, ratio_not_successful_cell)
                self.cell_type_count.append(cell_types_percent)
            ratio_immobile = count_immobile/self.total_trajectories*100
            ratio_confined = count_confined/self.total_trajectories*100
            ratio_free = count_free/self.total_trajectories*100
            ratio_not_successful = count_not_successful/self.total_trajectories*100
        else:
            ratio_immobile = 0
            ratio_confined = 0
            ratio_free = 0
            ratio_not_successful = 0
        print("%.2f %% are immobile" %(ratio_immobile))
        print("%.2f %% are confined" %(ratio_confined))
        print("%.2f %% are free" %(ratio_free)) 
        print("%.2f %% could not be analysed" %(ratio_not_successful)) 
        print("Total trajectories:", self.total_trajectories)

    def run_plot_diffusion_histogram(self, desired_bin_size):
        # bin size can only be something that can be converted to float (integer or float, comma separated)
        try:
            float(desired_bin_size)
        except:
            print("Insert a dot separated float or integer as bin size (e.g. 0.1).")
        # bin size can not be 0
        else:
            if float(desired_bin_size) != 0.0:
                self.clear_attributes()
                self.determine_max_min_diffusion()
                self.diffusions_log(float(desired_bin_size))
                self.calc_nonlogarithmic_diffusions()
                self.determine_mean_frequency()
                self.calc_mean_error()
                self.plot_bar_log_bins()
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
        self.diff_fig = plt.figure()
        plt.subplot(111, xscale="log")
        (_, caps, _) = plt.errorbar(self.hist_diffusion, self.mean_frequencies, yerr=self.mean_error, capsize=4, label="relative frequency")  # capsize length of cap
        for cap in caps:
            cap.set_markeredgewidth(1)  # markeredgewidth thickness of cap (vertically)
        plt.xlim(self.min_D, self.max_D)
        plt.legend()
        plt.title("Distribution of diffusion coefficients")
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
    
    # Save
    
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
        self.rossier_info = np.zeros([np.size(trajectories), 12])
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
                self.rossier_info[i,11] = 0
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
                self.rossier_info[i,10] = trajectory.dD_conf
                self.rossier_info[i,11] = trajectory.chi_MSD_fit
            # if analysis was not successful -> output is 0 by default
            elif not trajectory.analyse_successful and not trajectory.immobility:
                self.rossier_info[i,1] = trajectory.immobility
                self.rossier_info[i,2] = False
                self.rossier_info[i,3] = False
                self.rossier_info[i,4] = trajectory.analyse_successful
                self.rossier_info[i,5] = 0
                self.rossier_info[i,6] = 0
                self.rossier_info[i,7] = 0
                self.rossier_info[i,8] = 0
                self.rossier_info[i,9] = 0
                self.rossier_info[i,10] = 0
                self.rossier_info[i,11] = 0
            i += 1
 
    
def main():
    pass
    

if __name__ == "__main__":
    main()
    