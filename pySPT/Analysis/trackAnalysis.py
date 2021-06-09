"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Load tracked files, calculate diffusion coefficients, type analysis, save as *.h5.
"""

import numpy as np
import pandas as pd
import copy
import math
import matplotlib.pyplot as plt
import seaborn as sns


class TrackAnalysis():
    def __init__(self):
        self.seg_id = ""  # "seg id" or "track id"
        self.cell_trajectories = []  # [[],[]] contains list of cells, cells contain trajectories
        self.cell_trajectories_filtered = []  # deep copy of original cell trajectories
        self.cell_trajectories_index = []
        self.cell_trajectories_filtered_index = []  # deep copy of original cell trajectories index
        self.min_D = math.inf
        self.max_D = - math.inf
        self.total_trajectories = 0  # amount of trajectories in data set
        self.total_trajectories_cell = []  # amount of trajectories per cell
        self.cell_type_count = []  # tupel with percentage of types (immob, confined, free, notype %) per cell
        self.cell_sizes = []
        self.hist_log_Ds = []  # histograms (logD vs freq) from all cells as np arrays in this list
        self.diffusion_frequencies = []  # only freq (divided by cell size) of cells
        self.hist_diffusion = []  # diffusions from histogram calculation, transformed back -> 10^-(log10(D))
        self.mean_frequencies = []  # mean frequencies, size corrected
        self.mean_error = []  # standard error of mean value
        self.normalization_factor = 0.0  # 100/sum of all mean frequencies
        self.mean_D_cells = []  # average diffusion coefficient per type and cell
        self.mean_dD_cells = []  # SEM of average diffusion coefficient per type and cell
        self.mean_D = []  # average diffusion coefficient over all cells
        self.mean_length_cells = []  # average trajectory lengths per type and cells
        self.mean_dlength_cells = []  # SEM of average trajectory lengths per type and cell
        self.mean_length = []  # average trajectory lengths per type over all cells
        self.type_ratios = []  # average trajectory type in % over all cells
        # Save
        self.diffusion_info = []
        self.number_of_trajectories = 0
        self.rossier_info = []
        self.diff_plot = []
        self.rossier_plot = []
        self.diff_fig = []  # log diffusion plot
        self.diff_fig_types = []  # log diffusion plot types
        self.MSD_fig_types = []  # MSD mean values plot types
        self.trajectories_immob_cells = []
        self.trajectories_conf_cells = []
        self.trajectories_free_cells = []
        self.trajectories_notype_cells = []

    def run_statistics_no_filter(self):
        """
        If stared from JNB trackAnalysis, no filter are applied, slight deviation in function calls therefore.
        """
        self.get_index()
        self.create_init_filter_lst()
        self.type_percentage_pre()
        self.calc_mean_D_cells()
        self.calc_mean_length_cells()
        self.print_stats()
        
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
        Calculate percentage of immobile free and confined per cell and average over all cells.
        If no trajectory exists (in a type, total trajectories) percentages will be considered as zero.
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
            count_immobile, count_confined, count_free, count_not_successful = 0, 0, 0, 0
            for cell_index, cell in enumerate(self.cell_trajectories_filtered):
                count_immobile_cell, count_confined_cell, count_free_cell, count_not_successful_cell = 0, 0, 0, 0
                trajectories_immob, trajectories_conf, trajectories_free, trajectories_notype = [], [], [], []
                for trajectory in cell:
                    if trajectory.immobility and trajectory.confined and trajectory.analyse_successful:
                        count_immobile_cell += 1
                        count_immobile += 1
                        trajectories_immob.append(trajectory)
                    if trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                        count_confined_cell += 1
                        count_confined += 1
                        trajectories_conf.append(trajectory)
                    # has to be not confined AND not immobile (otherwise it will count the immobile particles as well)
                    if not trajectory.confined and not trajectory.immobility and trajectory.analyse_successful:
                        count_free_cell += 1
                        count_free += 1
                        trajectories_free.append(trajectory)
                    if not trajectory.analyse_successful:
                        count_not_successful_cell += 1
                        count_not_successful += 1
                        trajectories_notype.append(trajectory)
                self.trajectories_immob_cells.append(trajectories_immob)
                self.trajectories_conf_cells.append(trajectories_conf)
                self.trajectories_free_cells.append(trajectories_free)
                self.trajectories_notype_cells.append(trajectories_notype)
                ratio_immobile_cell = count_immobile_cell/self.total_trajectories_cell[cell_index]*100
                ratio_confined_cell = count_confined_cell/self.total_trajectories_cell[cell_index]*100
                ratio_free_cell = count_free_cell/self.total_trajectories_cell[cell_index]*100
                ratio_not_successful_cell = count_not_successful_cell/self.total_trajectories_cell[cell_index]*100
                cell_types_percent = (ratio_immobile_cell, ratio_confined_cell, ratio_free_cell, ratio_not_successful_cell)
                self.cell_type_count.append(cell_types_percent)
        self.type_ratios = np.nanmean(self.cell_type_count, axis=0)

    def calc_mean_D_cells(self):
        """
        Calculate the mean diffusion coefficient for four types per cell & the overall average of D.
        """
        for cell in self.cell_trajectories_filtered:
            immob = [trajectory.D for trajectory in cell if
                     trajectory.immobility and trajectory.confined and trajectory.analyse_successful]
            mean_D_immob = np.mean(immob)
            mean_dD_immob = np.std(immob, ddof=1) / math.sqrt(len(immob))
            conf = [trajectory.D for trajectory in cell if
                    not trajectory.immobility and trajectory.confined and trajectory.analyse_successful]
            mean_D_conf = np.mean(conf)
            mean_dD_conf = np.std(conf, ddof=1) / math.sqrt(len(conf))
            free = [trajectory.D for trajectory in cell if
                    not trajectory.immobility and not trajectory.confined and trajectory.analyse_successful]
            mean_D_free = np.mean(free)
            mean_dD_free = np.std(free, ddof=1) / math.sqrt(len(free))
            notype = [trajectory.D for trajectory in cell if not trajectory.analyse_successful]
            mean_D_notype = np.mean(notype)
            mean_dD_notype = np.std(notype, ddof=1) / math.sqrt(len(notype))
            mean_cell = (mean_D_immob, mean_D_conf, mean_D_free, mean_D_notype)
            mean_dD_cell = (mean_dD_immob, mean_dD_conf, mean_dD_free, mean_dD_notype)
            self.mean_D_cells.append(mean_cell)
            self.mean_dD_cells.append(mean_dD_cell)
        self.mean_D = np.nanmean(self.mean_D_cells, axis=0)

    def calc_mean_length_cells(self):
        """
        Calculate the mean trajectory length for four types per cell & the overall average length.
        """
        for cell in self.cell_trajectories_filtered:
            immob = [trajectory.length_trajectory for trajectory in cell if
                     trajectory.immobility and trajectory.confined and trajectory.analyse_successful]
            mean_length_immob = np.mean(immob)
            mean_dlength_immob = np.std(immob, ddof=1) / math.sqrt(len(immob))
            conf = [trajectory.length_trajectory for trajectory in cell if
                    not trajectory.immobility and trajectory.confined and trajectory.analyse_successful]
            mean_length_conf = np.mean(conf)
            mean_dlength_conf = np.std(conf, ddof=1) / math.sqrt(len(conf))
            free = [trajectory.length_trajectory for trajectory in cell if
                    not trajectory.immobility and not trajectory.confined and trajectory.analyse_successful]
            mean_length_free = np.mean(free)
            mean_dlength_free = np.std(free, ddof=1) / math.sqrt(len(free))
            notype = [trajectory.length_trajectory for trajectory in cell if not trajectory.analyse_successful]
            mean_length_notype = np.mean(notype)
            mean_dlength_notype = np.std(notype, ddof=1) / math.sqrt(len(notype))
            mean_cell = (mean_length_immob, mean_length_conf, mean_length_free, mean_length_notype)
            mean_dlength_cell = (mean_dlength_immob, mean_dlength_conf, mean_dlength_free, mean_dlength_notype)
            self.mean_length_cells.append(mean_cell)
            self.mean_dlength_cells.append(mean_dlength_cell)
        self.mean_length = np.nanmean(self.mean_length_cells, axis=0)

    def print_stats(self):
        print("%.2f %% are immobile, mean D = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(self.type_ratios[0], self.mean_D[0], self.mean_length[0]))
        print("%.2f %% are confined, mean D = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(self.type_ratios[1], self.mean_D[1], self.mean_length[1]))
        print("%.2f %% are free, mean D = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(self.type_ratios[2], self.mean_D[2], self.mean_length[2]))
        print("%.2f %% could not be analysed, mean D = %.5f \u03BCm\u00b2/s, mean length = %.0f frames" %(self.type_ratios[3], self.mean_D[3], self.mean_length[3]))
        print("Total trajectories:", self.total_trajectories)

    def run_plot_diffusion_histogram(self, desired_bin_size, MSD_delta_t_n, y_lim):
        try:
            float(desired_bin_size)
        except:
            print("Insert a dot separated float or integer as bin size (e.g. 0.1).")
        else:
            if float(desired_bin_size) != 0.0:
                # Histogram of log D
                self.clear_attributes()
                self.determine_max_min_diffusion()
                log_hist_immob, log_hist_conf, log_hist_free, log_hist_fit_fail = self.diffusions_log(
                    float(desired_bin_size))
                self.calc_nonlogarithmic_diffusions()
                self.determine_mean_frequency()
                self.calc_mean_error()
                self.plot_bar_log_bins()
                # Histogram of log D per type
                self.diffusion_hist_types(float(desired_bin_size))
                # MSD plot per type
                MSD_delta_t_n = None if MSD_delta_t_n == "None" else float(MSD_delta_t_n)
                y_lim = None if y_lim == "None" else float(y_lim)
                self.MSD_types(MSD_delta_t_n, y_lim)
                self.plot_sizes_ntracks()
                self.plot_type_distributions()
                self.plot_D_distributions()
                self.plot_length_distributions()
            else:
                print("Bin size can not be zero.")

    def plot_diffusion_hist_types(self, mean_log_hist_immob, error_log_hist_immob, mean_log_hist_conf,
                                  error_log_hist_conf, mean_log_hist_free, error_log_hist_free,
                                  mean_log_hist_fit_fail, error_log_hist_fit_fail):
        self.diff_fig_types = plt.figure()
        plt.subplot(111, xscale="log")
        (_, caps, _) = plt.errorbar(self.hist_diffusion, mean_log_hist_immob, yerr=error_log_hist_immob, capsize=4,
                                    label="relative frequency immobile", ecolor="#4169e1", color="#4169e1")
        for cap in caps:
            cap.set_markeredgewidth(1)

        (_, caps, _) = plt.errorbar(self.hist_diffusion, mean_log_hist_fit_fail, yerr=error_log_hist_fit_fail, capsize=4,
                                    label="relative frequency no type", ecolor="#8b008b", color="#8b008b")
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

        for c, i in enumerate(self.mean_frequencies):
            if i != 0 and c == 0:
                xlim_min = self.hist_diffusion[c]
                break
            elif i != 0:
                xlim_min = self.hist_diffusion[c-1]
                break
        for c, i in enumerate(np.flip(self.mean_frequencies)):
            if i != 0 and c == 0:
                xlim_max = self.hist_diffusion[len(self.hist_diffusion)-c]
                break
            elif i != 0:
                xlim_max = self.hist_diffusion[len(self.hist_diffusion)-c+1]
                break
        plt.xlim(xlim_min, xlim_max)
        plt.legend()
        plt.title("Distribution of diffusion coefficients per type")
        plt.ylabel("Normalized relative occurence [%]")
        plt.xlabel("D [\u03BCm\u00b2/s]")
        plt.show()

    def diffusion_hist_types(self, desired_bin_size):
        log_Ds_immob = []
        for cell in self.trajectories_immob_cells:
            log_Ds_immob_cell = [np.log10(trajectory.D) for trajectory in cell if trajectory.D > 0]
            log_Ds_immob.append(log_Ds_immob_cell)

        log_Ds_conf = []
        for cell in self.trajectories_conf_cells:
            log_Ds_conf_cell = [np.log10(trajectory.D) for trajectory in cell if trajectory.D > 0]
            log_Ds_conf.append(log_Ds_conf_cell)

        log_Ds_free = []
        for cell in self.trajectories_free_cells:
            log_Ds_free_cell = [np.log10(trajectory.D) for trajectory in cell if trajectory.D > 0]
            log_Ds_free.append(log_Ds_free_cell)

        log_Ds_notype = []
        for cell in self.trajectories_notype_cells:
            log_Ds_notype_cell = [np.log10(trajectory.D) for trajectory in cell if trajectory.D > 0]
            log_Ds_notype.append(log_Ds_notype_cell)

        hist_immob, hist_conf, hist_free, hist_notype = [], [], [], []
        for i, cell_size in enumerate(self.cell_sizes):
            hist_immob_cell = self.calc_diffusion_frequencies(log_Ds_immob[i], desired_bin_size, cell_size)
            hist_conf_cell = self.calc_diffusion_frequencies(log_Ds_conf[i], desired_bin_size, cell_size)
            hist_free_cell = self.calc_diffusion_frequencies(log_Ds_free[i], desired_bin_size, cell_size)
            hist_notype_cell = self.calc_diffusion_frequencies(log_Ds_notype[i], desired_bin_size, cell_size)
            hist_immob.append(hist_immob_cell[:, 1])  # only frequency counts
            hist_conf.append(hist_conf_cell[:, 1])
            hist_free.append(hist_free_cell[:, 1])
            hist_notype.append(hist_notype_cell[:, 1])

        mean_log_hist_immob = np.nanmean(hist_immob, axis=0) * self.normalization_factor
        error_log_hist_immob = np.nanstd(hist_immob, axis=0, ddof=1) / math.sqrt(len(hist_immob)) * self.normalization_factor
        mean_log_hist_conf = np.nanmean(hist_conf, axis=0) * self.normalization_factor
        error_log_hist_conf = np.nanstd(hist_conf, axis=0, ddof=1) / math.sqrt(len(hist_conf)) * self.normalization_factor
        mean_log_hist_free = np.nanmean(hist_free, axis=0) * self.normalization_factor
        error_log_hist_free = np.nanstd(hist_free, axis=0, ddof=1) / math.sqrt(len(hist_free)) * self.normalization_factor
        mean_log_hist_notype = np.nanmean(hist_notype, axis=0) * self.normalization_factor
        error_log_hist_notype = np.nanstd(hist_notype, axis=0, ddof=1) / math.sqrt(len(hist_notype)) * self.normalization_factor

        self.plot_diffusion_hist_types(mean_log_hist_immob, error_log_hist_immob, mean_log_hist_conf,
                                       error_log_hist_conf, mean_log_hist_free, error_log_hist_free,
                                       mean_log_hist_notype, error_log_hist_notype)

    def plot_MSD_types(self, mean_MSD_immob, mean_error_MSD_immob, mean_MSD_conf, mean_error_MSD_conf,
                       mean_MSD_free, mean_error_MSD_free, mean_MSD_fit_fail, mean_error_MSD_fit_fail,
                       MSD_delta_t_n, y_lim):
        """
        Plot the mean MSD curve per diffusion type
        """
        camera_time = self.cell_trajectories[0][0].dt  # camera integration time in s

        delta_t_immob = [camera_time * i for i in range(1, len(mean_MSD_immob)+1)]
        delta_t_conf = [camera_time * i for i in range(1, len(mean_MSD_conf)+1)]
        delta_t_free = [camera_time * i for i in range(1, len(mean_MSD_free)+1)]
        delta_t_fit_fail = [camera_time * i for i in range(1, len(mean_MSD_fit_fail)+1)]

        self.MSD_fig_types = plt.figure()
        (_, caps, _) = plt.errorbar(delta_t_immob, mean_MSD_immob, yerr=mean_error_MSD_immob, capsize=4,
                                    label="MSD immobile", ecolor="#4169e1", color="#4169e1")
        for cap in caps:
            cap.set_markeredgewidth(1)

        (_, caps, _) = plt.errorbar(delta_t_fit_fail, mean_MSD_fit_fail, yerr=mean_error_MSD_fit_fail, capsize=4,
                                    label="MSD no type", ecolor="#8b008b", color="#8b008b")
        for cap in caps:
            cap.set_markeredgewidth(1)

        (_, caps, _) = plt.errorbar(delta_t_conf, mean_MSD_conf, yerr=mean_error_MSD_conf, capsize=4,
                                    label="MSD confined", ecolor="#228b22", color="#228b22")
        for cap in caps:
            cap.set_markeredgewidth(1)

        (_, caps, _) = plt.errorbar(delta_t_free, mean_MSD_free, yerr=mean_error_MSD_free, capsize=4,
                                    label="MSD free", ecolor="#ff8c00", color="#ff8c00")
        for cap in caps:
            cap.set_markeredgewidth(1)

        x_lim_max = None if MSD_delta_t_n is None else MSD_delta_t_n
        plt.xlim(0, x_lim_max)
        plt.ylim(0, y_lim)
        plt.legend()
        plt.title("Mean MSD values per type")
        plt.ylabel("Mean MSD [\u03BCm\u00b2]")
        plt.xlabel("Time step [s]")
        plt.show()

    def plot_violine(self, lst, title, x_axis, y_axis, color="cornflowerblue"):
        df = pd.DataFrame(lst)
        fig = plt.figure(figsize=(3, 5))
        ax = sns.violinplot(data=df, color=color, showmeans=True,
                         meanprops={"marker": "s", "markerfacecolor": "white", "markeredgecolor": "0.25"})
        ax = sns.swarmplot(data=df, color="0.25")
        ax.set_ylabel(y_axis)
        ax.set_xlabel(x_axis)
        ax.set_title(title)
        ax.set(xticklabels=[])
        plt.show()


    def plot_violine_subplots(self, lsts, title, x_axiss, y_axis, colors):
        fig, axs = plt.subplots(1, len(lsts), sharey=True)
        for c, (lst, x_axis, color) in enumerate(zip(lsts, x_axiss, colors)):
            df = pd.DataFrame(lst)
            sns.violinplot(data=df, color=color, showmeans=True, meanprops={"marker": "s", "markerfacecolor": "white", "markeredgecolor": "0.25"}, ax=axs[c])
            sns.swarmplot(data=df, color="0.25", ax=axs[c])
            axs[c].set_xlabel(x_axis)
            axs[c].set(xticklabels=[])
        fig.text(0.04, 0.5, y_axis, va="center", rotation="vertical")
        plt.suptitle(title)
        plt.show()

    def plot_D_distributions(self):
        # average D per cell, global and per type
        all_Ds = []
        for cell in self.cell_trajectories:
             all_Ds.append(np.mean([trajectory.D for trajectory in cell]))
        # get mean D values per cell and type
        immob_Ds = []
        for cell in self.trajectories_immob_cells:
             immob_Ds.append(np.mean([trajectory.D for trajectory in cell]))
        conf_Ds = []
        for cell in self.trajectories_conf_cells:
            conf_Ds.append(np.mean([trajectory.D for trajectory in cell]))
        free_Ds = []
        for cell in self.trajectories_free_cells:
            free_Ds.append(np.mean([trajectory.D for trajectory in cell]))
        notype_Ds = []
        for cell in self.trajectories_notype_cells:
            notype_Ds.append(np.mean([trajectory.D for trajectory in cell]))
        lsts = [all_Ds, immob_Ds, conf_Ds, free_Ds, notype_Ds]
        title = "mean diffusion coefficients per cell"
        x_axiss = ["all", "immobile", "confined", "free", "notype"]
        y_axis = "diffusion coefficient [\u00B5m\u00B2/s]"
        colors = ["cornflowerblue", "#4169e1", "#228b22", "#ff8c00", "#8b008b"]
        self.plot_violine_subplots(lsts, title, x_axiss, y_axis, colors)

    def plot_sizes_ntracks(self):
        # cell sizes
        target = "segments" if self.seg_id == "seg id" else "trajectories"  # SEG OR TRACK
        self.plot_violine(self.cell_sizes, "cell sizes", "cells", "cell size [\u00B5m\u00B2]")
        self.plot_violine([len(tracks) for tracks in self.cell_trajectories], "number of " + target + " per cell", "cells", "number of " + target)

    def plot_length_distributions(self):
        # average lengths per cell, global and per type
        target = "segment" if self.seg_id == "seg id" else "trajectory"  # SEG OR TRACK
        all_lengths = []
        for cell in self.cell_trajectories:
             all_lengths.append(np.mean([trajectory.length_trajectory for trajectory in cell]))
        # get mean D values per cell and type
        immob_lengths = []
        for cell in self.trajectories_immob_cells:
             immob_lengths.append(np.mean([trajectory.length_trajectory for trajectory in cell]))
        conf_lengths = []
        for cell in self.trajectories_conf_cells:
            conf_lengths.append(np.mean([trajectory.length_trajectory for trajectory in cell]))
        free_lengths = []
        for cell in self.trajectories_free_cells:
            free_lengths.append(np.mean([trajectory.length_trajectory for trajectory in cell]))
        notype_lengths = []
        for cell in self.trajectories_notype_cells:
            notype_lengths.append(np.mean([trajectory.length_trajectory for trajectory in cell]))

        lsts = [all_lengths, immob_lengths, conf_lengths, free_lengths, notype_lengths]
        title = "mean " + target + " lengths per cell"
        x_axiss = ["all", "immobile", "confined", "free", "notype"]
        y_axis = "length [frames]"
        colors = ["cornflowerblue", "#4169e1", "#228b22", "#ff8c00", "#8b008b"]
        self.plot_violine_subplots(lsts, title, x_axiss, y_axis, colors)


    def plot_type_distributions(self):
        target = "segments" if self.seg_id == "seg id" else "trajectories"  # SEG OR TRACK
        immob_percentages = [cell[0] for cell in self.cell_type_count]
        conf_percentages = [cell[1] for cell in self.cell_type_count]
        free_percentages = [cell[2] for cell in self.cell_type_count]
        notype_percentages = [cell[3] for cell in self.cell_type_count]

        lsts = [immob_percentages, conf_percentages, free_percentages, notype_percentages]
        title = "mean type percentages per cell"
        x_axiss = ["immobile", "confined", "free", "notype"]
        y_axis = "type [\u0025]"
        colors = ["#4169e1", "#228b22", "#ff8c00", "#8b008b"]
        self.plot_violine_subplots(lsts, title, x_axiss, y_axis, colors)

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

    def MSD_types(self, MSD_delta_t_n, y_lim):
        average_MSDs_immob = []
        average_MSDs_immob_error = []
        for cell in self.trajectories_immob_cells:
            MSDs_cell = [trajectory.MSDs for trajectory in cell]
            mean_MSDs_cell, error_MSDs_cell = self.calc_mean_error_different_lengths(MSDs_cell)
            average_MSDs_immob.append(mean_MSDs_cell)
            average_MSDs_immob_error.append(error_MSDs_cell)
        mean_MSDs_immob, _ = self.calc_mean_error_different_lengths(average_MSDs_immob)
        mean_error_MSDs_immob, _ = self.calc_mean_error_different_lengths(average_MSDs_immob_error)

        average_MSDs_conf = []
        average_MSDs_conf_error = []
        for cell in self.trajectories_conf_cells:
            MSDs_cell = [trajectory.MSDs for trajectory in cell]
            mean_MSDs_cell, error_MSDs_cell = self.calc_mean_error_different_lengths(MSDs_cell)
            average_MSDs_conf.append(mean_MSDs_cell)
            average_MSDs_conf_error.append(error_MSDs_cell)
        mean_MSDs_conf, _ = self.calc_mean_error_different_lengths(average_MSDs_conf)
        mean_error_MSDs_conf, _ = self.calc_mean_error_different_lengths(average_MSDs_conf_error)

        average_MSDs_free = []
        average_MSDs_free_error = []
        for cell in self.trajectories_free_cells:
            MSDs_cell = [trajectory.MSDs for trajectory in cell]
            mean_MSDs_cell, error_MSDs_cell = self.calc_mean_error_different_lengths(MSDs_cell)
            average_MSDs_free.append(mean_MSDs_cell)
            average_MSDs_free_error.append(error_MSDs_cell)
        mean_MSDs_free, _ = self.calc_mean_error_different_lengths(average_MSDs_free)
        mean_error_MSDs_free, _ = self.calc_mean_error_different_lengths(average_MSDs_free_error)

        average_MSDs_notype = []
        average_MSDs_notype_error = []
        for cell in self.trajectories_notype_cells:
            MSDs_cell = [trajectory.MSDs for trajectory in cell]
            mean_MSDs_cell, error_MSDs_cell = self.calc_mean_error_different_lengths(MSDs_cell)
            average_MSDs_notype.append(mean_MSDs_cell)
            average_MSDs_notype_error.append(error_MSDs_cell)
        mean_MSDs_notype, _ = self.calc_mean_error_different_lengths(average_MSDs_notype)
        mean_error_MSDs_notype, _ = self.calc_mean_error_different_lengths(average_MSDs_notype_error)

        self.plot_MSD_types(mean_MSDs_immob, mean_error_MSDs_immob, mean_MSDs_conf, mean_error_MSDs_conf,
                       mean_MSDs_free, mean_error_MSDs_free, mean_MSDs_notype, mean_error_MSDs_notype,
                        MSD_delta_t_n, y_lim)


    def plot_bar_log_bins(self):
        self.diff_fig = plt.figure()
        plt.subplot(111, xscale="log")
        (_, caps, _) = plt.errorbar(self.hist_diffusion, self.mean_frequencies, yerr=self.mean_error, capsize=4,
                                    label="relative frequency")
        for cap in caps:
            cap.set_markeredgewidth(1)
        for c, i in enumerate(self.mean_frequencies):
            if i != 0 and c == 0:
                xlim_min = self.hist_diffusion[c]
                break
            elif i != 0:
                xlim_min = self.hist_diffusion[c-1]
                break
        for c, i in enumerate(np.flip(self.mean_frequencies)):
            if i != 0 and c == 0:
                xlim_max = self.hist_diffusion[len(self.hist_diffusion)-c]
                break
            elif i != 0:
                xlim_max = self.hist_diffusion[len(self.hist_diffusion)-c+1]
                break
        plt.xlim(xlim_min, xlim_max)
        plt.legend()
        plt.title("Distribution of diffusion coefficients")
        plt.ylabel("Normalized relative occurence [%]")
        plt.xlabel("D [\u03BCm\u00b2/s]")
        plt.show()

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
        log_hist_immob, log_hist_conf, log_hist_free, log_hist_fit_fail = [], [], [], []
        for cell_index in range(len(self.cell_trajectories_filtered)):
            log_Ds = np.zeros(len(self.cell_trajectories_filtered[cell_index]))
            trajectory_types = []
            cell_size = self.cell_sizes[cell_index]
            for trajectory_index in range(len(self.cell_trajectories_filtered[cell_index])):
                if self.cell_trajectories_filtered[cell_index][trajectory_index].D > 0:
                    log_Ds[trajectory_index] = np.log10(self.cell_trajectories_filtered[cell_index][trajectory_index].D)
                    trajectory_types.append((self.cell_trajectories_filtered[cell_index][trajectory_index].immobility,
                                             self.cell_trajectories_filtered[cell_index][trajectory_index].confined,
                                             self.cell_trajectories_filtered[cell_index][trajectory_index].analyse_successful))
            # log diffusion for all trajectories without D < 0
            log_Ds = np.delete(log_Ds, np.where(log_Ds == 0))
            log_diffusion_hist = self.calc_diffusion_frequencies(log_Ds, desired_bin_size, cell_size)
            self.hist_log_Ds.append(log_diffusion_hist)
            # log diffusion for trajectory types
            log_Ds_immob, log_Ds_conf, log_Ds_free, log_Ds_fit_fail = self.filter_types(log_Ds, trajectory_types)
            log_hist_immob_cell = self.calc_diffusion_frequencies(log_Ds_immob, desired_bin_size, cell_size)
            log_hist_conf_cell = self.calc_diffusion_frequencies(log_Ds_conf, desired_bin_size, cell_size)
            log_hist_free_cell = self.calc_diffusion_frequencies(log_Ds_free, desired_bin_size, cell_size)
            log_hist_fit_fail_cell = self.calc_diffusion_frequencies(log_Ds_fit_fail, desired_bin_size, cell_size)
            log_hist_immob.append(log_hist_immob_cell[:, 1])
            log_hist_conf.append(log_hist_conf_cell[:, 1])
            log_hist_free.append(log_hist_free_cell[:, 1])
            log_hist_fit_fail.append(log_hist_fit_fail_cell[:, 1])
        return log_hist_immob, log_hist_conf, log_hist_free, log_hist_fit_fail

    def filter_types(self, trajectory_info, trajectory_type):
        """
        Split an array of trajectory infos into 3 based on the 3 diffusion types.

        :param trajectory_info: Target info of trajectory.
        :param trajectory_type: Tuple of immobility and confined boolean, referring to the trajectory type.
        :return: 3 arrays of trajectory infos, separated based on diffusion type.
        """
        info_immob, info_conf, info_free, info_fit_fail = [], [], [], []
        for info, type in zip(trajectory_info, trajectory_type):
            if type == (1, 1, 1):
                info_immob.append(info)
            elif type == (0, 1, 1):
                info_conf.append(info)
            elif type == (0, 0, 1):
                info_free.append(info)
            elif type == (0, 1, 0):
                info_fit_fail.append(info)
            else:
                print("Error with type determination occurred, contact developer", info, type)
        return info_immob, info_conf, info_free, info_fit_fail

    def calc_diffusion_frequencies(self, log_diff, desired_bin, size):
        """
        :param log_diff: np array with log10(D) + frequency of one cell.
        :param size: cell size.
        :param desired_bin: bin size.
        """
        # min & max determined by diffusions_log_complete function
        min_bin = np.ceil(-6/desired_bin)*desired_bin
        max_bin = np.ceil(2/desired_bin)*desired_bin 
        bin_size = int(np.ceil((max_bin - min_bin)/desired_bin))
        hist = np.histogram(log_diff, range=(min_bin, max_bin), bins=bin_size)
        log_diffusion_hist = np.zeros([np.size(hist[0]), 2])
        log_diffusion_hist[:, 0] = hist[1][:-1]  # log(D)
        log_diffusion_hist[:, 1] = hist[0][:]  # freq
        log_diffusion_hist[:, 1] = log_diffusion_hist[:, 1] / size
        return log_diffusion_hist

    def calc_nonlogarithmic_diffusions(self):
        """
        Calculate the non logarithmic diffusion coefficients from log10(D) from histogram.
        """
        self.hist_diffusion = 10**self.hist_log_Ds[0][:,0]
        
    def determine_mean_frequency(self):
        """
        Mean frequency will be calculated based on all frequencies of cells.
        Normalize an array (sum of elements = 1) and represent is in percent (*100).
        """
        self.diffusion_frequencies = self.create_np_array(np.shape(self.hist_log_Ds)[1], len(self.cell_trajectories_filtered))
        for i in range (0, len(self.cell_trajectories_filtered)):
            self.diffusion_frequencies[:, i] = self.hist_log_Ds[i][:, 1]
        self.mean_frequencies = self.calc_mean_frequencies(self.diffusion_frequencies)
        
        self.normalization_factor = 100/np.sum(self.mean_frequencies)
        self.mean_frequencies = self.mean_frequencies * self.normalization_factor
        
    def calc_mean_error(self):
        """
        Standard deviation (N-1) divided by square root of number of elements.
        Normalize an array (sum of elements = 1) and represent is in percent (*100).
        """
        self.mean_error = np.std(self.diffusion_frequencies, axis=1, ddof=1)/(np.shape(self.diffusion_frequencies)[1])**(1/2)
        self.mean_error = self.mean_error * self.normalization_factor

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
        
    # Save

    def save_diff(self, trajectories):
        """
        diffusion info: col0 = id, col1= D, col2 = dD, col3 = MSD(0), col4 = chi2, col5 = length
        """
        self.diffusion_info = np.zeros([np.size(trajectories), 6])
        self.number_of_trajectories = np.shape(self.diffusion_info)[0]
        i = 0
        for trajectory in trajectories:
            self.diffusion_info[i, 0] = trajectory.trajectory_number
            self.diffusion_info[i, 1] = trajectory.D
            self.diffusion_info[i, 2] = trajectory.dD
            self.diffusion_info[i, 3] = trajectory.MSD_0
            self.diffusion_info[i, 4] = trajectory.chi_D
            self.diffusion_info[i, 5] = trajectory.length_trajectory
            i += 1
    
    def save_plots(self, trajectory):
        trajectory_number = trajectory.trajectory_number
        dt_D = trajectory.MSD_D[:, 0]
        MSD_D = trajectory.MSD_D[:, 1]
        fit_D = trajectory.MSD_D[:, 2]
        residues_D = trajectory.MSD_D[:, 3]
        dt_r = trajectory.MSD_fit[:, 0]
        MSD_r = trajectory.MSD_fit[:, 1]
        fit_r = trajectory.MSD_fit[:, 2]
        residues_r = trajectory.MSD_fit[:, 3]
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
            self.rossier_info[i, 0] = trajectory.trajectory_number
            # if trajectory is immobile the analysis was never made, all analysis output is 0 by default
            if trajectory.immobility:
                self.rossier_info[i, 1] = trajectory.immobility
                self.rossier_info[i, 2] = False
                self.rossier_info[i, 3] = False
                self.rossier_info[i, 4] = False  # no analyse was made -> False
                self.rossier_info[i, 5] = 0
                self.rossier_info[i, 6] = 0
                self.rossier_info[i, 7] = 0
                self.rossier_info[i, 8] = 0
                self.rossier_info[i, 9] = 0
                self.rossier_info[i, 10] = 0
                self.rossier_info[i, 11] = 0
            # if trajectory is confined it is not immobile and not free
            elif trajectory.confined:
                self.rossier_info[i, 1] = trajectory.immobility
                self.rossier_info[i, 2] = trajectory.confined
                self.rossier_info[i, 3] = not trajectory.confined
            # if trajectory is free it is not confined and not immobile
            elif not trajectory.confined:
                self.rossier_info[i, 1] = trajectory.immobility
                self.rossier_info[i, 2] = trajectory.confined
                self.rossier_info[i, 3] = not trajectory.confined
            # analysis is made for trajectories not immobile -> if analysis was successful, output gets values
            if trajectory.analyse_successful and not trajectory.immobility:
                self.rossier_info[i, 4] = trajectory.analyse_successful
                self.rossier_info[i, 5] = trajectory.tau
                self.rossier_info[i, 6] = trajectory.dtau
                self.rossier_info[i, 7] = trajectory.r
                self.rossier_info[i, 8] = trajectory.dr
                self.rossier_info[i, 9] = trajectory.D_conf
                self.rossier_info[i, 10] = trajectory.dD_conf
                self.rossier_info[i, 11] = trajectory.chi_MSD_fit
            # if analysis was not successful -> output is 0 by default
            elif not trajectory.analyse_successful and not trajectory.immobility:
                self.rossier_info[i, 1] = trajectory.immobility
                self.rossier_info[i, 2] = False
                self.rossier_info[i, 3] = False
                self.rossier_info[i, 4] = trajectory.analyse_successful
                self.rossier_info[i, 5] = 0
                self.rossier_info[i, 6] = 0
                self.rossier_info[i, 7] = 0
                self.rossier_info[i, 8] = 0
                self.rossier_info[i, 9] = 0
                self.rossier_info[i, 10] = 0
                self.rossier_info[i, 11] = 0
            i += 1
