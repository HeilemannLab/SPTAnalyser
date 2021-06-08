"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Determine the false positive rate of localizations per cell in %, based on background noise.
"""

import os
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class ExpNoiseRate():
    def __init__(self):
        self.cell_pd = []
        self.cell_names = []
        self.bg_names = []
        self.bg_pd = []
        self.roi_pd = []
        self.cell_locs = []
        self.bg_locs = []
        self.exp_noise_rates = []
        self.figures = []
        self.figure_names = ["density_cells", "density_background", "exp_noise_rate"]

    def get_files(self, path, suffix):
        """
        Read all *.csv files from a directory (localization files, ignore tracked files) and create pandas-frames.
        """
        files_pd, file_names = [], []
        for file in os.listdir(path):
            if file.endswith(suffix) and "tracked" not in file:
                file_pd = pd.read_csv(path + "\\" + file)
                files_pd.append(file_pd)
                file_names.append(os.path.splitext(os.path.split(path + "\\" + file)[1])[0])
        return files_pd, file_names

    def determine_number_locs(self):
        """
        Get number of localizations per frame and area.
        """
        cell_locs = []
        for i, cell_name in enumerate(self.cell_names):
            for j, cell_size in self.roi_pd.iterrows():
                if cell_size[0] in [cell_name, cell_name + ".csv"]:
                    cell_locs.append(max(self.cell_pd[i]["id"]) / (cell_size[1]*max(self.cell_pd[i]["frame"])))
        return cell_locs

    def determine_number_locs_bg(self, bg_size):
        """
        Get number of localizations per frame and area.
        """
        return [max(file["id"])/(bg_size*max(file["frame"])) for file in self.bg_pd]

    def calc_exp_noise_rates(self, mean_bg_loc):
        """
        Calculate exp_noise_rate: mean(localizations_background) / localizations_cell *  100.
        """
        return [mean_bg_loc / cell_loc * 100 for cell_loc in self.cell_locs]

    def plot_box(self, data_points, title, ylabel):
        df = pd.DataFrame(data_points, columns=[""])
        fig = plt.figure(figsize=(3, 5))
        ax = sns.violinplot(data=df, color="cornflowerblue", showmeans=True,
                         meanprops={"marker": "s", "markerfacecolor": "white", "markeredgecolor": "0.25"})
        ax = sns.swarmplot(data=df, color="0.25")
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        plt.show()
        self.figures.append(fig)

    def save_results(self, directory, folder_name, save_fig):
        os.mkdir(directory + "\\" + folder_name)
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        day = str(now.day)
        if len(month) == 1:
            month = str(0) + month
        if len(day) == 1:
            day = str(0) + day
        # save cell info & exp noise rate
        out_file_name = directory + "\\" + folder_name + "\\" + year + month + day + "_cells_exp_noise_rate.txt"
        header = "cell name\tcell density [localizations per px² frame]\texp_noise_rate [%]\t"
        max_name_length = max([len(i) for i in self.cell_names])
        data = np.zeros(np.array(self.cell_names).size, dtype=[("col1", "U"+str(max_name_length)),
                                                               ("col2", float), ("col3", float)])
        data["col1"] = np.array(self.cell_names)
        data["col2"] = np.array(self.cell_locs)
        data["col3"] = np.array(self.exp_noise_rates)
        np.savetxt(out_file_name, X=data, fmt=("%10s", "%.4e", "%.4e"), header=header)
        # save mean exp noise rate
        out_file_name = directory + "\\" + folder_name + "\\" + year + month + day + "_cells_exp_noise_rate_mean.txt"
        file = open(out_file_name, 'w+')
        if not file.closed:
            file.write("mean exp_noise_rate [%]\n")
            file.write("%.4f" %(np.mean(self.exp_noise_rates)))
            file.close()
        # save background info
        out_file_name = directory + "\\" + folder_name + "\\" + year + month + day + "_background.txt"
        header = "background name\tbackground density [localizations per px² frame]\t"
        max_name_length = max([len(i) for i in self.bg_names])
        data = np.zeros(np.array(self.bg_names).size, dtype=[("col1", "U"+str(max_name_length)), ("col2", float)])
        data["col1"] = np.array(self.bg_names)
        data["col2"] = np.array(self.bg_locs)
        np.savetxt(out_file_name, X=data, fmt=("%10s", "%.4e"), header=header)
        # save figures
        if save_fig:
            for name, fig in zip(self.figure_names, self.figures):
                fig.savefig(directory + "\\" + folder_name + "\\" + year + month + day + "_exp_noise_rate_" + name + ".pdf",
                               format="pdf", transparent=True, bbox_inches="tight")

    def run_exp_noise_rate(self, cell_dir, bg_dir, roi_path, suffix, bg_size):
        self.cell_pd, self.cell_names = self.get_files(cell_dir, suffix)
        self.bg_pd, self.bg_names = self.get_files(bg_dir, suffix)
        self.roi_pd = pd.read_csv(roi_path, skiprows=2)
        self.cell_locs = self.determine_number_locs()
        self.bg_locs = self.determine_number_locs_bg(bg_size)
        self.exp_noise_rates = self.calc_exp_noise_rates(np.mean(self.bg_locs))
        # check if rois for all cells are loaded
        if len(self.cell_names) != len(self.cell_locs):
            print("The number of localized cell files and provided areas in the log file has to be the same. "
                  "{} localized cell files and {} areas were loaded, please check your data!".format((len(self.cell_names)), len(self.cell_locs)))
        else:
            self.plot_box(self.cell_locs, "cells", "number of localizations per px² frame")
            self.plot_box(self.bg_locs, "background", "number of localizations per px² frame")
            self.plot_box(self.exp_noise_rates, "exp noise rate", "exp noise rate / %")
            print("name: locs per px² frame, exp noise rate %:")
            for name, loc, rate in zip(self.cell_names, self.cell_locs, self.exp_noise_rates):
                print(name + ":", str(loc) + ",", rate)
