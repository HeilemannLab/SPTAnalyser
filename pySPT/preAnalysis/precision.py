"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Extract the localization precision from a rapidSTORM (x and y) or thunderSTORM (x and y are the same -> 1 value) file
and build a mean value (log-transform values, average values, retransform average).
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
import datetime
import pandas as pd
import math
import os
from ..hmm import microscope


class Precision():
    def __init__(self):
        self.software = ""  # either rapidSTORM or thunderSTORM
        self.file_name = ""
        self.column_order = {}
        self.position_uncertainties = []  # uncertainties of x and y
        self.position_uncertainties_log = []  # log uncertainties of x and y
        self.position_uncertainties_hist_x = []  # x freq 
        self.position_uncertainties_hist_y = []  # x freq
        self.position_uncertainties_hist_log_x = []  # xbin xfreq xfit xres 
        self.position_uncertainties_hist_log_y = []  # ybin yfreq yfit yres
        self.mean_x = 0.
        self.mean_y = 0.
        self.init_a = 0.
        self.init_mu = 0.
        self.init_sigma = 0.
        self.mean_values = []  # precisions of all files in folder
        self.figures = []
        self.figure_names = []
        self.figure_box = []
        
    def run_precision(self):
        """
        Different functions for thunderSTORM and rapidSTORM because thunderSTORM has 1 uncertainty value in the
        x/y-plane, rapidSTORM has uncertainty values for x and y respectively.
        """
        self.load_localization_file()
        if self.software == "ThunderSTORM":
            self.ts_log_columns()
            self.hist_x()
            self.hist_x_log()
            self.gauss_fit()
            self.plot_hist(self.position_uncertainties_hist_x[:, 0], self.position_uncertainties_hist_x[:, 1], 0.5)
            self.plot_hist(self.position_uncertainties_hist_log_x[:, 0], self.position_uncertainties_hist_log_x[:, 1],
                           0.05, fit=True, fit_data=self.position_uncertainties_hist_log_x[:, 2], log=True)
        elif self.software == "rapidSTORM":
            self.rs_log_columns()
            self.hist_x()
            self.hist_x_log()
            self.hist_y()
            self.hist_y_log()
            self.gauss_fit()
            self.plot_hist(self.position_uncertainties_hist_x[:, 0], self.position_uncertainties_hist_x[:, 1], 0.5, direction="x")
            self.plot_hist(self.position_uncertainties_hist_y[:, 0], self.position_uncertainties_hist_y[:, 1], 0.5, direction="y")
            self.plot_hist(self.position_uncertainties_hist_log_x[:, 0], self.position_uncertainties_hist_log_x[:, 1],
                           0.05, fit=True, fit_data=self.position_uncertainties_hist_log_x[:, 2], log=True, direction="x")
            self.plot_hist(self.position_uncertainties_hist_log_y[:, 0], self.position_uncertainties_hist_log_y[:, 1],
                           0.05, fit=True, fit_data=self.position_uncertainties_hist_log_y[:, 2], log=True, direction="y")
        
    def load_localization_file(self):
        """
        Get position uncertainties of x (col0) and y (col1) from localization file.
        """
        if self.software == "ThunderSTORM":
            # different column names depending on TS version
            uncertainty_column_names = ["uncertainty_xy [nm]", "uncertainty [nm]"]
            file = pd.read_csv(self.file_name)
            for i in uncertainty_column_names:
                if i in list(file.columns.values):
                    column_name = i
            self.position_uncertainties = np.zeros([np.shape(file)[0], 2])
            self.position_uncertainties[:,0] = file[column_name]
        elif self.software == "rapidSTORM":
            x_uncertainty_index = list(self.column_order.keys())[list(self.column_order.values()).index('"Position-0-0-uncertainty"')]
            y_uncertainty_index = list(self.column_order.keys())[list(self.column_order.values()).index('"Position-1-0-uncertainty"')]
            # col0 = x uncertainty, col1 = y uncertainty
            self.position_uncertainties = np.loadtxt(self.file_name, usecols=(x_uncertainty_index, y_uncertainty_index))

    def rs_log_columns(self):
        """
        Natural logarithm of position uncertainties x and y.
        position_uncertainties_log col0 = ln(x)
        position_uncertainties_log col1 = ln(y)
        """
        self.position_uncertainties_log = np.zeros([np.size(self.position_uncertainties[:, 0]), 2])
        log_x = np.log(self.position_uncertainties[:, 0])
        log_y = np.log(self.position_uncertainties[:, 1])
        self.mean_x = math.exp(np.mean(log_x))
        self.mean_y = math.exp(np.mean(log_y))
        self.position_uncertainties_log[:, 0] = log_x
        self.position_uncertainties_log[:, 1] = log_y
        
    def ts_log_columns(self):
        """
        Natural logarithm of position uncertainties x.
        position_uncertainties_log col0 = ln(x)
        """
        self.position_uncertainties_log = np.zeros([np.size(self.position_uncertainties[:, 0]), 2])
        log = np.log(self.position_uncertainties[:, 0])
        self.mean_x = math.exp(np.mean(log))
        self.position_uncertainties_log[:, 0] = log

    def hist_x(self):
        """
        Create hist.
        position_uncertainties_hist_x col0 = bins
        position_uncertainties_hist_x col1 = normalized frequencies
        """
        max_bin = int(np.ceil(np.ceil(self.position_uncertainties[:, 0].max()/0.5)*0.5))
        bin_size = int(np.ceil(np.ceil(self.position_uncertainties[:, 0].max())/0.5))
        hist = np.histogram(self.position_uncertainties[:, 0], range=(0, max_bin), bins=bin_size, density=True)
        self.position_uncertainties_hist_x = np.zeros([np.size(hist[0]), 2])
        self.position_uncertainties_hist_x[:, 0] = hist[1][:-1]  # hist0 = freq, hist1 = bins, hist1>hist0
        self.position_uncertainties_hist_x[:, 1] = hist[0][:]
        self.position_uncertainties_hist_x[:, 1] = self.normalize_hist(self.position_uncertainties_hist_x[:, 1])

    def hist_x_log(self):
        max_bin = int(np.ceil(np.ceil(self.position_uncertainties_log[:, 0].max()/0.05)*0.05))
        bin_size = int(np.ceil(np.ceil(self.position_uncertainties_log[:, 0].max())/0.05))
        hist = np.histogram(self.position_uncertainties_log[:, 0], range=(0, max_bin), bins=bin_size, density=True)
        self.position_uncertainties_hist_log_x = np.zeros([np.size(hist[0]), 4])
        self.position_uncertainties_hist_log_x[:, 0] = hist[1][:-1]
        self.position_uncertainties_hist_log_x[:, 1] = hist[0][:]
        self.position_uncertainties_hist_log_x[:, 1] = self.normalize_hist(self.position_uncertainties_hist_log_x[:, 1])
        
    def hist_y(self):
        """
        Create hist.
        position_uncertainties_hist_y col0 = bins
        position_uncertainties_hist_y col1 = normalized frequencies
        """
        max_bin = int(np.ceil(np.ceil(self.position_uncertainties[:, 1].max()/0.5)*0.5))
        # divides the bin range in sizes -> desired bin = max_bin/bin_size
        bin_size = int(np.ceil(np.ceil(self.position_uncertainties[:, 1].max())/0.5))
        hist = np.histogram(self.position_uncertainties[:, 1], range=(0, max_bin), bins=bin_size, density=True)
        self.position_uncertainties_hist_y = np.zeros([np.size(hist[0]), 2])
        self.position_uncertainties_hist_y[:, 0] = hist[1][:-1]
        self.position_uncertainties_hist_y[:, 1] = hist[0][:]
        
    def hist_y_log(self):
        max_bin = int(np.ceil(np.ceil(self.position_uncertainties_log[:, 1].max()/0.05)*0.05))
        bin_size = int(np.ceil(np.ceil(self.position_uncertainties_log[:, 1].max())/0.05))
        hist = np.histogram(self.position_uncertainties_log[:, 1], range=(0, max_bin), bins=bin_size, density=True)
        self.position_uncertainties_hist_log_y = np.zeros([np.size(hist[0]), 4])
        self.position_uncertainties_hist_log_y[:, 0] = hist[1][:-1]
        self.position_uncertainties_hist_log_y[:, 1] = hist[0][:]
        self.position_uncertainties_hist_log_y[:, 1] = self.normalize_hist(self.position_uncertainties_hist_log_y[:, 1])
        
    def normalize_hist(self, normalized_col):
        """
        Normalize a column and return it.
        """
        normalized_col = normalized_col / np.sum(normalized_col)
        return normalized_col
      
    def gauss_func(self, x, A, mu, sigma):
        """
        Gausfunction not area normalized.
        :param p: tupel of amplitude, mu and sigma
        :param x: x
        :return: Gauss function
        """
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    
    def gauss_fit(self): 
        self.init_a = self.position_uncertainties_hist_log_x[:, 1].max()
        self.init_mu = np.median(self.position_uncertainties_hist_log_x[:, 0])
        self.init_sigma = np.std(self.position_uncertainties_hist_log_x[:, 0])
        coeff_x, var_matrix_x = curve_fit(self.gauss_func, self.position_uncertainties_hist_log_x[:, 0],
                                          self.position_uncertainties_hist_log_x[:, 1],
                                          p0=(self.init_a, self.init_mu, self.init_sigma))
        self.position_uncertainties_hist_log_x[:,2] = self.gauss_func(self.position_uncertainties_hist_log_x[:, 0],
                                            coeff_x[0], coeff_x[1], coeff_x[2])  # fit
        self.position_uncertainties_hist_log_x[:,3] = self.position_uncertainties_hist_log_x[:,1] - self.position_uncertainties_hist_log_x[:, 2]  # residues
        if self.software == "ThunderSTORM":
            print("The mean localization uncertainty is %.3f nm in the x/y-plane." %(self.mean_x))
        if self.software == "rapidSTORM":
            coeff_y, var_matrix_y = curve_fit(self.gauss_func, self.position_uncertainties_hist_log_y[:, 0],
                                              self.position_uncertainties_hist_log_y[:, 1],
                                              p0=(self.init_a, self.init_mu, self.init_sigma))
            self.position_uncertainties_hist_log_y[:, 2] = self.gauss_func(self.position_uncertainties_hist_log_y[:, 0],
                                                  coeff_y[0], coeff_y[1], coeff_y[2])  # fit
            self.position_uncertainties_hist_log_y[:, 3] = self.position_uncertainties_hist_log_y[:, 1] - self.position_uncertainties_hist_log_y[:, 2]  # residues
            print("The mean localization uncertainty is %.3f nm in x and %.3f nm in y direction." %(self.mean_x, self.mean_y))
        
    def plot_hist(self, x_axis, y_axis, width, fit=False, fit_data=[], colour="gray", fit_style="--c",
                  log=False, direction=False):
        fig = plt.figure()
        sp = fig.add_subplot(1, 1, 1)
        sp.bar(x_axis, y_axis, align="center", width=width, color=colour, label="fraction")
        if fit:
            sp.plot(x_axis, fit_data, fit_style, label = "gauss fit")
        sp.legend()
        sp.set_ylabel("Fraction")
        if log and direction:
            sp.set_title("Histogram of logarithmic localization uncertainties in {} direction".format(direction))
            sp.set_xlabel("ln(localization uncertainty)")
            title = "localization_uncertainty_ln_histogram_" + str(direction)
        if not log and direction:
            sp.set_title("Histogram of localization uncertainties in {} direction".format(direction))
            sp.set_xlabel("Localization uncertainty [nm]")
            title = "localization_uncertainty_histogram_" + str(direction)
        if log and not direction:
            sp.set_title("Histogram of logarithmic localization uncertainties")
            sp.set_xlabel("ln(localization uncertainty)")
            title = "localization_uncertainty_ln_histogram"
        if not log and not direction:
            sp.set_title("Histogram of localization uncertainties")
            sp.set_xlabel("Localization uncertainty [nm]")
            title = "localization_uncertainty_histogram"
        plt.show()  # print the graph
        self.figures.append(fig)
        self.figure_names.append(title)
        
    def save_x_hist(self, directory, base_name):
        """
        Output file with cols:
        col0 = Position uncertainty in nm
        col1 = Frequency
        """
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        day = str(now.day)
        if len(month) == 1:
            month = str(0) + month
        if len(day) == 1:
            day = str(0) + day
        if self.software == "ThunderSTORM":
            out_file_name = directory + "\\" + year + month + day + "_" + base_name + "_localization_uncertainty" + "_histogram.txt"
            header = "Localization uncertainty [nm]\tfraction\t"
        elif self.software == "rapidSTORM":
            out_file_name = directory + "\\" + year + month + day + "_" + base_name + "_localization_uncertainty" + "_x_histogram.txt"
            header = "Localization uncertainty x [nm]\tfraction\t"
        np.savetxt(out_file_name, X=self.position_uncertainties_hist_x, fmt=("%.1f", "%.4e"), header=header)

    def save_y_hist(self, directory, base_name):
        """
        Output file with cols:
        col0 = Position uncertainty in nm
        col1 = Frequency
        """
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        day = str(now.day)
        if len(month) == 1:
            month = str(0) + month
        if len(day) == 1:
            day = str(0) + day
        out_file_name = directory + "\\" + year + month + day + "_" + base_name + "_localization_uncertainty" + "_y_histogram.txt"
        header = "Localization uncertainty y [nm]\tfraction\t"
        np.savetxt(out_file_name, X=self.position_uncertainties_hist_y, fmt=("%.1f", "%.4e"), header=header)
        
    def save_x_hist_log(self, directory, base_name):
        """
        Output file with cols:
        col0 = Position uncertainty in nm
        col1 = Frequency
        """
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        day = str(now.day)
        if len(month) == 1:
            month = str(0) + month
        if len(day) == 1:
            day = str(0) + day
        if self.software == "ThunderSTORM":
            out_file_name = directory + "\\" + year + month + day + "_" + base_name + "_localization_uncertainty" + "_ln_histogram.txt"
            header = "ln(localization uncertainty)\tfraction\tgauss fit\tresidues\t"
        elif self.software == "rapidSTORM":
            out_file_name = directory + "\\" + year + month + day + "_" + base_name + "_localization_uncertainty" + "_ln_x_histogram.txt"
            header = "ln(localization uncertainty) x\tfraction\tgauss fit\tresidues\t"
        np.savetxt(out_file_name, X=self.position_uncertainties_hist_log_x, fmt=("%.2f", "%.4e", "%.4e", "%.4e"),
                   header=header)
        
    def save_y_hist_log(self, directory, base_name):
        """
        Output file with cols:
        col0 = Position uncertainty in nm
        col1 = Frequency
        """
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        day = str(now.day)
        if len(month) == 1:
            month = str(0) + month
        if len(day) == 1:
            day = str(0) + day
        out_file_name = directory + "\\" + year + month + day + "_" + base_name + "_localization_uncertainty" + "_ln_y_histogram.txt"
        header = "ln(localization uncertainty) y\tfraction\tgauss fit\tresidues\t "
        np.savetxt(out_file_name, X=self.position_uncertainties_hist_log_y, fmt=("%.2f", "%.4e", "%.4e", "%.4e"),
                   header=header)
        
    def save_fit_results(self, directory, base_name):
        """
        Output file with p_bleach, k, kv.
        """
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        day = str(now.day)
        if len(month) == 1:
            month = str(0) + month
        if len(day) == 1:
            day = str(0) + day
        out_file_name = directory + "\\" + year + month + day + "_" + base_name + "_localization_uncertainty.txt"
        file = open(out_file_name, 'w')
        if not file.closed:
            if self.software == "ThunderSTORM":
                file.write("localization uncertainty [nm]\n")
                file.write("%.4e\n" %(self.mean_x))
            elif self.software == "rapidSTORM":
                file.write("localization uncertainty x in nm\tlocalization uncertainty y in nm\n")
                file.write("%.4e\t%.4e\n" %(self.mean_x, self.mean_y))
            file.close()
        else:
            print("error: could not open file %s. Make sure the folder does exist" %(out_file_name))
        
    def save_precision(self, directory, base_name):
        if self.software == "ThunderSTORM":
            self.save_x_hist(directory, base_name)
            self.save_x_hist_log(directory, base_name)
        elif self.software == "rapidSTORM":
            self.save_x_hist(directory, base_name)
            self.save_y_hist(directory, base_name)
            self.save_x_hist_log(directory, base_name)
            self.save_y_hist_log(directory, base_name)
        self.save_fit_results(directory, base_name)
        print("Results are saved.")
        
    def save_hmm_microscope(self, directory, px_size, integration_time):
        """
        For the HMM-analysis a microscope.txt file for each cell is needed which contains the localization uncertainty
        (static), camera pixel size and integration time. This file will be saved in the SPTAnalysis/preAnalysis
        per cell optionally.
        """
        if self.software == "ThunderSTORM":
            loc_uncertainty = self.mean_x
        elif self.software == "rapidSTORM":
            loc_uncertainty = (self.mean_x + self.mean_y) / 2
        out_file_name = directory + "\\" + "microscope.txt"
        file = open(out_file_name, "w")
        if not file.closed:
            file.write("# SMLMS Microscope File \n")
            file.write("# pxl Size[nm] \n")
            file.write("# integration Time [s] \n")
            file.write("# localization precision [nm] \n")
            file.write("%.6e \n" %(float(px_size)))
            file.write("%.6e \n" %(float(integration_time)))
            file.write("%.6e \n" %(loc_uncertainty))
            
    def run_save_plots(self, directory, base_name):
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        day = str(now.day)
        if len(month) == 1:
            month = str(0) + month
        if len(day) == 1:
            day = str(0) + day
        for figure, name in zip(self.figures, self.figure_names):
            figure.savefig(directory + "\\" + year + month + day + "_" + base_name + "_" + name + ".pdf", format="pdf", transparent=True)

    # precision per folder

    def get_loc_files(self, target_dir):
        """
        Get all localization csv files and return pandas frames.
        """
        files, file_names = [], []
        for file in os.listdir(target_dir):
            if file.endswith("csv") and "tracked" not in file:
                files.append(target_dir + "\\" + file)
                file_names.append(file)
        files = [pd.read_csv(i) for i in files]
        return files, file_names

    def get_precisions(self, files, file_names):
        """
        Return mean uncertainty per localization file.
        """
        self.mean_values = []
        for file in files:
            log_uncertainties = [np.log10(i) for i in file["uncertainty_xy [nm]"]]
            mean_log = np.mean(log_uncertainties)
            self.mean_values.append(10**mean_log)
        for name, precision in zip(file_names, self.mean_values):
            print(os.path.splitext(name)[0] + ":", str(precision))
        print(self.mean_values)

    def plot_box(self):
        df = pd.DataFrame(self.mean_values)
        fig = plt.figure(figsize=(3, 5))
        ax = sns.boxplot(data=df, color="cornflowerblue", showmeans=True,
                         meanprops={"marker": "s", "markerfacecolor": "white", "markeredgecolor": "0.25"})
        ax = sns.swarmplot(data=df, color="0.25")
        ax.set_ylabel("precision / nm")
        plt.show()
        self.figure_box = fig

    def save_precision_list(self, path, precision_lst, save_fig):
        try:
            os.mkdir(path)
        except FileExistsError:
            pass
        file = open(path + "\\precisions_list.txt", "w+")
        file.write(str(precision_lst))
        file.close()
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        day = str(now.day)
        if len(month) == 1:
            month = str(0) + month
        if len(day) == 1:
            day = str(0) + day
        if save_fig:
            self.figure_box.savefig(path + "\\" + year + month + day + "_precision.pdf",
                           format="pdf", transparent=True, bbox_inches= "tight")
        print("Results successfully saved.")

    def save_microscope(self, file_dir, pixel_size, dt):
        _, file_names = self.get_loc_files(file_dir)
        for i, error in zip(file_names, self.mean_values):
            folder_name = "SPTAnalyser_" + os.path.splitext(i)[0]
            try:
                os.mkdir(file_dir + "\\" + folder_name)
            except FileExistsError:
                pass
            try:
                os.mkdir(file_dir + "\\" + folder_name + "\\preAnalysis")
            except FileExistsError:
                pass
            mic = microscope.Microscope(dt, pixel_size, error, file_dir + "\\" + folder_name + "\\preAnalysis", ym_to_nm=False)
            mic.save_hmm_microscope()
