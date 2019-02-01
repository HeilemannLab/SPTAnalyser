# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 09:37:49 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import datetime


class Precision():
    def __init__(self):
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
        
    def load_localization_file(self):
        """
        Insert position uncertainties of x and y.
        position_uncertainties col0 = x
        position_uncertainties col1 = y
        """
        x_uncertainty_index = list(self.column_order.keys())[list(self.column_order.values()).index('"Position-0-0-uncertainty"')]
        y_uncertainty_index = list(self.column_order.keys())[list(self.column_order.values()).index('"Position-1-0-uncertainty"')]
        self.position_uncertainties = np.loadtxt(self.file_name, usecols = (x_uncertainty_index, y_uncertainty_index)) # col0 = x uncertainty, col1 = y uncertainty
            
    def log_columns(self):
        """
        Natural logarithm of position uncertainties x and y.
        position_uncertainties_log col0 = ln(x)
        position_uncertainties_log col1 = ln(y)
        """
        self.position_uncertainties_log = np.zeros([np.size(self.position_uncertainties[:,0]),2])
        log_x = np.log(self.position_uncertainties[:,0])
        log_y = np.log(self.position_uncertainties[:,1])
        self.position_uncertainties_log[:,0] = log_x
        self.position_uncertainties_log[:,1] = log_y
        
    def hist_x(self):
        """
        Create hist.
        position_uncertainties_hist_x col0 = bins
        position_uncertainties_hist_x col1 = normalized frequencies
        """
        max_bin = int(np.ceil(np.ceil(self.position_uncertainties[:,0].max()/0.5)*0.5))  # max mjd ceiled and divisible by 20
        bin_size = int(np.ceil(np.ceil(self.position_uncertainties[:,0].max())/0.5))  # divides the bin range in sizes -> desired bin = max_bin/bin_size
        hist = np.histogram(self.position_uncertainties[:,0],
                                        range = (0, max_bin),
                                        bins = bin_size,
                                        density = True)
        #print("histx", hist[0])
        self.position_uncertainties_hist_x = np.zeros([np.size(hist[0]),2])
        self.position_uncertainties_hist_x[:,0] = hist[1][:-1]  #hist0 = freq, hist1 = bins, hist1>hist0
        self.position_uncertainties_hist_x[:,1] = hist[0][:]
        self.position_uncertainties_hist_x[:,1] = self.normalize_hist(self.position_uncertainties_hist_x[:,1])

    def hist_x_log(self):
        max_bin = int(np.ceil(np.ceil(self.position_uncertainties_log[:,0].max()/0.05)*0.05))
        min_bin = np.ceil(self.position_uncertainties_log[:,0].min())
        bin_size = int(np.ceil(np.ceil(self.position_uncertainties_log[:,0].max())/0.05))
        hist = np.histogram(self.position_uncertainties_log[:,0],
                            range = (0, max_bin),
                            bins = bin_size,
                            density = True)
        self.position_uncertainties_hist_log_x = np.zeros([np.size(hist[0]),4])
        self.position_uncertainties_hist_log_x[:,0] = hist[1][:-1]
        self.position_uncertainties_hist_log_x[:,1] = hist[0][:]
        self.position_uncertainties_hist_log_x[:,1] = self.normalize_hist(self.position_uncertainties_hist_log_x[:,1])
        #print("hist x log", hist)
        
    def hist_y(self):
        """
        Create hist.
        position_uncertainties_hist_y col0 = bins
        position_uncertainties_hist_y col1 = normalized frequencies
        """
        max_bin = int(np.ceil(np.ceil(self.position_uncertainties[:,1].max()/0.5)*0.5))  # max mjd ceiled and divisible by 20
        bin_size = int(np.ceil(np.ceil(self.position_uncertainties[:,1].max())/0.5))  # divides the bin range in sizes -> desired bin = max_bin/bin_size
        hist = np.histogram(self.position_uncertainties[:,1],
                                        range = (0, max_bin),
                                        bins = bin_size,
                                        density = True)
        self.position_uncertainties_hist_y = np.zeros([np.size(hist[0]),2])
        self.position_uncertainties_hist_y[:,0] = hist[1][:-1]  
        self.position_uncertainties_hist_y[:,1] = hist[0][:]
        
    def hist_y_log(self):
        max_bin = int(np.ceil(np.ceil(self.position_uncertainties_log[:,1].max()/0.05)*0.05))
        #min_bin = np.ceil(self.position_uncertainties_log[:,1].min())
        bin_size = int(np.ceil(np.ceil(self.position_uncertainties_log[:,1].max())/0.05))
        hist = np.histogram(self.position_uncertainties_log[:,1],
                            range = (0, max_bin),
                            bins = bin_size,
                            density = True)
        self.position_uncertainties_hist_log_y = np.zeros([np.size(hist[0]),4])
        self.position_uncertainties_hist_log_y[:,0] = hist[1][:-1]
        self.position_uncertainties_hist_log_y[:,1] = hist[0][:]
        self.position_uncertainties_hist_log_y[:,1] = self.normalize_hist(self.position_uncertainties_hist_log_y[:,1])
        #print("hist y log", hist)
        
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
        :return: gaus function
        """
        #A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/2.*sigma**2)
    
    def gauss_fit(self): 
        self.init_a = self.position_uncertainties_hist_log_x[:,1].max()
        self.init_mu = np.median(self.position_uncertainties_hist_log_x[:,0])
        self.init_sigma = np.std(self.position_uncertainties_hist_log_x[:,0])
        #print(self.init_a, self.init_mu, self.init_sigma)
        coeff_x, var_matrix_x = curve_fit(self.gauss_func, self.position_uncertainties_hist_log_x[:,0],
                                          self.position_uncertainties_hist_log_x[:,1],
                                          p0=(self.init_a, self.init_mu, self.init_sigma))
        self.position_uncertainties_hist_log_x[:,2] = self.gauss_func(self.position_uncertainties_hist_log_x[:,0],
                                            coeff_x[0], coeff_x[1], coeff_x[2])  # fit
        self.position_uncertainties_hist_log_x[:,3] = self.position_uncertainties_hist_log_x[:,1] - self.position_uncertainties_hist_log_x[:,2]  # residues
        
        coeff_y, var_matrix_y = curve_fit(self.gauss_func, self.position_uncertainties_hist_log_y[:,0],
                                          self.position_uncertainties_hist_log_y[:,1],
                                          p0=(self.init_a, self.init_mu, self.init_sigma))
        self.position_uncertainties_hist_log_y[:,2] = self.gauss_func(self.position_uncertainties_hist_log_y[:,0],
                                              coeff_y[0], coeff_y[1], coeff_y[2])  # fit
        self.position_uncertainties_hist_log_y[:,3] = self.position_uncertainties_hist_log_y[:,1] - self.position_uncertainties_hist_log_y[:,2]  # residues
        self.mean_x = np.exp(coeff_x[1])
        self.mean_y = np.exp(coeff_y[1])
        print("The mean position uncertainty is %.3f nm in x and %.3f nm in y direction." %(self.mean_x, self.mean_y))
        
    def plot_hist(self, x_axis, y_axis, width, fit=False, fit_data=[], colour="gray", fit_style="--"):
        fig = plt.figure()
        sp = fig.add_subplot(1, 1, 1)  # only 1 plot
        sp.bar(x_axis, y_axis, 
               align = "center",
               width = width,  # width = bin size
               color = colour)
        if fit:
            sp.plot(x_axis, fit_data, fit_style)  # "b--" change colour, line style "m-" ...
        sp.set_title("Histogram of position uncertainties in y direction")
        sp.set_xlabel("Position uncertainty [nm]")
        sp.set_ylabel("Fraction")
        #plt.savefig(out_file_name)
        plt.show()  # print the graph
        
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
        out_file_name = directory + "\ " + year + month + day + "_" + base_name + "_localization_uncertainty" + "_x_frequencies.txt" # System independent?
        header = "Position uncertainty [nm]\t fraction\t"
        np.savetxt(out_file_name, 
                   X=self.position_uncertainties_hist_x,
                   fmt = ("%.1f","%.4e"),
                   header = header)   

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
        out_file_name = directory + "\ " + year + month + day + "_" + base_name + "_localization_uncertainty" + "_y_frequencies.txt" # System independent?
        header = "Position uncertainty [nm]\t fraction\t"
        np.savetxt(out_file_name, 
                   X=self.position_uncertainties_hist_y,
                   fmt = ("%.1f","%.4e"),
                   header = header) 
        
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
        out_file_name = directory + "\ " + year + month + day + "_" + base_name + "_localization_uncertainty" + "_ln_x_frequencies.txt" # System independent?
        header = "ln(position uncertainty) \t fraction\t eponential fit\t residues\t"
        np.savetxt(out_file_name, 
                   X=self.position_uncertainties_hist_log_x,
                   fmt = ("%.2f","%.4e","%.4e","%.4e"),
                   header = header) 
        
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
        out_file_name = directory + "\ " + year + month + day + "_" + base_name + "_localization_uncertainty" + "_ln_y_frequencies.txt" # System independent?
        header = "ln(position uncertainty) \t fraction\t exponential fit\t residues\t "
        np.savetxt(out_file_name, 
                   X=self.position_uncertainties_hist_log_y,
                   fmt = ("%.2f","%.4e","%.4e","%.4e"),
                   header = header) 
        
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
        out_file_name = directory + "\ " + year + month + day + "_" + base_name + "_localization_uncertainty.txt"
        file = open(out_file_name, 'w')
        if not (file.closed):
            file.write("localization_uncertainty x in nm\t localization_uncertainty y in nm\n")
            file.write("%.4e\t%.4e\n" %(self.mean_x, self.mean_y))
            file.close()
        else:
            print("error: could not open file %s. Make sure the folder does exist" %(out_file_name))

    def run_precision(self):
        self.load_localization_file()
        self.log_columns()
        self.hist_x()
        self.hist_x_log()
        self.hist_y()
        self.hist_y_log()
        self.gauss_fit()
        self.plot_hist(self.position_uncertainties_hist_x[:,0], self.position_uncertainties_hist_x[:,1], 0.5)
        self.plot_hist(self.position_uncertainties_hist_y[:,0], self.position_uncertainties_hist_y[:,1], 0.5)
        self.plot_hist(self.position_uncertainties_hist_log_x[:,0], self.position_uncertainties_hist_log_x[:,1], 0.05, 
                       fit=True, fit_data = self.position_uncertainties_hist_log_x[:,2])     
        self.plot_hist(self.position_uncertainties_hist_log_y[:,0], self.position_uncertainties_hist_log_y[:,1], 0.05, 
                       fit=True, fit_data = self.position_uncertainties_hist_log_y[:,2])
        
    def save_precision(self, directory, base_name):
        #directory = directory
        #base_name = base_name
        self.save_x_hist(directory, base_name)
        self.save_y_hist(directory, base_name)
        self.save_x_hist_log(directory, base_name)
        self.save_y_hist_log(directory, base_name)
        self.save_fit_results(directory, base_name)


def main():
    file_name = "C:\\Users\\pcoffice37\\Documents\\rapidStorm_loc\\cell_17_MMStack_Pos0.ome.txt"  # testing file name
    column_order = {1: '"Position-0-0-uncertainty"', 3: '"Position-1-0-uncertainty"'}
    precision = Precision(file_name, column_order)
    precision.run_precision()

    
if __name__ == "__main__":
    main()
    