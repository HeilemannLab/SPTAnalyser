# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 15:40:08 2019

@author: Johanna Rahm, Research Group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt am Main.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

class PBleach():
    def __init__(self):
        self.file_name = ""
        self.mjds = []
        self.mjd_n_histogram = []
        self.p_bleach_results = []
        self.init_k = 0.01
        self.k = 0.01
        self.kcov = 0
        self.p_bleach = 0.0
        
    def load_seg_file(self, file_name):
        if not (file_name == ""):
            self.mjds = np.loadtxt(file_name, usecols = (1, 2)) # col0 = mjd, col1 = mjd_n
            #print(self.mjds)
        else:
            print("Insert a file name.")
    
    def count_mjd_n_frequencies(self):
        """
        Create histogram with bins = mjd_n and frequencies as np.ndarray
        """
        max_bin = self.mjds[:,1].max()  # max mjd_n value
        bin_size = int(max_bin)  # divides the bin range in sizes -> desired bin = max_bin/bin_size
        hist = np.histogram(self.mjds[:,1],
                            range = (0, max_bin),
                            bins = bin_size,
                            density = True)
        self.mjd_n_histogram = np.zeros([np.size(hist[0]),4])
        self.mjd_n_histogram [:,0] = hist[1][:-1]  # col0 = bins
        self.mjd_n_histogram [:,1] = hist[0][:]  # col1 = frequencies
        self.normalized_mjd_ns()  # normalize the histogram by the sum
        
    def normalized_mjd_ns(self):
        """
        Create normalized mjd_ns for histogram -> the sum equals 1.
        """
        self.mjd_n_histogram[:,1] = self.mjd_n_histogram[:,1]/np.sum(self.mjd_n_histogram[:,1])    
        
    def exp_decay_func(self, t, k):
        return k*np.exp(-t*k)
    
    def cum_exp_decay_func(self, t, k):
        """
        Die (kumulative) Verteilungsfunktion der Exponentialverteilung.
        Sie erlaubt die Berechnung der Wahrscheinlichkeit des Auftretens des nächsten Ereignisses im Intervall von 0 bis x, x = frames = 1.
        k = Ereignisrate.
        """
        return 1-np.exp(-t*k)
    
    def calc_k_bleach(self):
        """
        p_bleach = bleaching probability per particle and frame in the range of 0.0-1.0
        (1) The initial k value is needed to determine a k with its covarianz matrix.
        (2) Calculate the cumulative distribution function -> equals p_bleach.
        """
        self.k, self.kcov = curve_fit(self.exp_decay_func, self.mjd_n_histogram[:,0], self.mjd_n_histogram[:,1], p0 = self.init_k, method = "lm")
        self.p_bleach = self.cum_exp_decay_func(1, self.k)
        
    def calc_decay(self):
        """
        col2 = exp decay func, with bins and calculated k value.
        col3 = residues of fit -> (values - fit).
        """
        self.mjd_n_histogram [:,2] = self.exp_decay_func(self.mjd_n_histogram[:,0], self.k)
        self.mjd_n_histogram [:,3] = self.mjd_n_histogram [:,1] - self.mjd_n_histogram [:,2]
    
    def save_mjd_n_frequencies(self, directory, base_name):
        out_file_name = self.file_name[:-4] + "mjd_n_frequencies.txt"
        header = "mjd_n\t fraction\t exponential fit\t residues\t"
        np.savetxt(out_file_name, 
                   X=self.mjd_n_histogram,
                   fmt = ("%i","%.4e","%.4e","%.4e"),
                   header = header)    
    
    def save_fit_results(self):
        out_file_name = self.file_name[:-4] + "p_bleach.txt"
        self.p_bleach_results = np.zeros(3) # a np.array is like a column
        self.p_bleach_results[0], self.p_bleach_results[1], self.p_bleach_results[2] = self.p_bleach, self.k, self.kcov # fill it with values
        self.p_bleach_results = np.matrix(self.p_bleach_results) # convert it to a matrix to be able to plot results in a horizontal line
        header = "p_bleach\t k\t variance of k\t"
        np.savetxt(out_file_name,
                   X=self.p_bleach_results,
                   fmt = ("%.4e","%.4e","%.4e"),
                   header = header)
        
    def plot_mjd_frequencies(self):
        out_file_name = self.file_name[:-4] + "mjd_n_frequencies.pdf"
        fig = plt.figure()
        gridspec.GridSpec(4,4)  # set up supbplot grid 
        sp_1 = plt.subplot2grid((4,4), (0,0), colspan=4, rowspan=3)  # start left top = (0,0) = (row,column)
        sp_1.tick_params(axis='x',  # changes apply to the x-axis
                        which='both',  # both major and minor ticks are affected
                        bottom='off',  # ticks along the bottom edge are off
                        top='off',  # ticks along the top edge are off
                        labelbottom='off') # labels along the bottom edge are off
        #sp_1 = fig.add_subplot(2, 1, 1)  # (row, column, index)
        sp_1.bar(self.mjd_n_histogram [:,0], self.mjd_n_histogram [:,1], 
               align = "center",
               width = 1,
               color = "gray")  # (x, height of the bars, width of bars)
        sp_1.plot(self.mjd_n_histogram [:,0], self.mjd_n_histogram [:,2], "--")  # "b--" change colour, line style "m-" ...
        sp_1.set_title("PDF of number of data points used in MJD calculation")
        #sp_1.set_xlabel("Number of data points used in MJD calculation")
        sp_1.set_ylabel("Fraction")
        sp_2 = plt.subplot2grid((4,4), (3,0), colspan=4, rowspan=1)
        #sp_2 = fig.add_subplot(2, 1, 2) 
        sp_2.plot(self.mjd_n_histogram [:,0], self.mjd_n_histogram [:,3], "--")
        sp_2.set_ylabel("Residue")
        sp_2.set_xlabel("Number of MJDs per track")  # Number of data points used in MJD calculation
        plt.savefig(out_file_name)
        plt.show()
     
        
def main():
    p_bleach = PBleach()
    p_bleach.file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.tif.tracked.seg.txt"
    p_bleach.load_seg_file()
    p_bleach.count_mjd_n_frequencies()
    p_bleach.calc_k_bleach()
    p_bleach.calc_decay()
    p_bleach.save_mjd_n_frequencies()
    p_bleach.save_fit_results()
    p_bleach.plot_mjd_frequencies()
    
if __name__ == "__main__":
    main()
    
    