# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 09:25:47 2019

@author: pcoffice37
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class PBleach():
    def __init__(self):
        self.file_name = ""
        self.mjds = []
        self.mjd_n_histogram = []
        self.init_k = 0.01
        self.k = 0.01
        self.kcov = 0
        self.p_bleach = 0.0
        
    def load_seg_file(self):
        if not (self.file_name == ""):
            self.mjds = np.loadtxt(self.file_name, usecols = (1, 2)) # col1 = mjd, col2 = mjd_n
            print(self.mjds)
        else:
            print("Insert a file name.")
            
    def normalized_mjds(self):
        self.mjd_n_histogram [:,1] = self.mjd_n_histogram [:,1]/np.sum(self.mjd_n_histogram [:,1])
    
    def count_mjd_n_frequencies(self):
        """
        Create histogram with bins = mjd_n and frequencies as np.ndarray
        """
        
        max_bin = self.mjds[:,1].max()  # max mjd_n value
        bin_size = int(max_bin)  # divides the bin range in sizes -> desired bin = max_bin/bin_size
        print("bin size", bin_size)
        print("max_bin", max_bin)
        print(self.mjds[:,1].max())
        hist = np.histogram(self.mjds[:,1],
                                        range = (0, max_bin),
                                        bins = bin_size,
                                        density = True)
        self.mjd_n_histogram = np.zeros([np.size(hist[0]),4])
        self.mjd_n_histogram [:,0] = hist[1][:-1]  # bins
        self.mjd_n_histogram [:,1] = hist[0][:]  # frequencies
        self.normalized_mjds()
        
        print("hist:", hist)
        print("mjd histogram", self.mjd_n_histogram)  

    
    def exp_decay_func(self, t, k):
        return k*np.exp(-t*k)
    
    def cum_exp_decay_func(self, t, k):
        return 1-np.exp(-t*k)
    
    def calc_k_bleach(self):
        """
        p_bleach = bleaching probability per particle and frame in the range of 0.0-1.0
        """
        self.k, self.kcov = curve_fit(self.exp_decay_func, self.mjd_n_histogram[:,0], self.mjd_n_histogram[:,1], p0 = self.init_k, method = "lm")
        print("po, pcov", self.k, self.kcov)
        self.p_bleach = self.cum_exp_decay_func(1, self.k)
        print("p_bleach", self.p_bleach)
        
    def calc_decay(self):
         self.mjd_n_histogram [:,2] = self.exp_decay_func(self.mjd_n_histogram [:,0], self.k)
         self.mjd_n_histogram [:,3] = self.mjd_n_histogram [:,1] - self.mjd_n_histogram [:,2]
    
    def save_mjd_n_frequencies(self):
        out_file_name = self.file_name[:-4] + "mjd_n_frequencies.txt"
#        header = "The expected displacement is %i [nm].\nThe corresponding frequency is %.4e.\n" %(self.mjd_max, self.mjd_frequency_max)
#        header += "mjd [nm]\t fraction\t"
        header = "mjd_n\t fraction\t"
        np.savetxt(out_file_name, 
                   X=self.mjd_n_histogram,
                   fmt = ("%i","%.4e","%.4e","%.4e"),
                   header = header)    
        
    def plot_mjd_frequencies(self):
        out_file_name = self.file_name[:-4] + "mjd_n_frequencies.pdf"
        fig = plt.figure()
        sp_1 = fig.add_subplot(2, 1, 1)  # only 1 plot
        sp_1.bar(self.mjd_n_histogram [:,0], self.mjd_n_histogram [:,1], 
               align = "center",
               width = 1,
               color = "gray")
        sp_1.plot(self.mjd_n_histogram [:,0], self.mjd_n_histogram [:,2],
                "--",
                "b")
        sp_1.set_title("PDF of number of data points used in MJD calculation")
        sp_1.set_xlabel("Number of data points used in MJD calculation")
        sp_1.set_ylabel("Fraction")
        sp_2 = fig.add_subplot(2, 1, 2)  # only 1 plot
        sp_2.plot(self.mjd_n_histogram [:,0], self.mjd_n_histogram [:,3],
                "--",
                "b")
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
    p_bleach.plot_mjd_frequencies()
    
if __name__ == "__main__":
    main()