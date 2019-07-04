# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 15:40:08 2019

@author: Johanna Rahm

Research Group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt am Main.

P_bleach is the probability per particle and frame to bleach [0..1]. Calc k with exp decay funktion a*exp(-dt*k). With k, calc cumulative
distribution function 1-exp(-dt*k) = probability of event in the interval [0..1], 1 = 1 frame.
"""

#It can be estimated by taking
#the inverse of the total number of frames a particle is visible, not counting the frames in
#which it remains undetected (due to blinking)... from manual

import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import pandas as pd

class PBleach():
    def __init__(self):
        self.software = ""  # either thunderSTORM or rapidSTORM
        self.file_name = ""
        self.column_order = {}  # {0: '"track_id"', 4: '"mjd"', 6: '"mjd_n"'}
        self.mjds = []
        self.mjd_n_histogram = []  # mjd_n, frequencies, exp fit, resudie
        self.p_bleach_results = []  # p_bleach, k, kv
        self.dt = 0.02  # integration time in s
        self.init_k = 0.01
        self.p_bleach = 0.0
        self.a = 0.01
        self.k = 0.01
        self.kcov = 0
        self.ignore_points = 0  # number of points to ignore for the exp. fit
        self.figure = []
        
    def run_p_bleach(self):
        self.load_seg_file()
        self.count_mjd_n_frequencies()
        self.calc_k_bleach()
        self.calc_decay()
        self.plot_mjd_frequencies()

        
    def load_seg_file(self):
        mjd_n_index = list(self.column_order.keys())[list(self.column_order.values()).index('"seg.mjd_n"')]
        if self.software == "ThunderSTORM":
            df = pd.read_csv(self.file_name)
            df_mjd_n = df.iloc[:,mjd_n_index]
            self.mjds = np.zeros(np.shape(df)[0])
            self.mjds = df_mjd_n  # create numpy array with mjd_ns     
        elif self.software == "rapidSTORM":
            self.mjds = np.loadtxt(self.file_name, usecols = (mjd_n_index))  # col0 = mjd_n

    def count_mjd_n_frequencies(self):
        """
        Create histogram with bins = mjd_n and frequencies as np.ndarray.
        Multiply the bins with camera integration time -> [s].
        col(0) = time lag
        col(1) = time lag * dt
        col(2) = frequencies
        """
        max_bin = self.mjds.max()  # max mjd_n value
        bin_size = int(max_bin)  # divides the bin range in sizes -> desired bin = max_bin/bin_size
        hist = np.histogram(self.mjds,
                            range = (0, max_bin),
                            bins = bin_size,
                            density = True)
        self.mjd_n_histogram = np.zeros([np.size(hist[0]),5])
        self.mjd_n_histogram [:,0] = hist[1][:-1] # col0 = time lag
        np.multiply(self.mjd_n_histogram[:,0], float(self.dt), self.mjd_n_histogram [:,1])
        #self.mjd_n_histogram [:,1] = self.mjd_n_histogram[:,0]*self.dt  # col1 = s
        self.mjd_n_histogram [:,2] = hist[0][:]  # col2 = frequencies
        # The first frequency is the sum of all mjd_ns, then the mjd_n = 1 is substracted, new sum is build and mjd_n = 2 is substracted ... 
        decay_frequency = np.zeros(np.size(self.mjd_n_histogram[:,2]))
        frequency_sum = np.sum(self.mjd_n_histogram[:,2])
        decay_frequency[0] = frequency_sum
        for i in range(np.size(self.mjd_n_histogram[:,2])-1):
            frequency_sum -= self.mjd_n_histogram[i,2]
            decay_frequency[i+1] = frequency_sum
        self.mjd_n_histogram[:,2] = decay_frequency/np.sum(self.mjd_n_histogram[self.ignore_points:,2])
        #self.normalized_mjd_ns()  # normalize the histogram by the sum
        
    def normalized_mjd_ns(self):
        """
        Create normalized frequencies for histogram -> the sum equals 1.
        """
        self.mjd_n_histogram[:,2] = self.mjd_n_histogram[:,2]/np.sum(self.mjd_n_histogram[:,2])    
        
    def exp_decay_func(self, t, a, k):
        """
        Describing the exp decay for the histogram of track lengths.
        """
        return a*np.exp(-t*k)
    
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
        #init_a = self.mjd_n_histogram[0,1] #old version
        init_a = self.mjd_n_histogram[self.ignore_points:,2].max()  # always 1
        #  if one would like to neglect trajectories > valid_length for fitting a and k because most of them are 0 for one data set
        #[self.a, self.k], self.kcov = curve_fit(self.exp_decay_func, self.mjd_n_histogram[:self.valid_length,0], self.mjd_n_histogram[:self.valid_length,1], p0 = (init_a, init_k), method = "lm") #func, x, y, 
        [self.a, self.k], self.kcov = curve_fit(self.exp_decay_func, self.mjd_n_histogram[self.ignore_points:,1], self.mjd_n_histogram[self.ignore_points:,2], p0 = (init_a, self.init_k), method = "lm") #func, x, y, "lm"
        print("a&k", self.a, self.k)
        print(self.mjd_n_histogram)
        print(curve_fit(self.exp_decay_func, self.mjd_n_histogram[self.ignore_points:,1], self.mjd_n_histogram[self.ignore_points:,2], p0 = (init_a, self.init_k), method = "lm"))
        self.p_bleach = self.cum_exp_decay_func(self.dt, self.k)
        print("Results: p_bleach = %.3f, k = %.4e, kv = %.4e" %(self.p_bleach, self.k, self.kcov[1,1]))  # Output for Jupyter Notebook File
        
    def calc_decay(self):
        """
        col3 = exp decay func, with bins and calculated k value.
        col4 = residues of fit -> (values - fit).
        """
        self.mjd_n_histogram [self.ignore_points:,3] = self.exp_decay_func(self.mjd_n_histogram[self.ignore_points:,1], self.a, self.k)
        self.mjd_n_histogram [self.ignore_points:,4] = self.mjd_n_histogram [self.ignore_points:,2] - self.mjd_n_histogram [self.ignore_points:,3]
        
    def plot_mjd_frequencies(self):
        x1, x2 = 0, self.mjd_n_histogram[self.ignore_points:,1].max()  # x1 = min, x2 = max
        sp1_y1, sp1_y2 = 0, self.mjd_n_histogram[self.ignore_points:,2].max()
        sp2_y1, sp2_y2 = self.mjd_n_histogram[self.ignore_points:,4].min(), self.mjd_n_histogram[self.ignore_points:,4].max()
        fig = plt.figure()
        gridspec.GridSpec(4,4)  # set up supbplot grid 
        sp_1 = plt.subplot2grid((4,4), (0,0), colspan=4, rowspan=3)  # start left top = (0,0) = (row,column)
        sp_1.tick_params(axis='x',  # changes apply to the x-axis
                        which='both',  # both major and minor ticks are affected
                        bottom=False,  # ticks along the bottom edge are off
                        top=False,  # ticks along the top edge are off
                        labelbottom=False) # labels along the bottom edge are off
        #sp_1 = fig.add_subplot(2, 1, 1)  # (row, column, index)
        sp_1.bar(self.mjd_n_histogram [self.ignore_points:,1], self.mjd_n_histogram [self.ignore_points:,2], 
               align = "center",
               width = self.dt,
               color = "gray",
               label = "fraction")  # (x, height of the bars, width of bars)
        sp_1.plot(self.mjd_n_histogram [self.ignore_points:,1], self.mjd_n_histogram [self.ignore_points:,3], "--c", label = "exp fit")  # "b--" change colour, line style "m-" ...
        sp_1.legend()
        sp_1.set_title("Amount of tracks existing after time lag")
        #sp_1.set_xlabel("Number of data points used in MJD calculation")
        sp_1.set_ylabel("Fraction")
        sp_1.axis((x1, x2, sp1_y1, sp1_y2))
        sp_2 = plt.subplot2grid((4,4), (3,0), colspan=4, rowspan=1)
        #sp_2 = fig.add_subplot(2, 1, 2) 
        residue_line = np.zeros(len(self.mjd_n_histogram [self.ignore_points:,1]))
        sp_2.plot(self.mjd_n_histogram[self.ignore_points:,1], residue_line, ":", color = "0.75")        
        sp_2.plot(self.mjd_n_histogram [self.ignore_points:,1], self.mjd_n_histogram [self.ignore_points:,4], "*", color = "0.5", label= "residues")
        sp_2.legend()
        sp_2.set_ylabel("Residue")
        sp_2.set_xlabel("Time lag [s]")  # Number of data points used in MJD calculation
        sp_2.axis((x1, x2, sp2_y1, sp2_y2))
        plt.show() 
        self.figure = fig
        
    def save_mjd_n_frequencies(self, directory, base_name):
        """
        Output file with cols:
        col0 = mjd_n
        col1 = mjd_n*dt -> in s
        col2 = fraction
        col3 = exp decay func, with bins and calculated k value.
        col4 = residues of fit -> (values - fit).
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
        out_file_name = directory + "\ " + year + month + day + "_" + base_name + "_p_bleach" + "_histogram.txt" # System independent?
        header = "frames [count]\tduration [s]\tfraction\texponential fit\tresidues\t"
        np.savetxt(out_file_name, 
                   X=self.mjd_n_histogram,
                   fmt = ("%i","%.4e","%.4e","%.4e", "%.4e"),
                   header = header)   
        print("Results are saved.")    
        
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
        out_file_name = directory + "\ " + year + month + day + "_" + base_name + "_p_bleach.txt"
# old attempt
# =============================================================================
#         self.p_bleach_results = np.zeros(3) # a np.array is like a column
#         self.p_bleach_results[0], self.p_bleach_results[1], self.p_bleach_results[2] = self.p_bleach, self.k, self.kcov[1,1] # fill it with values
#         self.p_bleach_results = np.matrix(self.p_bleach_results) # convert it to a matrix to be able to plot results in a horizontal line
#         header = "p_bleach\t k\t variance of k\t"
#         np.savetxt(out_file_name,
#                    X=self.p_bleach_results,
#                    fmt = ("%.4e","%.4e","%.4e"),
#                    header = header)
# =============================================================================
        file = open(out_file_name, 'w')
        if not (file.closed):
            file.write("# p_bleach\tk [1/s]\tvariance of k [1/s\u00b2]\tnumber of points masked\n")
            file.write("%.4e\t%.4e\t%.4e\t%d\n" %(self.p_bleach, self.k, self.kcov[1,1], self.ignore_points))
            file.close()
        else:
            print("error: could not open file %s. Make sure the folder does exist" %(out_file_name))
            
    def save_plot(self, directory, base_name):
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        day = str(now.day)
        if len(month) == 1:
            month = str(0) + month
        if len(day) == 1:
            day = str(0) + day
        self.figure.savefig(directory + "\ " + year + month + day + "_" + base_name + "_p_bleach_histogram.pdf", format="pdf", transparent=True)

       
        
def main():
    file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.tif.tracked.seg.txt"
    p_bleach = PBleach()
    p_bleach.file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.tif.tracked.seg.txt"
    p_bleach.load_seg_file(file_name)
    p_bleach.count_mjd_n_frequencies()
    p_bleach.calc_k_bleach()
    p_bleach.calc_decay()
    p_bleach.save_mjd_n_frequencies(file_name)
    p_bleach.save_fit_results(file_name)
    p_bleach.plot_mjd_frequencies()
    p_bleach.save_plot_mjd_frequencies(file_name)
    
    
if __name__ == "__main__":
    main()
    