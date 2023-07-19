"""
@author: Johanna Rahm, Alexander Niedrig
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

P_bleach is the probability per particle and frame to bleach [0..1]. Calc k with exp decay function a*exp(-dt*k).
With k, calc cumulative distribution function 1-exp(-dt*k) = probability of event in the interval [0..1], 1 = 1 frame.
Draws no graphs to allow batch processing.
"""

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
        self.mjd_n_histogram = []  # mjd_n, frequencies, exp fit, residue
        self.p_bleach_results = []  # p_bleach, k, kv
        self.dt = 0.02  # integration time in s
        self.init_k = 0.5
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
        
    def load_seg_file(self):
        mjd_n_index = list(self.column_order.keys())[list(self.column_order.values()).index('"seg.mjd_n"')]
        if self.software == "ThunderSTORM":
            df = pd.read_csv(self.file_name)
            df_mjd_n = df.iloc[:,mjd_n_index]
            self.mjds = np.zeros(np.shape(df)[0])
            self.mjds = df_mjd_n  # create numpy array with mjd_ns     
        elif self.software == "rapidSTORM":
            self.mjds = np.loadtxt(self.file_name, usecols=(mjd_n_index))  # col0 = mjd_n

    def count_mjd_n_frequencies(self):
        """
        Create histogram with bins = mjd_n and frequencies as np.ndarray.
        Multiply the bins with camera integration time -> [s].
        col(0) = time lag, col(1) = time lag * dt, col(2) = frequencies.
        """
        max_bin = self.mjds.max()  # max mjd_n value
        bin_size = int(max_bin)  # divides the bin range in sizes -> desired bin = max_bin/bin_size
        hist = np.histogram(self.mjds, range=(0, max_bin), bins=bin_size, density=True)
        self.mjd_n_histogram = np.zeros([np.size(hist[0]), 5])
        self.mjd_n_histogram[:, 0] = hist[1][:-1]  # col0 = time lag
        np.multiply(self.mjd_n_histogram[:, 0], float(self.dt), self.mjd_n_histogram[:, 1])
        self.mjd_n_histogram[:, 2] = hist[0][:]  # col2 = frequencies
        # the first frequency is the sum of all mjd_ns, then the mjd_n = 1 is subtracted
        # new sum is build and mjd_n = 2 is subtracted ...
        decay_frequency = np.zeros(np.size(self.mjd_n_histogram[:, 2]))
        frequency_sum = np.sum(self.mjd_n_histogram[:, 2])
        decay_frequency[0] = frequency_sum
        for i in range(np.size(self.mjd_n_histogram[:, 2])-1):
            frequency_sum -= self.mjd_n_histogram[i, 2]
            decay_frequency[i+1] = frequency_sum
        self.mjd_n_histogram[:, 2] = decay_frequency/np.sum(self.mjd_n_histogram[self.ignore_points:,2])

    def normalized_mjd_ns(self):
        """
        Create normalized frequencies for histogram -> the integral equals 1.
        """
        self.mjd_n_histogram[:, 2] = self.mjd_n_histogram[:, 2]/np.sum(self.mjd_n_histogram[:, 2])
        
    def exp_decay_func(self, t, a, k):
        """
        Describing the exp decay for the histogram of track lengths.
        """
        return a*np.exp(-t*k)
    
    def cum_exp_decay_func(self, t, k):
        """
        Cumulative distribution function of the exp function. Corresponds to the probability of an event happening in
        the next time interval.
        :param t: time interval [s].
        :param k: event rate [1/s].
        """
        return 1-np.exp(-t*k)
    
    def calc_k_bleach(self):
        """
        Calculate p_bleach = bleaching probability per particle and frame [0-1].
        (1) The initial k value is needed to determine a k with its covariance matrix.
        (2) Calculate the cumulative distribution function -> equals p_bleach.
        """
        init_a = self.mjd_n_histogram[self.ignore_points:,2].max()  # always 1
        [self.a, self.k], self.kcov = curve_fit(self.exp_decay_func, self.mjd_n_histogram[self.ignore_points:, 1],
                                                self.mjd_n_histogram[self.ignore_points:, 2],
                                                p0=(init_a, self.init_k), method="lm")
        self.p_bleach = self.cum_exp_decay_func(self.dt, self.k)
        print("Results: p_bleach = %.3f, k = %.4e s\u207B\u00B9, kv = %.4e s\u207B\u00B2" %(self.p_bleach, self.k, self.kcov[1, 1]))
        
    def calc_decay(self):
        """
        col3 = exp decay func, with bins and calculated k value.
        col4 = residues of fit -> (values - fit).
        """
        self.mjd_n_histogram [self.ignore_points:, 3] = self.exp_decay_func(self.mjd_n_histogram[self.ignore_points:, 1], self.a, self.k)
        self.mjd_n_histogram [self.ignore_points:, 4] = self.mjd_n_histogram [self.ignore_points:, 2] - self.mjd_n_histogram [self.ignore_points:, 3]
        
    def plot_mjd_frequencies(self):
        x1, x2 = 0, self.mjd_n_histogram[self.ignore_points:, 1].max()  # x1 = min, x2 = max
        sp1_y1, sp1_y2 = 0, self.mjd_n_histogram[self.ignore_points:, 2].max()
        sp2_y1, sp2_y2 = self.mjd_n_histogram[self.ignore_points:, 4].min(), self.mjd_n_histogram[self.ignore_points:, 4].max()
        fig = plt.figure()
        gridspec.GridSpec(4, 4)
        sp_1 = plt.subplot2grid((4, 4), (0, 0), colspan=4, rowspan=3)  # start left top = (0,0) = (row,column)
        sp_1.tick_params(axis="x",  # changes apply to the x-axis
                        which="both",  # both major and minor ticks are affected
                        bottom=False,  # ticks along the bottom edge are off
                        top=False,  # ticks along the top edge are off
                        labelbottom=False)  # labels along the bottom edge are off
        sp_1.bar(self.mjd_n_histogram[self.ignore_points:, 1], self.mjd_n_histogram[self.ignore_points:, 2],
               align="center", width=self.dt, color="gray", label="fraction")  # (x, height of the bars, width of bars)
        sp_1.plot(self.mjd_n_histogram[self.ignore_points:, 1], self.mjd_n_histogram[self.ignore_points:, 3], "--c",
                  label="exp fit")
        sp_1.legend()
        sp_1.set_title("Amount of tracks existing after time lag")
        sp_1.set_ylabel("Fraction")
        sp_1.axis((x1, x2, sp1_y1, sp1_y2))
        sp_2 = plt.subplot2grid((4, 4), (3, 0), colspan=4, rowspan=1)
        residue_line = np.zeros(len(self.mjd_n_histogram[self.ignore_points:, 1]))
        sp_2.plot(self.mjd_n_histogram[self.ignore_points:, 1], residue_line, ":", color="0.75")
        sp_2.plot(self.mjd_n_histogram[self.ignore_points:, 1], self.mjd_n_histogram[self.ignore_points:, 4],
                  "*", color="0.5", label="residues")
        sp_2.legend()
        sp_2.set_ylabel("Residue")
        sp_2.set_xlabel("Time lag [s]")
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
        out_file_name = directory + "\\" + year + month + day + "_" + base_name + "_p_bleach" + "_histogram.txt"
        header = "frames [count]\ttime [s]\tfraction\texponential fit\tresidues\t"
        np.savetxt(out_file_name, X=self.mjd_n_histogram, fmt=("%i", "%.4e", "%.4e", "%.4e", "%.4e"), header=header)
        print("Results are saved at", directory)
        
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
        out_file_name = directory + "\\" + year + month + day + "_" + base_name + "_p_bleach.txt"
        file = open(out_file_name, 'w+')
        if not file.closed:
            file.write("# p_bleach\tk [1/s]\tvariance of k [1/s\u00b2]\tnumber of points masked\n")
            file.write("%.4e\t%.4e\t%.4e\t%d\n" %(self.p_bleach, self.k, self.kcov[1, 1], self.ignore_points))
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
        self.figure.savefig(directory + "\\" + year + month + day + "_" + base_name + "_p_bleach_histogram.pdf",
                            format="pdf", transparent=True)
