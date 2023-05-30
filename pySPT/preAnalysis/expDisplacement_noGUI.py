"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Determine how far a particle moves between frames in the x/y-plane [nm].
Fraction vs mjd [nm] weighted by mjd_n -> determine the average as initial exp_displacement value.
Based on segments.
"""

import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class ExpDisplacement():
    def __init__(self):
        self.software = ""  # either rapidSTORM or thunderSTORM
        self.file_name = ""
        self.column_order = {}  # {0: '"track_id"', 4: '"mjd"', 6: '"mjd_n"'}
        self.filter_immob = True  # if true ignore values from immobile labeled trajectories in swift.
        self.mjd = []  # mjd/mjd_n
        self.mjd_histogram = []
        self.average_mjd = 0
        self.max_mjd = 0
        self.fig = []

    def load_seg_file(self):
        """
        Create self.mjd (numpy.ndarray) col0 = mjd, col1 = mjd_n.
        """
        # get the key for a certain value
        mjd_index = list(self.column_order.keys())[list(self.column_order.values()).index('"seg.mjd"')]
        mjd_n_index = list(self.column_order.keys())[list(self.column_order.values()).index('"seg.mjd_n"')]
        if self.software == "ThunderSTORM":
            df = pd.read_csv(self.file_name)
            if self.filter_immob:
                df = df[df["seg.motion"] == "diffusion"]
            df_mjd = df.iloc[:, mjd_index]  # get csv columns by index
            df_mjd_n = df.iloc[:, mjd_n_index]
            self.mjd = np.zeros([np.shape(df)[0], 2])
            self.mjd[:, 0] = df_mjd  # create numpy array with col0 = mjd and col1 = mjd_n
            self.mjd[:, 1] = df_mjd_n
        elif self.software == "rapidSTORM":
            self.mjd = np.loadtxt(self.file_name, usecols=(mjd_index, mjd_n_index))  # col0 = mjd, col1 = mjd_n

    def count_mjd_frequencies(self):
        """
        Create histogram with col0 = bins = mjd and col1 = frequencies as np.ndarray.
        """
        max_bin = int(np.ceil(self.mjd[:, 0].max()/20)*20)  # max mjd ceiled and divisible by 20
        bin_size = int(np.ceil(self.mjd[:, 0].max()/20))  # -> desired bin = max_bin/bin_size
        hist = np.histogram(self.mjd[:, 0], range=(0, max_bin), bins=bin_size, density=True, weights=self.mjd[:, 1])
        self.mjd_histogram = np.zeros([np.size(hist[0]), 2])
        self.mjd_histogram[:, 0] = hist[1][:-1]  # mjd [nm]
        self.mjd_histogram[:, 1] = hist[0][:]  # frequencies
        
    def calc_exp_displacement(self):
        """
        Calculate the expected displacement (=average) and max displacement.
        """
        mjd_no_zeros = np.ma.masked_array(self.mjd[:, 0], self.mjd[:, 0] == 0)
        mjd_n_no_zeros = np.ma.masked_array(self.mjd[:, 1], self.mjd[:, 1] == 0)
        self.average_mjd = np.sum(np.multiply(mjd_no_zeros, mjd_n_no_zeros))/np.sum(mjd_n_no_zeros)
        self.max_mjd = self.mjd[:, 0].max()
        print("The expected displacement is %.3f nm." %(self.average_mjd)) 
        print("The max mean jump distance is %.3f nm." %(self.max_mjd))
        
    def plot_mjd_frequencies(self):
        self.fig = plt.figure()
        sp = self.fig.add_subplot(1, 1, 1)
        sp.bar(self.mjd_histogram[:, 0], self.mjd_histogram[:, 1], align="center", width=20, color="gray",
               edgecolor="black", label="fraction")
        sp.set_title("Histogram of Mean Jump Distance")
        sp.set_xlabel("Mean Jump Distance [nm]")
        sp.set_ylabel("Fraction")
        sp.legend()
        plt.show()

    def save_mjd_frequencies(self, directory, base_name):
        """
        Create a *.txt file with histogram information, col0 = mjd, col1 = frequencies, weighted by mjd_n.
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
        out_file_name = directory + "\\" + year + month + day + "_" + base_name + "_exp_displacement" + "_histogram.txt"
        header = "mjd [nm]\tfraction\t"
        np.savetxt(out_file_name, X=self.mjd_histogram, fmt=("%i", "%.4e"), header=header)
        
    def save_exp_disp(self, directory, base_name):
        """
        Create a *.txt file with exp_displacement and max_displacement.
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
        out_file_name = directory + "\\" + year + month + day + "_" + base_name + "_exp_displacement.txt"
        file = open(out_file_name, 'w+')
        file.write("exp_displacement [nm]\tmax_mjd [nm]\n")
        file.write("%.3f\t%.3f\n" %(self.average_mjd, self.max_mjd))
        file.close()
        
    def run_exp_displacement(self):
        self.load_seg_file()
        self.count_mjd_frequencies()
        self.calc_exp_displacement()
                
    def save_exp_displacement(self, directory, base_name, save_fig):
        self.save_exp_disp(directory, base_name)
        self.save_mjd_frequencies(directory, base_name)
        if save_fig:
            self.run_save_plot(directory, base_name)
        print("Results are saved at", directory)
        
    def run_save_plot(self, directory, base_name):
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        day = str(now.day)
        if len(month) == 1:
            month = str(0) + month
        if len(day) == 1:
            day = str(0) + day
        self.fig.savefig(directory + "\\" + year + month + day + "_" + base_name + "_exp_displacement_histogram.pdf",
                         format="pdf", transparent=True)
