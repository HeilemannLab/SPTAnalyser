# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:36:48 2019

@author: pcoffice37
"""
import numpy as np
import matplotlib.pyplot as plt


class TrajectoryLengths():
    
    def __init__(self):
        self.file_name = ""
        self.mjd = []
        self.trajectory_lengths_histogram = []
        self.trajectory_lengths = []
        
    def load_seg_file(self):
        """
        Load file if file_name is inserted
        """
        if not (self.file_name == ""):
            self.mjd = np.loadtxt(self.file_name, usecols = (1, 2))  
        else:
            print("Insert a file name.")
            
    def calc_trajectory_lengths(self):
        print(len(self.mjd[:,0]), len(self.mjd[:,1]))
        self.trajectory_lengths = np.multiply(self.mjd[:,0], self.mjd[:,1])
        print(type(self.trajectory_lengths))
        print(self.trajectory_lengths)
            
    def count_trajectory_length_frequencies(self):
        max_bin = np.ceil(self.trajectory_lengths.max())
        bin_size = int(np.ceil(self.trajectory_lengths.max()))
        print("bin size", bin_size)
        print("max bin", max_bin)
        max_bin_index = self.trajectory_lengths.argmax()
        hist = np.histogram(self.trajectory_lengths,
                                range = (0, max_bin),
                                bins = bin_size,
                                density = True)
        print("hist", hist)
        print("size", np.size(hist[0]), len(hist[1]))
        self.trajectory_lengths_histogram = np.zeros((np.size(hist[0]),2))
        self.trajectory_lengths_histogram[:,0] = hist[1][:-1]  # bins = length of trajectory, :-1 because hist[0] is 1 shorter
        self.trajectory_lengths_histogram[:,1] = hist[0][:]  # frequencies
        print(self.trajectory_lengths_histogram)
        
    def calc_p_bleach(self):
        pass
    
    def save_trajectory_frequencies(self):
        out_file_name = self.file_name[:-4] + "trajectory_frequencies.txt"
        header = "trajectory length [nm]\t fraction"
        np.savetxt(out_file_name, 
                   X=self.trajectory_lengths_histogram,
                   fmt = ("%i","%.4e"),
                   header = header)
        
    def plot_trajectory_lengths_frequencies(self):
        out_file_name = self.file_name[:-4] + "trajectory_frequencies.pdf"
        fig = plt.figure()
        sp = fig.add_subplot(1, 1, 1)  # only 1 plot
        sp.bar(self.trajectory_lengths_histogram [:,0], self.trajectory_lengths_histogram [:,1], 
               align = "center",
               width = 1,
               color = "gray",
               edgecolor = "black")
        sp.set_title("xxx")
        sp.set_xlabel("Trajectory Length [nm]")
        sp.set_ylabel("Fraction")
        plt.savefig(out_file_name)
        plt.show()
       # self.trajectory_lengths_histogram = np.zeros([np.size(hist[0]),2])
        
       #self.trajectory_length_histogram [:,0] = hist[1][:-1]
        #self.trajectory_length_histogram [:,1] = hist[0][:]


def main():
    trajectory_lengths = TrajectoryLengths()
    trajectory_lengths.file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.tif.tracked.seg.txt"
    trajectory_lengths.load_seg_file()
    trajectory_lengths.calc_trajectory_lengths()
    trajectory_lengths.count_trajectory_length_frequencies()
    trajectory_lengths.calc_p_bleach()
    trajectory_lengths.save_trajectory_frequencies()
    trajectory_lengths.plot_trajectory_lengths_frequencies()



if __name__ == "__main__":
    main()