# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 13:26:16 2019

@author: pcoffice37
"""
import numpy as np
import matplotlib.pyplot as plt


class ExpDisplacement():
    def __init__(self):
        self.mjd = []
        self.file_name = ""
        self.mjd_histogram = []
        #self.mjd_frequency_max = 0#
        #self.mjd_max_index = 0#
        #self.mjd_max = 0
        self.average_mjd = 0
         
    def load_seg_file(self):
        if not (self.file_name == ""):
            self.mjd = np.loadtxt(self.file_name, usecols = (1, 2)) # col1 = mjd, col2 = mjd_n
            print(self.mjd)
        else:
            print("Insert a file name.")
        
    def count_mjd_frequencies(self):
        """
        Create histogram with bins = mjd and frequencies as np.ndarray
        """
        max_bin = int(np.ceil(self.mjd[:,0].max()/20)*20)  # max mjd ceiled and divisible by 20
        bin_size = int(np.ceil(self.mjd[:,0].max()/20))  # divides the bin range in sizes -> desired bin = max_bin/bin_size
        print(self.mjd[:,0].max())
        hist = np.histogram(self.mjd[:,0],
                                        range = (0, max_bin),
                                        bins = bin_size,
                                        density = True,
                                        weights = self.mjd[:,1])
        self.mjd_histogram = np.zeros([np.size(hist[0]),2])
        self.mjd_histogram [:,0] = hist[1][:-1]  # mjd [nm]
        self.mjd_histogram [:,1] = hist[0][:]  # frequencies
        
    def calc_exp_displacement(self):
        """
        Exp_displacement = average displacement 
        """
        #self.mjd_frequency_max = self.mjd_histogram[:,1].max()
        #self.mjd_max_index = self.mjd_histogram[:,1].argmax()
        #self.mjd_max =  self.mjd_histogram[self.mjd_max_index, 0]
        #print("The expected displacement is %i [nm].\nThe corresponding frequency is %.4e." %(self.mjd_max, self.mjd_frequency_max))
        mjd_no_zeros = np.ma.masked_array(self.mjd[:,0], self.mjd[:,0] == 0)
        self.average_mjd = mjd_no_zeros.mean()

    def save_mjd_frequencies(self):
        out_file_name = self.file_name[:-4] + "mjd_frequencies.txt"
        #header = "The expected displacement is %i [nm].\nThe corresponding frequency is %.4e.\n" %(self.mjd_max, self.mjd_frequency_max)
        header = "The expected displacement is %.4e [nm].\n" %(self.average_mjd)
        header += "mjd [nm]\t fraction\t"
        np.savetxt(out_file_name, 
                   X=self.mjd_histogram,
                   fmt = ("%i","%.4e"),
                   header = header)
        
    def plot_mjd_frequencies(self):
        out_file_name = self.file_name[:-4] + "mjd_frequencies.pdf"
        fig = plt.figure()
        sp = fig.add_subplot(1, 1, 1)  # only 1 plot
        sp.bar(self.mjd_histogram [:,0], self.mjd_histogram [:,1], 
               align = "center",
               width = 20,
               color = "gray",
               edgecolor = "black")
        sp.set_title("PDF of Mean Jump Distance")
        sp.set_xlabel("Mean Jump Distance [nm]")
        sp.set_ylabel("Fraction")
        plt.savefig(out_file_name)
        plt.show()
        
    
def main():
    exp_displacement = ExpDisplacement()
    exp_displacement.file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.tif.tracked.seg.txt"
    exp_displacement.load_seg_file()
    exp_displacement.count_mjd_frequencies()
    exp_displacement.calc_exp_displacement()
    exp_displacement.save_mjd_frequencies()
    exp_displacement.plot_mjd_frequencies()
    
if __name__ == "__main__":
    main()
 
    