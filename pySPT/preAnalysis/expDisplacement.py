# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 13:26:16 2019

@author: Johanna Rahm, Research Group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt am Main.
"""

import datetime
import numpy as np
import matplotlib.pyplot as plt

class ExpDisplacement():
    def __init__(self):
        self.mjd = []
        self.mjd_histogram = []
        self.average_mjd = 0
        self.fig = []
        self.header = ""
        self.identifier = "identifier"  # identifier
        self.number_columns = 0
        self.significant_words = ["track_id", "mjd", "mjd_n"]
        self.sub_headers = []  # index in list = index of column in file
        self.column_order = {}
        
    def load_header(self, file_name):
        file = open(file_name)
        self.header = file.readline()  # get the header as first line
        
    def sub_header(self):
        cut_header = self.header
        while cut_header.find(self.identifier) != -1:
            # find returns the index of the first character of identifier
            slice_index = cut_header.find(self.identifier)
            sub_header = cut_header[:slice_index+len(self.identifier)]
            self.sub_headers.append(sub_header)
            cut_header = cut_header[slice_index+len(self.identifier):]
        # the last cut_header will not have idenfitier but the significant word in it, append it
        self.sub_headers.append(cut_header)
        # the first item will have only the idenfitier but not the sig word -> delete it
        self.sub_headers.pop(0)
        self.number_columns = len(self.sub_headers)
        print("sub_headers", self.sub_headers)
        
    def column_index(self):
        for sub_header in self.sub_headers:
            for word in self.significant_words:
                if word in sub_header and word not in self.column_order:
                    #append word (value) and index of sub_head (key) to dictionary
                    self.column_order[self.sub_headers.index(sub_header)] = word
        print("Dict", self.column_order)

    def load_seg_file(self, file_name):
        """
        If True Create self.mjd (numpy.ndarray) col0 = mjd, col1 = mjd_n else raise error.
        
        :param file_name: Name of the inserted file by widgetExpDisp. 
        """
        # get the key for a certain value
        if not (file_name == ""):
            mjd_index = list(self.column_order.keys())[list(self.column_order.values()).index("mjd")]
            mjd_n_index = list(self.column_order.keys())[list(self.column_order.values()).index("mjd_n")]
            print(mjd_index, mjd_n_index)
            self.mjd = np.loadtxt(file_name, usecols = (mjd_index, mjd_n_index)) # col0 = mjd, col1 = mjd_n
        else:
            print("Insert a file name.")
        
    def count_mjd_frequencies(self):
        """
        Create histogram with col0 = bins = mjd and col1 = frequencies as np.ndarray.
        
        """
        max_bin = int(np.ceil(self.mjd[:,0].max()/20)*20)  # max mjd ceiled and divisible by 20
        bin_size = int(np.ceil(self.mjd[:,0].max()/20))  # divides the bin range in sizes -> desired bin = max_bin/bin_size
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
        Exp_displacement = average displacement.
        """
        #self.mjd_frequency_max = self.mjd_histogram[:,1].max()
        #self.mjd_max_index = self.mjd_histogram[:,1].argmax()
        #self.mjd_max =  self.mjd_histogram[self.mjd_max_index, 0]
        #print("The expected displacement is %i [nm].\nThe corresponding frequency is %.4e." %(self.mjd_max, self.mjd_frequency_max))
        mjd_no_zeros = np.ma.masked_array(self.mjd[:,0], self.mjd[:,0] == 0)
        self.average_mjd = mjd_no_zeros.mean()
        print("The expected displacement is %.3f nm." %(self.average_mjd)) 

    def save_mjd_frequencies(self, directory, base_name):
        """
        Create a mjd_frequencies.txt file with exp displacement in header, col0 = mjd, col1 = frequencies.
        """
        now = datetime.datetime.now()
        year = str(now.year)
        year = year[2:]
        month = str(now.month)
        if len(month) == 1:
            month = str(0) + month
        day = str(now.day)
        out_file_name = directory + "\ " + year + month + day + "_" + base_name + "_mjd_frequencies.txt"  # Betriebssystemunabhängig?!?!?!
        #header = "The expected displacement is %i [nm].\nThe corresponding frequency is %.4e.\n" %(self.mjd_max, self.mjd_frequency_max)
        header = "The expected displacement is %.3f [nm].\n" %(self.average_mjd)
        header += "mjd [nm]\t fraction\t"
        np.savetxt(out_file_name, 
                   X=self.mjd_histogram,
                   fmt = ("%i","%.4e"),
                   header = header)
        
    def plot_mjd_frequencies(self):
        self.fig = plt.figure()
        sp = self.fig.add_subplot(1, 1, 1)  # only 1 plot
        sp.bar(self.mjd_histogram [:,0], self.mjd_histogram [:,1], 
               align = "center",
               width = 20,  # width = bin size
               color = "gray",
               edgecolor = "black")
        sp.set_title("PDF of Mean Jump Distance")
        sp.set_xlabel("Mean Jump Distance [nm]")
        sp.set_ylabel("Fraction")
        #plt.savefig(out_file_name)
        plt.show()  # print the graph
        #return fig


def main():
    exp_displacement = ExpDisplacement()
    #exp_displacement.file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.tif.tracked.seg.txt"
    exp_displacement.load_seg_file()
    exp_displacement.count_mjd_frequencies()
    exp_displacement.calc_exp_displacement()
    exp_displacement.save_mjd_frequencies()
    exp_displacement.plot_mjd_frequencies()
    
if __name__ == "__main__":
    main()
 
    