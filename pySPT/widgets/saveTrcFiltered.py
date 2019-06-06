# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:13:53 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

After filtering in the JNB trackStatistics, a filtered .trc file is saved for each cell.
"""

import datetime
import numpy as np

class SaveTrcFiltered():
    def __init__(self, trc_filtered, directory, folder, base_name):
        self.trc_filtered = trc_filtered  # col0 = track id, col1 = frame, col2 = x, col3 = y, col4 = placeholder, col5= intensity, col6 = seg id 
        self.directory = directory
        self.folder = folder
        self.base_name = base_name

    def continuous_numbering(self):
        """
        For HMM continuous indexing is needed. After HMM the results will be merged with the
        filtered analysis hdf5 file and the original idx are applied.
        """
        continuous_index_track = 1
        for i in range(len(self.trc_filtered)-1):
            if self.trc_filtered[i][0] == self.trc_filtered[i+1][0]:
                self.trc_filtered[i][0] = continuous_index_track
            else:
                self.trc_filtered[i][0] = continuous_index_track
                continuous_index_track += 1 
        # continuously index the seg id starting from 1
        continuous_index_seg = 1
        for i in range(len(self.trc_filtered)-1):
            if self.trc_filtered[i][6] == self.trc_filtered[i+1][6]:
                self.trc_filtered[i][6] = continuous_index_seg
            else:
                self.trc_filtered[i][6] = continuous_index_seg
                continuous_index_seg += 1 
        # the last entry can't compare the following entry because it is the last. It will get the index from the entry before
        self.trc_filtered[len(self.trc_filtered)-1][0] = self.trc_filtered[len(self.trc_filtered)-2][0]  # track id
        self.trc_filtered[len(self.trc_filtered)-1][6] = self.trc_filtered[len(self.trc_filtered)-2][6]  # seg id
                
    def adjust_base_name(self):
        """
        original base name: 190604_cell_1_MMStack_Pos0_trc_analysis
        target base name: xxx
        """
        left = len("DDMMYY_")
        right = len("_trc_analysis")
        self.base_name = self.base_name[left:]
        self.base_name = self.base_name[:-right]
        
    def save_filtered_trc_hmm(self):
        """
        Columns: track id, frame, x, y, placeholder, intensity.
        Target folder: Chosen directory for filtered dataset.
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
        #out_file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\pySPT_cell_1_MMStack_Pos0\\preAnalysis\\sorted.txt"
        out_file_name = self.directory + "\\" + self.folder + "\\" + year + month + day + "_" + self.base_name + "_trc_filtered_hmm.trc"
        header = "track_id\t frame\t x [pixel]\t y [pixel]\t placeholder\t intensity [photon]\t"
        trc_hmm = list(map(lambda row: list(row)[:6], self.trc_filtered))  # get rid of seg id
        np.savetxt(out_file_name, 
                   X=trc_hmm,
                   header=header,
                   fmt = ("%i","%i", "%.3f", "%.3f", "%i", "%.3f"))     
        
        