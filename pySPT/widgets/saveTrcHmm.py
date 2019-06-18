# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 11:24:25 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

In the trackAnalysis JNB a trc file for the hmm analysis with no tracks D < 0 and
continuous indexing is created.
"""

import datetime
import numpy as np 


class SaveTrcHmm():
    def __init__(self, trc_file, pixel_size, save_dir, raw_base_name):
        self.trc_file = trc_file
        self.pixel_size = pixel_size
        self.raw_base_name = raw_base_name
        self.save_dir = save_dir
        
    def run_save(self):
        self.ym_to_px()
        self.continuous_index()
        self.save_trc()
    
    def show_trc(self):
        """
        Testing purpose.
        """
        print(self.trc_file)
        
    def continuous_index(self):
        continuous_index_track = 1
        for i in range(len(self.trc_file)-1):
            if self.trc_file[i][0] == self.trc_file[i+1][0]:
                self.trc_file[i][0] = continuous_index_track
            else:
                self.trc_file[i][0] = continuous_index_track
                continuous_index_track += 1 
        self.trc_file[len(self.trc_file)-1][0] = self.trc_file[len(self.trc_file)-2][0]
        
    def ym_to_px(self):
        """
        The trc file has to be converted to pixel for hmm analysis.
        pixel_size [nm], x and y [ym]
        """
        self.trc_file[:,2] /= (self.pixel_size/1000)
        self.trc_file[:,3] /= (self.pixel_size/1000) 
        
    def save_trc(self):
        """
        Columns: track id, frame, x, y, placeholder, intensity.
        Target folder: pySPT\\hmm
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
        out_file_name = self.save_dir + "\\" + year + month + day + "_" + self.raw_base_name + "_trc_hmm.trc"
        header = "track_id\t frame\t x [pixel]\t y [pixel]\t placeholder\t intensity [photon]\t"
        #trc_hmm = list(map(lambda row: list(row)[:6], self.trc_filtered))  # get rid of seg id
        np.savetxt(out_file_name, 
                   X=self.trc_file,
                   header=header,
                   fmt = ("%i","%i", "%.3f", "%.3f", "%i", "%.3f"))     