# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 18:07:36 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

For the Hidden Markov Analysis a microscopy file per cell with camera integration time, 
pixel size, dynamic localization error is needed. The file will be saved in the
pySPT\hmm folder.
"""

class Microscope():
    def __init__(self, dt, pixel_size, sigma_dyn, save_dir):
        self.dt = dt
        self.pixel_size = pixel_size
        self.sigma_dyn = sigma_dyn * 1000  # sigma_dyn from cell analysis is in ym -> convert it to nm
        self.save_dir = save_dir

    def write_microscope_file(self, file_path, error=True):
        error = self.sigma_dyn if error else 0
        file = open(file_path, 'w')
        if not (file.closed):
            file.write("# SMLMS Microscope File \n")
            file.write("# pxl Size[nm] \n")
            file.write("# integration Time [s] \n")
            file.write("# localization precision [nm] \n")
            file.write("%.6e \n" %(float(self.pixel_size)))
            file.write("%.6e \n" %(float(self.dt)))
            file.write("%.6e \n" %(float(error)))

    def save_hmm_microscope(self, no_error=True):
        """
        For the HMM-analysis a microscope file for each cell is needed which contains the localization uncertainty,
        camera pixel size and integration time. This file will be saved in the pySPT preAnalysis folder per cell optionally.
        :param no_error: If true, two files are saved microscope with 0 localization precision and microscope_dyn_error
        with dynamic localization precision, if false microscope contains dynamic loc precision.
        """
        if no_error:
            out_file_name = self.save_dir + "\\" + "microscope.txt"
            self.write_microscope_file(out_file_name, error=False)
            out_file_name = self.save_dir + "\\" + "microscope_dyn_error.txt"
            self.write_microscope_file(out_file_name, error=True)

        else:
            out_file_name = self.save_dir + "\\" + "microscope.txt"
            self.write_microscope_file(out_file_name, error=True)
