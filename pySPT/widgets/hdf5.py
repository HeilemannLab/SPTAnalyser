# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 15:01:01 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Class for creating a .h5 file for one cell.
"""

import numpy as np
import h5py
import os

class Hdf5():
    # one cell
    def __init__(self):
        self.h5_file = []  # h5 file
        self.grp00 = []  # groups for structure
        self.grp01 = []
        self.grp02 = []
        self.grp03 = []
        self.grp04 = []
        self.trc_file_hdf5 = ""  # path of file with .h5 ending
        
    def create_h5(self, path):
        self.create_h5_name(path)
        self.create_h5_file()
        self.groups()
        
    def create_h5_file(self):
        self.h5_file = h5py.File(self.trc_file_hdf5, "w")  # w- or x = Create file, fail if exists
                        
    def create_h5_name(self, path):     
        # splitext -> tupel with path split from .* ending. It splits at the last dot in name.
        self.trc_file_hdf5 = os.path.splitext(path)[0] + ".h5"       

    def groups(self):
        self.grp00 = self.h5_file.create_group("trc")
        self.grp01 = self.h5_file.create_group("MSD")
        self.grp02 = self.h5_file.create_group("diffusion")
        self.subgrp02 = self.grp02.create_group("diffusionPlots")
        self.grp03 = self.h5_file.create_group("rossier")
        self.subgrp03 = self.grp03.create_group("rossierPlots")
        self.grp04 = self.h5_file.create_group("settings")
        
    def data_settings(self, dt, pixelsize, pixelamount, cell_size, tau_threshold, tau_min_length, fit_area, dof, dloc_dyn):
        dset = self.grp04.create_dataset("settings", (1,1), dtype = np.dtype([("dt [s]", float),
                                                      ("pixelsize [nm]", int),
                                                      ("pixel amount", int),
                                                      ("cell size [\u03BCm\u00b2]", float),
                                                      ("tau threshold [s]", float),
                                                      ("tau min trajectory length", float),
                                                      ("fit area", float),
                                                      ("dof", int),
                                                      ("\u0394 loc dyn [\u03BCm\u00b2/s]", float)]))
        dset["dt [s]"] = dt
        dset["pixelsize [nm]"] = pixelsize
        dset["pixel amount"] = pixelamount
        dset["cell size [\u03BCm\u00b2]"] = cell_size
        dset["tau threshold [s]"] = tau_threshold
        dset["tau min trajectory length"] = tau_min_length
        dset["fit area"] = fit_area
        dset["dof"] = dof
        dset["\u0394 loc dyn [\u03BCm\u00b2/s]"] = dloc_dyn
        
    def data_diffusion_info(self, number, trajectory_id, diffusion_coeff, ddiffusion_coeff, MSD_0, chi2, length):
        dset = self.grp02.create_dataset("diffusionInfos", (number,), dtype = np.dtype([("trajectory id", int),
                                                             ("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                             ("\u0394 diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                             ("MSD(0) [\u03BCm\u00b2]", float),
                                                             ("chi\u00b2 [\u03BCm\u2074]", float),
                                                             ("length [nm]", int)]))
        dset["trajectory id"] = trajectory_id
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion_coeff
        dset["\u0394 diffusion coefficient [\u03BCm\u00b2/s]"] = ddiffusion_coeff
        dset["MSD(0) [\u03BCm\u00b2]"] = MSD_0
        dset["chi\u00b2 [\u03BCm\u2074]"] = chi2
        dset["length [nm]"] = length
        
    def data_diffusion_plots(self, trajectory_number, dt, MSD, fit, residues):  #"diffusion plot {}".format(str(trajectory_number))
        # for correct order of numbering, fill the trajectory number with 000 -> 0001, 0002 ... 
        trajectory_number = str(int(trajectory_number))
        while len(trajectory_number) < 4:
            trajectory_number = "0" + trajectory_number        
        dset = self.subgrp02.create_dataset("diffusionPlot{}".format(trajectory_number), (4,), dtype = np.dtype([("dt [s]", float),
                                            ("MSD [\u03BCm\u00b2]", float),
                                            ("linear fit [\u03BCm\u00b2]", float),
                                            ("residues [\u03BCm\u00b2]", float)]))
        dset["dt [s]"] = dt
        dset["MSD [\u03BCm\u00b2]"] = MSD
        dset["linear fit [\u03BCm\u00b2]"] = fit
        dset["residues [\u03BCm\u00b2]"] = residues
        
    def data_rossier_plots(self, trajectory_number, dt, MSD, fit, residues):  #"diffusion plot {}".format(str(trajectory_number))
        # for correct order of numbering, fill the trajectory number with 000 -> 0001, 0002 ... 
        trajectory_number = str(int(trajectory_number))
        while len(trajectory_number) < 4:
            trajectory_number = "0" + trajectory_number        
        dset = self.subgrp03.create_dataset("rossierPlot{}".format(trajectory_number), (np.shape(dt)[0],), dtype = np.dtype([("dt [s]", float),
                                            ("MSD [\u03BCm\u00b2]", float),
                                            ("rossier fit [\u03BCm\u00b2]", float),
                                            ("residues [\u03BCm\u00b2]", float)]))
        dset["dt [s]"] = dt
        dset["MSD [\u03BCm\u00b2]"] = MSD
        dset["rossier fit [\u03BCm\u00b2]"] = fit
        dset["residues [\u03BCm\u00b2]"] = residues
        
    def data_rossier_info(self, number, trajectory_id, type_immobile, type_confined, type_free, analyse_success, tau, dtau, r, dr, dconfined, chi2):
        dset = self.grp03.create_dataset("rossierStatistics", (number,), dtype = np.dtype([("trajectory id", int),
                                                           ("type immobile", int),
                                                           ("type confined", int),("type free", int),
                                                           ("analyse succesful", int),
                                                           ("tau [s]", float),
                                                           ("\u0394 tau [s]", float),
                                                           ("r [\u03BCm]", float),
                                                           ("\u0394 r [\u03BCm]", float),
                                                           ("Diffusion confined [\u03BCm\u00b2/s]", float),
                                                           ("chi\u00b2 [\u03BCm\u2074]", float)]))
        dset["trajectory id"] = trajectory_id
        dset["type immobile"] = type_immobile
        dset["type confined"] = type_confined
        dset["type free"] = type_free
        dset["analyse succesful"] = analyse_success
        dset["tau [s]"] = tau
        dset["\u0394 tau [s]"] = dtau
        dset["r [\u03BCm]"] = r
        dset["\u0394 r [\u03BCm]"] = dr
        dset["Diffusion confined [\u03BCm\u00b2/s]"] = dconfined
        dset["chi\u00b2 [\u03BCm\u2074]"] = chi2
        self.h5_file.close()

        
    def trc(self, shape, trajectory_id, frame, x, y, placeholder, intensity):
        dset = self.grp00.create_dataset("trcFile", (shape[0],), dtype = np.dtype([("trajectory id", int),
                                                      ("frame", int),
                                                      ("x position [\u03BCm]", float),
                                                      ("y position [\u03BCm]", float),
                                                      ("placeholder", int),
                                                      ("intensity [A/D count]", float)]))
        dset["trajectory id"] = trajectory_id
        dset["frame"] = frame
        dset["x position [\u03BCm]"] = x
        dset["y position [\u03BCm]"] = y
        dset["placeholder"] = placeholder
        dset["intensity [A/D count]"] = intensity
        
    def msd(self, trajectory_number, dt, MSD):
        trajectory_number = str(int(trajectory_number))
        while len(trajectory_number) < 4:
            trajectory_number = "0" + trajectory_number     
        dset = self.grp01.create_dataset("Trajectory{}".format(trajectory_number), (np.shape(dt)[0],), dtype = np.dtype([("dt [s]", float),("MSD [\u03BCm\u00b2]", float)]))
        dset["dt [s]"] = dt
        dset["MSD [\u03BCm\u00b2]"] = MSD    
        
    def test(self, x, y):
        dset = self.grp00.create_dataset("test", (np.shape(x)[0],), dtype = np.dtype([("x", int), ("y", int)]))
        dset["x"] = x
        dset["y"] = y
        self.h5_file.close()
        
    
        
                       
def main():
    pass


if __name__ == "__main__":
    main()
