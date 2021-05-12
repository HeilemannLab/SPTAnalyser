"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Save filtered data set as h5 per cell in trackStatistics.ipynb.
"""

import numpy as np
import h5py
import os


class SaveFiltered():
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
        self.h5_file = h5py.File(self.trc_file_hdf5, "w")

    def create_h5_name(self, path):
        self.trc_file_hdf5 = os.path.splitext(path)[0] + ".h5"  

    def groups(self):
        self.grp00 = self.h5_file.create_group("trc")
        self.grp01 = self.h5_file.create_group("MSD")
        self.grp02 = self.h5_file.create_group("diffusion")
        self.subgrp02 = self.grp02.create_group("diffusionPlots")
        self.grp03 = self.h5_file.create_group("rossier")
        self.subgrp03 = self.grp03.create_group("rossierPlots")
        self.grp04 = self.h5_file.create_group("settings")
        self.grp05 = self.h5_file.create_group("filterInfo")
        self.grp06 = self.h5_file.create_group("statistics")
        
    def filter_info(self, filter_settings, filter_thresholds_values):
        dset = self.grp05.create_dataset("filters", (1, 1), dtype=np.dtype([("filter immobile", int),
                                                         ("filter confined", int),
                                                         ("filter free", int),
                                                         ("type determination not successful", int),
                                                         ("min trajectory length [frames]", int),
                                                        ("max trajectory length [frames]", int),
                                                        ("min diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                         ("max diffusion coefficient [\u03BCm\u00b2/s]", float)]))
        dset["filter immobile"] = filter_settings[0]
        dset["filter confined"] = filter_settings[1]
        dset["filter free"] = filter_settings[2]
        dset["type determination not successful"] = filter_settings[3]
        dset["min trajectory length [frames]"] = filter_thresholds_values[0]
        dset["max trajectory length [frames]"] = filter_thresholds_values[1]
        dset["min diffusion coefficient [\u03BCm\u00b2/s]"] = filter_thresholds_values[2]
        dset["max diffusion coefficient [\u03BCm\u00b2/s]"] = filter_thresholds_values[3]
        
    def statistics_4(self, type_percentages, trajectories_included, trajectories_excluded, mean_Ds, mean_dDs,
                     mean_lengths, mean_dlengths):
        dset = self.grp06.create_dataset("statistics_4", (1, 1), dtype=np.dtype([("immobile [%]", float),
                                                         ("confined [%]", float),
                                                         ("free [%]", float),
                                                         ("no type [%]", float),
                                                         ("trajectories included", int),
                                                         ("trajectories excluded", int),
                                                         ("mean D immobile [\u03BCm\u00b2/s]", float),
                                                         ("mean D confined [\u03BCm\u00b2/s]", float),
                                                         ("mean D free [\u03BCm\u00b2/s]", float),
                                                         ("mean D no type [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D immobile [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D confined [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D free [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D no type [\u03BCm\u00b2/s]", float),
                                                         ("mean length immobile [frames]", float),
                                                         ("mean length confined [frames]", float),
                                                         ("mean length free [frames]", float),
                                                         ("mean length no type [frames]", float),
                                                         ("\u0394 mean length immobile [frames]", float),
                                                         ("\u0394 mean length confined [frames]", float),
                                                         ("\u0394 mean length free [frames]", float),
                                                         ("\u0394 mean length no type [frames]", float)]))

        dset["immobile [%]"] = type_percentages[0]
        dset["confined [%]"] = type_percentages[1]
        dset["free [%]"] = type_percentages[2]
        dset["no type [%]"] = type_percentages[3]
        dset["trajectories included"] = trajectories_included
        dset["trajectories excluded"] = trajectories_excluded
        dset["mean D immobile [\u03BCm\u00b2/s]"] = mean_Ds[0]
        dset["mean D confined [\u03BCm\u00b2/s]"] = mean_Ds[1]
        dset["mean D free [\u03BCm\u00b2/s]"] = mean_Ds[2]
        dset["mean D no type [\u03BCm\u00b2/s]"] = mean_Ds[3]
        dset["\u0394 mean D immobile [\u03BCm\u00b2/s]"] = mean_dDs[0]
        dset["\u0394 mean D confined [\u03BCm\u00b2/s]"] = mean_dDs[1]
        dset["\u0394 mean D free [\u03BCm\u00b2/s]"] = mean_dDs[2]
        dset["\u0394 mean D no type [\u03BCm\u00b2/s]"] = mean_dDs[3]
        dset["mean length immobile [frames]"] = mean_lengths[0]
        dset["mean length confined [frames]"] = mean_lengths[1]
        dset["mean length free [frames]"] = mean_lengths[2]
        dset["mean length no type [frames]"] = mean_lengths[3]
        dset["\u0394 mean length immobile [frames]"] = mean_dlengths[0]
        dset["\u0394 mean length confined [frames]"] = mean_dlengths[1]
        dset["\u0394 mean length free [frames]"] = mean_dlengths[2]
        dset["\u0394 mean length no type [frames]"] = mean_dlengths[3]

    def statistics_3(self, type_percentages, trajectories_included, trajectories_excluded, mean_Ds, mean_dDs,
                     mean_lengths, mean_dlengths):
        dset = self.grp06.create_dataset("statistics_3", (1, 1), dtype=np.dtype([("immobile+notype [%]", float),
                                                         ("confined [%]", float),
                                                         ("free [%]", float),
                                                         ("trajectories included", int),
                                                         ("trajectories excluded", int),
                                                         ("mean D immobile+notype [\u03BCm\u00b2/s]", float),
                                                         ("mean D confined [\u03BCm\u00b2/s]", float),
                                                         ("mean D free [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D immobile+notype [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D confined [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D free [\u03BCm\u00b2/s]", float),
                                                         ("mean length immobile+notype [frames]", float),
                                                         ("mean length confined [frames]", float),
                                                         ("mean length free [frames]", float),
                                                         ("\u0394 mean length immobile+notype [frames]", float),
                                                         ("\u0394 mean length confined [frames]", float),
                                                         ("\u0394 mean length free [frames]", float)]))

        dset["immobile+notype [%]"] = type_percentages[4]
        dset["confined [%]"] = type_percentages[1]
        dset["free [%]"] = type_percentages[2]
        dset["trajectories included"] = trajectories_included
        dset["trajectories excluded"] = trajectories_excluded
        dset["mean D immobile+notype [\u03BCm\u00b2/s]"] = mean_Ds[4]
        dset["mean D confined [\u03BCm\u00b2/s]"] = mean_Ds[1]
        dset["mean D free [\u03BCm\u00b2/s]"] = mean_Ds[2]
        dset["\u0394 mean D immobile+notype [\u03BCm\u00b2/s]"] = mean_dDs[4]
        dset["\u0394 mean D confined [\u03BCm\u00b2/s]"] = mean_dDs[1]
        dset["\u0394 mean D free [\u03BCm\u00b2/s]"] = mean_dDs[2]
        dset["mean length immobile+notype [frames]"] = mean_lengths[4]
        dset["mean length confined [frames]"] = mean_lengths[1]
        dset["mean length free [frames]"] = mean_lengths[2]
        dset["\u0394 mean length immobile+notype [frames]"] = mean_dlengths[4]
        dset["\u0394 mean length confined [frames]"] = mean_dlengths[1]
        dset["\u0394 mean length free [frames]"] = mean_dlengths[2]
        
    def data_settings(self, dt, pixelsize, pixelamount, cell_size, tau_threshold, min_length_type, fit_area, dof, D_min,
                      dloc_dyn_type, seg_bool, min_length_hmm, dloc_dyn_hmm):
        dset = self.grp04.create_dataset("settings", (1, 1), dtype=np.dtype([("dt [s]", float),
                                                      ("pixelsize [nm]", int),
                                                      ("pixel amount", int),
                                                      ("cell size [\u03BCm\u00b2]", float),
                                                      ("tau threshold [s]", float),
                                                      ("min trajectory length type", float),
                                                      ("fit area", float),
                                                      ("dof", int),
                                                      ("D min [\u03BCm\u00b2/s]", float),
                                                      ("\u0394 loc dyn type [\u03BCm]", float),
                                                      ("track id", int),
                                                      ("seg id", int),
                                                      ("min trajectory length hmm", int),
                                                      ("\u0394 loc dyn hmm [\u03BCm]", float)]))
        dset["dt [s]"] = dt
        dset["pixelsize [nm]"] = pixelsize
        dset["pixel amount"] = pixelamount
        dset["cell size [\u03BCm\u00b2]"] = cell_size
        dset["tau threshold [s]"] = tau_threshold
        dset["min trajectory length type"] = min_length_type
        dset["fit area"] = fit_area
        dset["dof"] = dof
        dset["D min [\u03BCm\u00b2/s]"] = D_min
        dset["\u0394 loc dyn type [\u03BCm]"] = dloc_dyn_type
        dset["track id"] = not seg_bool
        dset["seg id"] = seg_bool
        dset["min trajectory length hmm"] = min_length_hmm
        dset["\u0394 loc dyn hmm [\u03BCm]"] = dloc_dyn_hmm
        
    def data_diffusion_info(self, number, trajectory_id, diffusion_coeff, ddiffusion_coeff, MSD_0, chi2, length):
        dset = self.grp02.create_dataset("diffusionInfos", (number,), dtype=np.dtype([("trajectory id", int),
                                                             ("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                             ("\u0394 diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                             ("MSD(0) [\u03BCm\u00b2]", float),
                                                             ("chi\u00b2 [\u03BCm\u2074]", float),
                                                             ("length [frames]", int)]))
        dset["trajectory id"] = trajectory_id
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion_coeff
        dset["\u0394 diffusion coefficient [\u03BCm\u00b2/s]"] = ddiffusion_coeff
        dset["MSD(0) [\u03BCm\u00b2]"] = MSD_0
        dset["chi\u00b2 [\u03BCm\u2074]"] = chi2
        dset["length [frames]"] = length
        
    def data_diffusion_plots(self, trajectory_number, dt, MSD, fit, residues, points_fit_D):
        # for correct order of numbering, fill the trajectory number with 000 -> 0001, 0002 ... 
        trajectory_number = str(int(trajectory_number))
        while len(trajectory_number) < 4:
            trajectory_number = "0" + trajectory_number        
        dset = self.subgrp02.create_dataset("diffusionPlot{}".format(trajectory_number), (points_fit_D,),
                                            dtype=np.dtype([("dt [s]", float),
                                            ("MSD [\u03BCm\u00b2]", float),
                                            ("linear fit [\u03BCm\u00b2]", float),
                                            ("residues [\u03BCm\u00b2]", float)]))
        dset["dt [s]"] = dt
        dset["MSD [\u03BCm\u00b2]"] = MSD
        dset["linear fit [\u03BCm\u00b2]"] = fit
        dset["residues [\u03BCm\u00b2]"] = residues
        
    def data_rossier_plots(self, trajectory_number, dt, MSD, fit, residues):
        # for correct order of numbering, fill the trajectory number with 000 -> 0001, 0002 ... 
        trajectory_number = str(int(trajectory_number))
        while len(trajectory_number) < 4:
            trajectory_number = "0" + trajectory_number        
        dset = self.subgrp03.create_dataset("rossierPlot{}".format(trajectory_number), (np.shape(dt)[0],),
                                            dtype=np.dtype([("dt [s]", float),
                                            ("MSD [\u03BCm\u00b2]", float),
                                            ("rossier fit [\u03BCm\u00b2]", float),
                                            ("residues [\u03BCm\u00b2]", float)]))
        dset["dt [s]"] = dt
        dset["MSD [\u03BCm\u00b2]"] = MSD
        dset["rossier fit [\u03BCm\u00b2]"] = fit
        dset["residues [\u03BCm\u00b2]"] = residues
        
    def data_rossier_info(self, number, trajectory_id, type_immobile, type_confined, type_free, analyse_success, tau,
                          dtau, r, dr, dconfined, ddconfined, chi2):
        dset = self.grp03.create_dataset("rossierStatistics", (number,), dtype=np.dtype([("trajectory id", int),
                                                           ("type immobile", int),
                                                           ("type confined", int),("type free", int),
                                                           ("analyse succesful", int),
                                                           ("tau [s]", float),
                                                           ("\u0394 tau [s]", float),
                                                           ("r [\u03BCm]", float),
                                                           ("\u0394 r [\u03BCm]", float),
                                                           ("D confined [\u03BCm\u00b2/s]", float),
                                                           ("\u0394 D confined [\u03BCm\u00b2/s]", float),
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
        dset["D confined [\u03BCm\u00b2/s]"] = dconfined
        dset["\u0394 D confined [\u03BCm\u00b2/s]"] = ddconfined
        dset["chi\u00b2 [\u03BCm\u2074]"] = chi2
        self.h5_file.close()

    def trc_type(self, shape, track_id, frame, x, y, placeholder, intensity, seg_id):
        dset = self.grp00.create_dataset("trcType", (shape[0],), dtype=np.dtype([("track id", int),
                                                      ("frame", int),
                                                      ("x position [\u03BCm]", float),
                                                      ("y position [\u03BCm]", float),
                                                      ("placeholder", int),
                                                      ("intensity [photon]", float),
                                                      ("seg id", int)]))
        dset["track id"] = track_id
        dset["frame"] = frame
        dset["x position [\u03BCm]"] = x
        dset["y position [\u03BCm]"] = y
        dset["placeholder"] = placeholder
        dset["intensity [photon]"] = intensity
        dset["seg id"] = seg_id
        
    def trc_hmm(self, shape, track_id, frame, x, y, placeholder, intensity):
        dset = self.grp00.create_dataset("trcHmm", (shape[0],), dtype=np.dtype([("track id", int),
                                                      ("frame", int),
                                                      ("x position [\u03BCm]", float),
                                                      ("y position [\u03BCm]", float),
                                                      ("placeholder", int),
                                                      ("intensity [photon]", float)]))
        dset["track id"] = track_id
        dset["frame"] = frame
        dset["x position [\u03BCm]"] = x
        dset["y position [\u03BCm]"] = y
        dset["placeholder"] = placeholder
        dset["intensity [photon]"] = intensity
        
    def msd(self, trajectory_number, dt, MSD):
        trajectory_number = str(int(trajectory_number))
        while len(trajectory_number) < 4:
            trajectory_number = "0" + trajectory_number     
        dset = self.grp01.create_dataset("Trajectory{}".format(trajectory_number), (np.shape(dt)[0],),
                                         dtype=np.dtype([("dt [s]", float), ("MSD [\u03BCm\u00b2]", float)]))
        dset["dt [s]"] = dt
        dset["MSD [\u03BCm\u00b2]"] = MSD    
