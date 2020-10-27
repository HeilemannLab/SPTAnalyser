"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Save statistics.h5 file in TrackStatistics.ipynb.
"""

import h5py
import os
import numpy as np


class SaveStatistics():
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
        self.grp00 = self.h5_file.create_group("cellInfo")
        self.grp01 = self.h5_file.create_group("cellCounts")
        self.grp04 = self.h5_file.create_group("filterInfo")
        self.grp05 = self.h5_file.create_group("diffusionHistogram")
        self.grp06 = self.h5_file.create_group("statistics")
        self.grp07 = self.h5_file.create_group("MSDPlots")

    def groups_bg(self):
        self.grp02 = self.h5_file.create_group("backgroundInfo")
        self.grp03 = self.h5_file.create_group("backgroundCounts")
        
    def cells(self, data):
        my_datatype = np.dtype([("cell name", h5py.special_dtype(vlen=str)),
                                ("cell size [\u03BCm\u00b2]", float)])
        dset = self.grp00.create_dataset("cells", (np.shape(data)[0],), dtype=my_datatype)
        data_array = np.array(data, dtype=my_datatype)
        dset[...] = data_array
        self.h5_file.close()
        
    def backgrounds(self, data):
        my_datatype = np.dtype([("background name", h5py.special_dtype(vlen=str)),
                                ("background size [\u03BCm\u00b2]", float)])
        dset = self.grp02.create_dataset("backgrounds", (np.shape(data)[0],), dtype=my_datatype)
        data_array = np.array(data, dtype=my_datatype)
        dset[...] = data_array

    def cell_counts(self, cell_name, number, diffusion_coeff, counts):
        dset = self.grp01.create_dataset(cell_name, (number,),
                                         dtype=np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                         ("counts / area [1/\u03BCm\u00b2]", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion_coeff
        dset["counts / area [1/\u03BCm\u00b2]"] = counts
        
    def bg_counts(self, bg_name, number, diffusion_coeff, counts):
        dset = self.grp03.create_dataset(bg_name, (number,),
                                         dtype=np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                         ("counts / area [1/\u03BCm\u00b2]", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion_coeff
        dset["counts / area [1/\u03BCm\u00b2]"] = counts
        
    def filter_info(self, filter_settings, filter_thresholds_values):
        dset = self.grp04.create_dataset("filters", (1, 1),
                                         dtype=np.dtype([("filter immobile", int),
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
        
    def statistics_4(self, type_percentages, dtype_percentages, trajectories_included, trajectories_excluded,
                     mean_Ds, mean_dDs, mean_lengths, mean_dlengths):
        dset = self.grp06.create_dataset("statistics_4", (1, 1),
                                         dtype=np.dtype([("immobile [%]", float),
                                                         ("confined [%]", float),
                                                         ("free [%]", float),
                                                         ("no type [%]", float),
                                                         ("\u0394 immobile [%]", float),
                                                         ("\u0394 confined [%]", float),
                                                         ("\u0394 free [%]", float),
                                                         ("\u0394 no type [%]", float),
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
        dset["\u0394 immobile [%]"] = dtype_percentages[0]
        dset["\u0394 confined [%]"] = dtype_percentages[1]
        dset["\u0394 free [%]"] = dtype_percentages[2]
        dset["\u0394 no type [%]"] = dtype_percentages[3]
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

    def statistics_3(self, type_percentages, dtype_percentages, trajectories_included, trajectories_excluded,
                     mean_Ds, mean_dDs, mean_lengths, mean_dlengths):
        dset = self.grp06.create_dataset("statistics_3", (1, 1),
                                         dtype=np.dtype([("immobile+no type [%]", float),
                                                         ("confined [%]", float),
                                                         ("free [%]", float),
                                                         ("\u0394 immobile+no type [%]", float),
                                                         ("\u0394 confined [%]", float),
                                                         ("\u0394 free [%]", float),
                                                         ("trajectories included", int),
                                                         ("trajectories excluded", int),
                                                         ("mean D immobile+no type [\u03BCm\u00b2/s]", float),
                                                         ("mean D confined [\u03BCm\u00b2/s]", float),
                                                         ("mean D free [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D immobile+no type [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D confined [\u03BCm\u00b2/s]", float),
                                                         ("\u0394 mean D free [\u03BCm\u00b2/s]", float),
                                                         ("mean length immobile+no type [frames]", float),
                                                         ("mean length confined [frames]", float),
                                                         ("mean length free [frames]", float),
                                                         ("\u0394 mean length immobile+no type [frames]", float),
                                                         ("\u0394 mean length confined [frames]", float),
                                                         ("\u0394 mean length free [frames]", float)]))
        dset["immobile+no type [%]"] = type_percentages[4]
        dset["confined [%]"] = type_percentages[1]
        dset["free [%]"] = type_percentages[2]
        dset["\u0394 immobile+no type [%]"] = dtype_percentages[4]
        dset["\u0394 confined [%]"] = dtype_percentages[1]
        dset["\u0394 free [%]"] = dtype_percentages[2]
        dset["trajectories included"] = trajectories_included
        dset["trajectories excluded"] = trajectories_excluded
        dset["mean D immobile+no type [\u03BCm\u00b2/s]"] = mean_Ds[4]
        dset["mean D confined [\u03BCm\u00b2/s]"] = mean_Ds[1]
        dset["mean D free [\u03BCm\u00b2/s]"] = mean_Ds[2]
        dset["\u0394 mean D immobile+no type [\u03BCm\u00b2/s]"] = mean_dDs[4]
        dset["\u0394 mean D confined [\u03BCm\u00b2/s]"] = mean_dDs[1]
        dset["\u0394 mean D free [\u03BCm\u00b2/s]"] = mean_dDs[2]
        dset["mean length immobile+no type [frames]"] = mean_lengths[4]
        dset["mean length confined [frames]"] = mean_lengths[1]
        dset["mean length free [frames]"] = mean_lengths[2]
        dset["\u0394 mean length immobile+no type [frames]"] = mean_dlengths[4]
        dset["\u0394 mean length confined [frames]"] = mean_dlengths[1]
        dset["\u0394 mean length free [frames]"] = mean_dlengths[2]
        
    def diffusion_plot(self, number, diffusion, mean_cell, dmean_cell):
        """
        Diffusion plot.
        """
        dset = self.grp05.create_dataset("histogram values", (number,),
                                         dtype=np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                         ("mean frequency cells", float),
                                                         ("\u0394 mean frequency cells", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion    
        dset["mean frequency cells"] = mean_cell
        dset["\u0394 mean frequency cells"] = dmean_cell
        
    def diffusion_plot_bg(self, number, diffusion, mean_cell, dmean_cell, mean_bg, dmean_bg, mean_cell_corr, dmean_cell_corr):
        """
        Diffusion plot, background correction, normalized.
        """
        dset = self.grp05.create_dataset("histogram values", (number,),
                                         dtype=np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                         ("mean frequency cells", float),
                                                         ("\u0394 mean frequency cells", float),
                                                         ("mean frequency background", float),
                                                         ("\u0394 mean frequency background", float),
                                                         ("mean frequency corrected", float),
                                                         ("\u0394 mean frequency corrected", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion
        dset["mean frequency cells"] = mean_cell
        dset["\u0394 mean frequency cells"] = dmean_cell
        dset["mean frequency background"] = mean_bg
        dset["\u0394 mean frequency background"] = dmean_bg
        dset["mean frequency corrected"] = mean_cell_corr
        dset["\u0394 mean frequency corrected"] = dmean_cell_corr
    
    def diffusion_plot_normalized(self, number, diffusion, mean_cell_percent, dmean_cell_percent):
        """
        Diffusion plot, normalized.
        """
        dset = self.grp05.create_dataset("histogram values normalized", (number,),
                                         dtype=np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                         ("mean frequency cells [%]", float),
                                                         ("\u0394 mean frequency cells [%]", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion    
        dset["mean frequency cells [%]"] = mean_cell_percent
        dset["\u0394 mean frequency cells [%]"] = dmean_cell_percent

    def diffusion_plot_normalized_types(self, number, diffusion, mean_cell_immob, dmean_cell_immob, mean_cell_conf,
                                        dmean_cell_conf, mean_cell_free, dmean_cell_free, mean_cell_notype,
                                        dmean_cell_notype, mean_cell_immob_notype, dmean_cell_immob_notype):
        """
        Diffusion plot, types, normalized.
        """
        dset = self.grp05.create_dataset("histogram values normalized types", (number,),
                                         dtype=np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                         ("mean frequency immob [%]", float),
                                                         ("\u0394 mean frequency immob [%]", float),
                                                         ("mean frequency conf [%]", float),
                                                         ("\u0394 mean frequency conf [%]", float),
                                                         ("mean frequency free [%]", float),
                                                         ("\u0394 mean frequency free [%]", float),
                                                         ("mean frequency notype [%]", float),
                                                         ("\u0394 mean frequency notype [%]", float),
                                                         ("mean frequency immob+notype [%]", float),
                                                         ("\u0394 mean frequency immob+notype [%]", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion
        dset["mean frequency immob [%]"] = mean_cell_immob
        dset["\u0394 mean frequency immob [%]"] = dmean_cell_immob
        dset["mean frequency conf [%]"] = mean_cell_conf
        dset["\u0394 mean frequency conf [%]"] = dmean_cell_conf
        dset["mean frequency free [%]"] = mean_cell_free
        dset["\u0394 mean frequency free [%]"] = dmean_cell_free
        dset["mean frequency notype [%]"] = mean_cell_notype
        dset["\u0394 mean frequency notype [%]"] = dmean_cell_notype
        dset["mean frequency immob+notype [%]"] = mean_cell_immob_notype
        dset["\u0394 mean frequency immob+notype [%]"] = dmean_cell_immob_notype

    def diffusion_plot_bg_normalized(self, number, diffusion, mean_cell_percent, dmean_cell_percent,
                                     mean_cell_corr_percent, dmean_cell_corr_percent):
        """
        Diffusion plot, background correction, normalized.
        """
        dset = self.grp05.create_dataset("histogram values normalized", (number,),
                                         dtype=np.dtype([("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                         ("mean frequency cells [%]", float),
                                                         ("\u0394 mean frequency cells [%]", float),
                                                         ("mean frequency corrected [%]", float),
                                                         ("\u0394 mean frequency corrected [%]", float)]))
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion
        dset["mean frequency cells [%]"] = mean_cell_percent
        dset["\u0394 mean frequency cells [%]"] = dmean_cell_percent
        dset["mean frequency corrected [%]"] = mean_cell_corr_percent
        dset["\u0394 mean frequency corrected [%]"] = dmean_cell_corr_percent
    
    def diffusion_bin_size(self, bin_size):
        dset = self.grp05.create_dataset("bin size", (1, 1), dtype=np.dtype([("bin size", float)]))
        dset["bin size"] = bin_size

    def average_MSD(self, dataset_name, number, delta_t, MSD_values, MSD_errors):
        dset = self.grp07.create_dataset(dataset_name, (number,),
                                         dtype=np.dtype([("delta_t [s]", float),
                                                         ("MSD [\u03BCm\u00b2]", float),
                                                         ("\u0394 MSD [\u03BCm\u00b2]", float)]))
        dset["delta_t [s]"] = delta_t
        dset["MSD [\u03BCm\u00b2]"] = MSD_values
        dset["\u0394 MSD [\u03BCm\u00b2]"] = MSD_errors
