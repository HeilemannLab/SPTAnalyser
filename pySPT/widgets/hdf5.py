# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 15:01:01 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import numpy as np
import h5py
import tkinter as tk
from tkinter import filedialog
import os


class Hdf5():
    # one cell
    def __init__(self):
        self.h5_file = []
        self.grp00 = []
        self.grp01 = []
        self.grp02 = []
        self.grp03 = []
        self.grp04 = []
        self.test_path = "C:\\Users\\pcoffice37\\Documents\\Python Scripts\\Hdf5\\test_files\\" + "dataset.h5"
        self.trc_file = ""
        self.trc_file_hdf5 = ""
        
        
    def create_h5_file(self):
        self.h5_file = h5py.File(self.test_path, "w")  # w- or x = Create file, fail if exists
        print("success")
    
    def get_trc_file(self):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.txt_file = filedialog.askopenfilename(filetypes=(("text files", "*.trc"), ("all files", "*.*")))
        root.update()
        root.destroy()
        self.trc_file = root.txt_file
    
    def load_trc(self):
        self.trc = np.loadtxt(self.trc_file, skiprows = 0, usecols=(0,1,2,3,4,5))
        print(self.trc)
                        
    def create_h5_name(self):        
        self.trc_file_hdf5 = os.path.splitext(self.trc_file)[0] + ".h5"       
        
    def safe_h5_file(self):
        shape = np.shape(self.trc)  # lxb
        file = h5py.File(self.trc_file_hdf5, "w")  # w = write, r = read ...
        dset = file.create_dataset("trc", (shape[0],),
                                   dtype = np.dtype([("id", int),
                                                      ("frame", int),
                                                      ("x [pixel]", float),
                                                      ("y [pixel]", float),
                                                      ("intensity", float)]))
        dset["frame"] = [self.trc[:,0]]
        dset["x [pixel]"] = [self.trc[:,1]]
        dset["y [pixel]"] = [self.trc[:,2]]
        dset["intensity"] = [self.trc[:,3]]
        # tidy
        file.close()
        
    def groups(self):
        self.grp00 = self.h5_file.create_group("trc")
        self.grp01 = self.h5_file.create_group("MSD")
        self.grp02 = self.h5_file.create_group("diffusion")
        self.grp03 = self.h5_file.create_group("rossier")
        self.grp04 = self.h5_file.create_group("settings")
        #self.h5_file.close()
        
    def data_settings(self, dt, pixelsize, tau_threshold, fit_area, dof, dloc_dyn):
        dset = self.grp04.create_dataset("settings", (1,1), dtype = np.dtype([("dt [s]", float), ("pixelsize [nm]", int), ("tau threshold [s]", float), ("fit area", float), ("dof", int), ("dloc dyn [\u03BCm\u00b2/s]", float)]))
        dset["dt [s]"] = dt
        dset["pixelsize [nm]"] = pixelsize
        dset["tau threshold [s]"] = tau_threshold
        dset["fit area"] = fit_area
        dset["dof"] = dof
        dset["dloc dyn [\u03BCm\u00b2/s]"] = dloc_dyn
        #self.h5_file.close()
        
    def data_diffusion_info(self, number, trajectory_id, diffusion_coeff, ddiffusion_coeff, MSD_0, chi2, length):
        dset = self.grp02.create_dataset("diffusion infos", (number,), dtype = np.dtype([("trajectory id", int),
                                                             ("diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                             ("d diffusion coefficient [\u03BCm\u00b2/s]", float),
                                                             ("MSD(0) [\u03BCm\u00b2]", float),
                                                             ("chi\u00b2 [\u03BCm\u2074]", float),
                                                             ("length [nm]", int)]))
        dset["trajectory id"] = trajectory_id
        dset["diffusion coefficient [\u03BCm\u00b2/s]"] = diffusion_coeff
        dset["d diffusion coefficient [\u03BCm\u00b2/s]"] = ddiffusion_coeff
        dset["MSD(0) [\u03BCm\u00b2]"] = MSD_0
        dset["chi\u00b2 [\u03BCm\u2074]"] = chi2
        dset["length [nm]"] = length
        #self.h5_file.close()

    def data_rossier_info(self, number, trajectory_id, type_immobile, type_confined, type_free, analyse_success, tau, dtau, r, dr, dconfined, chi2):
        dset = self.grp03.create_dataset("rossier infos", (number,), dtype = np.dtype([("trajectory id", int),("type immobile", int),("type confined", int),("type free", int),("analyse succesful", int),("tau [s]", float),("dtau [s]", float),("r [\u03BCm]", float),("dr [\u03BCm]", float),("Diffusion confined [\u03BCm\u00b2/s]", float), ("chi\u00b2 [\u03BCm\u2074]", float)]))
        dset["trajectory id"] = trajectory_id
        dset["type immobile"] = type_immobile
        dset["type confined"] = type_confined
        dset["type free"] = type_free
        dset["analyse succesful"] = analyse_success
        dset["tau [s]"] = tau
        dset["dtau [s]"] = dtau
        dset["r [\u03BCm]"] = r
        dset["dr [\u03BCm]"] = dr
        dset["Diffusion confined [\u03BCm\u00b2/s]"] = dconfined
        dset["chi\u00b2 [\u03BCm\u2074]"] = chi2
        self.h5_file.close()
        
    def data_rossier(self, dt, MSD, fit, residues):
        pass
    
    
        
# =============================================================================
#     def data_rossier_info(self, number, trajectory_id, tau, dtau, r, dr, dconfined, ddconfined, analyse_successful):
#         dset = self.grp03.create_dataset("rossier infos", (number,), dtype = np.dtype([("trajectory id", int),("tau [s]", float),("dtau [s]", float),("r [\u03BCm]", float),("dr [\u03BCm]", float),("Diffusion confined [\u03BCm\u00b2/s]", float),("dDiffusion confined [\u03BCm\u00b2/s]", float), ("analyse successful",str)]))
#         dset["trajectory id"] = trajectory_id
# # =============================================================================
# #         dset["type immobile"] = type_immobile
# #         dset["type confined"] = type_confined
# #         dset["type free"] = type_free
# # =============================================================================
# # =============================================================================
# #         dset["analyse succesful"] = analyse_success
# # =============================================================================
#         dset["tau [s]"] = tau
#         dset["dtau [s]"] = dtau
#         dset["r [\u03BCm]"] = r
#         dset["dr [\u03BCm]"] = dr
#         dset["Diffusion confined [\u03BCm\u00b2/s]"] = dconfined
#         dset["dDiffusion confined [\u03BCm\u00b2/s]"] = ddconfined
#         dset["analyse succesful"] = [np.string_(i) for i in analyse_successful]
#         self.h5_file.close()
# =============================================================================
        
        
        
        
                              
                              
def main():
    hdf5 = Hdf5()
    #hdf5.get_trc_file()
    print(hdf5.trc_file)
    #hdf5.trc_file = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.MIA\\tracking\\cell_1_MMStack_Pos0.ome_MIA-1-D.txt"
    #hdf5.get_trc_file()
    #hdf5.load_trc()
    #hdf5.create_h5_name()
    #hdf5.safe_h5_file()
    hdf5.create_h5_file()
    print(hdf5.h5_file.name)
    hdf5.groups()
    hdf5.group_settings()

if __name__ == "__main__":
    main()
        
    
    