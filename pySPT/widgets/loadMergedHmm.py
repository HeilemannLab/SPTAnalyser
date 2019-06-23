# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 17:37:57 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import h5py
import os


class LoadMergedHmm():
    def __init__(self, file_path):
        self.hmm_cell_hdf5 = h5py.File(file_path, "r")
        self.hmm_cell_name = os.path.splitext(os.path.split(file_path)[1])[0]  # raw base name
        # hmm
        self.transition_matrix = []  # numpy.ndarray with numpy.ndarrays for each state
        self.equilibrium_matrix = []
        self.observation_matrix = []
        self.observation_alphabet = []  # numpy.ndarray
        # hmm -> statistics
        self.number_of_states = 0  # floats/int
        self.symbols = 0
        self.log_likelohood = 0
        self.dof = 0
        self.bic = 0
        self.aic = 0
        # judi
        self.judi = []
        # physical model
        self.diffusion_coef = []  # numpy.ndarray with numpy.voids
        self.weight_coef = []
        # trc
        self.trc_hmm = []  # numpy.ndarray with numpy.voids
        self.trc_type = []
        #cell size
        self.cell_size = 0.0
        # groups
        self.group_msd = self.hmm_cell_hdf5["MSD"]
        self.group_diffusion = self.hmm_cell_hdf5["diffusion"]
        self.group_hmm = self.hmm_cell_hdf5["hmm"]
        self.group_judi = self.hmm_cell_hdf5["judi"]
        self.group_physical_model = self.hmm_cell_hdf5["physicalModel"]
        self.group_rossier = self.hmm_cell_hdf5["rossier"]
        self.group_settings = self.hmm_cell_hdf5["settings"]
        self.group_statistics = self.hmm_cell_hdf5["statistics"]
        self.group_trc = self.hmm_cell_hdf5["trc"]
        self.group_settings = self.hmm_cell_hdf5["settings"]

    def get_statistics_info(self):
        """
        Statistics dataset contains states, symbols, loglikelihood, dof, bic, aic.
        """
        dset_statistics = self.group_hmm["statistics"]
        self.number_of_states = dset_statistics["states"][0]
        self.symbols = dset_statistics["symbols"][0]
        self.log_likelohood = dset_statistics["logLikelihood"][0]
        self.dof = dset_statistics["dof"][0]
        self.bic = dset_statistics["bic"][0]
        self.aic = dset_statistics["aic"][0]

    def get_transition_matrix(self):
        dset_transition_matrix = self.group_hmm["transitionMatrix"]
        self.transition_matrix = dset_transition_matrix[:,:]

    def get_equilibrium_matrix(self):
        dset_equilibrium_matrix = self.group_hmm["equilibriumMatrix"]
        self.equilibrium_matrix = dset_equilibrium_matrix[0]

    def get_observation_alphabet(self):
        ##### typo
        dset_observation_alphabet = self.group_hmm["observationAlphabet"]
        self.observation_alphabet = dset_observation_alphabet[0]

    def get_observation_matrix(self):
        dset_observation_matrix = self.group_hmm["observationMatrix"]
        self.observation_matrix = dset_observation_matrix[:,:]

    def get_diffusion_coef(self):
        dset_diffusion_coef = self.group_physical_model["diffusionCoefficient"]
        self.diffusion_coef = dset_diffusion_coef[:]
        #print(self.diffusion_coef)

    def get_weight_coef(self):
        dset_weight_coef = self.group_physical_model["weightCoefficient"]
        self.weight_coef = dset_weight_coef[:]
        #print(self.weight_coef)

    def get_judi(self):
        #####
        dset_judi = self.group_judi["judiFile"]
        #self.judi = dset_judi
        #print(dset_judi[:,])

    def get_trcs(self):
        dset_trc_hmm = self.group_trc["trcHmm"]
        self.trc_hmm = dset_trc_hmm[:]
        dset_trc_type = self.group_trc["trcType"]
        self.trc_type = dset_trc_type[:]
        #print(self.trc)

    def show_hdf5_file(self):
        print(self.hmm_cell_hdf5)

    def close_hdf5_file(self):
        self.hmm_cell_hdf5.close()
        
    def get_cell_size(self):
        dset_settings = self.group_settings["settings"]
        self.cell_size = dset_settings[0][0][3]

    def run(self):
        self.get_statistics_info()
        self.get_transition_matrix()
        self.get_equilibrium_matrix()
        self.get_observation_matrix()
        self.get_observation_alphabet()
        self.get_weight_coef()
        self.get_diffusion_coef()
        self.get_trcs()
        self.get_cell_size()
        self.close_hdf5_file()
