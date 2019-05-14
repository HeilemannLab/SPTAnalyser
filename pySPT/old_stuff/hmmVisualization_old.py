# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:19:31 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Visualize the HMM analysis results.
"""

import h5py 
import os
import numpy as np
import matplotlib.pyplot as plt


class HMMCell():
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
        self.trc = []  # numpy.ndarray with numpy.voids
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
        dset_observation_alphabet = self.group_hmm["observationAlpabet"]
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
        
    def get_trc(self):
        dset_trc = self.group_trc["trcFile"]
        self.trc = dset_trc[:]
        print(self.trc)
        
    def show_hdf5_file(self):
        print(self.hmm_cell_hdf5)
        
    def close_hdf5_file(self):
        self.hmm_cell_hdf5.close()
        
    def run(self):
        self.get_statistics_info()
        self.get_transition_matrix()
        self.get_equilibrium_matrix()
        self.get_observation_matrix()
        self.get_observation_alphabet()
        self.get_weight_coef()
        self.get_diffusion_coef()
        self.get_trc()
        
        self.close_hdf5_file()
    
        
def main():
    cell01_path = "C:\\Users\\pcoffice37\\Documents\\testing_file_search\\testHdf5Merge\\CS5_Cell01Cell02\\190412_cell_1_MMStack_Pos0_trc_format_merged_3States.h5"
    cell02_path = "C:\\Users\\pcoffice37\\Documents\\testing_file_search\\testHdf5Merge\\CS5_Cell01Cell02\\190412_cell_2_MMStack_Pos0_trc_format_merged_3States.h5"
    hmm_cell01 = HMMCell(cell01_path)
    hmm_cell02 = HMMCell(cell02_path)
    hmm_cell01.run()
    hmm_cell02.run()
    
    # checkbox which cell is added to cs? 
    cover_slip = []
    cover_slip.append(hmm_cell01)
    cover_slip.append(hmm_cell02)
    get_cell_names(cover_slip)
    get_AIC(cover_slip)
    calc_mean_hmm_pi(cover_slip)
    calc_mean_physmod_pi(cover_slip)
    calc_mean_tp(cover_slip)
    calc_mean_states(cover_slip)
    plot_jd(hmm_cell01)
    plot_D(cover_slip)


# visualization

number_of_states = 3
colour_palett = ["orangered", "royalblue", "forestgreen", "darkmagenta", "orange"]

def get_cell_names(cover_slip):
    cell_names = []
    for cell in cover_slip:
        cell_names.append(cell.hmm_cell_name)
    print("Cell names: ", cell_names)

def get_AIC(cover_slip):
    aic_values = np.zeros(len(cover_slip))
    for cell in cover_slip:
        cell_index = cover_slip.index(cell)
        aic_values[cell_index] = cell.aic
    print("AIC values: ", aic_values)

def plot_D(cover_slip):
    mean_diff_coeff = np.zeros(number_of_states)
    diff_coeffs = np.zeros(number_of_states)
    print(diff_coeffs)
    print(mean_diff_coeff)
    np.zeros((number_of_states, len(cover_slip)))
    for state in range(number_of_states):
        diff_cell = np.zeros()
        
    # create matrix with rows = states, columns = cells
# =============================================================================
#     for state in range(number_of_states):
#         #diff_cell = 
#         for cell in cover_slip:
#             cell_index = cover_slip.index(cell)
#             diffusion_coef_values[state][cell_index] = cell.diffusion_coef[state][0]
#     print("HMM D: ", diffusion_coef_values)
# =============================================================================
    
    
    
def plot_jd(cell):
    pass

def plot_bar(bars, label_name, title_name, error=None, y_error=False):
    x = np.arange(len(bars))
    x_name = [i+1 for i in range(number_of_states)]
    fig, ax = plt.subplots()
    ax.set_ylim(0,1)
    if y_error:
        plt.bar(x, bars, yerr=error, capsize=5, edgecolor="black", color=colour_palett[:number_of_states])  # label=label_name
        plt.xticks(x, x_name)
    else:
        plt.bar(x, bars, capsize=5, edgecolor="black", color=colour_palett[:number_of_states])
    #plt.legend()
    plt.title(title_name)
    plt.show()    
    
def calc_mean_hmm_pi(cover_slip):
    number_of_cells = len(cover_slip)
    mean_hmm_pis = np.zeros(number_of_states)
    mean_hmm_pis_error = np.zeros(number_of_states)
    for i in range(number_of_states):
        mean_hmm_pi = np.zeros(number_of_cells)
        for cell in cover_slip:
            cell_index = cover_slip.index(cell)
            mean_hmm_pi[cell_index] = cell.equilibrium_matrix[i]
        mean_hmm_pis[i] = np.mean(mean_hmm_pi)
        mean_hmm_pis_error[i] = np.std(mean_hmm_pi, ddof=1)/(number_of_cells)**(1/2)
    plot_bar(mean_hmm_pis, "probability", "State distribution based on equilibrium matrix", error=mean_hmm_pis_error, y_error=True)
    print("mean values", mean_hmm_pis)
    
def calc_mean_physmod_pi(cover_slip):
    number_of_cells = len(cover_slip)
    mean_physmod_pis = np.zeros(number_of_states)
    mean_physmod_pis_error = np.zeros(number_of_states)
    for i in range(number_of_states):
        mean_physmod_pi = np.zeros(number_of_cells)
        for cell in cover_slip:
            cell_index = cover_slip.index(cell)
            mean_physmod_pi[cell_index] = cell.weight_coef[i][0]
        mean_physmod_pis[i] = np.mean(mean_physmod_pi)
        mean_physmod_pis_error[i] = np.std(mean_physmod_pi, ddof=1)/(number_of_cells)**(1/2)
    plot_bar(mean_physmod_pis, "probability", "State distribution based on weighted physical model", error=mean_physmod_pis_error, y_error=True)
    print("mean values", mean_physmod_pis)
    
def calc_mean_states(cover_slip):
    number_of_cells = len(cover_slip)
    mean_states = np.zeros(number_of_states)
    mean_states_error = np.zeros(number_of_states)
    cell_states = np.zeros(number_of_cells)
    for cell in cover_slip:
        cell_index = cover_slip.index(cell)
        state_counter = np.zeros(number_of_states)
        for i in range(len(cell.trc)):
            for state_number in range(number_of_states):
                if cell.trc[i][4] == state_number:
                    state_counter[state_number] += 1
        state_counter = np.divide(state_counter, np.sum(state_counter))
        ######## save that
    
    
def calc_mean_tp(cover_slip):
    number_of_cells = len(cover_slip)
    print("tp")
    mean_tps = np.zeros(np.shape(cover_slip[0].transition_matrix))
    mean_tps_error = np.zeros(np.shape(cover_slip[0].transition_matrix))
    for row in range(number_of_states):
        for column in range(number_of_states):
            mean_tp = np.zeros(number_of_cells)
            for cell in cover_slip:
                cell_index = cover_slip.index(cell)
                mean_tp[cell_index] = cell.transition_matrix[column][row]
            mean_tps[column][row] = np.mean(mean_tp)
            mean_tps_error[column][row]= np.std(mean_tp, ddof=1)/(number_of_cells)**(1/2)
    print("mean tp values", mean_tps)

    

    
if __name__ == "__main__":
    main()
