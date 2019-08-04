# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 10:53:42 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Create a cell object which creats trajectory objects. Cell contains informations about cell size and name.
"""

import numpy as np
import math
import copy
from . import trajectory
from ..hmm import microscope
#from multiprocessing import Pool
from tqdm import tqdm_notebook as tqdm
#from .analysis import trajectory
#from pySPT.analysis import trajectory

class Cell():
    def __init__(self):
        self.trc_file_type = []  # col0 = track id, col1 = frames, col2 = x, col3 = y, col4 = placeholder, col5 = intensity, col6 = seg id ~ trc file
        self.trc_file_hmm = []  # col0 = track id, col1 = track id continuously, col2 = frames, col3 = x, col4 = y, col5 = placeholder, col6 = intensity
        self.filtered_trc_file_hmm = []  # trc file after analysis with trajectories sigma D > 0
        self.converted_trc_file_type = []  # trc file as np arrays instead of lists in list
        self.trajectories = []  # contains trajectory objects
        self.analysed_trajectories = []  # contains analysed trajectory objects
        self.trajectories_hmm = []
        self.analysed_trajectories_hmm = []
        self.pixel_size = 0.0  # [nm] multiply with palmtracer in nm -> *10^-3 micrometer!
        self.pixel_amount = 0 # amount of pixels of detector (256x256) 
        self.size = 0.0  # size from roi [ym^2]
        self.name = ""  # name of cell file (raw base name)
        self.tau_threshold = 0.0  # hand down to trajectory
        self.dt = 0.0   # hand down to trajectory
        self.dof = 0.0  # hand down to trajectory
        self.rossier_fit_area = 0.0  # hand down to trajectory
        self.D_min = 0.0  # hand down to trajectory
        self.points_fit_D = 0  # hand down to trajectory
        self.min_track_length_type = 0.0 
        self.min_track_length_hmm = 0.0
        self.seg_id = True  # if True the seg id will be loaded as trajectory id, else the track id will be loaded
        self.sigma_dyn_type = 0.0  # dynamic localization error
        self.sigma_dyn_hmm = 0.0
        
        
    def run_analysis(self):
        # trc hmm
        self.cell_size()
        self.create_trajectories_hmm()  # based on the trc hmm file trajectories are created
        self.run_analysis_hmm()  # create hmm trajectory objects, if D > 0 -> trajectories are taken into account
        self.filter_trc_hmm()  # filter trc file for trajectories with D > 0
        #self.filter_trc_hmm_del(2)
        # trc type
        self.create_trajectories()
        self.analyse_trajectories()
        self.calc_sigma_dyn_hmm()  # take mean value of sigma dyn of filtered trajectories
        # self.calc_sigma_dyn_type()  # at this point calculating sigma dyn does not make sence, because D < 0 and not fitted trajectories exist.
        self.convert_trc_type()  # np arrays instead of lists in list
            
    def create_trajectories_hmm(self):
        trc_file = np.zeros([len(self.trc_file_hmm),6])
        trc_file[:,0] = list(map(lambda row: row[0], self.trc_file_hmm))  # col0 = track id
        trc_file[:,1] = list(map(lambda row: row[1], self.trc_file_hmm))  # frame
        trc_file[:,2] = list(map(lambda row: row[2]*int(self.pixel_size)*10**(-3), self.trc_file_hmm))  # x in ym
        trc_file[:,3] = list(map(lambda row: row[3]*int(self.pixel_size)*10**(-3), self.trc_file_hmm))  # y in ym
        trc_file[:,4] = list(map(lambda row: row[4], self.trc_file_hmm))  # placeholder
        trc_file[:,5] = list(map(lambda row: row[5], self.trc_file_hmm))  # intensity
        for trajectory_number in range(int(trc_file[:,0].min()), int(trc_file[:,0].max())+1):    
            idx = trc_file[:,0] == trajectory_number
            # get rid of first localizations
            #idx = self.localizations_del(2, idx)
            localizations = trc_file[idx,:]
            #print("Trajectory number len loc", trajectory_number, len(localizations))
            if not (localizations.size==0):
                self.trajectories_hmm.append(trajectory.Trajectory(localizations, self.tau_threshold, self.dt, self.dof, self.D_min, self.points_fit_D, self.rossier_fit_area))

    def run_analysis_hmm(self):
        """
        The diffusion coefficient and sigma dyn are calculated. If D < 0 the trajectory
        will be neglected.
        """
        for trajectory in self.trajectories_hmm:
            trajectory.analyse_particle()
            if trajectory.D > 0:
                self.analysed_trajectories_hmm.append(trajectory)
        #print("neg D", len([trajectory.D for trajectory in self.analysed_trajectories_hmm]))
        #print("neg D neg MSD0", len([(trajectory.D, trajectory.MSD_0) for trajectory in self.analysed_trajectories_hmm if trajectory.D < 0 and trajectory.MSD_0 < 0]))
# =============================================================================
#         count_D_MSD0_pos = [(trajectory.D, trajectory.MSD_0) for trajectory in self.analysed_trajectories_hmm if trajectory.D > 0 and trajectory.MSD_0 > 0]
#         count_D_MSD0_neg = [(trajectory.D, trajectory.MSD_0) for trajectory in self.analysed_trajectories_hmm if trajectory.D < 0 and trajectory.MSD_0 < 0]
#         count_Dpos_MSD0neg = [(trajectory.D, trajectory.MSD_0) for trajectory in self.analysed_trajectories_hmm if trajectory.D > 0 and trajectory.MSD_0 < 0]
#         count_Dneg_MSD0pos = [(trajectory.D, trajectory.MSD_0) for trajectory in self.analysed_trajectories_hmm if trajectory.D < 0 and trajectory.MSD_0 > 0]
#         print(count_D_MSD0_pos, count_D_MSD0_neg, count_Dpos_MSD0neg, count_Dneg_MSD0pos)
#         
#         summe = len(count_D_MSD0_pos) + len(count_D_MSD0_neg) + len(count_Dpos_MSD0neg) + len(count_Dneg_MSD0pos) 
#         print(len(count_D_MSD0_pos)/summe, len(count_D_MSD0_neg)/summe, len(count_Dpos_MSD0neg)/summe, len(count_Dneg_MSD0pos)/summe)
# =============================================================================
        self.analysed_trajectories_hmm = self.trajectories_hmm  
        
    def filter_trc_hmm(self):
        """
        The trc file is filtered so that it only contain trajectories with positive D and track length >= min track length hmm.
        """
        trajectory_idx = [trajectory.trajectory_number for trajectory in self.analysed_trajectories_hmm]
        trc_file = np.zeros([len(self.trc_file_hmm),6])
        trc_file[:,0] = list(map(lambda row: row[0], self.trc_file_hmm))  # col0 = track id
        trc_file[:,1] = list(map(lambda row: row[1], self.trc_file_hmm))  # frame
        trc_file[:,2] = list(map(lambda row: row[2]*int(self.pixel_size)*10**(-3), self.trc_file_hmm))  # x in ym
        trc_file[:,3] = list(map(lambda row: row[3]*int(self.pixel_size)*10**(-3), self.trc_file_hmm))  # y in ym
        trc_file[:,4] = list(map(lambda row: row[4], self.trc_file_hmm))  # placeholder
        trc_file[:,5] = list(map(lambda row: row[5], self.trc_file_hmm))  # intensity
        trc_idx = np.isin(trc_file[:,0], trajectory_idx)
        self.filtered_trc_file_hmm = trc_file[trc_idx,:]
        
    def localizations_del(self, locs_to_del, idx):
        """
        Delete the first localizations of a trajectory. Target: localizations based on which the trajectories are created.
        :param locs_to_del: Amount of localizations to delete.
        """
        # get rid of first two localizations
        count = 0
        idx_count = 0
        for i in idx:
            if i and count < 2:
                idx[idx_count] = False
                count += 1
            idx_count += 1
        return idx
                
    def filter_trc_hmm_del(self, locs_to_del):
        """
        Delete the first localizations of a trajectory. Target: saved filtered trc file hmm.
        :param locs_to_del: Amount of localizations to delete.
        """
        trajectory_idx = [trajectory.trajectory_number for trajectory in self.analysed_trajectories_hmm]
        for idx in trajectory_idx:
            entry_idx = 0
            count = 0
            for entry in self.filtered_trc_file_hmm[:,0]:    
                if idx == entry and count < locs_to_del:
                    self.filtered_trc_file_hmm = np.delete(self.filtered_trc_file_hmm, entry_idx, 0)
                    count += 1
                    entry_idx -= 1
                entry_idx += 1
    
    def calc_sigma_dyn_hmm(self):
        """
        Calculate the dynamic localization error, based on the mean D, mean MSD_0, dt and dof values.
        Only trajectories with D < 0 are taken into account.
        """
        MSD_0 = np.mean([trajectory.MSD_0 for trajectory in self.analysed_trajectories_hmm])
        D = np.mean([trajectory.D for trajectory in self.analysed_trajectories_hmm])
        #print("MSD0 & D", MSD_0, D, [trajectory.MSD_0 for trajectory in self.analysed_trajectories_hmm], [trajectory.D for trajectory in self.analysed_trajectories_hmm])
        self.sigma_dyn_hmm = math.sqrt((MSD_0+(4/3)*D*self.dt)/4)
    
    def calc_sigma_dyn_type(self):
        """
        Calculate the dynamic localization error, based on the mean D, mean MSD_0, dt and dof values.
        """
        MSD_0 = np.mean([trajectory.MSD_0 for trajectory in self.analysed_trajectories if trajectory.D > 0])
        D = np.mean([trajectory.D for trajectory in self.analysed_trajectories if trajectory.D > 0])
        self.sigma_dyn_type = math.sqrt((MSD_0+(4/3)*D*self.dt)/4)
        
    def create_trajectories(self):
        """
        Create list with trajectory objects from a trc file.
        Convert x & y positions from px to ym. 
        Choose which trajectory id will be taken for creating tracks.
        """
        trc_file = np.zeros([len(self.trc_file_type),7])
        if self.seg_id:
            trc_file[:,0] = list(map(lambda row: row[6], self.trc_file_type))  # col0 = seg id
        else: 
            trc_file[:,0] = list(map(lambda row: row[0], self.trc_file_type))  # col0 = track id
        trc_file[:,1] = list(map(lambda row: row[1], self.trc_file_type))  # frame
        trc_file[:,2] = list(map(lambda row: row[2]*int(self.pixel_size)*10**(-3), self.trc_file_type))  # x in ym
        trc_file[:,3] = list(map(lambda row: row[3]*int(self.pixel_size)*10**(-3), self.trc_file_type))  # y in ym
        trc_file[:,4] = list(map(lambda row: row[4], self.trc_file_type))  # placeholder
        trc_file[:,5] = list(map(lambda row: row[5], self.trc_file_type))  # intensity
        #print("trc", sorted(list(set(trc_file[:,0]))))
        for trajectory_number in range(int(trc_file[:,0].min()), int(trc_file[:,0].max())+1):    
            idx = trc_file[:,0] == trajectory_number
            localizations = trc_file[idx,:]
            if not (localizations.size==0):
                self.trajectories.append(trajectory.Trajectory(localizations, self.tau_threshold, self.dt, self.dof, self.D_min, self.points_fit_D, self.rossier_fit_area))
        
    def convert_trc_type(self):
        """
        To save the trc type file later, all information is stored.
        """
        trc_file = np.zeros([len(self.trc_file_type),7])
        trc_file[:,0] = list(map(lambda row: row[0], self.trc_file_type))  # col0 = track id
        trc_file[:,1] = list(map(lambda row: row[1], self.trc_file_type))  # frame
        trc_file[:,2] = list(map(lambda row: row[2]*int(self.pixel_size)*10**(-3), self.trc_file_type))  # x in ym
        trc_file[:,3] = list(map(lambda row: row[3]*int(self.pixel_size)*10**(-3), self.trc_file_type))  # y in ym
        trc_file[:,4] = list(map(lambda row: row[4], self.trc_file_type))  # placeholder
        trc_file[:,5] = list(map(lambda row: row[5], self.trc_file_type))  # intensity
        trc_file[:,6] = list(map(lambda row: row[6], self.trc_file_type))  # col0 = seg id
        self.converted_trc_file_type = trc_file
        
    def cell_size(self):
        """
        Convert cell size in pixel^2 to micrometer^2.
        If no Roi was loaded, cell size = amount of pixel^2 * pixelsize^2.
        """
        if self.size == 0:
            self.size = self.pixel_amount * (self.pixel_size/1000)**2  # in micrometer
        else:
            self.size = self.size * (self.pixel_size/1000)**2

    def analyse_trajectories(self):
        """
        Analyse trajectories without multiprocessing
        """
        for trajectory in tqdm(self.trajectories):
            trajectory.analyse_particle()
        self.analysed_trajectories = self.trajectories
        #print("immobile", len([(trajectory.D) for trajectory in self.analysed_trajectories if trajectory.D < self.D_min]), [(trajectory.D) for trajectory in self.analysed_trajectories if trajectory.D < self.D_min])
        #print("confined", len([(trajectory.tau, trajectory.D) for trajectory in self.analysed_trajectories if trajectory.tau < self.tau_threshold and not trajectory.immobility and trajectory.analyse_successful]), [(trajectory.tau, trajectory.D) for trajectory in self.analysed_trajectories if trajectory.tau < self.tau_threshold and not trajectory.immobility and trajectory.analyse_successful])
        #print("free", len([(trajectory.tau, trajectory.D) for trajectory in self.analysed_trajectories if trajectory.tau > self.tau_threshold and not trajectory.immobility and trajectory.analyse_successful]), [(trajectory.tau, trajectory.D) for trajectory in self.analysed_trajectories if trajectory.tau > self.tau_threshold and not trajectory.immobility and trajectory.analyse_successful])
        #print("not success", len([trajectory for trajectory in self.analysed_trajectories if not trajectory.analyse_successful]), [trajectory for trajectory in self.analysed_trajectories if not trajectory.analyse_successful])
        #print("total trajectories", len(self.analysed_trajectories))

        
    def plot_trajectory(self, trajectory_number):
        """
        Plot a defined trajectory by its number.
        """
        self.analysed_trajectories[trajectory_number-1].plot_particle()
        
    # Multiprocessing -> bugged ...
    @staticmethod
    def analyse_trajectory(self, trajectory):
        """
        Multiprocessing job: analyse 1 trajectory object.
        """
        trajectory.analyse_particle()
        return trajectory      
     
    def run_processes(self):
        """
        Target: analyse trajectory,
        arg: one trajectory out of the self.trajectories list
        """
        with Pool(None) as p:
            self.analysed_trajectories = p.map(self.analyse_trajectory, self.trajectories)
            return self.analysed_trajectories
        
        
def main():
    pass


if __name__ == "__main__":
    main()

