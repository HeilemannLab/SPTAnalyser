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
#from multiprocessing import Pool
from tqdm import tqdm_notebook as tqdm
#from .analysis import trajectory
#from pySPT.analysis import trajectory

class Cell():
    def __init__(self):
        self.trc_file_type = []  # col0 = track id, col1 = frames, col2 = x, col3 = y, col4 = placeholder, col5 = intensity, col6 = seg id ~ trc file
        self.trc_file_hmm = []  # col0 = track id, col1 = track id continuously, col2 = frames, col3 = x, col4 = y, col5 = placeholder, col6 = intensity
        self.filter_trc_file_hmm = []
        self.trajectories = []  # contains trajectory objects
        self.analysed_trajectories = []  # contains analysed trajectory objects
        self.trajectories_hmm = []
        self.analysed_trajectories_hmm = []
        self.pixel_size = 158  # [nm] multiply with palmtracer in nm -> *10^-3 micrometer!
        self.pixel_amount = 65536 # amount of pixels of detector (256x256) 
        self.size = 0.0  # size from roi [ym^2]
        self.name = ""  # name of cell file (raw base name)
        self.tau_threshold = 0.0  # hand down to trajectory
        self.dt = 0.0   # hand down to trajectory
        self.dof = 0.0  # hand down to trajectory
        self.D_min = 0.0  # hand down to trajectory
        self.points_fit_D = 4  # hand down to trajectory
        self.min_track_length_type = 0.0 
        self.min_track_length_hmm = 0.0
        self.seg_id = True  # if True the seg id will be loaded as trajectory id, else the track id will be loaded
        self.sigma_dyn_type = 0.0  # dynamic localization error
        self.sigma_dyn_hmm = 0.0
        
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
            localizations = trc_file[idx,:]
            if not (localizations.size==0):
                self.trajectories_hmm.append(trajectory.Trajectory(localizations, self.tau_threshold, self.dt, self.dof, self.D_min, self.points_fit_D))

    def run_analysis_hmm(self):
        """
        The diffusion coefficient and sigma dyn are calculated. If D < 0 the trajectory
        will be neglected.
        """
        for trajectory in self.trajectories_hmm:
            trajectory.analyse_particle()
            if trajectory.D > 0:
                self.analysed_trajectories_hmm.append(trajectory)
        #self.analysed_trajectories_hmm = self.trajectories_hmm  
        
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
    
    def calc_sigma_dyn_hmm(self):
        """
        Calculate the dynamic localization error, based on the mean D, mean MSD_0, dt and dof values.
        """
        self.sigma_dyn_hmm = np.mean([trajectory.sigma_dyn for trajectory in self.analysed_trajectories_hmm])
    
    def calc_sigma_dyn_type(self):
        """
        Calculate the dynamic localization error, based on the mean D, mean MSD_0, dt and dof values.
        """
        self.sigma_dyn_type = np.mean([trajectory.sigma_dyn for trajectory in self.analysed_trajectories if trajectory.D > 0])
# =============================================================================
#         sigma_dyn_type = np.mean([trajectory.sigma_dyn for trajectory in self.analysed_trajectories])
#         print("sigma with D > 0", self.sigma_dyn_type)
#         print("sigma", sigma_dyn_type)
# =============================================================================
        
    def create_trajectories(self):
        """
        Create list with trajectory objects from a trc file.
        Convert x & y positions from px to ym.
        """
        trc_file = np.zeros([len(self.trc_file_type),6])
        if self.seg_id:
            trc_file[:,0] = list(map(lambda row: row[6], self.trc_file_type))  # col0 = seg id
        else: trc_file[:,0] = list(map(lambda row: row[0], self.trc_file_type))  # col0 = track id
        trc_file[:,1] = list(map(lambda row: row[1], self.trc_file_type))  # frame
        trc_file[:,2] = list(map(lambda row: row[2]*int(self.pixel_size)*10**(-3), self.trc_file_type))  # x in ym
        trc_file[:,3] = list(map(lambda row: row[3]*int(self.pixel_size)*10**(-3), self.trc_file_type))  # y in ym
        trc_file[:,4] = list(map(lambda row: row[4], self.trc_file_type))  # placeholder
        trc_file[:,5] = list(map(lambda row: row[5], self.trc_file_type))  # intensity

        for trajectory_number in range(int(trc_file[:,0].min()), int(trc_file[:,0].max())+1):    
            idx = trc_file[:,0] == trajectory_number
            localizations = trc_file[idx,:]
            if not (localizations.size==0):
                self.trajectories.append(trajectory.Trajectory(localizations, self.tau_threshold, self.dt, self.dof, self.D_min, self.points_fit_D))
                
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
        
    def plot_trajectory(self, trajectory_number):
        """
        Plot a defined trajectory by its number.
        """
        self.analysed_trajectories[trajectory_number-1].plot_particle()
        
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
            print(self.analysed_trajectories)
            return self.analysed_trajectories
        
        
def main():
    pass


if __name__ == "__main__":
    main()

