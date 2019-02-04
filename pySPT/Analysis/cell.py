# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 10:53:42 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import numpy as np
from . import trajectory
from multiprocessing import Pool
#from .analysis import trajectory
#from pySPT.analysis import trajectory
import time

class Cell():
    def __init__(self):
        #self.all_trajectories = []  # col0 = trajectory, col1 = frames, col2 = x, col3 = y, col4 = intensity ~ trc file
        self.pixel_size = 158  # [nm] multiply with palmtracer in nm -> *10^-3 micrometer!
        self.trajectories = []  # contains trajectory objects
        self.analysed_trajectories = []  # contains analysed trajectory objects
        self.D = []
        self.states = np.zeros([3,1])
        
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

    def load_file(self):
        """
        Load file & create list with trajectory objects, xy localizations for each object are initialized. 
        """
        file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.MIA\\tracking\\cell_1_MMStack_Pos0.ome_MIA.trc"
        all_trajectories = np.loadtxt(file_name, usecols = (0, 1, 2, 3, 5)) # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = intensity
        all_trajectories[:,2] = np.multiply(all_trajectories[:,2], int(self.pixel_size)*10**(-3))
        all_trajectories[:,3] = np.multiply(all_trajectories[:,3], int(self.pixel_size)*10**(-3))
        print(all_trajectories[:,2])
        print(all_trajectories[:,3])
        for trajectory_number in range(int(all_trajectories[:,0].min()), int(all_trajectories[:,0].max())+1):    
            idx = all_trajectories[:,0] == trajectory_number
            localizations = all_trajectories[idx,:]
            if not (localizations.size==0):
                self.trajectories.append(trajectory.Trajectory(localizations))

    def analyse_trajectories(self):
        """
        Analyse trajectories without multiprocessing
        """
        list(map(lambda x: x.analyse_particle(), self.trajectories))
        self.analysed_trajectories = self.trajectories
        #print(self.analysed_trajectories)
        #print(self.trajectories)
        #for i in self.trajectories:
        #    print("Diff", i.D)

        
#    def plot_trajectorie(self, trajectory_number):
#        self.trajectories[trajectory_number-1].plot_particle()
        
# =============================================================================
#     def load_file(self):
#         file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.MIA\\tracking\\cell_1_MMStack_Pos0.ome_MIA.trc"
#         self.all_trajectories = np.loadtxt(file_name, usecols = (0, 1, 2, 3, 5)) # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = intensity
# 
#     def get_trajectories(self):
#         for trajectory_number in range(int(self.all_trajectories[:,0].min()), int(self.all_trajectories[:,0].max())+1):            
#             one_trajectory = trajectory.Trajectory(trajectory_number, len(self.all_trajectories[:,0]), self.all_trajectories)
#             self.trajectories.append(one_trajectory)
#             
#     def analyse_trajectories(self):
#         for i in self.trajectories:
#             i.analyse_particle()
#         #self.trajectories[:].analyse_particle()
#             
#     def plot_trajectorie(self, trajectory_number):
#         self.trajectories[trajectory_number-1].plot_particle()
# 
#         
#         
#     def get_localization(self):       
#         for i in range(0, self.len_all_trajectories):
#             if self.all_trajectories[i,0] == self.trajectory_number:
#                 self.localizations.append((self.all_trajectories[i,2], self.all_trajectories[i,3]))
# =============================================================================

            
            #localizations = all_trajectories[idx,:]
            #idx = np.searchsorted(all_trajectories[:,0], [trajectory_number, trajectory_number+1])
            #print(idx)
            #localizations = all_trajectories[idx[0]:idx[1],:]
            #print(localizations)
            #localization_tuple = (localizations[:,2], localizations[:,3])
            #print (localization_tuple)
        
# =============================================================================
#             for i in range(0, len(all_trajectories[:,0])):
#                 if all_trajectories[i,0] == trajectory_number:
#                     localizations.append((all_trajectories[i,2], all_trajectories[i,3]))
#             one_trajectory = trajectory.Trajectory()
#             one_trajectory.localizations = localizations
#             self.trajectories.append(one_trajectory)
# =============================================================================
        #print(self.trajectories)


# =============================================================================
# def main():
#     #start = time.time()
#     cell = Cell()
#     file_type = "palmtracer"
#     if file_type == "swift":
#         cell.pixel_size = 1
#     elif file_type == "palmtracer":
#         cell.pixel_size = 158
#     cell.load_file()
#     #cell.get_trajectories()
#     #cell.analyse_trajectories()
#     #end = time.time()
#     #print("Evalutation took {} seconds".format(end-start))
#     #cell.plot_trajectorie(1)
#     
# =============================================================================
# #     
# if __name__ == "__main__":
#     main()
# #     
# =============================================================================
# =============================================================================
