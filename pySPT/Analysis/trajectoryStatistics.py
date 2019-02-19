# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 11:45:55 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

from . import trajectory
from . import cell
import numpy as np
import copy
import math

class TrajectoryStatistics():
    def __init__(self):
        self.cell_trajectories = [] # [[],[]] contains list of cells, cells contain trajectories
        self.cell_trajectories_filtered = []
        self.cell_trajectories_index = []
        self.cell_trajectories_filtered_index = []
        self.cell_trajectories_name = []
        self.cell_trajectories_log = []
        self.roi_size = 0.0
        
    #def run_statistics(self, min_length, max_length, min_D, max_D):
    def run_statistics(self, min_length, max_length, min_D, max_D):
        self.get_index()
        self.create_init_filter_lst()
        if max_length == "max length":
            max_length = self.get_max_length()
            print("max length:", max_length)
        if min_length == "min length":
            min_length = self.get_min_length()
            print("min length:", min_length)
        if max_D == "max value":
            max_D = self.get_max_D()
            print("max diffusion coefficient:", max_D)
        if min_D == "min value":
            min_D = self.get_min_D()
            print("min diffusion coefficient:", min_D)
        self.filter_length(int(min_length), int(max_length))
        self.filter_D(float(min_D), float(max_D))
        #self.filter_type()

        
    def create_init_filter_lst(self):
        """
        Create copy of initial cell trajectories & index list.
        """
        self.cell_trajectories_filtered = copy.deepcopy(self.cell_trajectories)
        self.cell_trajectories_filtered_index = copy.deepcopy(self.cell_trajectories_index)

    def get_index(self):
        """
        Each cell is a list, create i lists for the cells. In the lists append the numbers of elements = trajectory numbers.
        Starting with 1.
        """
        for cell in range(0, len(self.cell_trajectories)):
            i = []
            for trajectory in range(0, len(self.cell_trajectories[cell])):
                i.append(trajectory+1)  # trajectory numbering starts with 1
            self.cell_trajectories_index.append(i)
        
        
    def plot_trajectory(self, cell, number):
        cell = int(cell) - 1
        number = int(number) -1
        self.cell_trajectories[cell][number].plot_particle()

    def get_max_length(self):
        """
        Max length of trajectory.
        :return: int of max length.
        """
        max_length = - math.inf
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if max_length < trajectory.length_trajectory:
                    max_length = trajectory.length_trajectory
        return int(max_length)

    def get_min_length(self):
        """
        Max length of trajectory.
        :return: int of max length.
        """
        min_length = math.inf
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if min_length > trajectory.length_trajectory:
                    min_length = trajectory.length_trajectory
        return int(min_length)
                
# =============================================================================
#     def filter_length(self, min_length=None, max_length=None): # max_length=default
#         """
#         Filter by length of trajectory (PALMTracer min_length = 20). Trajectories in the range between
#         min_length and max_length will be filtered -> min_length <= trajectory <= max_length.
#         :param min_length: minimal length of trajectories to be considered.
#         :param max_length: by default None -> None will be evaluated as the max length trajectory
#         """
#         if max_length == None or max_length == "":
#             max_length = self.get_max_length()
#             print("max_length", max_length)
#         if min_length == None or min_length == "":
#             min_length = self.get_min_length()
#             print("min length", min_length)
#         for cell in range(0, len(self.cell_trajectories_filtered)):
#             for trajectory in self.cell_trajectories_filtered[cell]:
#                 if trajectory.length_trajectory < int(min_length) or trajectory.length_trajectory > int(max_length):
#                     self.crop_lst(cell, trajectory)
#         # show selected trajectory objects & index
#         #print(self.cell_trajectories_filtered_index) 
#         # show selected trajectory lengths
#         for cell in range(0, len(self.cell_trajectories_filtered)):
#             for trajectory in self.cell_trajectories_filtered[cell]:
#                 pass
#                 #print(trajectory.length_trajectory)          
# =============================================================================
                
    def filter_length(self, min_length, max_length): # max_length=default
        """
        Filter by length of trajectory (PALMTracer min_length = 20). Trajectories in the range between
        min_length and max_length will be filtered -> min_length <= trajectory <= max_length.
        :param min_length: minimal length of trajectories to be considered.
        :param max_length: by default None -> None will be evaluated as the max length trajectory
        """
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if trajectory.length_trajectory < int(min_length) or trajectory.length_trajectory > int(max_length):
                    self.crop_lst(cell, trajectory)
        # show selected trajectory objects & index
        #print(self.cell_trajectories_filtered_index) 
        # show selected trajectory lengths
        #for cell in range(0, len(self.cell_trajectories_filtered)):
        #    for trajectory in self.cell_trajectories_filtered[cell]:
        #        pass
                #print(trajectory.length_trajectory)     

    def crop_lst(self, cell, trajectory):
        """
        Delete trajectory and its index.
        """
        trajectory_index = self.cell_trajectories_filtered[cell].index(trajectory)
        self.cell_trajectories_filtered[cell] = self.cell_trajectories_filtered[cell][0:trajectory_index] + \
        self.cell_trajectories_filtered[cell][trajectory_index+1:]
        self.cell_trajectories_filtered_index[cell] = self.cell_trajectories_filtered_index[cell][0:trajectory_index] + \
        self.cell_trajectories_filtered_index[cell][trajectory_index+1:] 
        
    def filter_type(self, trajectory_types):
        """
        Filter by trajectory type:
        Immobile, confined = True -> Confined, confined = False -> Free.
        :param trajectory_type: list with types ["immobile", "free"] would filter immob & free and neglect confined.
        """
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if "immobile" not in trajectory_types:
                    # immobile: immobility = True, confined = True
                    if trajectory.immobility and trajectory.confined:
                        self.crop_lst(cell, trajectory)
                if "confined" not in trajectory_types:
                    # confined: immobility = False, confinded = True
                    if trajectory.confined and not trajectory.immobility:
                        self.crop_lst(cell, trajectory)
                    # free: immobility = False, confined = False
                if "free" not in trajectory_types:
                    if not trajectory.confined and not trajectory.immobility:
                        self.crop_lst(cell, trajectory)
        # show selected trajectory objects & index
        print(len(self.cell_trajectories_filtered[0]))
        #print(self.cell_trajectories_filtered)
        print(self.cell_trajectories_filtered_index) 
        # show selected trajectory lengths
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                print("immob & confined", trajectory.immobility, trajectory.confined)
                print("length", trajectory.length_trajectory)
    
    def get_max_D(self):
        max_D = - math.inf
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if max_D < trajectory.D:
                    max_D = trajectory.D
        return float(max_D)
    
    def get_min_D(self):
        min_D = math.inf
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if min_D > trajectory.D:
                    min_D = trajectory.D
        return float(min_D)
        
    def filter_D(self, min_D, max_D):
        for cell in range(0, len(self.cell_trajectories_filtered)):
            for trajectory in self.cell_trajectories_filtered[cell]:
                if trajectory.D < min_D or trajectory.D > max_D:
                    self.crop_lst(cell, trajectory)
        # show selected trajectory objects & index
        #print(self.cell_trajectories_filtered_index) 
        # show selected trajectory lengths
        #for cell in range(0, len(self.cell_trajectories_filtered)):
        #    for trajectory in self.cell_trajectories_filtered[cell]:
        #        print(trajectory.D, trajectory.length_trajectory)
        
    def background_frequencies(self):
        """
        Create list with log10(D) for each trajectory in each cell of list.
        """
        self.cell_trajectories_log = []
        print(type(self.cell_trajectories_log))
        for cell in range(0, len(self.cell_trajectories_filtered)):
            i = []
            for trajectory in range(0, len(self.cell_trajectories_filtered[cell])):
                i.append(np.log10(self.cell_trajectories_filtered[cell][trajectory].D))
            self.cell_trajectories_log.append(i)
        #print(self.cell_trajectories_log)
        
    def frequency_count(self):
        max_bin = int(np.ceil(self.cell_trajectories_log.max()/0.05)*0.05)
        bin_size = int(np.ceil(self.cell_trajectories_log.max()/0.05))
        hist = np.histogram(self.cell_trajectories_log,
                            range = (0, max_bin),
                            bins = bin_size,
                            density = True)
        self.log_histogram = np.zeros([np.size(hist[0]),2])
        self.log_histogram[:,0] = hist[1][:-1]  # log(D)
        self.log_histogram[:,1] = hist[0][:]  # frequencies
        self.log_histogram[:,1] = self.normalize_hist(self.log_histogram[:,1])
        
    def normalize_hist(self, normalized_col):
        """
        Normalize a column and return it.
        """
        normalized_col = normalized_col / np.sum(normalized_col)
        return normalized_col
        
 
    
def main():
    pass
        
    
if __name__ == "__main__":
    main()
        