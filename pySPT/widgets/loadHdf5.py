# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 14:38:40 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import h5py
import numpy as np
import pandas as pd
import os


class LoadHdf5():
    def __init__(self):
        self.file_names = []  # coverslip.cell_files, list with file names in loading order
        self.names = []  # cell.name, raw base name of file, has to be distributed to indiv cells
        self.hdf5 = []  # list of xxx objects
        self.trajectory_numbers = []  # amount of trajectories per cell
        self.cell_numbers = []  # amount of cells
        self.sizes = []  # list of cell sizes -> cell.size
        self.pixel_sizes = []  # list of pixel sizes -> cell.pixel_size
        self.pixel_amounts = []  # list of amount of pixel of detector -> cell.pixel_amount
        self.dts = []  # integration times -> trajectory.dt
        self.tau_thresholds = []  # -> trajectory.tau_threshold
        self.fit_areas = []  # -> trajectory.fit_area
        self.dofs = []  # -> trajectory.dof
        self.D_mins = []  # -> trajectory.D_min
        
        self.cells_lengths_trajectories = []  # [[xxx],[yyy]] contains list objects = cells, values = trajectories -> trajectory.length_trajectory
        self.cells_lengths_MSDs = []  # -> trajectory.length_MSD
        
        self.cells = []  # list of cell objects ->cover_slip.cells -> needs entire list
        
    def test_create_file_names(self):
        file_name01 = "C:\\Users\\pcoffice37\\Documents\\testing_file_search\\cells\\cell_1_MMStack_Pos0.ome_MIA.h5"
        file_name02 = "C:\\Users\\pcoffice37\\Documents\\testing_file_search\\cells\\subdirectory\\cell_2_MMStack_Pos0.ome_MIA.h5"
        self.file_names.append(file_name01)
        self.file_names.append(file_name02)
        
    def read_h5(self):
        for file_name in self.file_names:
            h5_file = h5py.File(file_name, "r")
            self.hdf5.append(h5_file)
            
    def count_trajectory_numbers(self):
        """
        Create list with numbers of trajectories for each cell [int, int ...].
        """
        for h5 in self.hdf5:
            diffusion_group = h5["diffusion"]
            diffusion_infos_data = diffusion_group["diffusionInfos"].value
            self.trajectory_numbers.append(np.shape(diffusion_infos_data)[0])
        print(self.trajectory_numbers)

    def count_cells(self):
        self.cell_numbers = len(self.hdf5)
        
    def get_cell_name(self):
        """
        From a path create the raw base name.
        """
        for cell in self.file_names:
            base_name = os.path.basename(cell)
            raw_base_name = ""
            for i in base_name:
                if i == ".":
                    break
                else:
                    raw_base_name += i   
            self.names.append(raw_base_name)
        print("names", self.names)
        
    def get_cell_size(self):
        """
        Create list with cell sizes. 
        """
        for h5 in self.hdf5:
            size_group = h5["size"]
            size_data = size_group["size"].value  # sizedata = [(396.80278,)]
            size = size_data[0][0]
            self.sizes.append(size)
        print("size", self.sizes)        
    
    def settings(self):
        """
        Handling pixel size & amount.
        """
        for h5 in self.hdf5:
            settings_group = h5["settings"]
            settings_data = settings_group["settings"].value  # [[(0.02, 158, 65536, 0.12, 0.6, 4, 0.0065)]]
            self.dts.append(settings_data[0][0][0])
            self.pixel_sizes.append(settings_data[0][0][1])
            self.pixel_amounts.append(settings_data[0][0][2])
            self.tau_thresholds.append(settings_data[0][0][3])
            self.fit_areas.append(settings_data[0][0][4])
            self.dofs.append(settings_data[0][0][5])
            self.D_mins.append(settings_data[0][0][6])
        print("settings", self.dts, self.pixel_sizes, self.pixel_amounts, self.tau_thresholds, self.fit_areas, self.dofs, self.D_mins)
        
    def create_np_array(self, length, columns=1):
        """
        :param length: number of np arrays.
        :param columns: amount of entries in a numpy array.
        :return: return a numpy array.
        """
        np_array = np.zeros((length,columns))
        return np_array
        
    def get_lengths_trajectories_MSDs(self):         
        for h5 in self.hdf5:
            h5_index = self.hdf5.index(h5)
            max_trajectory_index = self.trajectory_numbers[h5_index]
            diffusion_group = h5["diffusion"]
            diffusion_infos_data = diffusion_group["diffusionInfos"].value
            length_trajectories = self.create_np_array(max_trajectory_index)
            length_MSDs = self.create_np_array(max_trajectory_index)
            for trajectory_index in range(0, max_trajectory_index):
                length_trajectory = diffusion_infos_data[trajectory_index][5]
                length_trajectories[trajectory_index] = length_trajectory
                length_MSDs[trajectory_index] = length_trajectory - 1
            self.cells_lengths_trajectories.append(length_trajectories)
            self.cells_lengths_MSDs.append(length_MSDs)

# =============================================================================
#     def get_lengths_trajectories_MSDs(self):
#         """
#         Handling length_MSD, length_trajectory for trajectory objects.
#         """
#         for h5 in self.hdf5:
#             h5_index = self.hdf5.index(h5)
#             max_trajectory_index = self.trajectory_numbers[h5_index]
#             lengths_trajectories = []
#             lengths_MSDs = []
#             diffusion_group = h5["diffusion"]
#             diffusion_infos_data = diffusion_group["diffusionInfos"].value
#             for trajectory_index in range(0, max_trajectory_index):
#                 length_trajectory = diffusion_infos_data[trajectory_index][5]
#                 lengths_trajectories.append(length_trajectory)
#                 length_MSD = length_trajectory - 1
#                 lengths_MSDs.append(length_MSD)
#             self.lengths_trajectories.append(lengths_trajectories)
#             self.lengths_MSDs.append(lengths_MSDs)
# =============================================================================


            
            
            
            
            
            
            
            
# =============================================================================
#             data_frame_size = h5.get("size/size")
#             self.sizes.append(data_frame_size.values[0][0])  # cell size is stored as [[size]] -> acess float by [0][0]
#         print(self.sizes)
# =============================================================================
        
        
            #print(h5.get("/size/size"))
            #print(type(h5.get("/size/size")))


            #df = h5.get("/size/size")
            #print(df.columns)
            
            #df = h5.get("/diffusion/diffusion plots/diffusion plot 0001")
            #print(df)
            #print(df.keys())
            #print(type(df.values))
            
            #data_frame_size = h5.get("size/size")
            
            #df.frame["size [\u03BCm\u00b2]"]
    

        #print(self.hdf5[0].keys())
# =============================================================================
#         data_frame_trc = self.hdf5[0].get("trc/trc file")
#         print(data_frame_trc)
# =============================================================================
        #for h5 in self.hdf5:
        #    print(h5.keys())
            #data_frame_pixelsize = h5.get("settings/settings")
            #print(data_frame_pixelsize)
# =============================================================================
#             self.pixel_sizes.append(data_frame_pixelsize)
#         print(self.pixel_sizes)
# =============================================================================
    

            

        
# =============================================================================
#         
#         file = h5py.File(self.file_names[0], "r")
#         print(file.keys())
#         keys = [key for key in file.keys()]
#         print(keys)
#         group = file["rossier"]
#         print(group)
#         for key in group.keys():
#             print(key)
#             
#             
#         data = group["rossierStatistics"].value
#         print(len(data))
#         print(type(data))
#         print(type(data[0]))
#         print(data[0][10])
#         #print(data)
# =============================================================================

        #data = list(file[key])

        
        
        
    

def main():
    load_h5 = LoadHdf5()
    load_h5.test_create_file_names()
    load_h5.read_h5()
    load_h5.get_cell_name()
    load_h5.count_trajectory_numbers()
    load_h5.count_cells()
    load_h5.get_cell_size()
    load_h5.settings()
    #load_h5.create_np_array(10)
    load_h5.get_lengths_trajectories_MSDs()
# =============================================================================
#     load_hdf5.load_pd()
#     load_hdf5.get_cell_name()
#     load_hdf5.get_cell_size()
#     load_hdf5.settings_cell()
# =============================================================================

if __name__ == "__main__":
    main()
    