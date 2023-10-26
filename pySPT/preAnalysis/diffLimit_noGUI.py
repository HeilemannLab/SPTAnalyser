"""
@author: Johanna Rahm, Alexander Niedrig
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Determine the minimal distance localizations can be told apart based on nearest neighbor analysis.
Draws no graphs to allow batch processing.
"""


import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm_notebook as tqdm


class GridNNSearch:
    """
    The computing time of the nearest neighbor analysis is accelerated by a grid, which lies over the space
    of the localizations, thereby creating subregions. Only neighbours in adjacent subregions are accepted as
    potential nearest neighbor candidates for a center.
    """

    def __init__(self, centers, neighbors, number_subregions, number_pixels, pixel_size, max_search_area):
        self.number_pixels = number_pixels  # 256
        self.pixel_size = pixel_size  # 158 nm
        self.number_subregions = number_subregions  # e.g. ~ number of neighbor points
        self.border_length = number_pixels * pixel_size / number_subregions  # length of subregion border
        self.neighbors = neighbors  # list of list of neighbors xy-locations
        self.centers = centers  # xy-locations of centers
        self.grid_neighbors = self.create_grid(number_subregions)  # subregions linked to their points (xy-locs)
        self.center_is_neighbor = self.center_is_neighbor()  # list of booleans if center list equals a neighbor list
        self.max_search_area = max_search_area  # stop criterium to find nearest neighbor for sparse regions

    def center_is_neighbor(self):
        """
        Check if the center list is equal with each of the neighbor lists.
        :return: List of booleans if center list equals neighbor lists
        """
        center_is_neighbor = []
        for neighbor_list in self.neighbors:
            is_same = True if np.array_equal(self.centers, neighbor_list) else False
            center_is_neighbor.append(is_same)
        return center_is_neighbor

    def create_grid(self, number_subregions):
        """
        Divide the space into subregions that store the localizations of neighbors in their subspace.
        :param number_subregions: Total number of subregions.
        :return: Dict of key = subregion index & value = list of neighbors, their original indices and reference
        to neighbor list.
        """
        # dict with key = XY and value = []
        grid_neighbors = {}
        for x in range(number_subregions):
            for y in range(number_subregions):
                grid_neighbors.update({(x, y): []})
        # store neighbor xy-locations in corresponding subregion
        for c, neighbor_list in enumerate(self.neighbors):
            for i, neighbor in enumerate(neighbor_list):
                subregion = self.point_to_subregion(neighbor)
                # i = idx in list, c = idx of neighbor list (FGFR1 -> 0, FGFR2 -> 1 ...)
                grid_neighbors[subregion].append((neighbor, i, c))
        
        return grid_neighbors

    def get_nn_distances(self):
        """
        Calculate the distance to the nearest neighbor and get its type for each center.
        :return: Nearest neighbor distances and types.
        """
        nn_distances = []
        nn_types = []
        for center_idx, center in enumerate(self.centers):
            subregion = self.point_to_subregion(center)
            search_area = 0
            found_neighbor = False
            while not found_neighbor:
                # check if neighbors are in search_area range & if the neighbor is not the particle itself
                # if no neighbors are found, increase the search_area
                # for each valid sub region get the xy-localizations
                neighbors = list(
                    map(lambda x: self.grid_neighbors[x], self.get_valid_sub_regions(subregion, search_area)))
                # merge elements of sublists to one list
                neighbors_conc = [j for i in neighbors for j in i]
                # delete the element in neighbors that equals center
                for c, neighbor in enumerate(neighbors_conc):
                    if np.array_equal(neighbor[0], center) and neighbor[1] == center_idx:
                        neighbors_conc.pop(c)
                # set appropriate search area (area that found neighbor +1)
                if neighbors_conc:
                    search_area += 1
                    found_neighbor = True
                elif search_area > self.max_search_area:
                    break
                else:
                    search_area += 1
            valid_sub_regions = self.get_valid_sub_regions(subregion, search_area)
            # get all potential nearest neighbor candidates from the valid subregions
            neighbors_in_area = list(map(lambda x: self.grid_neighbors[x], valid_sub_regions))
            neighbors_in_area_conc = [j for i in neighbors_in_area for j in i]
            for c, neighbor in enumerate(neighbors_in_area_conc):
                if np.array_equal(neighbor[0], center) and neighbor[1] == center_idx:
                    neighbors_in_area_conc.pop(c)
            # calc nearest neighbor distance and type (type refers to the index of lists of self.neighbors)
            min_distance, nn_type = self.calc_min_distance(center, neighbors_in_area_conc)
            nn_distances.append(min_distance)
            nn_types.append(nn_type)
        return nn_distances, nn_types

    def calc_min_distance(self, center, neighbor_candidates):
        """
        Calculate the euclidean distance of candidates in a list and return the nearest neighbor distance and type.
        :param center: Target center.
        :param neighbor_candidates: Potential nearest neighbors.
        :return: Nearest neighbor distance and type (type = index of list in self.neighbors).
        """
        min_distance = math.inf
        nn_type = ""
        for idx in range(len(neighbor_candidates)):
            distance = np.linalg.norm(center - neighbor_candidates[idx][0])
            if distance < min_distance:
                min_distance = distance
                nn_type = neighbor_candidates[idx][2]
        return min_distance, nn_type

    def point_to_subregion(self, point):
        """
        Sort a xy-localization to a subregion.
        :param point: xy-localization.
        :return: (x, y) refers to the subregion indices.
        """
        x = int(np.floor(point[0] / self.border_length))
        y = int(np.floor(point[1] / self.border_length))
        # this only happens if the localization is directly at the border of the measurement space.
        if x == self.number_subregions:
            x -= 1
        if y == self.number_subregions:
            y -= 1
        return x, y

    def get_valid_sub_regions(self, center, search_area):
        """
        From a center subregion get all valid subregions within the appropriate search area.
        :param center: Subregion indices.
        :param search_area: Valid deviations of center subregion-id.
        :return: List of valid subregion-ids.
        """
        valid_sub_regions = []
        for x in range(center[0] - search_area, center[0] + search_area + 1):
            for y in range(center[1] - search_area, center[1] + search_area + 1):
                if x in range(self.number_subregions) and y in range(self.number_subregions):
                    valid_sub_regions.append((x, y))
        return valid_sub_regions

class DiffLimit():
    def __init__(self):
        self.px_size = 0  # pixel size in nm
        self.n_px = 0  # number of pixels in a row
        self.nn_distances = []  # nearest neighbors per frame of localization files
        self.min_nn_distances = []  # minimal nearest neighbor per localization file
        self.files = []  #  pd localization files
        self.file_paths = []  # paths to localization files
        self.file_names = []  # names of localization files
        self.figure = []  # boxplot of nearest neighbor distances
        self.max_search_area = 0  # max search area (stop criterium for sparse regions)

    def clear_object(self):
        self.__init__()

    def get_files(self, dir, file_ending):
        file_path = dir
        files, file_paths, file_names = [], [], []
        files = []
        if file_ending == ".txt":
            localization_header = "# <localizations insequence="
            for i in os.listdir(file_path):
                if i.endswith(file_ending) and "tracked" not in i:
                    file_path = dir + "\\" + i
                    pd_file = pd.read_csv(file_path, sep=r"\t")
                    if localization_header in str(pd_file.columns):
                        file_names.append(i)
                        file_paths.append(file_path)
                        files.append(pd.read_csv(file_path, skiprows=1, header=None, sep=" "))
        elif file_ending == ".csv":
            for i in os.listdir(file_path):
                if i.endswith(file_ending) and "tracked" not in i:
                    file_path = dir + "\\" + i
                    pd_file = pd.read_csv(file_path)
                    if list(pd_file.columns) == ['id', 'frame', 'x [nm]', 'y [nm]', 'sigma [nm]', 'intensity [photon]',
                                                 'offset [photon]', 'bkgstd [photon]', 'uncertainty_xy [nm]']:
                        file_names.append(i)
                        file_paths.append(file_path)
                        files.append(pd_file)
        return files, file_paths, file_names

    def xy_stack(self, x, y):
        """
        Calculate the xy-coordinates in nm and stack them together.

        :param x: x coordinates [nm].
        :param y: y coordinates [nm].
        :return: xy-coordinates with shape (len, 2).
        :rtype: Numpy ndarray.
        """
        xy_stack = np.column_stack((x, y))
        return xy_stack

    def plot_min_nn_distances(self):
        fig = plt.figure(figsize=(3, 5))
        ax = sns.boxplot(data=self.min_nn_distances, color="cornflowerblue", showmeans=True,
                         meanprops={"marker": "s", "markerfacecolor": "white", "markeredgecolor": "0.25"})
        ax = sns.swarmplot(data=self.min_nn_distances, color="0.25")
        ax.set_ylabel("min nearest neighbor distance per file [nm]")
        self.figure = fig

    def run_diff_limit(self, path, file_ending):
        self.files, self.file_paths, self.file_names = self.get_files(path, file_ending)
        if len(self.files) == 0:
            print("No matching files found in ", path)
        else:
            with tqdm(total=len(self.files), desc="Localization file") as pbar:
                for file in self.files:
                    file_distances = []
                    if file_ending == ".csv":
                        for frame in range(min(file["frame"]), max(file["frame"])+1):
                            frame_file = file[file["frame"] == frame]
                            xy_positions = self.xy_stack(frame_file["x [nm]"].to_numpy(), frame_file["y [nm]"].to_numpy())
                            neighbor_positions = []
                            neighbor_positions.append(xy_positions)
                            if len(xy_positions) > 0:
                                grid_nn_search = GridNNSearch(xy_positions, list(neighbor_positions),
                                                              int(np.floor(math.sqrt(len(xy_positions)))), self.n_px,
                                                              self.px_size, self.max_search_area)
                                nn_distances, _ = grid_nn_search.get_nn_distances()
                                file_distances.extend(nn_distances)
                    elif file_ending == ".txt":
                        for frame in range(int(min(file.iloc[:,4])), int(max(file.iloc[:,4]))):
                            frame_file = file[file.iloc[:,4] == frame]
                            xy_positions = self.xy_stack(frame_file.iloc[:,0].to_numpy(), frame_file.iloc[:,2].to_numpy())
                            neighbor_positions = []
                            neighbor_positions.append(xy_positions)
                            if len(xy_positions) > 0:
                                grid_nn_search = GridNNSearch(xy_positions, list(neighbor_positions),
                                                              int(np.floor(math.sqrt(len(xy_positions)))), self.n_px,
                                                              self.px_size, self.max_search_area)
                                nn_distances, _ = grid_nn_search.get_nn_distances()
                                file_distances.extend(nn_distances)
                    self.nn_distances.append(file_distances)
                    self.min_nn_distances.append(min(file_distances))
                    pbar.update(1)
            self.plot_min_nn_distances()

    def save(self, dir_path, folder_name, file_header, save_plot):
        os.mkdir(dir_path + "\\" + folder_name)
        nn_distances_pd = pd.DataFrame(self.nn_distances)
        nn_distances_pd = nn_distances_pd.transpose()
        nn_distances_pd.columns = file_header
        nn_distances_pd.to_csv(dir_path + "\\" + folder_name + "\\nearest_neighbor_distances.csv", index=False)
        min_nn_distances_pd = pd.DataFrame(self.min_nn_distances)
        min_nn_distances_pd = min_nn_distances_pd.transpose()
        min_nn_distances_pd.columns = file_header
        min_nn_distances_pd.to_csv(dir_path + "\\" + folder_name + "\\min_nearest_neighbor_distances.csv", index=False)
        if save_plot:
            self.figure.savefig(dir_path + "\\" + folder_name + "\\min_nearest_neighbor_boxplot.pdf",
                format="pdf", transparent=True, bbox_inches="tight")
