"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Handle data loading for transitionCounts.
Data: h5 files after trackStatistics loaded from directory, "statistics.h5" is ignored.
"""


import os
import pandas as pd
import h5py


class DataLoader():
    def __init__(self, dir_paths):
        self.h5_names, self.h5_files = self.get_h5_files(dir_paths)
        self.tracked_names, self.tracked_files = self.get_tracked_files(dir_paths)
        self.got_files = True if self.h5_files and self.tracked_files else False

    def get_h5_files(self, dirs):
        """
        Read all h5 files of multiple directories, ignore files with "statistics" in name.
        :param dirs: List of directories.
        :return: List with file names and h5 objects.
        """
        all_file_names = []
        all_files = []
        for dir in dirs:
            file_names = []
            files = []
            for file in os.listdir(dir):
                if file.endswith(".h5") and "statistics" not in file:
                    file_names.append(file)
                    files.append(h5py.File(dir + "\\" + file, "r"))
            all_file_names.extend(file_names)
            all_files.extend(files)
        return all_file_names, all_files

    def get_tracked_files(self, dirs):
        """
        Read all tracked.csv files of multiple directories, ignore csv files without "tracked" in name.
        :param dirs: List of directories.
        :return: List with file names and pandas objects.
        """
        all_file_names = []
        all_files = []
        for dir in dirs:
            files = []
            file_names = []
            for file in os.listdir(dir):
                if file.endswith("csv") and "tracked" in file:
                    files.append(pd.read_csv(dir + "\\" + file))
                    file_names.append(file)
            all_file_names.extend(file_names)
            all_files.extend(files)
        return all_file_names, all_files