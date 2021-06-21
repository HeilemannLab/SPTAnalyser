"""
@author: Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Extract hdf5 file information.
"""


import os
import sys
import configparser
import numpy as np
import h5py
import math
import csv
from itertools import zip_longest


class IncorrectConfigException(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)


def get_files(dir_path, mask_names):
    files = []
    file_names = []
    for file in os.listdir(dir_path):
        print(file, mask_names)
        if file.endswith(".h5") and file not in mask_names:
            files.append(dir_path + "\\" + file)
            file_names.append(file)
    return files, file_names


def read_files(file_paths):
    h5_files = []
    for file in file_paths:
        h5_file = h5py.File(file, "r")
        h5_files.append(h5_file)
    return h5_files


def extract_values(files, folder, frame, idx):
    column_name = list(files[0][folder][frame].dtype.fields.keys())[idx]
    print("Target column: " + folder + " \u2192 " + frame + " \u2192 " + column_name)
    for file in files:
        rows = file[folder][frame][()]
        if np.shape(rows) == (1, 1):
            return get_values_single(files, folder, frame, idx)
        else:
            return get_values_column(files, folder, frame, idx)


def get_values_column(files, folder, frame, idx):
    all_values, mean_values, mean_errors_STD, mean_errors_SEM = [], [], [], []
    for file in files:
        values = []
        rows = file[folder][frame][()]
        for row in rows:
            values.append(row[idx])
        error_STD = np.std(values, ddof=1)
        error_SEM = np.std(values, ddof=1)/math.sqrt(len(values))
        mean = np.mean(values)
        all_values.append(values)
        mean_values.append(mean)
        mean_errors_STD.append(error_STD)
        mean_errors_SEM.append(error_SEM)
    return all_values, mean_values, mean_errors_STD, mean_errors_SEM


def get_values_single(files, folder, frame, idx):
    all_values, mean_values, mean_errors_STD, mean_errors_SEM = [], [], [], []
    for file in files:
        rows = file[folder][frame][()]
        for row in rows:
            mean_values.append(row[0][idx])
    mean_errors_STD = np.std(mean_values, ddof=1)
    mean_errors_SEM = np.std(mean_values, ddof=1)/math.sqrt(len(mean_values))
    return all_values, mean_values, mean_errors_STD, mean_errors_SEM


def save_mean_results(path, file_name, target_names, mean, STD, SEM):
    out_file_name = path+"\\"+file_name+"_mean.csv"
    header = "target name\tmean\tstandard deviation\tstandard error"
    max_name_length = max([len(i) for i in target_names])
    data = np.zeros(np.array(target_names).size, dtype=[("col1", "U" + str(max_name_length)),
                                                           ("col2", float), ("col3", float), ("col4", float)])
    data["col1"] = np.array(target_names)
    data["col2"] = np.array(mean)
    data["col3"] = np.array(STD)
    data["col4"] = np.array(SEM)
    np.savetxt(out_file_name, X=data, fmt=("%10s", "%.4e", "%.4e", "%.4e"), header=header)


def save_all_results(path, file_name, target_names, values):
    out_file_name = path+"\\"+file_name+"_values.csv"
    header = ",".join(target_names)+"\n"
    with open(out_file_name, "w+", newline="") as f:
        writer = csv.writer(f)
        f.write(header)
        for value in zip_longest(*values):
            writer.writerow(value)


def save_single_results(path, file_name, target_names, mean, STD, SEM):
    out_file_name = path+"\\"+file_name+"_single_values.csv"
    header = "target names, values\n"
    with open(out_file_name, "w+", newline="") as f:
        f.write(header)
        for name, value in zip(target_names, mean):
            f.write(name + ", " + str(value)+"\n")
        f.write("mean: " + str(np.mean(mean)) + "\n")
        f.write("standard deviation: " + str(STD) + "\n")
        f.write("standard error: " + str(SEM))


def main(config_path):
    config = configparser.ConfigParser()
    config.sections()
    config.read(config_path)

    try:
        path = config["INPUT_DIR"]["dir"]
    except KeyError:
        raise IncorrectConfigException("No input directory defined in config.")

    try:
        folder_name = config["PARAMETERS"]["folder_name"]
    except KeyError:
        raise IncorrectConfigException("Parameter folder_name missing in config.")

    try:
        file_name = config["PARAMETERS"]["file_name"]
    except KeyError:
        raise IncorrectConfigException("Parameter file_name missing in config.")

    try:
        column_idx = int(config["PARAMETERS"]["column_idx"])
    except KeyError:
        raise IncorrectConfigException("Parameter column_idx missing in config.")

    try:
        mask_names = config["PARAMETERS"]["mask_files"].strip("][").replace(" ", "").split(",")
    except KeyError:
        raise IncorrectConfigException("Parameter mask_files missing in config.")

    try:
        save_path = config["SAVE"]["save_dir"]
    except KeyError:
        raise IncorrectConfigException("Parameter save_dir missing in config.")

    try:
        save_name = config["SAVE"]["save_name"]
    except KeyError:
        raise IncorrectConfigException("Parameter save_name missing in config.")

    # path = r"T:\Chemie_phd\SPTAnalyser\test_data\trackStats"
    # save_path = r"T:\Chemie_phd\SPTAnalyser\test_data\trackStats"
    # save_name = "results"
    # folder_name = "statistics"  # "rossier" / statistics
    # file_name = "statistics_3"  # "rossierStatistics" statistics_3
    # column_idx = 8  # starting from 1
    # mask_names = ["Fab_CS2_cell05", "Fab_CS2_cell04.h5"]

    mask_names = [name + ".h5" if name[-3:] != ".h5" else name for name in mask_names]
    print(mask_names)
    file_paths, file_names = get_files(path, mask_names)
    h5_files = read_files(file_paths)
    all_values, mean_values, mean_errors_STD, mean_errors_SEM = extract_values(h5_files, folder_name, file_name, column_idx-1)
    if all_values:
        save_all_results(save_path, save_name, file_names, all_values)
        save_mean_results(save_path, save_name, file_names, mean_values, mean_errors_STD, mean_errors_SEM)
    else:
        save_single_results(save_path, save_name, file_names, mean_values, mean_errors_STD, mean_errors_SEM)


if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except IndexError:
        print("Usage: python write_TS_macro.py your_config_file.ini")
