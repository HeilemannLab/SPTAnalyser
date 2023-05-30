"""
@author: Alexander Niedrig, Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Determines cell sizes, expected noise rate, precision and diffraction limit.
"""

import os
import sys
import shutil
import configparser
from pySPT.notebooks import precision
from pySPT.notebooks import diffractionLimit
import pandas as pd
from pySPT.notebooks import exp_noise_rate


class IncorrectConfigException(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)

def get_matching_files(directory, target, exclusion_String):
    """
    Search in directory and its subdirectories for files with target in their name and specific file ending.

    :param directory: all files in this directory & subdirectories are checked
    :param target: displacement or p_Bleach
    :param ending: File ending
    :return: List of matching file paths
    """
    matching_files = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if target in name.lower():
                if exclusion_String not in name.lower():
                    matching_files.append(os.path.join(path, name))
    return matching_files

#functions to write the cell sizes

def read_areas(directory):
    files = os.listdir(directory)
    files = [file for file in files if file.endswith(".csv")]
    file_areas = [int(pd.read_csv(directory + "\\" + file)["Area"]) for file in files]
    return files, file_areas


def write_log(files, file_areas, save_log):
    f = open(save_log, "w+")
    f.write("Calibration:\n")
    f.write('"X:", 1, "Y:", 1, "Units:", "", "Gray Units:", ""\n')
    f.write('"Image Name", "Area"\n')

    for file, file_area in zip(files, file_areas):
        if file[-9:-4] == "_size":
            file = file[:-9] + ".csv"
        f.write('"' + file + '"' + ", " + str(file_area) + '\n')

    f.close()
    print(save_log + " is saved.")


def main(config_path):
    dirs = []
    config = configparser.ConfigParser()
    config.sections()
    config.read(config_path)
    try:
        if len([key for key in config["INPUT_DIRS"]]):
            for key in config["INPUT_DIRS"]:
                dirs.append(config["INPUT_DIRS"][key])
        else:
            raise IncorrectConfigException("No input directory defined in config.")
    except KeyError:
        raise IncorrectConfigException("Section INPUT_DIRS missing in config.")
    try:
        software = config["SOFTWARE"]["software"]
        if software == ("ThunderSTORM" or "rapidSTORM"):
            pass
        else:
            raise IncorrectConfigException("Invalid Software")
    except KeyError:
        raise IncorrectConfigException("Parameter software missing in config.")
    try:
        pixel_size = config["CAMERA"]["pixel_size"]
    except KeyError:
        raise IncorrectConfigException("Parameter pixel_size missing in config.")
    try:
        pixel_per_row = config["CAMERA"]["pixel_per_row"]
    except KeyError:
        raise IncorrectConfigException("Parameter pixel_per_row missing in config.")
    try:
        camera_integration_time = config["CAMERA"]["camera_integration_time"]
    except KeyError:
        raise IncorrectConfigException("Parameter camera_integration_time missing in config.")
    try:
        background_size = config["CAMERA"]["background_size"]
    except KeyError:
        raise IncorrectConfigException("Parameter background_size missing in config.")


    for dir in dirs:
        if os.path.exists(dir + "\\PreAnalysisParameters"): #resets previous calculation of PreAnalysisParameters and recreates foldes
            shutil.rmtree(dir + "\\PreAnalysisParameters")
        os.mkdir(dir + "\\PreAnalysisParameters")
        os.mkdir(dir + "\\PreAnalysisParameters\\parameters")
        files, file_areas = read_areas(dir + "\\cells\\rois")   #writes cell_sizes.log for directories
        write_log(files, file_areas, dir + "\\cells\\rois\\cell_sizes.LOG")
        #runs virtual precision notebooks
        precison_nb= precision.precisionNotebook(dir, software,pixel_size,camera_integration_time)
        precison_nb.run_analysis()
        precison_nb.save_analysis()
        #runs virtual exp_noise notebook
        noiseRate_nb = exp_noise_rate.noiseRateNotebook(dir,software,background_size)
        noiseRate_nb.run_analysis()
        noiseRate_nb.save_analysis()
        #runs virtual diffraction_limit notebook
        diffLimit_nb = diffractionLimit.diffLimitNotebook(dir,software,pixel_size,pixel_per_row)
        diffLimit_nb.run_analysis()
        diffLimit_nb.save_analysis()
        #changes dir for post calculation data sorting, copies the important files to the parameters subfolder and renames the exp_noise_rate
        dir = dir + "\\PreAnalysisParameters"
        shutil.copy(dir + "\\precision\\precisions.txt",dir +"\\parameters")
        noiseRate = get_matching_files(dir + "\\exp_noise_rate", "exp_noise_rate","mean")
        shutil.copy(noiseRate[0], dir+"\\parameters")
        os.rename(get_matching_files(dir+'\\parameters', "exp_noise_rate","mean")[0],dir+"\\parameters\\exp_noise_rate.txt")
        shutil.copy(dir + "\\diff_limit_nn\\min_nearest_neighbor_distances.csv", dir + "\\parameters")



if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except IndexError:
        print("Usage: python BatchPreAnalysis.py your_config_file.ini")
