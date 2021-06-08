"""
@author: Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Collect cell sizes per coverslip in .LOG file for SPTAnalyser.
"""


import os
import pandas as pd
import sys
import configparser


class IncorrectConfigException(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)


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

    for dir in dirs: 
        files, file_areas = read_areas(dir)
        write_log(files, file_areas, dir + "\\cell_sizes.LOG")


if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except IndexError:
        print("Usage: python write_TS_macro.py your_config_file.ini")
