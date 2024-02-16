"""
@author: Alexander Niedrig, Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
Analyzes tracks, then determines the dynamic localization error and calculates the minimal detectable diffusion coefficient with it.
"""
import configparser
import os
import shutil
import sys
import time
import warnings

import pandas as pd


class IncorrectConfigException(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)


def get_matching_files(directory, target, exclusion_string):
    """
    Search in directory and its subdirectories for files with target in their name but missing the exclusion string

    :param directory: all files in this directory & subdirectories are checked
    :param target: substring of filename
    :param exclusion_string: excludes files with this string in the name
    :return: List of matching file paths
    """
    matching_files = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if target.lower() in name.lower():
                if exclusion_string not in name.lower():
                    matching_files.append(os.path.join(path, name))
    return matching_files


def main(config_path):
    start_time = time.time()
    config = configparser.ConfigParser()
    config.sections()
    config.read(config_path)
    directories = []
    
    try:
        spt_analyser_dir = config["SPT_ANALYSER_DIR"]["spt_analyser_dir"]
        sys.path.append(spt_analyser_dir)  # Add SPTAnalyser directory to Python path
        from pySPT.notebookspy import trackAnalysis_noGUI as trackAnalysis
        from pySPT.notebookspy import trackStatistics_noGUI as sigma_dyn
        
    except KeyError:
        raise IncorrectConfigException("Parameter SPT_ANALYSER_DIR missing in config.")



    try:
        software = config["SOFTWARE"]["software"]
    except KeyError:
        raise IncorrectConfigException("Parameter software missing in config.")

    try:
        if len([key for key in config["INPUT_DIRS"]]):
            for key in config["INPUT_DIRS"]:
                directories.append(config["INPUT_DIRS"][key])
        else:
            raise IncorrectConfigException("No condition directory defined in config.")
    except KeyError:
        raise IncorrectConfigException("Section INPUT_DIRS missing in config.")

    try:
        mask_words = config["MASK_WORDS"]["mask"]
    except KeyError:
        raise IncorrectConfigException("Parameter mask missing in config.")

    try:
        if config["BACKGROUND"]["background"].lower() == "true":
            background = []
            for key in config["BACKGROUND_DIRS"]:
                background.append(config["BACKGROUND_DIRS"][key])
        elif config["BACKGROUND"]["background"].lower() == "false":
            background = ""
        else:
            raise IncorrectConfigException("Parameter background faulty.")
    except KeyError:
        raise IncorrectConfigException("Parameter background missing in config.")
    try:
        pixel_size = config["CAMERA"]["pixel_size"]
    except KeyError:
        raise IncorrectConfigException("Parameter pixel_size missing in config.")
    try:
        camera_integration_time = config["CAMERA"]["camera_integration_time"]
    except KeyError:
        raise IncorrectConfigException("Parameter camera_integration_time missing in config.")
    try:
        background_size = config["CAMERA"]["background_size"]
    except KeyError:
        raise IncorrectConfigException("Parameter background_size missing in config.")
    try:
        n_points = config["DIFFUSION_TYPE_PARAM"]["number_of_points"]
    except KeyError:
        raise IncorrectConfigException("Parameter number_of_points missing in config.")
    try:
        MSD_f_area = config["DIFFUSION_TYPE_PARAM"]["MSD_fit_area"]
    except KeyError:
        raise IncorrectConfigException("Parameter MSD_fit_area missing in config.")
    try:
        dof = config["DIFFUSION_TYPE_PARAM"]["degree_of_freedom"]
    except KeyError:
        raise IncorrectConfigException("Parameter degree_of_freedom missing in config.")
    # minimum detectable D should be zero to calculate this value
    min_dynamic_D = "0"
    try:
        min_tra_len = config["DIFFUSION_TYPE_PARAM"]["min_track_length"]
    except KeyError:
        raise IncorrectConfigException("Parameter min_track_length missing in config.")

    try:
        id_type = config["ID_TYPE"]["id_type"]
    except KeyError:
        raise IncorrectConfigException("Parameter id_type missing in config.")
    try:
        filter_min_length = config["FILTER_PARAMETERS"]["min_trajectory_len"]
    except KeyError:
        raise IncorrectConfigException("Parameter min_trajectory_len missing in config.")
    try:
        filter_max_length = config["FILTER_PARAMETERS"]["max_trajectory_len"]
    except KeyError:
        raise IncorrectConfigException("Parameter max_trajectory_len missing in config.")
    filter_min_D = "0" # must be set to 0 to determine Dmin
    #try:
    #    filter_min_D = config["FILTER_PARAMETERS"]["min_D"]
    #except KeyError:
    #    raise IncorrectConfigException("Parameter min_D missing in config.")
    try:
        filter_max_D = config["FILTER_PARAMETERS"]["max_D"]
    except KeyError:
        raise IncorrectConfigException("Parameter max_D missing in config.")
    
    # all data neccessary for sigma dyn determination
    filter_immobile = 'true'
    filter_confined = 'true'
    filter_free = 'true'
    filter_noType = 'true'
    #try:
    #    filter_immobile = config["FILTER_PARAMETERS"]["filter_for_immobile"]
    #except KeyError:
    #    raise IncorrectConfigException("Parameter filter_for_immobile missing in config.")
    #try:
    #    filter_confined = config["FILTER_PARAMETERS"]["filter_for_confined"]
    #except KeyError:
    #    raise IncorrectConfigException("Parameter filter_for_confined missing in config.")
    #try:
    #    filter_free = config["FILTER_PARAMETERS"]["filter_for_free"]
    #except KeyError:
    #    raise IncorrectConfigException("Parameter filter_for_free missing in config.")
    #try:
    #    filter_noType = config["FILTER_PARAMETERS"]["filter_for_noType"]
    #except KeyError:
    #    raise IncorrectConfigException("Parameter filter_for_noType missing in config.")
    
    try:
        save_dir = config["SAVE_DIR"]["save"]
    except KeyError:
        raise IncorrectConfigException("Parameter save missing in config.")
    # resets tracking directories
    try:
        os.mkdir(save_dir)
    except FileExistsError:
        pass
    try:
        os.mkdir(save_dir + "\\trackAnalysis")
        if len(background) != 0:
            os.mkdir(save_dir + "\\trackAnalysis\\backgrounds")
        os.mkdir(save_dir + "\\trackAnalysis\\cells")
    except FileExistsError:
        print("directory " + save_dir + "\\trackAnalysis exists. Will be overwritten")
        shutil.rmtree(save_dir + "\\trackAnalysis")
        os.mkdir(save_dir + "\\trackAnalysis")
        if len(background) != 0:
            os.mkdir(save_dir + "\\trackAnalysis\\backgrounds")
        os.mkdir(save_dir + "\\trackAnalysis\\cells")

    # runs analysis for each given directory

    for dir in directories:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cslog = get_matching_files(dir, ".log", "=")
            if len(cslog) == 1:
                cslog = cslog[0]
            elif len(cslog) == 0:
                raise IncorrectConfigException("cell sizes .log file missing in " + dir)
            elif len(cslog) > 1:
                raise IncorrectConfigException("multiple .log files in " + dir)
            else:
                print("This should not happen.")
            print("Analysing " + str(directories.index(dir) + 1) + "/" + str(len(directories)) + "  " + dir)
            notebook_analysis = trackAnalysis.analysisNotebook(software, mask_words, dir,
                                                               cslog,
                                                               pixel_size, background_size,
                                                               camera_integration_time, n_points, MSD_f_area, dof,
                                                               min_dynamic_D,
                                                               min_tra_len,
                                                               id_type)
            notebook_analysis.run_analysis()
            notebook_analysis.save_analysis()
            for file in get_matching_files(dir, ".h5", "background"):
                try:
                    shutil.move(file, save_dir + "\\trackAnalysis\\cells")
                except shutil.Error:
                    print(
                        "File " + file + " already exists in the trackAnalysis\\cells directory. Please check for duplicate files. Skipping...")
    # runs analysis for given backgrounds
    for bg in background:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            print("Analysing " + str(background.index(bg) + 1) + "/" + str(len(background)) + "  " + bg)
            notebook_analysis = trackAnalysis.analysisNotebook(software, mask_words, bg, "", pixel_size,
                                                               background_size,
                                                               camera_integration_time, n_points, MSD_f_area, dof,
                                                               min_dynamic_D,
                                                               min_tra_len,
                                                               id_type)
            notebook_analysis.run_analysis()
            notebook_analysis.save_analysis()
            for file in get_matching_files(bg, ".h5", "cell"):
                try:
                    shutil.move(file, save_dir + "\\trackAnalysis\\backgrounds")
                except shutil.Error:
                    print(
                        "File " + file + " already exists in the trackAnalysis\\background directory. Please check for duplicate files. Skipping...")
                    pass
    # creates the statistics notebook depending on whether backgrounds are given or not
    if len(background) == 0:
        notebook_statistics = sigma_dyn.statisticsNotebook(save_dir + "\\trackAnalysis\\cells", "",
                                                           filter_min_length, filter_max_length, filter_min_D,
                                                           filter_max_D, filter_immobile,
                                                           filter_confined, filter_free, filter_noType)
    else:
        notebook_statistics = sigma_dyn.statisticsNotebook(save_dir + "\\trackAnalysis\\cells",
                                                           save_dir + "\\trackAnalysis\\backgrounds",
                                                           filter_min_length, filter_max_length, filter_min_D,
                                                           filter_max_D, filter_immobile,
                                                           filter_confined, filter_free, filter_noType)
    # gets two lists, one with the cell names one with the sigma_dyn and writes them into a dataframe
    cells, sigma_dyns = notebook_statistics.calc()
    sigma_dyn_frame = pd.DataFrame(data={"cell": cells, "Sigma Dyn [um]": sigma_dyns})
    # extracts sigma_dyn
    sigma_dyn_value = sigma_dyn_frame.iloc[:, 1].quantile([0.75]).item()
    print("Sigma = " + str(sigma_dyn_value))
    # calculates d_min
    d_min = sigma_dyn_value ** 2 / (float(dof) * float(camera_integration_time) * 4)
    print("D_min = " + str(d_min))
    # output as file
    sigma_dyn_frame.to_csv(save_dir + "\\sigma_dyn_per_cell.csv", index=False)
    out = pd.Series([sigma_dyn_value, d_min], ["Sigma Dyn [um]", "D_min [um^2s^-1]"])
    out.to_csv(save_dir + "\\D_min.csv", header=False)
    shutil.rmtree(save_dir + "\\trackAnalysis")
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except IndexError:
        print("Usage: python Dmin_autocalculator.py your_config_file.ini")
