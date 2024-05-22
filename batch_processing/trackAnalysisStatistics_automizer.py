"""
@author: Alexander Niedrig, Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
Analyzes tracks and filters data, runs background correction, as well as rearranes it like the trackStatistics notebook does
"""
import configparser
import os
import shutil
import sys
import time
import warnings

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

def delete_empty_directories(directories):
    for directory in directories:
        for dirpath, dirnames, filenames in os.walk(directory, topdown=False):
            for dirname in dirnames:
                folder_path = os.path.join(dirpath, dirname)
                if not os.listdir(folder_path):
                    os.rmdir(folder_path)

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
        from pySPT.notebookspy import trackStatistics_noGUI as trackStatistics
        
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
    try:
        min_dynamic_D = config["DIFFUSION_TYPE_PARAM"]["min_detectable_D"]
    except KeyError:
        raise IncorrectConfigException("Parameter min_D missing in config.")
    try:
        min_tra_len = config["DIFFUSION_TYPE_PARAM"]["min_track_length"]
    except KeyError:
        raise IncorrectConfigException("Parameter min_track_length missing in config.")

    try:
        id_type = config["ID_TYPE"]["id_type"]
    except KeyError:
        raise IncorrectConfigException("Parameter id_type missing in config.")
    try:
        filter_min_traj = config["FILTER_PARAMETERS"]["min_trajectory_len"]
    except KeyError:
        raise IncorrectConfigException("Parameter min_trajectory_len missing in config.")
    try:
        filter_max_traj = config["FILTER_PARAMETERS"]["max_trajectory_len"]
    except KeyError:
        raise IncorrectConfigException("Parameter max_trajectory_len missing in config.")
    try:
        filter_min_D = config["FILTER_PARAMETERS"]["min_D"]
    except KeyError:
        raise IncorrectConfigException("Parameter min_D missing in config.")
    try:
        max_D = config["FILTER_PARAMETERS"]["max_D"]
    except KeyError:
        raise IncorrectConfigException("Parameter max_D missing in config.")
    try:
        filter_immobile = config["FILTER_PARAMETERS"]["filter_for_immobile"]
    except KeyError:
        raise IncorrectConfigException("Parameter filter_for_immobile missing in config.")
    try:
        filter_confined = config["FILTER_PARAMETERS"]["filter_for_confined"]
    except KeyError:
        raise IncorrectConfigException("Parameter filter_for_confined missing in config.")
    try:
        filter_free = config["FILTER_PARAMETERS"]["filter_for_free"]
    except KeyError:
        raise IncorrectConfigException("Parameter filter_for_free missing in config.")
    try:
        filter_noType = config["FILTER_PARAMETERS"]["filter_for_noType"]
    except KeyError:
        raise IncorrectConfigException("Parameter filter_for_noType missing in config.")
    try:
        save_dir = config["SAVE_DIR"]["save_dir"]
    except KeyError:
        raise IncorrectConfigException("Parameter save missing in config.")
        
    # resets tracking directories
    try:
        os.mkdir(save_dir)
    except FileExistsError:
        raise IncorrectConfigException("Saving directory already exists, aborting... Please state a new directory.")
    try:
        os.mkdir(save_dir + '\\trackAnalysis')
        if len(background) != 0:
            os.mkdir(save_dir + '\\trackAnalysis\\backgrounds')
        os.mkdir(save_dir + '\\trackAnalysis\\cells')
    except FileExistsError:
        print("directory " + save_dir + '\\trackAnalysis exists and will be overwritten')
        shutil.rmtree(save_dir + '\\trackAnalysis')
        os.mkdir(save_dir + '\\trackAnalysis')
        if len(background) != 0:
            os.mkdir(save_dir + '\\trackAnalysis\\backgrounds')
        os.mkdir(save_dir + '\\trackAnalysis\\cells')

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
            print("Analysing " + str(directories.index(dir) + 1) + "/" + str(len(directories)) + " : " + dir)
            notebook_analysis = trackAnalysis.analysisNotebook(software, mask_words, dir,
                                                               cslog, pixel_size,
                                                               background_size,
                                                               camera_integration_time, n_points, MSD_f_area, dof,
                                                               min_dynamic_D,
                                                               min_tra_len,
                                                               id_type)
            notebook_analysis.run_analysis()
            notebook_analysis.save_analysis()
        for file in get_matching_files(dir, '.h5', 'background'):
            try:
                shutil.move(file, save_dir + '\\trackAnalysis\\cells') 
            except shutil.Error:
                print(
                    'File ' + file + ' already exists in the trackAnalysis directory. Please check for duplicate files. Skipping...')

    # runs analysis for given backgrounds
    if len(background) > 0:
        for bg in background:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                print("Analysing " + str(background.index(bg) + 1) + "/" + str(len(background)) + "  " + bg)
                notebook_analysis = trackAnalysis.analysisNotebook(software, mask_words, bg, '', pixel_size,
                                                                   background_size,
                                                                   camera_integration_time, n_points, MSD_f_area, dof,
                                                                   min_dynamic_D,
                                                                   min_tra_len,
                                                                   id_type)
                notebook_analysis.run_analysis()
                notebook_analysis.save_analysis()
            for file in get_matching_files(bg, '.h5', 'cell'):
                try:
                    shutil.move(file, save_dir + '\\trackAnalysis\\backgrounds')
                except shutil.Error:
                    print(
                        'File ' + file + ' already exists in the trackAnalysis\\background directory. Please check for duplicate files. Skipping...')
                    pass
    print("running statistic filtering")
    if len(background) == 0:
        notebook_statistics = trackStatistics.statisticsNotebook(save_dir + '\\trackAnalysis\\cells',
                                                                 '', filter_min_traj, filter_max_traj,
                                                                 filter_min_D, max_D, filter_immobile,
                                                                 filter_confined, filter_free, filter_noType)

    else:
        notebook_statistics = trackStatistics.statisticsNotebook(save_dir + '\\trackAnalysis\\cells',
                                                                 save_dir + '\\trackAnalysis\\backgrounds',
                                                                 filter_min_traj,
                                                                 filter_max_traj,
                                                                 filter_min_D, max_D, filter_immobile,
                                                                 filter_confined, filter_free, filter_noType)
    notebook_statistics.save(save_dir)
    delete_empty_directories(directories)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except KeyError:
        print("Usage: python trackAnalysisStatistics.py your_config_file.ini")
