"""
@author: Alexander Niedrig, Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Creates a swift batch, runs it, and determines p_Bleach and exp_displacement for both,
then repeats with the new values until both have converged.
"""
import configparser
import os
import sys
import subprocess
import csv
from pySPT.notebookspy import expDisplacement_noGUI as expDisplacement
from pySPT.notebookspy import pBleach_noGUI as pBleach
import shutil
import numpy
import time
import pandas as pd


class IncorrectConfigException(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)


def get_loc_files(dir_path):
    """get localization files from directory"""
    files = []
    for file in os.listdir(dir_path):
        if file.endswith(".csv") and "tracked" not in file:
            files.append(dir_path + "\\" + file)
    return files


def write_statistics(output_path, dataset):
    """
    writes the contents of a statistics dictionary to a file
    """
    output = open(output_path, 'w', newline='')
    for i in dataset.keys():
        output.write(i + ',')
    output.write('\n')
    writer = csv.writer(output, delimiter=',')
    writer.writerows(zip(*dataset.values()))
    output.close()


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
            if target in name.lower():
                if exclusion_string not in name.lower():
                    matching_files.append(os.path.join(path, name))
    return matching_files


def define_swft_command(file_path, tau, exp_displacement, exp_noise_rate, max_displacement, max_displacement_pp,
                        p_bleach, p_switch, precision, diff_limit):
    """swft command per file"""
    file_command = 'swft.exe ' + file_path + ' --tau ' + '"' + tau + '" ' + '--out_values "all" '
    file_command += '--exp_displacement ' + '"' + exp_displacement + '"' + ' '
    file_command += '--exp_noise_rate ' + '"' + str(exp_noise_rate) + '"' + ' '
    file_command += '--max_displacement ' + '"' + max_displacement + '"' + ' '
    file_command += '--max_displacement_pp ' + '"' + max_displacement_pp + '"' + ' '
    file_command += '--p_bleach ' + '"' + p_bleach + '"' + ' '
    file_command += '--p_switch ' + '"' + p_switch + '"' + ' '
    file_command += '--precision ' + '"' + str(precision) + '"' + ' '
    file_command += '--diffraction_limit ' + '"' + diff_limit + '"' + ' '
    return file_command


def gather_Data(filepaths):
    """
    Collects the displacement and pBleach data of all notebook outputs in the given list of directories
    """
    #gathers a list of all relevant files
    displacement_files = []
    bleach_files = []
    for working_directory in filepaths:
        displacement_files += get_matching_files(working_directory, 'displacement', 'histogram')
        bleach_files += get_matching_files(working_directory, 'p_bleach', 'histogram')
    displacement_data = []
    bleach_data = []
    # reads the values and puts them into one list, then returns it
    for file in displacement_files:
        disp = float(open(file, 'r').readlines()[1].split('	')[0])
        displacement_data.append(disp)
    for file in bleach_files:
        bleach = float(open(file, 'r').readlines()[1].split('	')[0])
        bleach_data.append(bleach)
    return displacement_data, bleach_data


def write_swift(config_path, exp_displacement, p_bleach):
    """
    Effectively does what the main method of the write_swift_batch.py does,
    however the parameters exp_displacement and p_bleach are no longer read from file but given as arguments
    """
    file_paths = []
    n_file_paths = []
    precision_vals = []
    exp_noise_rate_vals = []
    config = configparser.ConfigParser()
    config.sections()
    config.read(config_path)
    """
    The script expects the files precision and exp_noise_rate to be in the folder
    \PreAnalysisParameter\parameters, unlike the write_swift_batch.py which allows for a directory to be specified.
    However this is where BatchPreAnalysis places them.
    """
    try:
        if len([key for key in config["INPUT_DIRS"]]):
            for key in config["INPUT_DIRS"]:
                file_paths_dir = get_loc_files(config["INPUT_DIRS"][key] + "\\cells\\tracks")
                n_file_paths.append(len(file_paths_dir))
                file_paths.extend(file_paths_dir)
                precision_vals.append(
                    config["INPUT_DIRS"][key] + '\\PreAnalysisParameters\\parameters\\precisions.txt')
                assert (len(n_file_paths) == len(precision_vals))
                exp_noise_rate_vals.append(
                    config["INPUT_DIRS"][
                        key] + '\\PreAnalysisParameters\\parameters\\exp_noise_rate.txt')
                assert (len(n_file_paths) == len(exp_noise_rate_vals))
        else:
            raise IncorrectConfigException("No input directory defined in config.")
    except KeyError:
        raise IncorrectConfigException("Section INPUT_DIRS missing in config.")

    try:
        tau = config["PARAMETERS_GLOBAL"]["tau"]
    except KeyError:
        raise IncorrectConfigException("Parameter tau missing in config.")

    try:
        max_displacement = config["PARAMETERS_GLOBAL"]["max_displacement"]
        max_displacement = str(float(max_displacement) * float(exp_displacement))
    except KeyError:
        raise IncorrectConfigException("Parameter max_displacement missing in config.")

    try:
        max_displacement_pp = config["PARAMETERS_GLOBAL"]["max_displacement_pp"]
        max_displacement_pp = str(float(max_displacement_pp) * float(exp_displacement))
    except KeyError:
        raise IncorrectConfigException("Parameter max_displacement_pp missing in config.")

    try:
        p_switch = config["PARAMETERS_GLOBAL"]["p_switch"]
    except KeyError:
        raise IncorrectConfigException("Parameter p_switch missing in config.")

    try:
        diffraction_limit = config["PARAMETERS_GLOBAL"]["diffraction_limit"]
    except KeyError:
        raise IncorrectConfigException("Parameter diffraction_limit missing in config.")

    try:
        batch_path = config["SAVE_DIRECTORY"]["save_path"] + '\\swift.bat'
    except KeyError:
        raise IncorrectConfigException("No save directory defined in config.")

    precision = []
    for i in range(len(n_file_paths)):
        if os.path.exists(precision_vals[i]):
            precision.extend(pd.read_csv(precision_vals[i], sep=" ", encoding="latin").iloc[:, 1].to_list())
        else:
            precision += [precision_vals[i]] * n_file_paths[i]

    exp_noise_rate = []
    for i in range(len(n_file_paths)):
        if os.path.exists(exp_noise_rate_vals[i]):
            exp_noise_rate.extend(pd.read_csv(exp_noise_rate_vals[i], sep=" ", encoding="latin").iloc[:, 2].to_list())
        else:
            exp_noise_rate += [exp_noise_rate_vals[i]] * n_file_paths[i]

    my_bat = open(batch_path, "w+")
    for c, file in enumerate(file_paths):
        my_bat.write(define_swft_command(file, tau, str(exp_displacement), exp_noise_rate[c],
                                         max_displacement, max_displacement_pp, str(p_bleach), p_switch,
                                         precision[c], diffraction_limit))
        my_bat.write('\n')
    my_bat.write('ECHO All swift jobs are done.')
    my_bat.close()
    print("Batch file successfully saved at", batch_path)


def converge(save_directory, config_path, software, values, exp_displacement, p_bleach, directories, iterationNumber,
             ini_k, cam_dt, n_points):
    """
    Recursively runs the steps to determine exp_displacement and p_bleach until they are converged
    """

    # removes tracked files; as '=' is not permitted in filepaths, it returns all files
    for working_directory in directories:
        for file in get_matching_files(working_directory, '.tracked', '='):
            os.remove(file)

    # Runs the creation of the .bat file and executes it

    write_swift(config_path, exp_displacement, p_bleach)
    batch_execution = subprocess.Popen(save_directory + '\\swift.bat', shell=True, stdin=subprocess.PIPE, text=True)
    while batch_execution.returncode is None:
        batch_execution.communicate('\n')  # ensures the powershell doesn't pause during execution
    batch_execution.wait()

    # Moves the new tracked files into it
    for working_directory in directories:
        for file in get_matching_files(working_directory, '.tracked', '='):
            shutil.move(file, working_directory + '\\converger')

    # creates the virtual notebooks for each of the tracked files and executes them
    for working_directory in directories:
        tracked = get_matching_files(working_directory, '.tracked', 'meta')
        for cell in tracked:
            notebook_expD = expDisplacement.displacementNotebook(cell=cell, directory=working_directory,
                                                                 software=software)
            notebook_expD.run_analysis()
            notebook_expD.save_analysis()

            notebook_pBleach = pBleach.bleachNotebook(cell=cell, directory=working_directory, software=software,
                                                      inital_k=ini_k, camera_integration_time=cam_dt,
                                                      number_of_points=n_points)
            notebook_pBleach.run_analysis()
            notebook_pBleach.save_analysis()

    # sets the values as the old ones and gathers the new output files
    exp_displacement_old = exp_displacement
    p_bleach_old = p_bleach
    displacement_data, bleach_data = gather_Data(directories)

    #determins the new values and writes them to the output list
    exp_displacement = numpy.mean(displacement_data)
    displacement_data += [exp_displacement, numpy.std(displacement_data),
                          (numpy.std(displacement_data) / numpy.sqrt(len(displacement_data)))]
    values['exp_displacement ' + str(iterationNumber + 1)] = displacement_data
    p_bleach = numpy.mean(bleach_data)
    bleach_data += [p_bleach, numpy.std(bleach_data),
                    (numpy.std(bleach_data) / numpy.sqrt(len(bleach_data)))]
    values['p_bleach ' + str(iterationNumber + 1)] = bleach_data

    # checks whether the values have converged, recursively executes itself if not
    exp_displacement = float(exp_displacement)
    p_bleach = float(p_bleach)
    if abs(exp_displacement - exp_displacement_old) < 0.001 * exp_displacement_old and abs(
            p_bleach - p_bleach_old) < 0.001 * p_bleach_old and iterationNumber < 10:
        print("values converged after " + str(iterationNumber + 1) + " iterations: exp_displacement=" + str(
            exp_displacement) + " p_bleach = " + str(p_bleach))
    elif iterationNumber >= 10:
        print("10 iterations reached")
    else:
        print("values not converged after " + str(iterationNumber + 1) + " iterations: exp_displacement=" + str(
            exp_displacement) + " p_bleach = " + str(
            p_bleach) + "\n Old Values: exp_displacement=" + str(exp_displacement_old) + " p_bleach = " + str(
            p_bleach_old))
        converge(save_directory, config_path, software, values, exp_displacement, p_bleach, directories,
                 iterationNumber + 1, ini_k, cam_dt, n_points)


def main(config_path):
    start_time = time.time()
    config = configparser.ConfigParser()
    config.sections()
    config.read(config_path)
    directories = []
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
        exp_displacement = float(config["PARAMETERS_GLOBAL"]["exp_displacement"])
    except KeyError:
        raise IncorrectConfigException("Parameter exp_displacement missing in config.")

    try:
        p_bleach = float(config["PARAMETERS_GLOBAL"]["p_bleach"])
    except KeyError:
        raise IncorrectConfigException("Parameter p_bleach missing in config.")

    try:
        initial_k = config["P_BLEACH_PARAMETERS"]["Initial_k"]
    except KeyError:
        raise IncorrectConfigException("Parameter Initial_k missing in config.")

    try:
        camera_integration_time = config["P_BLEACH_PARAMETERS"]["camera_integration_time"]
    except KeyError:
        raise IncorrectConfigException("Parameter camera_integration_time missing in config.")

    try:
        number_of_points = config["P_BLEACH_PARAMETERS"]["number_of_points"]
    except KeyError:
        raise IncorrectConfigException("Parameter number_of_points missing in config.")
    try:
        save_directory = config["SAVE_DIRECTORY"]["save_path"]
    except KeyError:
        raise IncorrectConfigException("Parameter save_path missing in config.")

    # prepares folders and output, resets working directory from old runs

    for coverslip in directories:
        try:
            shutil.rmtree(coverslip + "\\converger")
        except FileNotFoundError:
            pass
        os.mkdir(coverslip + "\\converger")

    #sets up the value per cell output

    values = {'Cell': []}
    for working_directory in directories:
        for cell in get_matching_files(working_directory + '\\cells\\tracks', 'cell', 'protocol'):
            name = cell.split('\\')[-1]
            name = name[:-4]
            values['Cell'].append(name)
    values['Cell'] += ['mean', 'SD', 'SEM']

    # starts up the recursive convergence process
    converge(save_directory, config_path, software, values, exp_displacement, p_bleach, directories, 0, initial_k,
             camera_integration_time, number_of_points)

    # writes output per cell
    write_statistics(save_directory + '\\value_history.csv', values)

    # prepares the output per coverslip

    values_per_cs = {'Coverslip': ['MEAN', 'SD', 'SEM']}
    for coverslip in directories:
        data_disp, data_bleach = gather_Data([coverslip])
        values_per_cs[coverslip.split('\\')[-1] + ' exp_displacement'] = [numpy.mean(data_disp), numpy.std(data_disp), (
                    numpy.std(data_disp) / numpy.sqrt(len(data_disp)))]
        values_per_cs[coverslip.split('\\')[-1] + ' p_Bleach'] = [numpy.mean(data_bleach), numpy.std(data_bleach),
                                                                  (numpy.std(data_bleach) / numpy.sqrt(
                                                                      len(data_bleach)))]
    write_statistics(save_directory + '\\mean_per_cs.csv', values_per_cs)

    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except IndexError:
        print("Usage: python DPB_converger.py your_config_file.ini")
