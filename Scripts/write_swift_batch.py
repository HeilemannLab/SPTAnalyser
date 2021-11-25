"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Write batch file for swift with set parameters for multiple files.
"""


import os
import sys
import configparser
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


def main(config_path):
    file_paths = []
    n_file_paths = []
    precision_vals = []
    exp_noise_rate_vals = []
    config = configparser.ConfigParser()
    config.sections()
    config.read(config_path)

    try:
        if len([key for key in config["INPUT_DIRS"]]):
            for key in config["INPUT_DIRS"]:
                file_paths_dir = get_loc_files(config["INPUT_DIRS"][key])
                n_file_paths.append(len(file_paths_dir))
                file_paths.extend(file_paths_dir)
        else:
            raise IncorrectConfigException("No input directory defined in config.")
    except KeyError:
        raise IncorrectConfigException("Section INPUT_DIRS missing in config.")

    try:
        exp_displacement = config["PARAMETERS_GLOBAL"]["exp_displacement"]
    except KeyError:
        raise IncorrectConfigException("Parameter exp_displacement missing in config.")

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
        p_bleach = config["PARAMETERS_GLOBAL"]["p_bleach"]
    except KeyError:
        raise IncorrectConfigException("Parameter p_bleach missing in config.")

    try:
        p_switch = config["PARAMETERS_GLOBAL"]["p_switch"]
    except KeyError:
        raise IncorrectConfigException("Parameter p_switch missing in config.")

    try:
        diffraction_limit = config["PARAMETERS_GLOBAL"]["diffraction_limit"]
    except KeyError:
        raise IncorrectConfigException("Parameter diffraction_limit missing in config.")

    try:
        for key in config["PRECISION_INDIVIDUAL"]:
            precision_vals.append(config["PRECISION_INDIVIDUAL"][key])
        assert(len(n_file_paths) == len(precision_vals))
    except KeyError:
        raise IncorrectConfigException("No precision defined in config.")
    except AssertionError:
        raise IncorrectConfigException("Number of input directories must equal number of precision entries.")

    try:
        for key in config["EXP_NOISE_RATE_INDIVIDUAL"]:
            exp_noise_rate_vals.append(config["EXP_NOISE_RATE_INDIVIDUAL"][key])
        assert(len(n_file_paths) == len(exp_noise_rate_vals))
    except KeyError:
        raise IncorrectConfigException("No exp_noise_rate defined in config.")
    except AssertionError:
        raise IncorrectConfigException("Number of input directories must equal number of exp noise rate entries.")

    try:
        batch_path = config["SAVE_BATCH"]["batch_path"]
    except KeyError:
        raise IncorrectConfigException("No save directory defined in config.")

    # precision = []
    # for i in range(len(n_file_paths)):
    #     if os.path.exists(precision_vals[i]):
    #         with open(precision_vals[i], "r") as f:
    #             first_line = f.readline()
    #             precision += first_line.strip('][\n').split(', ')
    #     else:
    #         precision += [precision_vals[i]] * n_file_paths[i]

    # exp_noise_rate = []
    # for i in range(len(n_file_paths)):
    #     if os.path.exists(exp_noise_rate_vals[i]):
    #         with open(exp_noise_rate_vals[i], "r") as f:
    #             first_line = f.readline()
    #             exp_noise_rate += first_line.strip('][\n').split(', ')
    #     else:
    #         exp_noise_rate += [exp_noise_rate_vals[i]] * n_file_paths[i]

    precision = []
    for i in range(len(n_file_paths)):
        if os.path.exists(precision_vals[i]):
            precision = pd.read_csv(precision_vals[i], sep=" ", encoding="latin").iloc[:, 1].to_list()
        else:
            precision += [precision_vals[i]] * n_file_paths[i]

    exp_noise_rate = []
    for i in range(len(n_file_paths)):
        if os.path.exists(exp_noise_rate_vals[i]):
            exp_noise_rate = pd.read_csv(exp_noise_rate_vals[i], sep=" ", encoding="latin").iloc[:, 2].to_list()
        else:
            exp_noise_rate += [exp_noise_rate_vals[i]] * n_file_paths[i]

    my_bat = open(batch_path, "w+")
    for c, file in enumerate(file_paths):
        my_bat.write(define_swft_command(file, tau, exp_displacement, exp_noise_rate[c],
                                         max_displacement, max_displacement_pp, p_bleach, p_switch,
                                         precision[c], diffraction_limit))
        my_bat.write('\n')
    my_bat.write('ECHO All swift jobs are done.')
    my_bat.write('\n')
    my_bat.write('PAUSE')
    my_bat.close()
    print("Batch file successfully saved at", batch_path)


if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except IndexError:
        print("Usage: python write_swift_batch.py your_config_file.ini")
