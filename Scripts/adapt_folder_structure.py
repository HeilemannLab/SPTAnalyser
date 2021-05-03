"""
@author: Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Create a directory and matching subdirectories for SPTData analysis:
root -> cells -> tifs
root -> cells -> rois
root -> cells -> tracks
root -> background -> tifs
root -> background -> tracks
"""


import os
import sys
import shutil
import configparser


class IncorrectConfigException(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)


def create_folders(directory):
    main_folder_name = os.path.basename(directory)
    try:
        os.mkdir(directory + "\\cells")
        os.mkdir(directory + "\\background")
        os.mkdir(directory + "\\cells\\tracks")
        os.mkdir(directory + "\\cells\\rois")
        os.mkdir(directory + "\\cells\\tifs")
        os.mkdir(directory + "\\background\\tracks")
        os.mkdir(directory + "\\background\\tifs")
    except FileExistsError:
        print("These folders already exist. If you want to create them again, you have to delete the 'cells' and"
              " 'background' folders in", directory)


def get_matching_files(directory, target, ending=".tif", dl=False):
    """
    Search in directory and its subdirectories for files with target in their name and specific file ending.

    :param directory: all files in this directory & subdirectories are checked
    :param target: Either cell or background
    :param ending: File ending
    :param dl: if True get dl.tif images (transmitted light) else get .tif movies
    :return: List of matching file paths
    """
    matching_files = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if target in name and os.path.splitext(name)[1] == ending:  # "dl not in name"
                if dl:
                    if "dl" in name.lower():
                        matching_files.append(os.path.join(path, name))
                else:
                    if "dl" not in name.lower():
                        matching_files.append(os.path.join(path, name))
    return matching_files


def sort_files_target_folder(files, folder):
    """
    Move all files to another folder.

    :param files: List of files
    :param folder: Target folder
    """
    for f in files:
        shutil.move(f, folder)


def rename_files(folder, remove_str="_MMStack.ome.tif"):
    """
    Rename files in a folder cell_1_MMStack.ome.tif -> cell_1.tif
    """
    for f in os.listdir(folder):
        if remove_str in f:
            os.rename(folder + "\\" + f, folder + "\\" + f[:-len(remove_str)] + ".tif")


def remove_empty_folders(dir):
    """
    Go through items in directory, if empty folders, delete them.
    :param dir: Directory
    """
    for f in os.listdir(dir):
        if os.path.isdir(dir + "\\" + f) and len(os.listdir(dir + "\\" + f)) == 0:
            os.rmdir(dir + "\\" + f)


def main(config_path):
    config = configparser.ConfigParser()
    config.sections()
    config.read(config_path)

    measurement_directories = []
    try:
        if len([key for key in config["INPUT_DIRS"]]):
            for key in config["INPUT_DIRS"]:
                measurement_directories.append(config["INPUT_DIRS"][key])
        else:
            raise IncorrectConfigException("No input directory defined in config.")
    except KeyError:
        raise IncorrectConfigException("Section INPUT_DIRS missing in config.")
    try:
        remove_str = config["PARAMETERS"]["remove_str"]
    except KeyError:
        raise IncorrectConfigException("Parameter remove_str missing in config.")

    for measurement_directory in measurement_directories:
        # create subfolder structure
        create_folders(measurement_directory)
        # get file paths based on conditions
        files_cell_tif = get_matching_files(measurement_directory, target="cell", ending=".tif", dl=False)
        files_cell_dl = get_matching_files(measurement_directory, target="cell", ending=".tif", dl=True)
        files_cell_txt = get_matching_files(measurement_directory, target="cell", ending=".txt", dl=False)
        files_background_tif = get_matching_files(measurement_directory, target="background", ending=".tif", dl=False)
        files_background_txt = get_matching_files(measurement_directory, target="background", ending=".txt", dl=False)
        # move files to fitting subfolder
        sort_files_target_folder(files_cell_tif, measurement_directory + "\\cells\\tifs")
        sort_files_target_folder(files_cell_txt, measurement_directory + "\\cells\\tifs")
        sort_files_target_folder(files_cell_dl, measurement_directory + "\\cells\\tifs")
        sort_files_target_folder(files_background_tif, measurement_directory + "\\background\\tifs")
        sort_files_target_folder(files_background_txt, measurement_directory + "\\background\\tifs")
        # rename tif files (cell_MMStack.ome.tif -> cell.tif)
        rename_files(measurement_directory + "\\cells\\tifs", remove_str)
        rename_files(measurement_directory + "\\background\\tifs", remove_str)
        # delete empty folders
        remove_empty_folders(measurement_directory)
        print("Folder structures successfully adapted.")


if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except IndexError:
        print("Usage: python adapt_folder_structure.py your_config_file.ini")

