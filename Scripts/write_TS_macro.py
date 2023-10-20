"""
@author: Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Write fiji-ThunderSTORM macro:
- define multiple paths of CS directories
- define path to background region of interest (= camera chip area)
- define photon intensity range for multi-emitter fit
- define save directory for macro
"""


import os
import sys
from datetime import datetime
import configparser


class IncorrectConfigException(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)


def write_macro(save_path, all_tifs, all_tif_names, all_rois, all_save_files, intensity_range):
    f = open(save_path, "w+")
    f.write('requires("1.45s");\n')
    f.write('setOption("ExpandableArrays", true);\n\n')
    f.write("target_cells = newArray;\n")

    for c, tif in enumerate(all_tifs):
        f.write('target_cells[' + str(c) + '] = "' + tif + '"\n')
        
    f.write("\n")    
    f.write("cell_tif_names = newArray;\n")
    
    for c, names in enumerate(all_tif_names):
        f.write('cell_tif_names[' + str(c) + '] = "' + names + '"\n')
    
    f.write("\n")
    f.write("target_rois = newArray;\n")
    
    for c, rois in enumerate(all_rois):
        f.write('target_rois[' + str(c) + '] = "' + rois + '"\n')
            
    f.write("\n")
    f.write("save_paths = newArray;\n")
    
    for c, saves in enumerate(all_save_files):
        f.write('save_paths[' + str(c) + '] = "' + saves + '"\n')
    
    f.write("\n")
    f.write("for (i=0; i<target_cells.length ;i++){\n")
    f.write("\topen(target_cells[i]);\n")
    f.write('\troiManager("Open", target_rois[i]);\n')
    f.write('\troiManager("Select", i);\n')
    f.write('\trun("Run analysis", "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector=[Local maximum] connectivity=8-neighbourhood threshold=std(Wave.F1) estimator=[PSF: Integrated Gaussian] sigma=1.6 fitradius=3 method=[Maximum likelihood] full_image_fitting=false mfaenabled=true keep_same_intensity=false nmax=3 fixed_intensity=true expected_intensity='+ intensity_range +' pvalue=1.0E-6 renderer=[Averaged shifted histograms] magnification=5.0 colorize=false threed=false shifts=2 repaint=50");\n')
    f.write('\trun("Show results table", "action=duplicates distformula=uncertainty_xy");\n')
    f.write('\trun("Export results", "floatprecision=1 filepath=" + save_paths[i] + " fileformat=[CSV (comma separated)] sigma=true intensity=true chi2=true offset=true saveprotocol=true x=true y=true bkgstd=true id=true uncertainty_xy=true frame=true");\n')
    f.write("\tselectWindow(cell_tif_names[i]);\n")
    f.write("\tclose();\n")
    f.write("}")


def get_files(dir, dir_extension, file_ending, ignore_str):
    file_path = dir + dir_extension
    file_names = []
    for i in os.listdir(file_path):
        if i.endswith(file_ending) and not any([x in i.lower() for x in ignore_str]):
            file_names.append(i)
    files = [dir + dir_extension + i for i in file_names]
    return files, file_names 


def convert_to_fiji_path(path):
    return path.replace("\\", "\\\\")


def main(config_path):
    dirs = []
    config = configparser.ConfigParser()
    config.sections()
    config.read(config_path)

    try:
        if len([key for key in config["INPUT_DIRS"]]):
            for key in config["INPUT_DIRS"]:
                dirs.append(convert_to_fiji_path(config["INPUT_DIRS"][key]))
        else:
            raise IncorrectConfigException("No input directory defined in config.")
    except KeyError:
        raise IncorrectConfigException("Section INPUT_DIRS missing in config.")
    try:
        background_roi = convert_to_fiji_path(config["PARAMETERS"]["background_roi"])
    except KeyError:
        raise IncorrectConfigException("Parameter background_roi missing in config.")
    try:
        intensity_range = config["PARAMETERS"]["intensity_range"]
    except KeyError:
        raise IncorrectConfigException("Parameter intensity_range missing in config.")
    try:
        ignore_str = config["PARAMETERS"]["ignore_str"].strip("][").replace(" ", "").split(",")
    except KeyError:
        raise IncorrectConfigException("Parameter ignore_str missing in config.")
    try:
        macro_directory = config["SAVE_DIR"]["macro_directory"]
    except KeyError:
        raise IncorrectConfigException("No save directory defined in config.")

    all_tifs, all_tif_names, all_rois, all_save_files = [], [], [], []
    for dir in dirs:
        if ignore_str != [""] and ignore_str != []:
            ignore_strs_tl = ["_dl", "_tl"] + ignore_str
            ignore_strs_size = ["_size"] + ignore_str
        else:
            ignore_strs_tl = ["_dl", "_tl"]
            ignore_strs_size = ["_size"]
            
        tif_files, tif_names = get_files(dir, r"\\cells\\tifs\\", ".tif", ignore_strs_tl)  # tif files with "_dl" will be ignored
        roi_files, _ = get_files(dir, r"\\cells\\rois\\", ".roi", ignore_strs_size)
        save_files = [dir + r"\\cells\\tracks\\" + os.path.splitext(i)[0] + ".csv" for i in tif_names]
        background_files, background_names = get_files(dir, r"\\background\\tifs\\", ".tif", ignore_strs_tl)
        save_background = [dir + r"\\background\\tracks\\" + os.path.splitext(i)[0] + ".csv" for i in background_names]
        all_tifs.extend(tif_files)
        all_tifs.extend(background_files)
        all_tif_names.extend(tif_names)
        all_tif_names.extend(background_names)
        all_rois.extend(roi_files)
        all_rois.extend([background_roi for i in range(len(background_files))])
        all_save_files.extend(save_files)
        all_save_files.extend(save_background)

    macro_name = datetime.today().strftime('%Y-%m-%d_%Hh-%Mm') + "_TS-macro.ijm"
    write_macro(macro_directory + "\\" + macro_name, all_tifs, all_tif_names, all_rois, all_save_files, intensity_range)
    print("Macro saved successfully at ", macro_directory + "\\" + macro_name)


if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except IndexError:
        print("Usage: python write_TS_macro.py your_config_file.ini")
