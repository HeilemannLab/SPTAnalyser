"""
@author: Johanna Rahm, Claudia Catapano, Annabelle Brauel
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Write fiji macro to speed up ROI selection:
- define multiple paths of CS directories
- define save directory for macro
- creates arrays for all *.tif files, dl files (and optionally files containing a specific string)
- automatically selects the polygon selection tool, z-projects the raw data, opens the ROI manager and B&C settings
- after accepting a ROI it is renamed, the area measured in px² and both automatically saved in the rois folder of the respective CS directories
"""


import os
import sys
from datetime import datetime
import configparser
import re


class IncorrectConfigException(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)


def write_macro_str(save_path, all_tifs, all_tif_names, all_dls, all_dl_names, all_strs, all_str_names, all_save_rois, all_save_areas, zoom_factor):
    f = open(save_path, "w+")
    f.write('requires("1.54f");\n')
    f.write('run("Fresh Start");')
    f.write('run("Brightness/Contrast...");')
    f.write('setOption("ExpandableArrays", true);\n\n')
    f.write('run("Close All");\n')
    f.write('roiManager("reset");\n')
    f.write('run("ROI Manager...");\n')
    f.write('run("Set Measurements...", "area redirect=None decimal=9");\n\n')
    f.write("target_cells = newArray;\n")

    for c, tif in enumerate(all_tifs):
        f.write('target_cells[' + str(c) + '] = "' + tif + '"\n')
        
    f.write("\n")    
    f.write("cell_tif_names = newArray;\n")
    
    for c, names in enumerate(all_tif_names):
        f.write('cell_tif_names[' + str(c) + '] = "' + names + '"\n')
    
    f.write("\n")
    f.write("target_dls = newArray;\n")
    
    for c, dlnames in enumerate(all_dls):
        f.write('target_dls[' + str(c) + '] = "' + dlnames + '"\n')
        
    f.write("\n")
    f.write("cell_dl_names = newArray;\n")
    
    for c, dls in enumerate(all_dl_names):
        f.write('cell_dl_names[' + str(c) + '] = "' + dls + '"\n')
           
    f.write("\n")
    f.write("target_strs = newArray;\n")
    
    for c, strs in enumerate(all_strs):
        f.write('target_strs[' + str(c) + '] = "' + strs + '"\n')
        
    f.write("\n")
    f.write("cell_str_names = newArray;\n")
    
    for c, strnames in enumerate(all_str_names):
        f.write('cell_str_names[' + str(c) + '] = "' + strnames + '"\n')
    
    f.write("\n")
    f.write("save_rois = newArray;\n")
    
    for c, roisaves in enumerate(all_save_rois):
        f.write('save_rois[' + str(c) + '] = "' + roisaves + '"\n')
        
    f.write("\n")
    f.write("save_areas = newArray;\n")
    
    for c, areasaves in enumerate(all_save_areas):
        f.write('save_areas[' + str(c) + '] = "' + areasaves + '"\n')
    
    f.write("\n")
    f.write("for (i=0; i<target_cells.length ;i++){\n")
    f.write("\tif (File.exists(save_rois[i])) {\n")
    f.write('\t\tprint("ROI file exists for " + cell_tif_names[i] + ". Skipped processing.");\n')
    f.write("\t\tcontinue;\n")
    f.write("\t}\n")
    f.write("\topen(target_cells[i]);\n")
    for _ in range(zoom_factor):
        f.write('\trun("In [+]");\n')
    f.write('\trun("Z Project...", "projection=[Sum Slices]");\n')
    f.write('\trun("Enhance Contrast", "saturated=0.35");\n')
    f.write('\trun("Enhance Contrast", "saturated=0.35");\n')
    for _ in range(zoom_factor):
        f.write('\trun("In [+]");\n')   
    f.write("\topen(target_dls[i]);\n")
    for _ in range(zoom_factor):
        f.write('\trun("In [+]");\n')
    f.write("\topen(target_strs[i]);\n")
    for _ in range(zoom_factor):
        f.write('\trun("In [+]");\n')    
    f.write('\trun("Tile");\n')
    f.write('\trun("Set Scale...", "distance=0 known=0 unit=pixel global");\n')
    f.write('\tsetTool("polygon");\n')
    f.write('\twaitForUser("Draw ROI","Press \'t\' to save ROI and click OK to continue. \\n \\nOnly the last ROI in the ROI manager will be saved.");\n')
    f.write('\tcount = roiManager("count");\n')
    f.write("\tif (count-1 >= 0) {\n")
    f.write('\t\troiManager("select", count-1);\n')
    f.write('\t\troiManager("Rename", File.getName(save_rois[i]));\n')
    f.write('\t\troiManager("Save", save_rois[i]);\n')
    f.write("\t\t} else {\n")
    f.write("\t\t\tdo {\n")
    f.write('\t\t\t\twaitForUser("Draw ROI","No ROI was added to the ROI manager. Please draw a ROI before continuing. \\n \\nPress \'t\' to save ROI and click OK to continue");\n')
    f.write('\t\t\t\tcount = roiManager("count");\n')
    f.write("\t\t\t\tif (count - 1 >= 0) {\n")
    f.write('\t\t\t\t\troiManager("select", count - 1);\n')
    f.write('\t\t\t\t\troiManager("Rename", File.getName(save_rois[i]));\n')
    f.write('\t\t\t\t\troiManager("Save", save_rois[i]);\n')
    f.write("\t\t\t\t} else {\n")
    f.write("\t\t\t\t}\n")
    f.write("\t\t\t} while (count - 1 < 0);\n")
    f.write("\t\t}\n")
    f.write('\troiManager("Select", count-1);\n')
    f.write('\trun("Measure");\n')
    f.write('\tsaveAs("Results", save_areas[i]);\n')
    f.write('\trun("Clear Results");\n')
    f.write('\trun("Close All");\n')
    f.write('\troiManager("reset");\n')
    f.write("}\n")
    f.write("cleanUp();\n")
    f.write("function cleanUp() {\n")
    f.write('\tif (isOpen("Results")) {\n')
    f.write('\t\tselectWindow("Results"); \n')
    f.write('\t\trun("Close" );\n')
    f.write("\t}\n")
    f.write('\tif (isOpen("Log")) {\n')
    f.write('\t\tselectWindow("Log");\n')
    f.write('\t\trun("Close" );\n')
    f.write("\t}\n")
    f.write('\tif (isOpen("ROI Manager")) {\n')
    f.write('\t\tselectWindow("ROI Manager");\n')
    f.write('\t\trun("Close" );\n')
    f.write("\t}\n")
    f.write("\twhile (nImages()>0) {\n")
    f.write("\t\tselectImage(nImages());\n")
    f.write('\t\trun("Close");\n')
    f.write("\t}\n")
    f.write("}\n")
    f.write('waitForUser("Congrats, you\'re done!", "All cells in the defined directories are processed.");')
    
def write_macro(save_path, all_tifs, all_tif_names, all_dls, all_dl_names, all_save_rois, all_save_areas, zoom_factor):
    f = open(save_path, "w+")
    f.write('requires("1.54f");\n')
    f.write('run("Fresh Start");')
    f.write('run("Brightness/Contrast...");')
    f.write('setOption("ExpandableArrays", true);\n\n')
    f.write('run("Close All");\n')
    f.write('roiManager("reset");\n')
    f.write('run("ROI Manager...");\n')
    f.write('run("Set Measurements...", "area redirect=None decimal=9");\n\n')
    f.write("target_cells = newArray;\n")

    for c, tif in enumerate(all_tifs):
        f.write('target_cells[' + str(c) + '] = "' + tif + '"\n')
        
    f.write("\n")    
    f.write("cell_tif_names = newArray;\n")
    
    for c, names in enumerate(all_tif_names):
        f.write('cell_tif_names[' + str(c) + '] = "' + names + '"\n')
    
    f.write("\n")
    f.write("target_dls = newArray;\n")
    
    for c, dlnames in enumerate(all_dls):
        f.write('target_dls[' + str(c) + '] = "' + dlnames + '"\n')
        
    f.write("\n")
    f.write("cell_dl_names = newArray;\n")
    
    for c, dls in enumerate(all_dl_names):
        f.write('cell_dl_names[' + str(c) + '] = "' + dls + '"\n')
               
    f.write("\n")
    f.write("save_rois = newArray;\n")
    
    for c, roisaves in enumerate(all_save_rois):
        f.write('save_rois[' + str(c) + '] = "' + roisaves + '"\n')
        
    f.write("\n")
    f.write("save_areas = newArray;\n")
    
    for c, areasaves in enumerate(all_save_areas):
        f.write('save_areas[' + str(c) + '] = "' + areasaves + '"\n')
    
    f.write("\n")
    f.write("for (i=0; i<target_cells.length ;i++){\n")
    f.write("\tif (File.exists(save_rois[i])) {\n")
    f.write('\t\tprint("ROI file exists for " + cell_tif_names[i] + ". Skipped processing.");\n')
    f.write("\t\tcontinue;\n")
    f.write("\t}\n")
    f.write("\topen(target_cells[i]);\n")
    for _ in range(zoom_factor):
        f.write('\trun("In [+]");\n')
    f.write('\trun("Z Project...", "projection=[Sum Slices]");\n')
    f.write('\trun("Enhance Contrast", "saturated=0.35");\n')
    f.write('\trun("Enhance Contrast", "saturated=0.35");\n')
    for _ in range(zoom_factor):
        f.write('\trun("In [+]");\n') 
    f.write("\topen(target_dls[i]);\n")
    for _ in range(zoom_factor):
        f.write('\trun("In [+]");\n')
    f.write('\trun("Tile");\n')
    f.write('\trun("Set Scale...", "distance=0 known=0 unit=pixel global");\n')
    f.write('\tsetTool("polygon");\n')
    f.write('\twaitForUser("Draw ROI","Press \'t\' to save ROI and click OK to continue. \\n \\nOnly the last ROI in the ROI manager will be saved.");\n')
    f.write('\tcount = roiManager("count");\n')
    f.write("\tif (count-1 >= 0) {\n")
    f.write('\t\troiManager("select", count-1);\n')
    f.write('\t\troiManager("Rename", File.getName(save_rois[i]));\n')
    f.write('\t\troiManager("Save", save_rois[i]);\n')
    f.write("\t\t} else {\n")
    f.write("\t\t\tdo {\n")
    f.write('\t\t\t\twaitForUser("Draw ROI","No ROI was added to the ROI manager. \\n \\nPlease draw a ROI before continuing. \\n \\nPress \'t\' to save ROI and click OK to continue");\n')
    f.write('\t\t\t\tcount = roiManager("count");\n')
    f.write("\t\t\t\tif (count - 1 >= 0) {\n")
    f.write('\t\t\t\t\troiManager("select", count - 1);\n')
    f.write('\t\t\t\t\troiManager("Rename", File.getName(save_rois[i]));\n')
    f.write('\t\t\t\t\troiManager("Save", save_rois[i]);\n')
    f.write("\t\t\t\t} else {\n")
    f.write("\t\t\t\t}\n")
    f.write("\t\t\t} while (count - 1 < 0);\n")
    f.write("\t\t}\n")
    f.write('\troiManager("Select", count-1);\n')
    f.write('\trun("Measure");\n')
    f.write('\tsaveAs("Results", save_areas[i]);\n')
    f.write('\trun("Clear Results");\n')
    f.write('\trun("Close All");\n')
    f.write('\troiManager("reset");\n')
    f.write("}\n")
    f.write("cleanUp();\n")
    f.write("function cleanUp() {\n")
    f.write('\tif (isOpen("Results")) {\n')
    f.write('\t\tselectWindow("Results"); \n')
    f.write('\t\trun("Close" );\n')
    f.write("\t}\n")
    f.write('\tif (isOpen("Log")) {\n')
    f.write('\t\tselectWindow("Log");\n')
    f.write('\t\trun("Close" );\n')
    f.write("\t}\n")
    f.write('\tif (isOpen("ROI Manager")) {\n')
    f.write('\t\tselectWindow("ROI Manager");\n')
    f.write('\t\trun("Close" );\n')
    f.write("\t}\n")
    f.write("\twhile (nImages()>0) {\n")
    f.write("\t\tselectImage(nImages());\n")
    f.write('\t\trun("Close");\n')
    f.write("\t}\n")
    f.write("}\n")
    f.write('waitForUser("Congrats, you\'re done!", "All cells in the defined directories are processed.");')

def get_files(dir, dir_extension, file_ending, ignore_str):
    file_path = dir + dir_extension
    file_names = []
    for i in os.listdir(file_path):
        if i.endswith(file_ending) and not any([x in i.lower() for x in ignore_str]):
            file_names.append(i)
    file_names.sort(key=lambda x: [int(c) if c.isdigit() else c.lower() for c in re.split('(\d+)', x)])
    files = [dir + dir_extension + i for i in file_names]
    return files, file_names 
    
def get_dls(dir, dir_extension, dl_ending, ignore_str):
    dl_path = dir + dir_extension
    dl_names = []
    for i in os.listdir(dl_path):
        if i.endswith(dl_ending) and not any([x in i.lower() for x in ignore_str]):
            dl_names.append(i)
    dl_names.sort(key=lambda x: [int(c) if c.isdigit() else c.lower() for c in re.split('(\d+)', x)])
    dls = [dir + dir_extension + i for i in dl_names]
    return dls, dl_names 
    
def get_strs(dir, dir_extension, str_ending, ignore_str, add_str):
    str_path = dir + dir_extension
    str_names = []
    for i in os.listdir(str_path):
        if i.endswith(str_ending) and not any([x in i.lower() for x in ignore_str]):
            if any([x in i.lower() for x in add_str]):
                str_names.append(i)
    str_names.sort(key=lambda x: [int(c) if c.isdigit() else c.lower() for c in re.split('(\d+)', x)])
    strs = [dir + dir_extension + i for i in str_names]
    return strs, str_names 

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
        zoom_factor = int(config["PARAMETERS"]["zoom_factor"])
    except KeyError:
        zoom_factor = 0  # Default zoom factor if not specified in the config file
    try:
        ignore_str = config["PARAMETERS"]["ignore_str"].strip("][").replace(" ", "").split(",")
    except KeyError:
        raise IncorrectConfigException("Parameter ignore_str missing in config.")
    try:
        activate_str = config["PARAMETERS"].getboolean("activate_str")
    except KeyError:
        activate_str = False  # Set default value if activate_str is not found in the config file
    try:
        add_str = config["PARAMETERS"]["add_str"].strip("][").replace(" ", "").split(",")
    except KeyError:
        raise IncorrectConfigException("Parameter add_str missing in config.")
    try:
        macro_directory = config["SAVE_DIR"]["macro_directory"]
    except KeyError:
        raise IncorrectConfigException("No save directory defined in config.")
    try:
        macro_name = config["SAVE_DIR"]["macro_name"]
    except KeyError:
        macro_name = datetime.today().strftime('%Y-%m-%d_%Hh-%Mm') + "_ROI-macro"
    
    all_tifs, all_tif_names, all_dls, all_dl_names, all_strs, all_str_names, all_save_rois, all_save_areas = [], [], [], [], [], [], [], []
    
    for dir in dirs:
        if ignore_str != [""] and ignore_str != []:
            ignore_strs_tl = ["_dl", "_tl"] + ignore_str + add_str
            ignore_strs_cells = ["_cell_"] + ignore_str + add_str
        else:
            ignore_strs_tl = ["_dl", "_tl"]
            ignore_strs_cells = ["_cell_"]
            
        tif_files, tif_names = get_files(dir, r"\\cells\\tifs\\", ".tif", ignore_strs_tl)  # tif files with "_dl" will be ignored
        dl_files, dl_names = get_dls(dir, r"\\cells\\tifs\\", ".tif", ignore_strs_cells)
        str_files, str_names = get_strs(dir, r"\\cells\\tifs\\", ".tif", ignore_str, add_str) 
        save_rois = [dir + r"\\cells\\rois\\" + os.path.splitext(i)[0] + ".roi" for i in tif_names]
        save_areas = [dir + r"\\cells\\rois\\" + os.path.splitext(i)[0] + ".csv" for i in tif_names]
        all_tifs.extend(tif_files)
        all_tif_names.extend(tif_names)
        all_dls.extend(dl_files)
        all_dl_names.extend(dl_names)
        all_strs.extend(str_files)
        all_str_names.extend(str_names)
        all_save_rois.extend(save_rois)
        all_save_areas.extend(save_areas)
        
    final_macro_name = "ROI-macro_" + macro_name + "_" + datetime.today().strftime('%Y-%m-%d_%Hh-%Mm') + ".ijm"
    
    if activate_str:
        write_macro_str(macro_directory + "\\" + final_macro_name, all_tifs, all_tif_names, all_dls, all_dl_names, all_strs, all_str_names, all_save_rois, all_save_areas, zoom_factor)
    else:
        write_macro(macro_directory + "\\" + final_macro_name, all_tifs, all_tif_names, all_dls, all_dl_names, all_save_rois, all_save_areas, zoom_factor)
     
    print("Macro saved successfully at ", macro_directory + "\\" + final_macro_name)


if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except IndexError:
        print("Usage: python write_ROI_macro.py your_config_file.ini")
