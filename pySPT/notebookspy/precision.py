"""
@author: Alexander Niedrig, Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

virtual python only version of the precision.ipynb notebook for batch processing analysis
"""
import os
from pySPT.widgets import widgetDirectoryStructure
from pySPT.widgets import widgetColumnSort
from pySPT.widgets import widgetPrecision
from pySPT.preAnalysis import precision_noGUI as precision

class precisionNotebook():
    def __init__(self,directory,software,pixel_size,camera_dt):
        self.widget_dir_structure = widgetDirectoryStructure.WidgetDirStructure()
        self.widget_precision = widgetPrecision.WidgetPrecision(pixel_size=str(pixel_size), camera_dt=str(camera_dt))  # adjust the default parameters
        self.widget_precision.dir_box.value = directory + "\\cells\\tracks"
        self.widget_precision.software_button.value = software
        self.widget_precision.dir_box_save.value = directory + "\\PreAnalysisParameters"
        self.precision = precision.Precision()
        self.widget_precision.save_fig_checkbox.value = False


    def run_analysis(self):
        self.widget_precision.create_clear_output()
        files, file_names = self.precision.get_loc_files(self.widget_precision.dir_box.value)
        self.precision.get_precisions(files, file_names)
        self.precision.analysis_executed = True

    # In[5]:

    def save_analysis(self):
        self.widget_precision.create_clear_output()
        if self.precision.analysis_executed:
            _, file_names = self.precision.get_loc_files(self.widget_precision.dir_box.value)
            self.precision.save_precision_list(self.widget_precision.dir_box_save.value + "\\" + self.widget_precision.box_foldername.value, self.precision.mean_values, self.widget_precision.save_fig_checkbox.value, file_names)
        else:
            print("Please run the analysis first, by clicking at the 'run' button.")

