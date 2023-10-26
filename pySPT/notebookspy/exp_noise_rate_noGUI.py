"""
@author: Alexander Niedrig, Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

virtual python only version of the exp_noise_rate.ipynb notebook for batch processing analysis
"""
import os
from pySPT.widgets import widgetExpNoiseRate
from pySPT.widgets import widgetDirectoryStructure
from pySPT.widgets import widgetColumnSort
from pySPT.preAnalysis import expNoiseRate_noGUI as expNoiseRate

class noiseRateNotebook():

    def __init__(self,directory,software,background_size):

        self.widget_exp_noise_rate = widgetExpNoiseRate.WidgetExpNoiseRate(area_camera=background_size)  # adjust the default parameters
        self.widget_exp_noise_rate.software_button.value = software
        self.widget_dir_structure = widgetDirectoryStructure.WidgetDirStructure()
        self.exp_noise_rate = expNoiseRate.ExpNoiseRate()


        '''File handling'''
        self.widget_exp_noise_rate.dir_name = directory + "\\cells\\tracks"
        self.widget_exp_noise_rate.dir_name_bg = directory + "\\background\\tracks"
        self.widget_exp_noise_rate.roi_path = directory + "\\cells\\rois\\cell_sizes.LOG"
        self.widget_exp_noise_rate.dir_box_save.value = directory + "\\swift_analysis_parameter"
        self.widget_exp_noise_rate.box_foldername.value = "exp_noise_rate"

    def run_analysis(self):
        self.widget_exp_noise_rate.create_clear_output()
        self.exp_noise_rate.run_exp_noise_rate(self.widget_exp_noise_rate.dir_name, self.widget_exp_noise_rate.dir_name_bg, self.widget_exp_noise_rate.roi_path, self.widget_exp_noise_rate.determine_suffix(), int(self.widget_exp_noise_rate.background_size_box.value))

    def save_analysis(self):
        self.widget_exp_noise_rate.create_clear_output()
        self.exp_noise_rate.save_results(self.widget_exp_noise_rate.dir_box_save.value, self.widget_exp_noise_rate.box_foldername.value, False)
