"""
@author: Alexander Niedrig, Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

virtual python only version of the diffraction_limit.ipynb notebook for batch processing analysis
"""
from pySPT.widgets import widgetDiffLimit
from pySPT.widgets import widgetColumnSort
from pySPT.preAnalysis import diffLimit_noGUI as diffLimit

class diffLimitNotebook():

    def __init__(self,directory,software,pixel_size,pixel_per_row):

        self.widget_diff_limit = widgetDiffLimit.WidgetDiffLimit(pixel_size, pixel_per_row)  # adjust the default parameters
        self.diff_limit = diffLimit.DiffLimit()
        self.widget_diff_limit.software_button.value = software
        self.widget_diff_limit.dir_name = directory + "\\cells\\tracks"
        self.widget_diff_limit.dir_box_save.value = directory + "\\swift_analysis_parameter"
        self.widget_diff_limit.box_foldername.value="diff_limit_nn"

    def run_analysis(self):
        self.widget_diff_limit.create_clear_output()
        self.diff_limit.clear_object()
        self.diff_limit.px_size, self.diff_limit.n_px, self.diff_limit.max_search_area = int(self.widget_diff_limit.px_size_box.value), int(self.widget_diff_limit.n_px_box.value), 100
        self.diff_limit.run_diff_limit(self.widget_diff_limit.dir_name, self.widget_diff_limit.determine_suffix())


    def save_analysis(self):
        self.widget_diff_limit.create_clear_output()
        self.diff_limit.save(self.widget_diff_limit.dir_box_save.value, self.widget_diff_limit.box_foldername.value, self.diff_limit.file_names, False)
        print("Results are saved at", self.widget_diff_limit.dir_box_save.value + "\\"  + self.widget_diff_limit.box_foldername.value)


