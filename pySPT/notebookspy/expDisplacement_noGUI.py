"""
@author: Alexander Niedrig, Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

virtual python only version of the expDisplacement.ipynb notebook for batch processing analysis
"""
import sys
import os
from pySPT.widgets import widgetExpDisp
from pySPT.widgets import widgetDirectoryStructure
from pySPT.widgets import widgetColumnSort
from pySPT.preAnalysis import expDisplacement_noGUI as expDisplacement

class displacementNotebook():
    def __init__(self,directory,cell,software,):
        self.widget_exp_disp = widgetExpDisp.WidgetExpDisp()

        self.widget_exp_disp.software_button.value = software

        self.widget_exp_disp.file_text_box.value = cell
        self.exp_displacement = expDisplacement.ExpDisplacement()
        self.widget_dir_structure = widgetDirectoryStructure.WidgetDirStructure()
        self.widget_exp_disp.save_fig_checkbox.value = False

    def run_analysis(self):
        self.widget_exp_disp.create_clear_output()
        if self.widget_exp_disp.is_file(self.widget_exp_disp.file_text_box.value):
            widget_column_sort = widgetColumnSort.WidgetColumnSort(self.widget_exp_disp.file_text_box.value, self.widget_exp_disp.software_button.value, [('"seg.id"',), ('"seg.mjd"',), ('"seg.mjd_n"',), ('"seg.motion"',)])
            widget_column_sort.check_header()
            if widget_column_sort.correct_header:
                widget_column_sort.run_column_sort()
                self.exp_displacement.file_name = self.widget_exp_disp.file_text_box.value
                self.exp_displacement.software = self.widget_exp_disp.software_button.value
                self.exp_displacement.filter_immob = self.widget_exp_disp.filter_immobile_checkbox.value
                self.exp_displacement.column_order = widget_column_sort.column_order
                self.exp_displacement.run_exp_displacement()
            else:
                self.widget_exp_disp.warning_wrong_file()
        else:
            self.widget_exp_disp.warning_wrong_file_path()

    def save_analysis(self):
        self.widget_exp_disp.create_clear_output()
        self.widget_dir_structure.name_handling(self.widget_exp_disp.file_text_box.value)
        self.widget_dir_structure.create_raw_base_name()
        self.widget_dir_structure.sub_folder = "\\preAnalysis"
        self.widget_dir_structure.create_folder()
        self.exp_displacement.save_exp_displacement(self.widget_dir_structure.sub_folder_dir, self.widget_dir_structure.raw_base_name, False)

