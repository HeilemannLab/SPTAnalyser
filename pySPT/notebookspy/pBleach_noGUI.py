"""
@author: Alexander Niedrig, Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

virtual python only version of the p_Bleach.ipynb notebook for batch processing analysis
"""
import sys
import os
from pySPT.widgets import widgetDirectoryStructure
from pySPT.widgets import widgetColumnSort
from pySPT.widgets import widgetPBleach
from pySPT.preAnalysis import pBleach_noGUI as pBleach

class bleachNotebook():


    def __init__(self,directory,cell,software,inital_k,camera_integration_time,number_of_points,):
        self.widget_p_bleach = widgetPBleach.WidgetPBleach(k=inital_k, camera_dt=camera_integration_time, n_points=number_of_points)
        self.widget_p_bleach.software_button.value = software

        self.p_bleach = pBleach.PBleach()
        self.widget_dir_structure = widgetDirectoryStructure.WidgetDirStructure()
        self.widget_p_bleach.file_button.on_click(self.widget_p_bleach.open_file)
        self.widget_p_bleach.file_text_box.observe(self.widget_p_bleach.change_file_box)


        self.widget_dir_structure.name_handling(self.widget_p_bleach.file_name)
        self.widget_dir_structure.create_raw_base_name()

        self.widget_p_bleach.file_text_box.value = cell
    def run_analysis(self):
        self.widget_p_bleach.create_clear_output()
        if self.widget_p_bleach.got_file_name:
            widget_column_sort = widgetColumnSort.WidgetColumnSort(self.widget_p_bleach.file_text_box.value, self.widget_p_bleach.software_button.value, [('"seg.id"',), ('"seg.mjd_n"',)])
            widget_column_sort.check_header()
            if widget_column_sort.correct_header:
                widget_column_sort.run_column_sort()
                self.p_bleach.ignore_points = int(self.widget_p_bleach.ignore_points.value)
                self.p_bleach.file_name = self.widget_p_bleach.file_text_box.value
                self.p_bleach.software = self.widget_p_bleach.software_button.value
                self.p_bleach.column_order = widget_column_sort.column_order
                self.p_bleach.dt = float(self.widget_p_bleach.integration_time.value)
                self.p_bleach.init_k = float(self.widget_p_bleach.init_k.value)
                self.p_bleach.file_name = self.widget_p_bleach.file_text_box.value
                self.p_bleach.run_p_bleach()
            else:
                self.widget_p_bleach.warning_wrong_file()
        else:
            self.widget_p_bleach.warning_wrong_file_path()

    def save_analysis(self):
        self.widget_p_bleach.create_clear_output()
        self.widget_dir_structure.name_handling(self.widget_p_bleach.file_text_box.value)
        self.widget_dir_structure.create_raw_base_name()
        self.widget_dir_structure.sub_folder = "\\preAnalysis"
        self.widget_dir_structure.create_folder()
        self.p_bleach.save_fit_results(self.widget_dir_structure.sub_folder_dir, self.widget_dir_structure.raw_base_name)
        self.p_bleach.save_mjd_n_frequencies(self.widget_dir_structure.sub_folder_dir, self.widget_dir_structure.raw_base_name)

