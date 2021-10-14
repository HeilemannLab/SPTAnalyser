"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Handle widgets of transitionCounts jupyter notebook
"""

# TODO be able to load multiple directories
# TODO Add run widgets
# TODO Add parameter widgets
# TODO Add plot widgets
# TODO Add save widgets

import tkinter as tk
import os
from ipywidgets import widgets
import tkinter.filedialog as fd
from IPython.display import clear_output


class DefineDirectory():
    def __init__(self, description, value=""):
        """
        Define *.hdf5 files for analysis as center or neighbor, give short name to file.
        :param description: Short description of target file.
        :param title: Title of opened filedialog.
        :param filetype: File ending.
        :param idx: Idx of file.
        """
        self.got_dir = False
        self.dir_name = ""
        self.dir_box = self.create_dir_box(description, value)
        self.dir_button = self.create_dir_button()

    def create_dir_button(self):
        button = widgets.Button(
            description="browse",
            disabled=False,
            button_style="",
            tooltip="browse for directory")
        return button

    def open_dir(self, b):
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        root.lift()
        root.update()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title="Please select a directory")
        root.update()
        root.destroy()
        self.dir_name = root.name
        self.dir_box.value = self.dir_name
        self.got_dir = True

    def create_dir_box(self, description, val=""):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="directory to be searched in", description=str(description),
                            disabled=False, style=style)
        return text

    def change_dir_box(self, change):
        self.dir_name = self.dir_box.value
        self.got_dir = True


class RunAnalysis():
    def __init__(self):
        self.run_analysis_button = self.create_run_analysis_button()
        self.run_plot_button = self.create_run_plot_button()

    def create_run_analysis_button(self):
        button = widgets.Button(
            description="run & save",
            disabled=False,
            button_style="",
            tooltip="run analysis")
        return button

    def create_run_plot_button(self):
        button = widgets.Button(
            description="plot",
            disabled=False,
            button_style="",
            tooltip="display & save")
        return button

    def create_clear_output(self):
        clear_output()


class Parameter():
    def __init__(self, n_types, mask, counts_path, trajectory_path, save_dir, save_folder):
        self.n_diff_states_box = self.create_n_diff_states_box(n_types)
        self.mask_box = self.create_mask_box(mask)
        self.counts_file_box = self.create_counts_file_box(counts_path)
        self.trajectory_file_box = self.create_trajectory_file_box(trajectory_path)
        self.save_dir_box = self.create_save_dir_box(save_dir)
        self.save_folder_box = self.create_save_folder_box(save_folder)


    def create_n_diff_states_box(self, val=""):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="number of diffusion states", description="Diffusion states",
                            disabled=False, style=style)
        return text

    def create_mask_box(self, val=""):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="directory to save results", description="Mask value",
                            disabled=False, style=style)
        return text

    def create_counts_file_box(self, val=""):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="path to counts file", description="Counts file",
                            disabled=False, style=style)
        return text

    def create_trajectory_file_box(self, val=""):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="path to trajectory table file", description="Trajectory file",
                            disabled=False, style=style)
        return text

    def create_save_dir_box(self, val=""):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="directory to save results", description="Save dir",
                            disabled=False, style=style)
        return text

    def create_save_folder_box(self, val=""):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="name of results folder", description="Save folder",
                            disabled=False, style=style)
        return text