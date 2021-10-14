"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Widget handling for diffractionimit.ipynb.
"""

import tkinter as tk
import os
import tkinter.filedialog as fd
from ipywidgets import widgets
from IPython.display import clear_output


class WidgetDiffLimit():
    def __init__(self, px_size, n_px):
        self.software_button = self.create_software_button()
        self.got_dir = False
        self.dir_name = ""
        self.dir_button = self.create_dir_button()
        self.dir_box = self.create_dir_box()
        self.got_dir_bg = False
        self.dir_name_bg = ""
        self.dir_button_bg = self.create_dir_button()
        self.dir_box_bg = self.create_dir_box()
        self.roi_path = ""
        self.got_roi_path = False
        self.roi_text_box = self.create_roi_box()
        self.px_size_box = self.create_px_size_box(px_size)
        self.n_px_box = self.create_n_px_box(n_px)
        self.roi_button = self.create_roi_button()
        self.run_button = self.create_run_button()
        self.save_button = self.create_save_button()
        self.save_fig_checkbox = self.create_save_fig_checkbox()
        self.got_dir_save = False
        self.box_foldername = self.create_box_foldername()
        self.dir_button_save = self.create_dir_button()
        self.dir_box_save = self.create_dir_box(placeholder="Directory to save")

    def create_software_button(self):
        button = widgets.RadioButtons(
                options=["ThunderSTORM", "rapidSTORM"],
                disabled=False)
        return button

    def create_dir_button(self):
        button = widgets.Button(
            description="browse",
            disabled=False,
            button_style="",
            tooltip="browse for directory")
        return button

    def create_dir_box(self, val="", desc="Directory", placeholder="directory to be searched in"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder=placeholder, description=str(desc), disabled=False, style=style)
        return text

    def open_dir(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.attributes("-topmost", True)
        root.lift()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title="Please select a directory")
        root.update()
        root.destroy()
        self.dir_name = root.name
        self.dir_box.value = self.dir_name
        self.got_dir = True

    def change_dir_box(self, change):
        self.dir_name = self.dir_box.value
        self.got_dir = True

    def open_dir_bg(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.attributes("-topmost", True)
        root.lift()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title="Please select a directory")
        root.update()
        root.destroy()
        self.dir_name_bg = root.name
        self.dir_box_bg.value = self.dir_name_bg
        self.got_dir_bg = True

    def change_dir_box_bg(self, change):
        self.dir_name_bg = self.dir_box_bg.value
        self.got_dir_bg = True

    def determine_suffix(self):
        """
        Depending on the chosen software, the file ending of the target files differs.
        """
        if self.software_button.value == "ThunderSTORM":
            return ".csv"
        elif self.software_button.value == "rapidSTORM":
            return ".txt"

    def create_roi_button(self):
        button = widgets.Button(
            description="browse",
            disabled=False,
            button_style="",
            tooltip="browse for file")
        return button

    def open_file(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.attributes("-topmost", True)
        root.lift()
        root.name = fd.askopenfilename(title="Import .roi file",
                                       filetypes=(("roi logs", "*.log"), ("all files", "*.*")))
        self.roi_path = root.name
        root.update()
        root.destroy()
        self.roi_text_box.value = self.roi_path
        if os.path.isfile(self.roi_path):
            self.got_roi_name = True

    def create_roi_box(self, val="", desc="Complete path"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="insert path of roi.log", description=str(desc), disabled=False,
                            style=style)
        return text

    def create_px_size_box(self, val, desc="Pixel size [nm]"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="pixel size in nm", description=str(desc), disabled=False,
                            style=style)
        return text

    def create_n_px_box(self, val, desc="Number of pixels per row"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="nmber of pixels per row", description=str(desc), disabled=False,
                            style=style)
        return text

    def change_roi_box(self, change):
        self.roi_path = self.roi_text_box.value
        if os.path.isfile(self.roi_path):
            self.got_roi_path = True

    def create_background_size_box(self, val="65536", desc="Area of camera in px²"):
        """
        Box for inserting the amount of pixel on the camera.
        """
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="insert area", description=str(desc), disabled=False, style=style)
        return text

    def create_run_button(self):
        button = widgets.Button(
                description="run",
                disabled=False,
                button_style="",
                tooltip="run the analysis")
        return button

    def create_clear_output(self):
        clear_output()

    def create_save_button(self):
        button = widgets.Button(
            description="save",
            disabled=False,
            button_style="",
            tooltip="save the results")
        return button

    def create_save_fig_checkbox(self):
        check_box = widgets.Checkbox(
                value=True,
                description="Save plot",
                disabled=False)
        return check_box

    def open_dir_save(self, b):
        root = tk.Tk()
        root.withdraw()
        root.update()
        root.attributes("-topmost", True)
        root.lift()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title="Please select a directory")
        root.update()
        root.destroy()
        self.dir_name_save = root.name
        self.dir_box_save.value = self.dir_name_save
        self.got_dir_save = True

    def change_dir_box_save(self, change):
        self.dir_name_save = self.dir_box_save.value
        self.got_dir_save = True

    def create_save_box(self, val="", desc="Complete path"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="insert path", description=str(desc), disabled=False,
                            style=style)
        return text

    def create_box_foldername(self, val="diff_limit_nn", desc="Foldername"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="name of folder", description=str(desc),
                            disabled=False, style=style)
        return text
