"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Widget handling for hiddenMarkovModeling.ipynb
"""

import os
from ipywidgets import widgets
from IPython.display import clear_output
import tkinter as tk
import tkinter.filedialog as fd


class WidgetInitHMM():
    def __init__(self, dt=0.02, init_n_start=1, init_n_end=6, x_axis=300):
        self.got_dir = False
        self.dir_name = ""
        self.dir_button = self.create_dir_button()
        self.dir_box = self.create_dir_box()
        self.dt_box = self.create_dt_box(dt)
        self.init_n_start = self.create_init_n_start(init_n_start)
        self.init_n_end = self.create_init_n_end(init_n_end)
        self.run_button = self.create_run_button()
        self.n_states_box = self.create_n_states_box(init_n_start)
        self.x_axis_box = self.create_x_axis_box(x_axis)
        self.cell_options = []
        self.drop_down_cells = self.create_drop_down_cells()
        self.plot_button = self.create_plot_button()
        self.box_foldername = self.create_box_foldername()
        self.figure_type_box = self.create_figure_type_box()
        self.dir_button_save = self.create_dir_button()
        self.dir_box_save = self.create_dir_box(placeholder="Directory to save")
        self.save_button = self.create_save_button()

    def create_clear_output(self):
        clear_output()

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
        root.attributes("-topmost", True)
        root.lift()
        root.update()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title="Please select a directory")
        root.update()
        root.destroy()
        self.dir_name = root.name
        self.dir_box.value = self.dir_name
        self.got_dir = True

    def change_dir_box(self, change):
        self.dir_name = self.dir_box.value
        self.got_dir = True

    def create_drop_down_cells(self):
        drop_down_cells = widgets.Dropdown(
            options=self.cell_options,
            description="Cell:",
            disabled=False)
        return drop_down_cells

    def create_dt_box(self, val, desc="Camera integration time [s]"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Define time between two acquired frames",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_init_n_start(self, val, desc="Min number of states"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Define an integer", description=str(desc), disabled=False,
                            style=style)
        return text

    def create_init_n_end(self, val, desc="Max number of states"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Define an integer", description=str(desc), disabled=False,
                            style=style)
        return text

    def create_run_button(self):
        button = widgets.Button(
                description="run",
                disabled=False,
                button_style="",
                tooltip="run the analysis")
        return button

    def create_n_states_box(self, val, desc="Number of states"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="e.g. with best scores",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_x_axis_box(self, val, desc="x axis range [nm]"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="max range of x axis in nm",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_plot_button(self):
        button = widgets.Button(
                description="plot",
                disabled=False,
                button_style="",
                tooltip="plot results for chosen number of states")
        return button

    def create_save_button(self):
        button = widgets.Button(
            description="save",
            disabled=False,
            button_style="",
            tooltip="save the results")
        return button

    def open_dir_save(self, b):
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        root.lift()
        root.update()
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

    def create_box_foldername(self, val="hidden_markov_modeling", desc="Foldername"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="name of folder", description=str(desc),
                            disabled=False, style=style)
        return text

    def create_figure_type_box(self, desc="Figure format"):
        style = {"description_width": "initial"}
        text = widgets.Text(value="pdf", placeholder="png, pdf, svg ...",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_save_button(self):
        button = widgets.Button(
            description="save",
            disabled=False,
            button_style="",
            tooltip="save the results")
        return button


class WidgetHMM():
    def __init__(self, dt=0.02, n_states=2, stability=0.9, diffs="1002.8, 50800.10", weights="0.3, 0.7", min_len=4,
                 immob_threshold=25, epsilon=20, x_axis=300, graphviz_bin=r"C:\Program Files (x86)\Graphviz2.38\bin"):
        self.got_dir = False
        self.dir_name = ""
        self.dir_button = self.create_dir_button()
        self.dir_box = self.create_dir_box()
        self.got_graphviz_dir = False
        self.dir_got_graphviz_name = ""
        self.dir_graphviz_button = self.create_dir_button()
        self.dir_graphviz_box = self.create_dir_graphviz_box(graphviz_bin)
        self.dt_box = self.create_dt_box(dt)
        self.run_button = self.create_run_button()
        self.n_states_box = self.create_n_states_box(n_states)
        self.stability_box = self.create_stability_box(stability)
        self.min_len_box = self.create_min_len_box(min_len)
        self.D_box = self.create_D_box(diffs)
        self.w_box = self.create_w_box(weights)
        self.train_ws_checkbox = self.create_train_ws_checkbox()
        self.train_ds_checkbox = self.create_train_ds_checkbox()
        self.train_tps_checkbox = self.create_train_tps_checkbox()
        self.choose_weights_button = self.create_choose_weights_button()
        self.immob_threshold_box = self.create_immob_threshold_box(immob_threshold)
        self.immob_check_button = self.create_immob_check_button()
        self.epsilon_box = self.create_epsilon_box(epsilon)
        self.create_d_corr_button = self.create_d_corr_button()
        self.x_axis_box = self.create_x_axis_box(x_axis)
        self.cell_options = []
        self.drop_down_cells = self.create_drop_down_cells()
        self.plot_button = self.create_plot_button()
        self.box_foldername = self.create_box_foldername()
        self.figure_type_box = self.create_figure_type_box()
        self.dir_button_save = self.create_dir_button()
        self.dir_box_save = self.create_dir_box(placeholder="Directory to save")
        self.save_button = self.create_save_button()

    def create_clear_output(self):
        clear_output()

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
        root.attributes("-topmost", True)
        root.lift()
        root.update()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title="Please select a directory")
        root.update()
        root.destroy()
        self.dir_name = root.name
        self.dir_box.value = self.dir_name
        self.got_dir = True

    def change_dir_box(self, change):
        self.dir_name = self.dir_box.value
        self.got_dir = True

    def create_dt_box(self, val, desc="Camera integration time [s]"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Define time between two acquired frames",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_run_button(self):
        button = widgets.Button(
                description="run",
                disabled=False,
                button_style="",
                tooltip="run the analysis")
        return button

    def create_dir_graphviz_button(self):
        button = widgets.Button(
            description="browse",
            disabled=False,
            button_style="",
            tooltip="browse for bin folder")
        return button

    def change_dir_graphviz_box(self, change):
        self.dir_name = self.dir_box.value
        self.got_dir = True

    def open_graphviz_dir(self, b):
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        root.lift()
        root.update()
        root.name = fd.askdirectory(initialdir=os.getcwd(), title="Please select a directory")
        root.update()
        root.destroy()
        self.dir_graphviz_name = root.name
        self.dir_graphviz_box.value = self.dir_graphviz_name
        self.got_graphviz_dir = True

    def create_dir_graphviz_box(self, val, desc="Graphviz bin"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="define path to Graphviz bin folder",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_n_states_box(self, val, desc="Number of states"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="e.g. with best scores",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_stability_box(self, val, desc="Stability"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Probability of staying in the state per timestep",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_min_len_box(self, val, desc="Min trajectory length"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="trajectories < min length will be filtered out",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_D_box(self, val, desc="Initial diffusion coefficients [\u00B5m²/s]"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Insert values comma separated",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_w_box(self, val, desc="Initial weights"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Insert values comma separated",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_train_ws_checkbox(self, desc="Train weights"):
        text = widgets.Checkbox(value=True, description=str(desc), disabled=False, indent=False)
        return text

    def create_train_ds_checkbox(self, desc="Train diffusion coefficients"):
        text = widgets.Checkbox(value=True, description=str(desc), disabled=False, indent=False)
        return text

    def create_train_tps_checkbox(self, desc="Train transition probabilities"):
        text = widgets.Checkbox(value=True, description=str(desc), disabled=False, indent=False)
        return text

    def create_choose_weights_button(self, desc="Node size"):
        text = widgets.RadioButtons(options=["jdd fit weights", "occurrence", "state probabilities"], value="jdd fit weights",
                                    layout={"width": "max-content"}, description=desc, disabled=False)
        return text

    def create_immob_threshold_box(self, val, desc="Immobile threshold [nm]"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Insert values comma separated",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_immob_check_button(self):
        button = widgets.Button(
                description="run",
                disabled=False,
                button_style="",
                tooltip="check if apparent diffusion states by HMM are immobile")
        return button

    def create_epsilon_box(self, val, desc="Static error [nm]"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="Insert values comma separated",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_d_corr_button(self):
        button = widgets.Button(
                description="run",
                disabled=False,
                button_style="",
                tooltip="correct diffusion coefficients by dynamic and static errors")
        return button

    def create_x_axis_box(self, val, desc="x axis range [nm]"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="max range of x axis in nm",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_drop_down_cells(self):
        drop_down_cells = widgets.Dropdown(
            options=self.cell_options,
            description="Cell:",
            disabled=False)
        return drop_down_cells

    def create_plot_button(self):
        button = widgets.Button(
                description="plot",
                disabled=False,
                button_style="",
                tooltip="plot results for chosen number of states")
        return button

    def create_save_button(self):
        button = widgets.Button(
            description="save",
            disabled=False,
            button_style="",
            tooltip="save the results")
        return button

    def open_dir_save(self, b):
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        root.lift()
        root.update()
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

    def create_box_foldername(self, val="hidden_markov_modeling", desc="Foldername"):
        style = {"description_width": "initial"}
        text = widgets.Text(value=str(val), placeholder="name of folder", description=str(desc),
                            disabled=False, style=style)
        return text

    def create_figure_type_box(self, desc="Figure format"):
        style = {"description_width": "initial"}
        text = widgets.Text(value="pdf", placeholder="png, pdf, svg ...",
                            description=str(desc), disabled=False, style=style)
        return text

    def create_save_button(self):
        button = widgets.Button(
            description="save",
            disabled=False,
            button_style="",
            tooltip="save the results")
        return button
