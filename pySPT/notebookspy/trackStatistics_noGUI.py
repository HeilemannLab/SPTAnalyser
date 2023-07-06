"""
@author: Alexander Niedrig, Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

virtual python only version of the trackStatistics.ipynb notebook for batch processing analysis
"""
import os
import warnings

from tqdm import tqdm_notebook as tqdm

from pySPT.Analysis import coverSlip
from pySPT.Analysis import trajectoryStatistics_noGUI as trajectoryStatistics
from pySPT.widgets import loadHdf5
from pySPT.widgets import saveStatistics
from pySPT.widgets import widgetLoadHdf5
from pySPT.widgets.widgetNotebooks import init_filter_notebook
from pySPT.widgets.widgetNotebooks import init_save_filtered_analysis
from pySPT.widgets.widgetNotebooks import init_save_filtered_trc
from pySPT.widgets.widgetNotebooks import init_save_track_stats


class statisticsNotebook():

    def __init__(self, directory, background, traj_min, traj_max, d_min, d_max, immobile, confined, free,
                 determination):
        print(directory)
        warnings.filterwarnings("ignore")
        self.widget_load_hdf5 = widgetLoadHdf5.WidgetLoadHdf5(bin_size="0.1", x_range="2",
                                                              y_range="0.5")  # adjust the default parameters
        self.load_hdf5 = loadHdf5.LoadHdf5()
        self.cover_slip = coverSlip.CoverSlip()
        self.track_stats = trajectoryStatistics.TrajectoryStatistics()

        # initializes files:
        self.widget_load_hdf5.dir_box.value = directory
        self.widget_load_hdf5.change_dir_box(directory)

        init_filter_notebook(self.cover_slip, self.widget_load_hdf5, self.load_hdf5, is_cell=True)
        self.track_stats.cells = self.cover_slip.cells
        self.track_stats.cell_trajectories = self.cover_slip.cell_trajectories
        self.track_stats.cell_sizes = [cell.size for cell in self.cover_slip.cells]
        self.track_stats.create_filtered_framework()

        # initializes background if given:

        self.widget_load_hdf5.dir_box_bg.value = background
        self.widget_load_hdf5.change_dir_box_bg(background)

        init_filter_notebook(self.cover_slip, self.widget_load_hdf5, self.load_hdf5, is_cell=False)
        self.track_stats.backgrounds = self.cover_slip.backgrounds
        self.track_stats.background_trajectories = self.cover_slip.background_trajectories
        self.track_stats.bg_sizes = [background.size for background in self.cover_slip.backgrounds]

        # filter options
        self.widget_load_hdf5.create_clear_output()
        self.widget_load_hdf5.min_length_box.value = traj_min
        self.widget_load_hdf5.max_length_box.value = traj_max
        self.widget_load_hdf5.min_D_box.value = d_min
        self.widget_load_hdf5.max_D_box.value = d_max
        if immobile == 'false':
            self.widget_load_hdf5.immob_type_check_box.value = False
        else:
            self.widget_load_hdf5.immob_type_check_box.value = True
        if confined == 'false':
            self.widget_load_hdf5.confined_type_check_box.value = False
        else:
            self.widget_load_hdf5.confined_type_check_box.value = True
        if free == 'false':
            self.widget_load_hdf5.free_type_check_box.value = False
        else:
            self.widget_load_hdf5.free_type_check_box.value = True
        if determination == 'false':
            self.widget_load_hdf5.analyse_not_successful_check_box.value = False
        else:
            self.widget_load_hdf5.analyse_not_successful_check_box.value = True

        self.widget_load_hdf5.bin_size_box.value = '0.1'
        self.widget_load_hdf5.MSD_y_lim.value = '0.5'

        self.track_stats.run_statistics(self.widget_load_hdf5.min_length_box.value,
                                        self.widget_load_hdf5.max_length_box.value,
                                        self.widget_load_hdf5.min_D_box.value, self.widget_load_hdf5.max_D_box.value,
                                        self.widget_load_hdf5.immob_type_check_box.value,
                                        self.widget_load_hdf5.confined_type_check_box.value,
                                        self.widget_load_hdf5.free_type_check_box.value,
                                        self.widget_load_hdf5.analyse_not_successful_check_box.value)
        self.track_stats.run_diffusion_histogram(self.widget_load_hdf5.bin_size_box.value,
                                                 self.widget_load_hdf5.MSD_delta_t_n.value,
                                                 self.widget_load_hdf5.MSD_y_lim.value, 0)

    def calc(self):
        self.widget_load_hdf5.create_clear_output()
        cells = []
        sigmas = []
        for cell_idx in range(len(self.track_stats.cell_trajectories_filtered)):
            print("{}: {} \u03BCm".format(self.track_stats.cells[cell_idx].name, self.track_stats.sigma_dyns[cell_idx]))
            cells.append(self.track_stats.cells[cell_idx].name)
            sigmas.append(self.track_stats.sigma_dyns[cell_idx])
        return cells, sigmas

    def save(self, save_folder):
        # sets up save location
        self.widget_load_hdf5.save_dir_box.value = save_folder
        self.widget_load_hdf5.change_save_dir_box(save_folder)
        # saves
        self.widget_load_hdf5.create_clear_output()
        h5_stats = saveStatistics.SaveStatistics()
        names = [i.name for i in self.track_stats.cells]
        if len(names) != len(set(names)):
            print(
                "Error: A target file name existed multiple times. Please make sure to name each file individually to ensure proper analysis and saving. Currently loaded file names:",
                names)
        else:
            if os.path.exists(
                    self.widget_load_hdf5.save_dir_box.value + "\\" + self.widget_load_hdf5.save_folder_name_box.value):
                print("Directory already exsists. Please choose another directory or folder name.")
            else:
                os.makedirs(
                    self.widget_load_hdf5.save_dir_box.value + "\\" + self.widget_load_hdf5.save_folder_name_box.value)
                if self.widget_load_hdf5.save_dir_box.value:
                    init_save_track_stats(h5_stats, self.track_stats, self.widget_load_hdf5.save_dir_box.value,
                                          self.widget_load_hdf5.save_folder_name_box.value,
                                          self.widget_load_hdf5.save_name_box.value)
                    if self.widget_load_hdf5.filtered_dataset_checkbox.value:
                        self.track_stats.filter_cell_trc()
                        for cell_index in tqdm(range(len(self.track_stats.cells))):
                            if self.track_stats.cell_trajectories_filtered[cell_index]:
                                init_save_filtered_analysis(self.cover_slip, cell_index, self.track_stats,
                                                            self.widget_load_hdf5.save_dir_box.value,
                                                            self.widget_load_hdf5.save_folder_name_box.value)
                        print("The filtered dataset is saved.")
                    if self.widget_load_hdf5.hmm_trc_checkbox.value:
                        print("The trc files of the filtered dataset will be saved.")
                        self.track_stats.filter_cell_trc()
                        for cell_index in tqdm(range(len(self.track_stats.cells))):
                            if self.track_stats.cell_trajectories_filtered[cell_index]:
                                init_save_filtered_trc(self.track_stats, self.widget_load_hdf5.save_dir_box.value,
                                                       self.widget_load_hdf5.save_folder_name_box.value)
                    if self.widget_load_hdf5.Dplot_checkbox.value and self.track_stats.diff_fig:
                        self.track_stats.save_diff_fig(self.widget_load_hdf5.save_dir_box.value,
                                                       self.widget_load_hdf5.save_folder_name_box.value)
