"""
@author: Alexander Niedrig, Johanna Rahm, Claudia Catapano
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

virtual python only version of the TrackAnalysis.ipynb notebook for batch processing analysis
"""
from pySPT.Analysis import coverSlip
from pySPT.Analysis import trackAnalysis
from pySPT.widgets import widgetColumnSort
from pySPT.widgets import widgetDirectoryStructure
from pySPT.widgets import widgetTrackAnalysis
from pySPT.widgets.widgetNotebooks import init_save_track_analysis


class analysisNotebook():

    def __init__(self, software, mask_words, directory, roi, pixel_size, total_area, camera_integration_time,
                 number_of_points, fit_area, dof, min_D, min_track_length, id_type):

        self.widget_track_analysis = widgetTrackAnalysis.WidgetTrackAnalysis(pixel_size=pixel_size,
                                                                             area_camera=total_area,
                                                                             camera_dt=camera_integration_time,
                                                                             n_points_D=number_of_points,
                                                                             fit_area_MSD=fit_area, dof_D=dof,
                                                                             min_D=min_D, min_length=min_track_length,
                                                                             hmm_min_length="20", hmm_float="10",
                                                                             bin_size="0.1", x_range="2",
                                                                             y_range="0.5")  # adjust the default parameters
        self.track_analysis = trackAnalysis.TrackAnalysis()
        self.cover_slip = coverSlip.CoverSlip()
        self.widget_track_analysis.software_button.value = software
        self.widget_track_analysis.ignore_words_box.value = mask_words
        self.widget_track_analysis.dir_box.value = directory
        self.widget_track_analysis.roi_box.value = roi
        self.widget_track_analysis.trajectory_id_button.value = id_type
        self.widget_track_analysis.camera_pixel_size_box.value = pixel_size
        self.widget_track_analysis.camera_pixel_amount_box.value = total_area
        self.widget_track_analysis.camera_integration_time_box.value = camera_integration_time
        self.widget_track_analysis.change_dir_box('')
        self.widget_track_analysis.change_roi_box('')
        self.widget_track_analysis.points_D_fit_box.value = number_of_points
        self.widget_track_analysis.rossier_fit_area_box.value = fit_area
        self.widget_track_analysis.dof_box.value = dof
        self.widget_track_analysis.D_min_box.value = min_D
        self.widget_track_analysis.min_track_length_box.value = min_track_length
        self.widget_dir_structure = widgetDirectoryStructure.WidgetDirStructure()

    def run_analysis(self):
        self.widget_track_analysis.create_clear_output()
        self.widget_track_analysis.searchSubFolders(self.widget_track_analysis.dir_name)
        if self.widget_track_analysis.got_dir:
            self.cover_slip.rossier_fit_area, self.cover_slip.software, self.cover_slip.min_track_length_type, self.cover_slip.min_track_length_hmm, self.cover_slip.dt, self.cover_slip.pixel_size, self.cover_slip.pixel_amount, self.cover_slip.dof, self.cover_slip.D_min, self.cover_slip.roi_file, self.cover_slip.cell_files, self.cover_slip.points_fit_D, self.cover_slip.seg_id = self.widget_track_analysis.rossier_fit_area_box.value, self.widget_track_analysis.software_button.value, self.widget_track_analysis.min_track_length_box.value, self.widget_track_analysis.min_track_length_hmm_box.value, self.widget_track_analysis.camera_integration_time_box.value, self.widget_track_analysis.camera_pixel_size_box.value, self.widget_track_analysis.camera_pixel_amount_box.value, self.widget_track_analysis.dof_box.value, self.widget_track_analysis.D_min_box.value, self.widget_track_analysis.roi_name, self.widget_track_analysis.file_names, self.widget_track_analysis.points_D_fit_box.value, self.widget_track_analysis.trajectory_id_button.value
            for cell_idx in range(len(self.cover_slip.cell_files)):
                if self.widget_track_analysis.software_button.value != "PALMTracer":
                    if self.widget_track_analysis.software_button.value == "ThunderSTORM":
                        self.widget_column_sort = widgetColumnSort.WidgetColumnSort(
                            self.cover_slip.cell_files[cell_idx], self.widget_track_analysis.software_button.value,
                            [('"track.id"',), ('"x [nm]"',), ('"y [nm]"',), ('"frame"',), ('"intensity [photon]"',),
                             ('"seg.id"',)])
                    elif self.widget_track_analysis.software_button.value == "rapidSTORM":
                        self.widget_column_sort = widgetColumnSort.WidgetColumnSort(
                            self.cover_slip.cell_files[cell_idx], self.widget_track_analysis.software_button.value,
                            [('"track.id"',), ('"Position-0-0"',), ('"Position-1-0"',), ('"Amplitude-0-0"',),
                             ('"ImageNumber-0-0"',), ('"seg.id"',)])
                    self.widget_column_sort.check_header()
                    if self.widget_column_sort.correct_header:
                        self.widget_column_sort.run_column_sort()
                        self.cover_slip.column_orders.append(self.widget_column_sort.column_order)
            self.cover_slip.create_cells()
            self.track_analysis.cell_sizes = [cell.size for cell in self.cover_slip.cells]
            self.track_analysis.cell_trajectories = self.cover_slip.cell_trajectories
            self.track_analysis.run_statistics_no_filter()
        else:
            self.widget_track_analysis.warning_trc_file()
        self.widget_track_analysis.cells = self.cover_slip.cells

    def plot_diffusions(self):
        self.widget_track_analysis.create_clear_output()
        self.track_analysis.run_plot_diffusion_histogram(self.widget_track_analysis.bin_size_box.value,
                                                         self.widget_track_analysis.MSD_delta_t_n.value,
                                                         self.widget_track_analysis.MSD_y_lim.value)

    def save_analysis(self):
        for cell_index in range(0, len(self.cover_slip.cells)):
            init_save_track_analysis(self.cover_slip, cell_index, self.track_analysis, self.widget_track_analysis)
