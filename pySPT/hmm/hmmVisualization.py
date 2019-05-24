# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:19:31 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Visualize the HMM analysis results.
"""

import h5py
import os
import numpy as np
import matplotlib.pyplot as plt
from graphviz import Digraph
import os

os.environ["PATH"] += os.pathsep + 'C:\\Program Files (x86)\\Graphviz2.38\\bin'


class HMMCell():
    def __init__(self, file_path):
        self.hmm_cell_hdf5 = h5py.File(file_path, "r")
        self.hmm_cell_name = os.path.splitext(os.path.split(file_path)[1])[0]  # raw base name
        # hmm
        self.transition_matrix = []  # numpy.ndarray with numpy.ndarrays for each state
        self.equilibrium_matrix = []
        self.observation_matrix = []
        self.observation_alphabet = []  # numpy.ndarray
        # hmm -> statistics
        self.number_of_states = 0  # floats/int
        self.symbols = 0
        self.log_likelohood = 0
        self.dof = 0
        self.bic = 0
        self.aic = 0
        # judi
        self.judi = []
        # physical model
        self.diffusion_coef = []  # numpy.ndarray with numpy.voids
        self.weight_coef = []
        # trc
        self.trc = []  # numpy.ndarray with numpy.voids
        #cell size
        self.cell_size = 0.0
        # groups
        self.group_msd = self.hmm_cell_hdf5["MSD"]
        self.group_diffusion = self.hmm_cell_hdf5["diffusion"]
        self.group_hmm = self.hmm_cell_hdf5["hmm"]
        self.group_judi = self.hmm_cell_hdf5["judi"]
        self.group_physical_model = self.hmm_cell_hdf5["physicalModel"]
        self.group_rossier = self.hmm_cell_hdf5["rossier"]
        self.group_settings = self.hmm_cell_hdf5["settings"]
        self.group_statistics = self.hmm_cell_hdf5["statistics"]
        self.group_trc = self.hmm_cell_hdf5["trc"]
        self.group_settings = self.hmm_cell_hdf5["settings"]


    def get_statistics_info(self):
        """
        Statistics dataset contains states, symbols, loglikelihood, dof, bic, aic.
        """
        dset_statistics = self.group_hmm["statistics"]
        self.number_of_states = dset_statistics["states"][0]
        self.symbols = dset_statistics["symbols"][0]
        self.log_likelohood = dset_statistics["logLikelihood"][0]
        self.dof = dset_statistics["dof"][0]
        self.bic = dset_statistics["bic"][0]
        self.aic = dset_statistics["aic"][0]

    def get_transition_matrix(self):
        dset_transition_matrix = self.group_hmm["transitionMatrix"]
        self.transition_matrix = dset_transition_matrix[:,:]

    def get_equilibrium_matrix(self):
        dset_equilibrium_matrix = self.group_hmm["equilibriumMatrix"]
        self.equilibrium_matrix = dset_equilibrium_matrix[0]

    def get_observation_alphabet(self):
        ##### typo
        dset_observation_alphabet = self.group_hmm["observationAlphabet"]
        self.observation_alphabet = dset_observation_alphabet[0]

    def get_observation_matrix(self):
        dset_observation_matrix = self.group_hmm["observationMatrix"]
        self.observation_matrix = dset_observation_matrix[:,:]

    def get_diffusion_coef(self):
        dset_diffusion_coef = self.group_physical_model["diffusionCoefficient"]
        self.diffusion_coef = dset_diffusion_coef[:]
        #print(self.diffusion_coef)

    def get_weight_coef(self):
        dset_weight_coef = self.group_physical_model["weightCoefficient"]
        self.weight_coef = dset_weight_coef[:]
        #print(self.weight_coef)

    def get_judi(self):
        #####
        dset_judi = self.group_judi["judiFile"]
        #self.judi = dset_judi
        #print(dset_judi[:,])

    def get_trc(self):
        dset_trc = self.group_trc["trcFile"]
        self.trc = dset_trc[:]
        #print(self.trc)

    def show_hdf5_file(self):
        print(self.hmm_cell_hdf5)

    def close_hdf5_file(self):
        self.hmm_cell_hdf5.close()
        
    def get_cell_size(self):
        dset_settings = self.group_settings["settings"]
        self.cell_size = dset_settings[0][0][3]

    def run(self):
        self.get_statistics_info()
        self.get_transition_matrix()
        self.get_equilibrium_matrix()
        self.get_observation_matrix()
        self.get_observation_alphabet()
        self.get_weight_coef()
        self.get_diffusion_coef()
        self.get_trc()
        self.get_cell_size()
        self.close_hdf5_file()


def main():
    # home paths
    #cell01_path = "C:\\Users\\Minakon\\Downloads\\CS5_Cell01Cell02\\CS5_Cell01Cell02\\190412_cell_1_MMStack_Pos0_trc_format_merged_3States.h5"
    #cell02_path = "C:\\Users\\Minakon\\Downloads\\CS5_Cell01Cell02\\CS5_Cell01Cell02\\190412_cell_2_MMStack_Pos0_trc_format_merged_3States.h5"
    # test paths
    #cell01_path = "C:\\Users\\pcoffice37\\Documents\\testing_file_search\\testHdf5Merge\\CS5_Cell01Cell02\\190412_cell_1_MMStack_Pos0_trc_format_merged_3States.h5"  
    #cell02_path = "C:\\Users\\pcoffice37\\Documents\\testing_file_search\\testHdf5Merge\\CS5_Cell01Cell02\\190412_cell_2_MMStack_Pos0_trc_format_merged_3States.h5"  
    # CS5 Fab 
    cover_slip = []
    Fab_CS5_paths = []
# =============================================================================
#     for i in range(1, 11):
# 
#         cell_path = "C:\\Users\\pcoffice37\\Documents\\Datasets_MET_ML\\Fab_MET_data\\160404_coverslip5_Fab-Atto647N_025nM_in_IM\\160404_CS5_restingMET_pBleach_new\\HMM\\160404_CS5_restingMET\\cell" + str(i) +"_3States.h5"
#         Fab_CS5_paths.append(cell_path)    
# =============================================================================
    
    for i in range(1, 21):
        cell_path = "C:\\Users\\pcoffice37\\Documents\\Datasets_MET_ML\\Fab_MET_data\\160404_coverslip5_Fab-Atto647N_025nM_in_IM\\160404_CS5_restingMET_pBleach_new\\HMM\\160404_CS5_restingMET_Dmin0\\cell" + str(i) + "_3States_Dmin0.h5"
        Fab_CS5_paths.append(cell_path)
    print(Fab_CS5_paths)
    
    #cell_path_10 = "C:\\Users\\pcoffice37\\Documents\\Datasets_MET_ML\\Fab_MET_data\\160404_coverslip5_Fab-Atto647N_025nM_in_IM\\160404_CS5_restingMET_pBleach_new\\HMM\\160404_CS5_restingMET\\cell10_3States.h5"
    #Fab_CS5_paths.append(cell_path_10)
    for cell_path in Fab_CS5_paths:
        hmm_cell = HMMCell(cell_path)
        hmm_cell.run()
        cover_slip.append(hmm_cell)
    # calc & plot
    get_cell_names(cover_slip)
    get_AIC(cover_slip)
# =============================================================================
#     calc_mean_hmm_pi(cover_slip)
#     calc_mean_physmod_pi(cover_slip)
# =============================================================================
# =============================================================================
#     calc_mean_states(cover_slip)
#     calc_mean_tp(cover_slip)
#     #plot_jd(hmm_cell01)
#     calc_D(cover_slip)
# =============================================================================
    mean_pi = calc_mean_states(cover_slip)[0]
    dmean_pi = calc_mean_states(cover_slip)[1]
    mean_tp = calc_mean_tp(cover_slip)[0]
    dmean_tp = calc_mean_tp(cover_slip)[1]
    diffusions = calc_D(cover_slip)[0]
    ddiffusions = calc_D(cover_slip)[1]
    single_diffusions = calc_D(cover_slip)[2]
    plot_D_boxplot(single_diffusions)
    plot_D(single_diffusions, cover_slip)
    loc_density = calc_loc_density(cover_slip)
    plot_loc_density(loc_density, cover_slip)
    #state_transition_diagram(mean_pi, mean_tp, dmean_pi, dmean_tp, diffusions, ddiffusions)
    print("mean pi based on frequencies : ", mean_pi)
    plot_bar(mean_pi, diffusions,
             "State distribution based on frequency of states",
             error=dmean_pi, y_error=True)
    
    plot_pie(mean_pi, "State distribution based on frequency of states", diffusions, mean_pi)
    print("mean tp values: ", mean_tp)
    state_transition_diagram(mean_pi, mean_tp, dmean_pi, dmean_tp, diffusions, ddiffusions)
    print("mean diff", diffusions)



# visualization

number_of_states = 3
number_of_cells = 12
#np.set_printoptions(precision=3)
colour_palett = ["royalblue", "forestgreen", "darkorange", "darkmagenta", "orangered"]
colour_palett_hex = ["#4169e1", "#228b22", "#ff8c00", "#8b008b", "#ff4500"]
colour_palett_rgba = [i+"80" for i in colour_palett_hex] # "#%2x%2x%2x%2x"; alpha channel hex opacity values: https://medium.com/@magdamiu/android-transparent-colors-a2d55a9b4e66

def get_cell_names(cover_slip):
    cell_names = []
    for cell in cover_slip:
        cell_names.append(cell.hmm_cell_name)
    print("Cell names: ", cell_names)

def get_AIC(cover_slip):
    aic_values = np.zeros(len(cover_slip))
    for cell in cover_slip:
        cell_index = cover_slip.index(cell)
        aic_values[cell_index] = cell.aic
    print("AIC values: ", aic_values)

def calc_D(cover_slip):
    """
    return: mean diff coeff, error, single diff coeffs
    """
    # columns = states, rows = cells
    mean_diff_coeff = np.zeros(number_of_states)
    mean_diff_coeff_error = np.zeros(number_of_states)
    #diff_coeffs = np.zeros([len(cover_slip), number_of_states])
    diff_coeffs = np.zeros([len(cover_slip), number_of_states])

    for cell in cover_slip:
        cell_diffusion = np.zeros(number_of_states)
        cell_index = cover_slip.index(cell)
        for state in range(number_of_states):
            cell_diffusion[state] = cell.diffusion_coef[state][0]
        diff_coeffs[cell_index] = cell_diffusion
    mean_diff_coeff = np.mean(diff_coeffs, 0)
    mean_diff_coeff_error = np.std(diff_coeffs, 0, ddof=1) / (len(cover_slip))**(1/2)
    return mean_diff_coeff, mean_diff_coeff_error, diff_coeffs

def plot_D(diff_coeffs, cover_slip):
    """
    Plot each diffusion coefficient vs the number of cells.
    """
    x = [i+1 for i in range(len(cover_slip))]
    fig = plt.figure()
    plt.title("Diffusion coefficients of single cells")
    ax = fig.add_subplot(111)
    ax.set_axisbelow(True)
    ax.grid(linestyle=':', alpha=0.5)
    plt.xticks(np.arange(1, len(cover_slip)+1, step=1))
    for state in range(number_of_states):
        plt.plot(x, diff_coeffs[:, state], "o", color=colour_palett[state])
        #plt.plot(x, diff_coeffs[:, state], "-", alpha=0.5, color=colour_palett[state])
    plt.ylabel("Diffusion coefficient [\u00B5m\u00B2/s]")
    plt.xlabel("Cell number")
    plt.show()

def plot_D_boxplot(diff_coeffs):
    """
    boxplot with matplotlib: Box = first to third quartile, lower to higher quartile,
    the lower quartile splits off 25 % from 75 % of the data
    the higher quartile splits off the highest 25 % of the data
    middle line: median, which splits the dataset in half
    whiskers show the range of the data.
    """
    data_to_plot = diff_coeffs
    print(data_to_plot)
    # Create a figure instance
    fig = plt.figure()
    plt.title("Diffusion coefficients of HMM states")
    # Create an axes instance
    ax = fig.add_subplot(111)
    ax.set_axisbelow(True)
    ax.grid(linestyle=':', alpha=0.5)
    plt.ylabel("Diffusion coefficient [\u00B5m\u00B2/s]")
    plt.xlabel("Number of state")
    # Create the boxplot
    bp = ax.boxplot(data_to_plot, patch_artist=True)
    ## change outline color, fill color and linewidth of the boxes
    for box in bp['boxes']:
        # change outline color
        #box.set(color='#7570b3', linewidth=2)
        # change fill color
        box.set(facecolor=colour_palett_hex[bp["boxes"].index(box)])
    # ## change color and linewidth of the whiskers
    # for whisker in bp['whiskers']:
    #     whisker.set(color='#7570b3', linewidth=2)
    # ## change color and linewidth of the caps
    # for cap in bp['caps']:
    #     cap.set(color='#7570b3', linewidth=2)
    # ## change color and linewidth of the medians
    # for median in bp['medians']:
    #     median.set(color='#b2df8a', linewidth=2)
    # ## change the style of fliers and their fill
    # for flier in bp['fliers']:
    #     flier.set(marker='o', color='#e7298a', alpha=0.5)
    plt.show()

def calc_loc_density(cover_slip):
    """
    Number of localizations included for trajectory building (not quite the density, but gives an impression).
    """
    loc_density = np.zeros(len(cover_slip))
    for cell in cover_slip:
        loc_density[cover_slip.index(cell)] = np.shape(cell.trc)[0]/cell.cell_size
    return loc_density
    
def plot_loc_density(loc_density, cover_slip):
    x = [i+1 for i in range(len(cover_slip))]
    fig = plt.figure()
    plt.title("Number of localizations per cell and \u00B5m\u00B2")
    ax = fig.add_subplot(111)
    ax.set_axisbelow(True)
    ax.grid(linestyle=':', alpha=0.5)
    plt.xticks(np.arange(1, len(cover_slip)+1, step=1))
    plt.plot(x, loc_density, "o", color="darkslategray")
    #plt.plot(x, number_locs, "-", alpha=0.5, color="darkslategray")
    plt.ylabel("Localization density [1/\u00B5m\u00B2]")
    plt.xlabel("Cell number")
    plt.show()

def plot_jd(cell):
    pass

def plot_pie(values, title_name, label, mean_pi):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    mean_pi = list(map(lambda x: x*100, mean_pi))
    label = list(map(lambda x: str(x)[:5]+" \u00B5m\u00B2/s", label))
    wedges, texts, autotexts = plt.pie(values, labels=label, colors=colour_palett[:number_of_states], autopct="%0.2f%%")
    #plt.setp(autotexts, size=8, weight="bold")
    plt.title(title_name)
    plt.show()
    
def plot_bar(bars, label_name, title_name, error=None, y_error=False):
    x = np.arange(len(bars))
    x_name = [i+1 for i in range(number_of_states)]
    label_name = list(map(lambda x: str(x)[:5]+" \u00B5m\u00B2/s", label_name))
    fig, ax = plt.subplots()
    ax.set_xticklabels(label_name)
    ax.set_axisbelow(True)
    ax.grid(linestyle=':', alpha=0.5)
    ax.set_ylim(0,1)
    if y_error:
        plt.bar(x, bars, yerr=error, capsize=5, edgecolor="black", color=colour_palett[:number_of_states], label=label_name)  # label=label_name
        plt.xticks(x, x_name)
    else:
        plt.bar(x, bars, capsize=5, edgecolor="black", color=colour_palett[:number_of_states])
    plt.legend()
    plt.title(title_name)
    plt.show()

# =============================================================================
# def calc_mean_hmm_pi(cover_slip):
#     number_of_cells = len(cover_slip)
#     mean_hmm_pis = np.zeros(number_of_states)
#     mean_hmm_pis_error = np.zeros(number_of_states)
#     for i in range(number_of_states):
#         mean_hmm_pi = np.zeros(number_of_cells)
#         for cell in cover_slip:
#             cell_index = cover_slip.index(cell)
#             mean_hmm_pi[cell_index] = cell.equilibrium_matrix[i]
#         mean_hmm_pis[i] = np.mean(mean_hmm_pi)
#         mean_hmm_pis_error[i] = np.std(mean_hmm_pi, ddof=1)/(number_of_cells)**(1/2)
#     plot_bar(mean_hmm_pis, "probability", "State distribution based on equilibrium matrix", error=mean_hmm_pis_error, y_error=True)
#     print("hmm pi mean", mean_hmm_pis)    
# 
# def calc_mean_physmod_pi(cover_slip):
#     number_of_cells = len(cover_slip)
#     mean_physmod_pis = np.zeros(number_of_states)
#     mean_physmod_pis_error = np.zeros(number_of_states)
#     for i in range(number_of_states):
#         mean_physmod_pi = np.zeros(number_of_cells)
#         for cell in cover_slip:
#             cell_index = cover_slip.index(cell)
#             mean_physmod_pi[cell_index] = cell.weight_coef[i][0]
#         mean_physmod_pis[i] = np.mean(mean_physmod_pi)
#         mean_physmod_pis_error[i] = np.std(mean_physmod_pi, ddof=1)/number_of_cells**(1/2)
#     plot_bar(mean_physmod_pis, "probability", "State distribution based on weighted physical model", error=mean_physmod_pis_error, y_error=True)
#     print("physmod pi mean", mean_physmod_pis)
# =============================================================================
    
def calc_mean_states(cover_slip):
    number_of_cells = len(cover_slip)
    cell_states = np.zeros([len(cover_slip), number_of_states])
    for cell in cover_slip:
        cell_index = cover_slip.index(cell)
        state_counter = np.zeros(number_of_states)
        for i in range(len(cell.trc)):
            for state_number in range(number_of_states):
                if cell.trc[i][4] == state_number:
                    state_counter[state_number] += 1
        state_counter = np.divide(state_counter, np.sum(state_counter))
        cell_states[cell_index] = state_counter
    mean_states = np.mean(cell_states, 0)
    mean_states_error = np.std(cell_states, 0,  ddof=1) / (number_of_cells) **(1/2)
    return mean_states, mean_states_error

def calc_mean_tp(cover_slip):
    number_of_cells = len(cover_slip)
    mean_tps = np.zeros(np.shape(cover_slip[0].transition_matrix))
    mean_tps_error = np.zeros(np.shape(cover_slip[0].transition_matrix))
    for row in range(number_of_states):
        for column in range(number_of_states):
            mean_tp = np.zeros(number_of_cells)
            for cell in cover_slip:
                cell_index = cover_slip.index(cell)
                mean_tp[cell_index] = cell.transition_matrix[column][row]
            mean_tps[column][row] = np.mean(mean_tp)
            mean_tps_error[column][row]= np.std(mean_tp, ddof=1)/(number_of_cells)**(1/2)
    
    return mean_tps, mean_tps_error


def state_transition_diagram2(mean_diff, mean_tp, dmean_diff, dmean_tp, diffusions, ddiffusions):
    # round significantly ...
    add_to_node_size = 0.5  # label are too large for node -> add value to which pi_diffusion is added
    transparency_correction = 0.2  # tp values can be very small -> add value to transparency to make them visible
    transparency_correction_threshold = 0.05  # if tp < threshold -> transparency correction is applied
    dashed_threshold = 0.002  # tp near 0 are displayed as dashed
    edge_fontsize = "10"
    error = False
    dot = Digraph(comment="State Transition Diagram")
    float_precision = "%.3f"
    float_precision_np = int(float_precision[2])+2
    mean_tp_str = mean_tp.astype("str")
    dmean_tp_str = dmean_tp.astype("str")
    mean_diff = list(map(lambda x: str(x+add_to_node_size), mean_diff))
    diffusions = list(map(lambda x: str(float_precision % x), diffusions))
    ddiffusions = list(map(lambda x: str(float_precision % x), ddiffusions))

    for i in range(len(diffusions)):
        if error:
            node_name = diffusions[i] + " \u00B1 " + ddiffusions[i]
        else:
            node_name = diffusions[i]
        dot.node(str(i+1), node_name,
                 color="black", fillcolor=colour_palett_hex[i], fontcolor="black",
                 # colour_palett_rgba[i]
                 style="filled", shape="circle", fixedsize="shape", width=mean_diff[i], pos="0,0!")

    for row in range(np.shape(mean_tp)[0]):
        for column in range(np.shape(mean_tp)[1]):
            if error:
                label_name = " "+mean_tp_str[column][row][:float_precision_np]+ " \u00B1 " +dmean_tp_str[column][row][:float_precision_np]
            else:
                label_name = " " + mean_tp_str[column][row][:float_precision_np]
            tp = mean_tp[row][column]
            tp_transparency = hex(int(np.round(tp * 255)))[2:]
            # correct for transparency (most tp values are very small -> nearly invisible)
            if tp < transparency_correction_threshold and tp > dashed_threshold:
                tp += transparency_correction # correction factor for transparency, otherwise arrow not visible at all
                tp_transparency = hex(int(np.round(tp * 255)))[2:] # 255 * tp results in float -> round (0.5 = 0, 0.6 = 1 ...) and convert to int for hex convertion
                dot.edge(str(column + 1), str(row + 1), label=label_name,
                         color=colour_palett_hex[row] + tp_transparency + ":" + colour_palett_hex[column] + tp_transparency,
                         fontsize=edge_fontsize)
            # to display difference between very small and nearly 0 tp values, draw nearly 0 with dashed lines
            elif tp < dashed_threshold:
                tp += transparency_correction
                tp_transparency = hex(int(np.round(tp * 255)))[2:]
                dot.edge(str(column + 1), str(row + 1), label=label_name,
                         color=colour_palett_hex[row] + tp_transparency + ":" + colour_palett_hex[column] + tp_transparency,
                         style="dashed", fontsize=edge_fontsize)
            else:
                dot.edge(str(column+1), str(row+1), label=label_name,
                         color=colour_palett_hex[row]+tp_transparency+":"+colour_palett_hex[column]+tp_transparency, fontsize=edge_fontsize)

    dot.render('test-output/round-table.gv', view=True)

def tp_percentage_rounded(tp):
    exponent = 4 
    while True:
        if tp >= 0.5*10**(-exponent):
            tp *= 10**exponent
            tp = int(np.round(tp))
            tp = tp / 10**(exponent-2)
            return tp
        exponent += 1
        
def tp_px_mapping(tp, min_log=-5, min_px=0.5, max_px=5):
    log_tp = np.log10(tp)
    if log_tp <= min_log:
        return min_px
    log_tp -= min_log  # positive range
    log_tp /= -min_log  # 0-1
    return max_px*log_tp + min_px*(1-log_tp)
    
def state_transition_diagram(mean_diff, mean_tp, dmean_diff, dmean_tp, diffusions, ddiffusions):
 
    min_node_size = 0.5
    mult_with_node_size = min_node_size / min(mean_diff) # label are too large for node    
    
    edge_fontsize = "10"
    dot = Digraph(comment="State Transition Diagram")
    float_precision = "%.3f"
    float_precision_np = int(float_precision[2])+2
    
    var_width = True
    colored_edges = True
    mean_diff_rounded = [str(tp_percentage_rounded(x)) for x in mean_diff]
    mean_diff_size = list(map(lambda x: str(float_precision % (float(x)*mult_with_node_size)**(0.5)), mean_diff))
    diffusions = list(map(lambda x: str(float_precision % x), diffusions))

    for i in range(len(diffusions)):
        dot.node(str(i+1), mean_diff_rounded[i]+"%",
                 color="black", fillcolor=colour_palett_hex[i], fontcolor="black",
                 # colour_palett_rgba[i]
                 style="filled", shape="circle", fixedsize="shape", width=mean_diff_size[i], pos="0,0!")

    for row in range(np.shape(mean_tp)[0]):
        for column in range(np.shape(mean_tp)[1]):
            #label_name = " " + mean_tp_str[column][row][:float_precision_np]
            label_name = str(tp_percentage_rounded(mean_tp[column][row]))+"%"
            tp = mean_tp[row][column]
            
            dot.edge(str(column+1), str(row+1), label=" "+label_name,
                     color=(gv_edge_color_gradient(colour_palett_hex[column], colour_palett_hex[row], 25) if colored_edges else "black"),
                     fontsize=edge_fontsize, style="filled", penwidth=(str(tp_px_mapping(tp)) if var_width else "1"))
    dot.render('test-output/Dmin.gv', view=True)
    

def gv_edge_color_gradient(c1, c2, res=50):
    colorList = ""
    c1 = Color.fromHex(c1)
    c2 = Color.fromHex(c2)
    cDiff = c2 - c1
    weight = str(1.0 / (res + 1.0))
    colorList += c1.toHex() + ";" + weight
    for i in range(res):
        colorList += ":"
        colorList += (c1 + cDiff.scalarMult(i / res)).toHex()
        colorList += ";" + weight
    return colorList

class Color:
    def __init__(self, r, g, b):
        self.r = r
        self.g = g
        self.b = b
        
    @staticmethod
    def fromHex(hexcode):
        hexcode = hexcode[1:]
        r = int(hexcode[:2], 16) / 255.0
        g = int(hexcode[2:4], 16) / 255.0
        b = int(hexcode[4:], 16) / 255.0
        return Color(r, g, b)
        
    def toHex(self):
        self.clamp()
        hexStr = "#"
        for channel in (self.r, self.g, self.b):
            hexStr += str(hex(int(channel * 255.0)))[2:]
        return hexStr
        
    def __add__(self, other):
        r = self.r + other.r
        g = self.g + other.g
        b = self.b + other.b
        return Color(r, g, b)

    def __sub__(self, other):
        r = self.r - other.r
        g = self.g - other.g
        b = self.b - other.b
        return Color(r, g, b)
        
    def scalarMult(self, scalar):
        r = self.r * scalar
        g = self.g * scalar
        b = self.b * scalar
        return Color(r, g, b)
    
    def clamp(self):
        self.r = max(0.0, min(1.0, self.r))
        self.g = max(0.0, min(1.0, self.g))
        self.b = max(0.0, min(1.0, self.b))


if __name__ == "__main__":
    main()
