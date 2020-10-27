"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Visualize the HMM analysis results.
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from graphviz import Digraph
import os
import seaborn as sns
import pandas as pd


class HmmVisualization():
    def __init__(self):
        self.state = ""
        self.bin_path = ""
        self.tmp_path = ""
        self.cells = []  # list of loadMergedHmm cell objects
        self.number_of_cells = 0
        self.number_of_states = 0  # number of hidden states = diffusion coefficients 
        self.cell_names = []
        self.aic_values = []
        self.bic_values = []
        self.log_likelihoods = []
        self.single_states_percentages = []
        self.states_percentages = []  # population of states in %
        self.states_percentages_error = []  # statistic error of mean
        self.single_tps = []
        self.mean_tps = []  # mean value of the transition probability tp matrix entries
        self.mean_tps_error = []
        self.single_Ds = []  # single values of diffusion coefficients per cell 
        self.mean_D = []  # mean values of diffusion coefficients
        self.mean_D_error = []
        self.loc_density = []  # localizations in the trc hmm file / cell size
        self.mean_node_size = []  # state population mapped on circle size
        self.mean_edge_size = []  # tp log mapped on arrow size
        # each mean diffusion coefficient has a colour index, the position in the colour index list refers to D
        # the value in the colour index refers to the colour.
        self.colour_index = []
        self.colour_palett = ["royalblue", "forestgreen", "darkorange", "darkmagenta", "orangered"]
        self.colour_palett_hex = ["#4169e1", "#228b22", "#ff8c00", "#8b008b", "#ff4500"]
        # "#%2x%2x%2x%2x"; alpha channel hex opacity values
        self.colour_palett_rgba = [i+"80" for i in self.colour_palett_hex]
        self.figures = []  # list with figure objects
        self.figure_names = []  # "figure_name"
        self.pixel_sizes = []
        # saving
        self.save_plots = False  # if true -> plots will be saved
        self.save_dir = ""
        self.save_folder_name = ""  

    def choose_state(self):
        if self.state == "physical model":
            self.calc_states_percentages_phys_mod()
        elif self.state == "equilibrium matrix":
            self.calc_states_percentages_eq_matrix()
        elif self.state == "state occurence":
            self.calc_states_percentages_occurence()

    def run(self):
        self.set_path()
        self.get_cell_names()
        self.get_number_of_cells()
        self.get_number_of_states()
        self.get_information_values()
        self.get_single_tps()
        self.get_pixel_sizes()
        self.choose_state()
        self.calc_mean_tp()
        self.calc_mean_D()
        self.get_colour_index()
        self.shuffle_colour_list()
        self.calc_loc_density()
        self.plot_D_boxplot()
        self.plot_D()
        self.plot_loc_density()
        self.plot_bar_state_percentages()
        self.plot_pie_state_percentages()
        self.plot_box_state_percentages()
        self.state_transition_diagram()

    def clear(self):
        self.__init__()

    def set_path(self):
        os.environ["PATH"] += os.pathsep + self.bin_path  # "C:\\Program Files (x86)\\Graphviz2.38\\bin"
        
    def run_save_plots(self):
        try:
            os.mkdir(self.save_dir + "\\" + self.save_folder_name)
        except OSError:
            print("Folder already exists.")
        else:
            for figure, name in zip(self.figures, self.figure_names):
                figure.savefig(self.save_dir + "\\" + self.save_folder_name + "\\" + name + ".pdf", format="pdf", transparent=True)
            self.state_transition_diagram()
        
    def get_number_of_states(self):
        self.number_of_states = np.shape(self.cells[0].diffusion_coef)[0]
    
    def get_number_of_cells(self):
        self.number_of_cells = len(self.cells)
        
    def get_pixel_sizes(self):
        for cell in self.cells:
            self.pixel_sizes.append(cell.pixel_size)

    def get_cell_names(self):
        """
        Get raw base name of loadMergedHmm cell objects
        """
        for cell in self.cells:
            self.cell_names.append(cell.hmm_cell_name)
        
    def get_information_values(self):
        self.aic_values = np.zeros(len(self.cells))
        self.bic_values = np.zeros(len(self.cells))
        self.log_likelihoods = np.zeros(len(self.cells))
        for cell in self.cells:
            cell_idx = self.cells.index(cell)
            self.aic_values[cell_idx] = cell.aic
            self.bic_values[cell_idx] = cell.bic
            self.log_likelihoods[cell_idx] = cell.log_likelihood
        
    def get_single_tps(self):
        for cell in self.cells:
            self.single_tps.append(cell.transition_matrix)
        
    def calc_states_percentages_eq_matrix(self):
        """
        Population of states in % based on the equilibrium matrix.
        """
        self.single_states_percentages = np.zeros([self.number_of_cells, self.number_of_states])
        self.state_percentages = np.zeros(self.number_of_states)
        self.states_percentages_error = np.zeros(self.number_of_states)
        eq_matrix = []
        for cell in self.cells:
            cell_idx = self.cells.index(cell)
            eq_matrix.append(cell.equilibrium_matrix)
            self.single_states_percentages[cell_idx] = cell.equilibrium_matrix 
        self.states_percentages = np.mean(eq_matrix, 0) 
        self.states_percentages_error = np.std(eq_matrix, 0,  ddof=1) / (self.number_of_cells)**(1/2)
        print("Mean population of states:", end=" ")
        for i in self.states_percentages:
            print("%.5f"%i, end=" ")
        print(" ")
    
    def calc_states_percentages_phys_mod(self):
        """
        Population of states in % based on the physical model.
        """
        self.single_states_percentages = np.zeros([self.number_of_cells, self.number_of_states])
        self.state_percentages = np.zeros(self.number_of_states)
        self.states_percentages_error = np.zeros(self.number_of_states)
        weight_coefs = []
        for cell in self.cells:
            cell_idx = self.cells.index(cell)
            cell_weights = []
            for state in range(self.number_of_states):
                cell_weights.append(cell.weight_coef[state][0])
            self.single_states_percentages[cell_idx] = cell_weights
            weight_coefs.append(cell_weights)
        self.states_percentages = np.mean(weight_coefs, 0)
        self.states_percentages_error = np.std(weight_coefs, 0,  ddof=1) / (self.number_of_cells)**(1/2)
        print("Mean population of states:", end=" ")
        for i in self.states_percentages:
            print("%.5f"%i, end=" ")
        print(" ")
 
    def calc_states_percentages_occurence(self):
        """
        Population of states in % based on the occurrence of states.
        """
        cell_states = np.zeros([self.number_of_cells, self.number_of_states])
        self.single_states_percentages = np.zeros([self.number_of_cells, self.number_of_states])
        for cell in self.cells:
            cell_idx = self.cells.index(cell)
            state_counter = np.zeros(self.number_of_states)
            for i in range(len(cell.trc_hmm)):
                for state_number in range(self.number_of_states):
                    if cell.trc_hmm[i][4] == state_number:
                        state_counter[state_number] += 1
            state_counter = np.divide(state_counter, np.sum(state_counter))
            self.single_states_percentages[cell_idx] = state_counter
            cell_states[cell_idx] = state_counter
        self.states_percentages = np.mean(cell_states, 0)
        self.states_percentages_error = np.std(cell_states, 0,  ddof=1) / (self.number_of_cells)**(1/2)
        print("Mean population of states:", end=" ")
        for i in self.states_percentages:
            print("%.5f"%i, end=" ")
        print(" ")
        
    def calc_mean_tp(self):
        """
        Mean value of the transition probabilities.
        """
        self.mean_tps = np.zeros(np.shape(self.cells[0].transition_matrix))
        self.mean_tps_error = np.zeros(np.shape(self.cells[0].transition_matrix))
        for row in range(self.number_of_states):
            for column in range(self.number_of_states):
                mean_tp = np.zeros(self.number_of_cells)
                for cell in self.cells:
                    cell_idx = self.cells.index(cell)
                    mean_tp[cell_idx] = cell.transition_matrix[column][row]
                self.mean_tps[column][row] = np.mean(mean_tp)
                self.mean_tps_error[column][row]= np.std(mean_tp, ddof=1) / (self.number_of_cells)**(1/2)
        print("Mean transition probabilities:") 
        for row in self.mean_tps:
            count = 0
            for i in row:
                count += 1
                if count != self.number_of_states:
                    print("%.5f"%i, end=" ")
                else:
                    print("%.5f"%i)
                
    def calc_mean_D(self):
        """
        Calc mean diffusion coefficient, error, extract diffusion coefficients per cell.
        """
        # columns = states, rows = cells
        self.mean_D = np.zeros(self.number_of_states)
        self.mean_D_error = np.zeros(self.number_of_states)
        self.single_Ds = np.zeros([self.number_of_cells, self.number_of_states])
        for cell in self.cells:
            cell_diffusion = np.zeros(self.number_of_states)
            cell_idx = self.cells.index(cell)
            for state in range(self.number_of_states):
                cell_diffusion[state] = cell.diffusion_coef[state][0]
            self.single_Ds[cell_idx] = cell_diffusion
        self.mean_D = np.mean(self.single_Ds, 0)
        self.mean_D_error = np.std(self.single_Ds, 0, ddof=1) / (self.number_of_cells)**(1/2) 
        print("Mean diffusion coefficients [\u00B5m\u00B2/s]:", end=" ")
        for i in self.mean_D:
            print("%.5f"%i, end=" ")
        print(" ")
        
    def get_colour_index(self):
        mean_D_sort = sorted(self.mean_D)
        for i in self.mean_D:
            idx = mean_D_sort.index(i)
            self.colour_index.append(idx)
            mean_D_sort[idx] = ""
        
    def shuffle_colour_list(self):
        self.colour_palett = ["royalblue", "forestgreen", "darkorange", "darkmagenta", "orangered"]
        self.colour_palett_hex = ["#4169e1", "#228b22", "#ff8c00", "#8b008b", "#ff4500"]
        self.colour_palett = self.colour_palett[:self.number_of_states]
        self.colour_palett_hex = self.colour_palett_hex[:self.number_of_states]
        shuffled_list = ["" for _ in self.colour_palett]
        shuffled_list_hex = ["" for _ in self.colour_palett]
        for j, i in enumerate(self.colour_index):
            shuffled_list[j] = self.colour_palett[i]
            shuffled_list_hex[j] = self.colour_palett_hex[i]
            j += 1
        self.colour_palett = shuffled_list
        self.colour_palett_hex = shuffled_list_hex       
            
    def calc_loc_density(self):
        """
        Number of localizations in the hmm trc file / cell size.
        """
        self.loc_density = np.zeros(self.number_of_cells)
        for cell in self.cells:
            self.loc_density[self.cells.index(cell)] = np.shape(cell.trc_hmm)[0]/cell.cell_size

    def plot_D_boxplot(self):
        sns.set_palette(sns.color_palette(self.colour_palett_hex))
        df = pd.DataFrame(self.single_Ds)
        fig = plt.figure(figsize=(3, 5))
        ax = sns.boxplot(data=df, showmeans=True,
                         meanprops={"marker": "s", "markerfacecolor": "white", "markeredgecolor": "0.25"})
        ax = sns.swarmplot(data=df, color="0.25")
        ax.set_title("Diffusion coefficients of states")
        ax.set_ylabel("Diffusion coefficient [\u00B5m\u00B2/s]")
        plt.show()
        self.figures.append(fig)
        self.figure_names.append("Diffusion_coefficients_boxplot")

    def plot_D(self):
        """
        Plot each diffusion coefficient vs the number of cells.
        """
        x = [i+1 for i in range(self.number_of_cells)]
        fig = plt.figure()
        plt.title("Diffusion coefficients of single cells")
        ax = fig.add_subplot(111)
        ax.set_axisbelow(True)
        ax.grid(linestyle=":", alpha=0.5)
        plt.xticks(np.arange(1, self.number_of_cells+1, step=1))
        for state in range(self.number_of_states):
            plt.plot(x, self.single_Ds[:, state], "o", color=self.colour_palett[state])
        plt.ylabel("Diffusion coefficient [\u00B5m\u00B2/s]")
        plt.xlabel("Cell number")
        plt.show()
        self.figures.append(fig)
        self.figure_names.append("Diffusion_coefficients_cells")
        
    def plot_loc_density(self):
        """
        Localization density vs cell number.
        """
        x = [i+1 for i in range(self.number_of_cells)]
        fig = plt.figure()
        plt.title("Number of localizations per cell and \u00B5m\u00B2")
        ax = fig.add_subplot(111)
        ax.set_axisbelow(True)
        ax.grid(linestyle=':', alpha=0.5)
        plt.xticks(np.arange(1, self.number_of_cells+1, step=1))
        plt.plot(x, self.loc_density, "o", color="darkslategray")
        plt.ylabel("Localization density [1/\u00B5m\u00B2]")
        plt.xlabel("Cell number")
        plt.show()
        self.figures.append(fig)
        self.figure_names.append("Localization_density")
        
    def plot_bar_state_percentages(self):
        bars = self.states_percentages
        label_name = self.mean_D
        title_name = "State distribution"
        error = self.states_percentages_error
        y_error = True
        x = np.arange(len(bars))
        x_name = [i+1 for i in range(self.number_of_states)]
        label_name = list(map(lambda x: self.D_rounded(x) + " \u00B5m\u00B2/s", label_name))
        fig, ax = plt.subplots()
        ax.set_xticklabels(label_name)
        ax.set_axisbelow(True)
        ax.grid(linestyle=":", alpha=0.5)
        ax.set_ylim(0, 1)
        plt.ylabel("Population")
        if y_error:
            lines = plt.bar(x, bars, yerr=error, capsize=5, edgecolor="black",
                            color=self.colour_palett[:self.number_of_states], label=label_name)
            plt.xticks(x, x_name)
        else:
            lines = plt.bar(x, bars, capsize=5, edgecolor="black", color=self.colour_palett[:self.number_of_states])
        plt.legend((lines),(label_name))
        plt.title(title_name)
        plt.show()
        self.figures.append(fig)
        self.figure_names.append("State_distrubution_bar_plot")

    def plot_box_state_percentages(self):
        sns.set_palette(sns.color_palette(self.colour_palett_hex))
        df = pd.DataFrame(self.single_states_percentages)
        fig = plt.figure(figsize=(3, 5))
        ax = sns.boxplot(data=df, showmeans=True,
                         meanprops={"marker": "s", "markerfacecolor": "white", "markeredgecolor": "0.25"})
        ax = sns.swarmplot(data=df, color="0.25")
        ax.set_title("State distribution")
        plt.show()
        self.figures.append(fig)
        self.figure_names.append("State_distrubution_box_plot")
 
    def plot_pie_state_percentages(self):
        values = self.states_percentages
        title_name = "State distribution"
        label = self.mean_D
        mean_pi = self.states_percentages
        fig = plt.figure()
        ax = fig.add_subplot(111)
        mean_pi = list(map(lambda x: x*100, mean_pi))
        label = list(map(lambda x: self.D_rounded(x) + " \u00B5m\u00B2/s", label))
        wedges, texts, autotexts = plt.pie(values, labels=label, colors=self.colour_palett[:self.number_of_states],
                                           autopct="%0.2f%%")
        plt.title(title_name)
        plt.show()
        self.figures.append(fig)
        self.figure_names.append("State_distrubution_pie_plot")

    def plot_trajectory(self, cell_name, trajectory_idx):
        """
        Plot a trajectory.
        :param cell_name: name of cell, index will be created to cell name.
        :param trajectory: number of trajectory -> index-1.
        """
        for cell in self.cells:
            if cell.hmm_cell_name == cell_name:
                cell_idx = self.cells.index(cell)
                trc_file = np.zeros([len(cell.trc_hmm), 6])
                trc_file[:, 0] = list(map(lambda row: row[0], cell.trc_hmm))  # col0 = track id
                trc_file[:, 1] = list(map(lambda row: row[1], cell.trc_hmm))  # frame
                trc_file[:, 2] = list(map(lambda row: row[2]*int(self.pixel_sizes[cell_idx])*10**(-3), cell.trc_hmm))  # x in ym
                trc_file[:, 3] = list(map(lambda row: row[3]*int(self.pixel_sizes[cell_idx])*10**(-3), cell.trc_hmm))  # y in ym
                trc_file[:, 4] = list(map(lambda row: row[4], cell.trc_hmm))  # state
                trc_file[:, 5] = list(map(lambda row: row[5], cell.trc_hmm))  # intensity
                idx = trc_file[:, 0] == trajectory_idx
                localizations = trc_file[idx, :]
                self.show_trajectory(localizations)
        
    def show_trajectory(self, localizations):
        plt.plot(localizations[:,2], localizations[:,3], linestyle="--",  label="localization", color="lightgray")
        for localization in localizations:
            state = int(localization[4])
            plt.plot(localization[2], localization[3], linestyle="--", marker="o", color=self.colour_palett[state])
        plt.title("Trajectory of one particle")
        plt.xlabel(r"$\mu$" + "m in x direction")
        plt.ylabel(r"$\mu$" + "m in y direction")
        plt.show()
    
    def D_rounded(self, D):
        exponent = 5  # floating point precision
        D *= 10**exponent
        D = int(np.round(D))
        D /= 10**exponent  # / instead of * -exp -> no floating point artifacts
        return str(D)

    def state_transition_diagram(self):
        min_node_size = 1.5
        # label are too large for node, the smallest % has to fit in the node
        mult_with_node_size = min_node_size / min(self.states_percentages)
        edge_fontsize = "10"
        dot = Digraph(comment="State Transition Diagram")
        float_precision = "%.3f"
        var_width = True
        colored_edges = True
        # return the state % between 0-100 %
        mean_diff_rounded = [str(self.tp_percentage_rounded(x)) for x in self.states_percentages]
        # A = pi * r^2 -> r = (A/pi)**0.5
        self.mean_node_size = list(map(lambda x: float(x)*mult_with_node_size/math.pi**(0.5), self.states_percentages))
        mean_node_size = list(map(lambda x: float_precision % (float(x)*mult_with_node_size/math.pi)**(0.5), self.states_percentages))
        diffusions = self.mean_D
        # represent D in ym^2/s with certain float precision
        diffusions = list(map(lambda x: str(float_precision % x), diffusions))
    
        for i in range(len(diffusions)):
            dot.node(str(i+1), mean_diff_rounded[i]+"%",
                     color="black", fillcolor=self.colour_palett_hex[i], fontcolor="black",
                     # colour_palett_rgba[i]
                     style="filled", shape="circle", fixedsize="shape", width=mean_node_size[i], pos="0,0!")
        self.mean_edge_size = np.zeros(np.shape(self.cells[0].transition_matrix))
        for row in range(np.shape(self.mean_tps)[0]):
            for column in range(np.shape(self.mean_tps)[1]):
                label_name = " " + str(self.tp_percentage_rounded(self.mean_tps[column][row]))+"%"
                tp = self.mean_tps[row][column]
                dot.edge(str(column+1), str(row+1), label=" "+label_name,
                         color=(self.gv_edge_color_gradient(self.colour_palett_hex[column], self.colour_palett_hex[row], 25) if colored_edges else "black"),
                         fontsize=edge_fontsize, style="filled", penwidth=(str(self.tp_px_mapping(tp, row, column)) if var_width else "1"))

        if not self.save_plots:
            dot.render(self.tmp_path + "/State_transition_diagram.svg", view=True)
        else:
            dot.render(self.save_dir + "\\" + self.save_folder_name + "\\" + "State_transition_diagram.svg", view=False)
            dot.format = "svg"
        
    def tp_percentage_rounded(self, tp):
        """
        Example: 0.4056357 -> 40.56
        """
        exponent = 4 
        while True:
            if tp >= 0.5*10**(-exponent):
                tp *= 10**exponent
                tp = int(np.round(tp))
                tp = tp / 10**(exponent-2)
                return tp
            exponent += 1
        
    def tp_px_mapping(self, tp, row, column, min_log=-5, min_px=0.5, max_px=5):
        log_tp = np.log10(tp)
        if log_tp <= min_log:
            return min_px
        log_tp -= min_log  # positive range
        log_tp /= -min_log  # 0-1
        self.mean_edge_size[row][column] = max_px*log_tp + min_px*(1-log_tp)
        return max_px*log_tp + min_px*(1-log_tp)  # sth between max & min px

    def gv_edge_color_gradient(self, c1, c2, res=50):
        color_list = ""
        c1 = Color.fromHex(c1)
        c2 = Color.fromHex(c2)
        cDiff = c2 - c1
        weight = str(1.0 / (res + 1.0))
        color_list += c1.toHex() + ";" + weight
        for i in range(res):
            color_list += ":"
            color_list += (c1 + cDiff.scalarMult(i / res)).toHex()
            color_list += ";" + weight
        return color_list
    

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
