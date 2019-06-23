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
import math
import matplotlib.pyplot as plt
from graphviz import Digraph
import os

os.environ["PATH"] += os.pathsep + 'C:\\Program Files (x86)\\Graphviz2.38\\bin'


class HmmVisualization():
    def __init__(self):
        self.cells = []  # list of loadMergedHmm cell objects
        self.number_of_cells = 0
        self.number_of_states = 3  # fill out function!!!!
        self.cell_names = []
        self.aic_values = []
        self.mean_aic_value = 0
        self.states_percentages = []  # population of states in %
        self.states_percentages_error = []  # statistic error of mean
        self.mean_tps = []  # mean value of the transition probability tp matrix entries
        self.mean_tps_error = []
        self.single_Ds = []  # single values of diffusion coefficients per cell 
        self.mean_D = []  # mean values of diffusion coefficients
        self.mean_D_error = []
        self.loc_density = []  # localizations in the trc hmm file / cell size
        self.colour_palett = ["royalblue", "forestgreen", "darkorange", "darkmagenta", "orangered"]
        self.colour_palett_hex = ["#4169e1", "#228b22", "#ff8c00", "#8b008b", "#ff4500"]
        self.colour_palett_rgba = [i+"80" for i in self.colour_palett_hex] # "#%2x%2x%2x%2x"; alpha channel hex opacity values: https://medium.com/@magdamiu/android-transparent-colors-a2d55a9b4e66
        
    def run(self):
        self.get_cell_names()
        self.get_number_of_cells()
        self.get_number_of_states()
        self.get_aic_values()
        self.calc_states_percentages()
        self.calc_mean_tp()
        self.calc_mean_D()
        self.calc_loc_density()
        self.plot_D_boxplot()
        self.plot_D()
        self.plot_loc_density()
        self.plot_bar_state_percentages()
        self.plot_pie_state_percentages()
        self.state_transition_diagram()
        
    def get_number_of_states(self):
        pass##############################################
    ### fix round error for D very small
    
    def get_number_of_cells(self):
        self.number_of_cells = len(self.cells)
        
    def get_cell_names(self):
        """
        Get raw base name of loadMergedHmm cell objects
        """
        for cell in self.cells:
            self.cell_names.append(cell.hmm_cell_name)
        #print("Cell names: ", self.cell_names)
        
    def get_aic_values(self):
        self.aic_values = np.zeros(len(self.cells))
        for cell in self.cells:
            cell_idx = self.cells.index(cell)
            self.aic_values[cell_idx] = cell.aic
        self.mean_aic_value = np.mean(self.aic_values)
        #print("AIC values: ", self.aic_values)
        print("Mean AIC value: %.3f"% self.mean_aic_value)
        
    def calc_states_percentages(self):
        """
        Population of states in %.
        """
        cell_states = np.zeros([self.number_of_cells, self.number_of_states])
        for cell in self.cells:
            cell_idx = self.cells.index(cell)
            state_counter = np.zeros(self.number_of_states)
            for i in range(len(cell.trc_hmm)):
                for state_number in range(self.number_of_states):
                    if cell.trc_hmm[i][4] == state_number:
                        state_counter[state_number] += 1
            state_counter = np.divide(state_counter, np.sum(state_counter))
            cell_states[cell_idx] = state_counter
        self.states_percentages = np.mean(cell_states, 0)
        self.states_percentages_error = np.std(cell_states, 0,  ddof=1) / (self.number_of_cells) **(1/2)
        print("Mean population of states:", end=' ')
        for i in self.states_percentages:
            print("%.5f"%i, end=' ')
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
                self.mean_tps_error[column][row]= np.std(mean_tp, ddof=1)/(self.number_of_cells)**(1/2)
        print("Mean transition probabilities:") 
        for row in self.mean_tps:
            count = 0
            for i in row:
                count += 1
                if count != self.number_of_states:
                    print("%.5f"%i, end=' ')
                else:
                    print("%.5f"%i)
                
    def calc_mean_D(self):
        """
        return: mean diff coeff, error, single diff coeffs
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
        print("Mean diffusion coefficients [\u00B5m\u00B2/s]", end=" ")
        for i in self.mean_D:
            print("%.5f"%i, end=' ')
        print(" ")
        #print("Single Ds", self.single_Ds)
        
    def calc_loc_density(self):
        """
        Number of localizations in the hmm trc file / cell size.
        """
        self.loc_density = np.zeros(self.number_of_cells)
        for cell in self.cells:
            self.loc_density[self.cells.index(cell)] = np.shape(cell.trc_hmm)[0]/cell.cell_size
        #print("Density", self.loc_density)

    def plot_D_boxplot(self):
        """
        boxplot with matplotlib: Box = first to third quartile, lower to higher quartile,
        the lower quartile splits off 25 % from 75 % of the data
        the higher quartile splits off the highest 25 % of the data
        middle line: median, which splits the dataset in half
        whiskers show the range of the data.
        """
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
        bp = ax.boxplot(self.single_Ds, patch_artist=True)
        ## change outline color, fill color and linewidth of the boxes
        for box in bp['boxes']:
            # change outline color
            #box.set(color='#7570b3', linewidth=2)
            # change fill color
            box.set(facecolor=self.colour_palett_hex[bp["boxes"].index(box)])
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
        
    def plot_D(self):
        """
        Plot each diffusion coefficient vs the number of cells.
        """
        x = [i+1 for i in range(self.number_of_cells)]
        fig = plt.figure()
        plt.title("Diffusion coefficients of single cells")
        ax = fig.add_subplot(111)
        ax.set_axisbelow(True)
        ax.grid(linestyle=':', alpha=0.5)
        plt.xticks(np.arange(1, self.number_of_cells+1, step=1))
        for state in range(self.number_of_states):
            plt.plot(x, self.single_Ds[:, state], "o", color=self.colour_palett[state])
            #plt.plot(x, diff_coeffs[:, state], "-", alpha=0.5, color=colour_palett[state])
        plt.ylabel("Diffusion coefficient [\u00B5m\u00B2/s]")
        plt.xlabel("Cell number")
        plt.show()
        
    def plot_loc_density(self):
        """
        Localization density vs cell number
        """
        x = [i+1 for i in range(self.number_of_cells)]
        fig = plt.figure()
        plt.title("Number of localizations per cell and \u00B5m\u00B2")
        ax = fig.add_subplot(111)
        ax.set_axisbelow(True)
        ax.grid(linestyle=':', alpha=0.5)
        plt.xticks(np.arange(1, self.number_of_cells+1, step=1))
        plt.plot(x, self.loc_density, "o", color="darkslategray")
        #plt.plot(x, number_locs, "-", alpha=0.5, color="darkslategray")
        plt.ylabel("Localization density [1/\u00B5m\u00B2]")
        plt.xlabel("Cell number")
        plt.show()
        
    def plot_bar_state_percentages(self):
        bars = self.states_percentages
        label_name = self.mean_D
        title_name = "State distribution based on frequency of states"
        error = self.states_percentages_error
        y_error = True
        
        x = np.arange(len(bars))
        x_name = [i+1 for i in range(self.number_of_states)]
        label_name = list(map(lambda x: str(x)[:5]+" \u00B5m\u00B2/s", label_name))  # nearly zero isnt handeled right !!!!!!!!
        fig, ax = plt.subplots()
        ax.set_xticklabels(label_name)
        ax.set_axisbelow(True)
        ax.grid(linestyle=':', alpha=0.5)
        ax.set_ylim(0,1)
        if y_error:
            plt.bar(x, bars, yerr=error, capsize=5, edgecolor="black", color=self.colour_palett[:self.number_of_states], label=label_name)  # label=label_name
            plt.xticks(x, x_name)
        else:
            plt.bar(x, bars, capsize=5, edgecolor="black", color=self.colour_palett[:self.number_of_states])
        plt.legend()
        plt.title(title_name)
        plt.show()
 
    def plot_pie_state_percentages(self):
        values = self.states_percentages
        title_name = "State distribution based on frequency of states"
        label = self.mean_D
        mean_pi = self.states_percentages
        fig = plt.figure()
        ax = fig.add_subplot(111)
        mean_pi = list(map(lambda x: x*100, mean_pi))
        label = list(map(lambda x: str(x)[:5]+" \u00B5m\u00B2/s", label))
        wedges, texts, autotexts = plt.pie(values, labels=label, colors=self.colour_palett[:self.number_of_states], autopct="%0.2f%%")
        #plt.setp(autotexts, size=8, weight="bold")
        plt.title(title_name)
        plt.show()
    # state_transition_diagram(mean_pi, mean_tp, dmean_pi, dmean_tp, diffusions, ddiffusions)

    def state_transition_diagram(self):
        
        min_node_size = 1.5
        mult_with_node_size = min_node_size / min(self.states_percentages) # label are too large for node, the smallest % has to fit in the node
        edge_fontsize = "10"
        dot = Digraph(comment="State Transition Diagram")
        float_precision = "%.3f"

        var_width = True
        colored_edges = True
        mean_diff_rounded = [str(self.tp_percentage_rounded(x)) for x in self.states_percentages]  # return the state % between 0-100 % 
        # A = pi * r^2 -> r = (A/pi)**0.5
        mean_diff_size = list(map(lambda x: float_precision % (float(x)*mult_with_node_size/math.pi)**(0.5), self.states_percentages))

        diffusions = self.mean_D
        diffusions = list(map(lambda x: str(float_precision % x), diffusions))  # represent D in ym^2/s with certain float precision
    
        for i in range(len(diffusions)):
            dot.node(str(i+1), mean_diff_rounded[i]+"%",
                     color="black", fillcolor=self.colour_palett_hex[i], fontcolor="black",
                     # colour_palett_rgba[i]
                     style="filled", shape="circle", fixedsize="shape", width=mean_diff_size[i], pos="0,0!")
        for row in range(np.shape(self.mean_tps)[0]):
            for column in range(np.shape(self.mean_tps)[1]):
                #label_name = " " + mean_tp_str[column][row][:float_precision_np]
                label_name = " " + str(self.tp_percentage_rounded(self.mean_tps[column][row]))+"%"
                tp = self.mean_tps[row][column]
                
                dot.edge(str(column+1), str(row+1), label=" "+label_name,
                         color=(self.gv_edge_color_gradient(self.colour_palett_hex[column], self.colour_palett_hex[row], 25) if colored_edges else "black"),
                         fontsize=edge_fontsize, style="filled", penwidth=(str(self.tp_px_mapping(tp)) if var_width else "1"))

        dot.render('test-output/Dmin.gv', view=True)

    def tp_percentage_rounded(self, tp):
        """
        0.4056357 -> 40.56
        """
        exponent = 4 
        while True:
            if tp >= 0.5*10**(-exponent):
                tp *= 10**exponent
                tp = int(np.round(tp))
                tp = tp / 10**(exponent-2)
                return tp
            exponent += 1
        
    def tp_px_mapping(self, tp, min_log=-5, min_px=0.5, max_px=5):
        log_tp = np.log10(tp)
        if log_tp <= min_log:
            return min_px
        log_tp -= min_log  # positive range
        log_tp /= -min_log  # 0-1
        return max_px*log_tp + min_px*(1-log_tp)  # sth between max & min px

    def gv_edge_color_gradient(self, c1, c2, res=50):
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
        
        