"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Initial parameter estimation and analysis of hidden markov modeling.
"""

import os
import math
import copy
import ermine as em
import numpy as np
from numpy.typing import ArrayLike
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from graphviz import Digraph
from tqdm import tqdm_notebook as tqdm


class Color:
    def __init__(self, r, g, b):
        self.r = r
        self.g = g
        self.b = b

    @staticmethod
    def from_hex(hexcode):
        hexcode = hexcode[1:]
        r = int(hexcode[:2], 16) / 255.0
        g = int(hexcode[2:4], 16) / 255.0
        b = int(hexcode[4:], 16) / 255.0
        return Color(r, g, b)

    def to_hex(self):
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

    def scalar_mult(self, scalar):
        r = self.r * scalar
        g = self.g * scalar
        b = self.b * scalar
        return Color(r, g, b)

    def clamp(self):
        self.r = max(0.0, min(1.0, self.r))
        self.g = max(0.0, min(1.0, self.g))
        self.b = max(0.0, min(1.0, self.b))


class InitHMM():
    """Get inital parameters for HMM by fitting jump distance distributions"""
    def __init__(self):
        self.dir_path, self.dt, self.n_start, self.n_end = "", 0.0, 0, 0
        self.names, self.jump_dfs = [], []
        self.model_dfs, self.information_criteria_plot = [], []
        self.n_states, self.n_BIC, self.n_AIC, self.n_AICc = [], [], [], []
        self.models_diffusion_coefficients, self.models_weights = [], []
        self.mean_diffs, self.mean_weights = [], []

    def clear_init(self):
        self.__init__()

    def read_files(self):
        """
        Get tracked.csv files from directory and return target names and jump distance information.
        """
        names, jump_dfs = [], []
        for file in os.listdir(self.dir_path):
            if file.endswith("csv") and "tracked" in file:
                file_path = self.dir_path + "\\" + file
                data_df = pd.read_csv(filepath_or_buffer=file_path)
                jump_dfs.append(em.preprocess_swift_data(data_df))
                names.append(file)
        return names, jump_dfs

    def calc_scores(self):
        """
        Calculate information criteria for all n states and return frame with cell names and results.
        """
        model_dfs = []
        self.models_diffusion_coefficients = [[] for _ in range(self.n_start, self.n_end+1)]
        self.models_weights = [[] for _ in range(self.n_start, self.n_end+1)]
        with tqdm(total=len(self.names), desc="Cells") as pbar:
            for jump_df in self.jump_dfs:
                count = 0
                x = jump_df["jump_distance"].values
                jmm = em.JumpDistanceMixtureModel(n_components=self.n_start,
                                                  degrees_of_freedom=4,
                                                  tau=self.dt,
                                                  init_params="wm",
                                                  params="wm")
                jmm.fit(x)
                model_df = pd.DataFrame(jmm.evaluate(x))
                self.models_diffusion_coefficients[count].append(jmm.diffusion_coefficients())
                self.models_weights[count].append(jmm._weights)
                count += 1
                for i in range(self.n_start+1, self.n_end+1):
                    jmm = em.JumpDistanceMixtureModel(n_components=i,
                                                      degrees_of_freedom=4,
                                                      tau=self.dt,
                                                      init_params="wm",
                                                      params="wm")
                    jmm.fit(x)
                    model_df = model_df.append(pd.DataFrame(jmm.evaluate(x)))
                    self.models_diffusion_coefficients[count].append(jmm.diffusion_coefficients())
                    self.models_weights[count].append(jmm._weights)
                    count += 1
                model_dfs.append(model_df)
                pbar.update(1)

        for name, model in zip(self.names, model_dfs):
            model.insert(0, "name", [name for _ in range(self.n_start, self.n_end+1)])
        model_dfs = pd.concat(model_dfs, ignore_index=True)

        return model_dfs

    def get_best_model(self):
        """
        Return the best n state model suggestion based on information criteria.
        """
        best_model_BIC = [0 for _ in range(self.n_start, self.n_end+1)]
        best_model_AIC = [0 for _ in range(self.n_start, self.n_end+1)]
        best_model_AICc = [0 for _ in range(self.n_start, self.n_end+1)]
        best_model_n = [i for i in range(self.n_start, self.n_end+1)]
        for name in self.names:
            cell_model = self.model_dfs[self.model_dfs["name"] == name]
            min_BIC_idx = cell_model["BIC"].idxmin()
            min_AIC_idx = cell_model["AIC"].idxmin()
            min_AICc_idx = cell_model["AICc"].idxmin()
            best_model_BIC[best_model_n.index(cell_model["classes"][min_BIC_idx])] += 1
            best_model_AIC[best_model_n.index(cell_model["classes"][min_AIC_idx])] += 1
            best_model_AICc[best_model_n.index(cell_model["classes"][min_AICc_idx])] += 1
        return best_model_n, best_model_BIC, best_model_AIC, best_model_AICc

    def plot_best_model(self):
        """
        Count the best scores per n of states and display as barplot.
        """
        fig = plt.figure()
        plt.subplot(1, 3, 1)
        plt.bar(self.n_states, self.n_BIC / np.sum(self.n_BIC) * 100)
        plt.xlabel("n states")
        plt.ylabel("best score [%]")
        plt.xticks(range(self.n_start, self.n_end + 1))
        plt.title("BIC")

        plt.subplot(1, 3, 2)
        plt.bar(self.n_states, self.n_AIC / np.sum(self.n_AIC) * 100)
        plt.xlabel("n states")
        plt.ylabel("best score [%]")
        plt.xticks(range(self.n_start, self.n_end + 1))
        plt.title("AIC")

        plt.subplot(1, 3, 3)
        plt.bar(self.n_states, self.n_AICc / np.sum(self.n_AICc) * 100)
        plt.xlabel("n states")
        plt.ylabel("best score [%]")
        plt.xticks(range(self.n_start, self.n_end + 1))
        plt.title("AICc")

        plt.tight_layout()
        plt.suptitle("Recommended number of states", y=1.03)
        plt.show()
        self.information_criteria_plot = fig

    def run_scores(self):
        """
        Main function to get initial parameters for defined n range
        """
        self.names, self.jump_dfs = self.read_files()
        self.model_dfs = self.calc_scores()
        self.n_states, self.n_BIC, self.n_AIC, self.n_AICc = self.get_best_model()
        self.plot_best_model()

    def get_average_params(self):
        """
        Get mean values of weights and diffusion coefficients per n states.
        """
        for i , n in enumerate(range(self.n_start, self.n_end + 1)):
            self.mean_diffs.append(np.mean(self.models_diffusion_coefficients[i], axis=0))

            self.mean_weights.append(np.mean(self.models_weights[i], axis=0))
            print("n states =", n)
            for c, (diff, w) in enumerate(zip(self.mean_diffs[i], self.mean_weights[i]), 1):
                print("averaged w"+str(c)+":", str(w)+";", "averaged D"+str(c)+" [um²/s]:", str(diff[0]/1000000))
            print("\n")
        np.divide(self.mean_diffs, 1000000)

    def plot_cell_results(self, n_states, cell, x_range):
        """
        Plot the jump distribution fit per cell (individual Ds and weights).
        """
        target_idx = self.names.index(cell)
        state_idx = range(self.n_start, self.n_end + 1).index(n_states)
        model_df = pd.DataFrame({"r": np.arange(0, x_range, 1), "superposition": np.zeros(x_range)})
        for i in range(n_states):
            unimodal_judi_model = em.JumpDistanceModel(
                diffusion_coefficient=self.models_diffusion_coefficients[state_idx][target_idx][i],
                degrees_of_freedom=4,
                tau=self.dt)
            component_string = str(" $\\omega$ = %.2f, D = %.2e [um²/s]" % (
            self.models_weights[state_idx][target_idx][i],
            self.models_diffusion_coefficients[state_idx][target_idx][i]/1000000))
            model_df[component_string] = self.models_weights[state_idx][target_idx][i] * unimodal_judi_model.pdf(
                distance=model_df["r"])
            model_df["superposition"] = model_df["superposition"] + model_df[component_string]

        sns.kdeplot(data=self.jump_dfs[target_idx], x="jump_distance", fill=True, bw_adjust=0.1, clip=[0, x_range])
        sns.lineplot(data=model_df.melt(id_vars=['r']), x="r", y="value", color="black", style="variable")
        plt.title(cell)
        plt.xlabel("jump distance [nm]")
        plt.ylabel("density")
        plt.show()


class HMM():
    def __init__(self):
        self.dir_path, self.dt, self.n_states, self.min_len, self.color_palette = "", 0, 0, 0, []
        self.names, self.jump_dfs = [], []
        self.init_Ds, self.init_ws, self.init_tps = [], [], []
        self.model_Ds, self.model_Ds_corr, self.model_stateprobs, self.model_tps, self.model_ws = [], [], [], [], []
        self.mean_model_Ds, self.mean_model_Ds_corr, self.mean_model_stateprobs, self.mean_model_tps, self.mean_model_ws = [], [], [], [], []
        self.std_model_Ds, self.std_model_Ds_corr, self.std_model_stateprobs, self.std_model_tps, self.std_model_ws = [], [], [], [], []
        self.sem_model_Ds, self.sem_model_Ds_corr, self.sem_model_stateprobs, self.sem_model_tps, self.sem_model_ws = [], [], [], [], []
        self.model_scores = []
        self.state_sequences = []
        self.mean_viterbi_counts, self.viterbi_counts = [], []
        self.trainable_params = ""
        self.color_palette = []
        self.figures = []

    def clear_init(self):
        self.__init__()

    def read_files(self):
        """
        Get tracked.csv files from directory and return target names and jump distance information.
        """
        names, jump_dfs = [], []
        for file in os.listdir(self.dir_path):
            if file.endswith("csv") and "tracked" in file:
                file_path = self.dir_path + "\\" + file
                data_df = pd.read_csv(filepath_or_buffer=file_path)
                data_df = data_df[data_df["track.lifetime"]+1 >= self.min_len]
                jump_dfs.append(em.preprocess_swift_data(data_df, min_track_length=self.min_len-1))
                names.append(file)
        return names, jump_dfs

    @staticmethod
    def create_observation_sequence(judi_df: pd.DataFrame) -> (ArrayLike, ArrayLike):
        """
        Creates an observation sequence of single particle tracking trajectories for the analysis with ErmineHMM.
        Parameters
        ----------
        judi_df : pd.DataFrame
            Pandas DataFrame object that comprises information on "jump_distance" and "track.id_departure".
        Returns
        -------
        (ArrayLike, ArrayLike)
            x: Jump distance observation sequence.
            lengths: Lengths of the individual sequences in x.
        """
        unique_track_id_vec = np.unique(judi_df["track.id_departure"].values)
        track_number = np.shape(unique_track_id_vec)[0]
        lengths = []
        x = np.ndarray([0, 1])
        judi_vec = judi_df["jump_distance"].values
        track_id_vec = judi_df["track.id_departure"].values
        frame_dep = judi_df["frame_departure"].values
        x_dep = judi_df["x [nm]_departure"].values
        y_dep = judi_df["y [nm]_departure"].values
        frame_dest = judi_df["frame_destination"].values
        x_dest = judi_df["x [nm]_destination"].values
        y_dest = judi_df["y [nm]_destination"].values
        track_id_dest = judi_df["track.id_destination"].values

        for i in np.arange(0, track_number, 1):
            idx = track_id_vec == unique_track_id_vec[i]
            track_x = np.expand_dims(judi_vec[idx], axis=0)
            lengths.append(np.shape(track_x)[1])
            x = np.concatenate([x, track_x.T])
        return (x, lengths, unique_track_id_vec, frame_dep, x_dep, y_dep, frame_dest, x_dest, y_dest, track_id_vec, track_id_dest)

    def fit_hmm(self):
        with tqdm(total=len(self.names), desc="Cells") as pbar:
            c = 0
            for jump_df in self.jump_dfs:
                hmm = em.ErmineHMM(n_components=self.n_states,
                                   diffusion_degrees_of_freedom=4,
                                   tau=self.dt,
                                   init_params="",
                                   params=self.trainable_params,
                                   n_iter=1000,
                                   tol=1e-5)
                hmm.startprob_ = self.init_ws
                hmm.transmat_ = self.init_tps
                hmm.diffusion_coefficients_ = self.init_Ds
                # x_hmm = jumps, lengths = trajectory length in frames
                x_hmm, lengths, unique_track_id_vec, frame_dep, x_dep, y_dep, frame_dest, x_dest, y_dest, track_id_vec, track_id_dest = self.create_observation_sequence(jump_df)
                hmm.fit(x_hmm, lengths)
                idx_sequence = []
                for i, item in zip(lengths, unique_track_id_vec):
                    idx_sequence += i * [item]
                self.model_stateprobs.append(hmm.startprob_)
                self.model_tps.append(hmm.transmat_)
                self.model_Ds.append(hmm.diffusion_coefficients_)
                self.model_scores.append(pd.DataFrame(hmm.evaluate(x_hmm, lengths)))

                state_sequence = hmm.predict(x_hmm, lengths)
                viterbi_df = pd.DataFrame({"idx": idx_sequence, "track.id_departure":track_id_vec, "frame_departure":frame_dep, "x [nm]_departure":x_dep, "y [nm]_departure":y_dep,
                                           "track.id_destination": track_id_dest, "frame_destination": frame_dest, "x [nm]_destination": x_dest, "y [nm]_destination": y_dest,
                                           "jump": x_hmm[:, 0],
                                           "state": state_sequence+1,
                                           "D": np.take(hmm.diffusion_coefficients_, state_sequence)/1000000})
                self.state_sequences.append(viterbi_df)
                pbar.update(1)
                # get viterbi state occurrence
                counts = viterbi_df["state"].value_counts(normalize=True)
                count_keys = list(counts.keys())
                count_keys.sort()
                zipped_lists = sorted(zip(list(counts), list(counts.keys())))
                counts_sorted = [percentage for percentage, count in zipped_lists]
                self.viterbi_counts.append(counts_sorted)
        self.mean_viterbi_counts = np.mean(self.viterbi_counts, axis=0)

    def fit_hmm_(self):
        with tqdm(total=len(self.names), desc="Cells") as pbar:
            for jump_df in self.jump_dfs:
                hmm = em.ErmineHMM(n_components=self.n_states,
                                   diffusion_degrees_of_freedom=4,
                                   tau=self.dt,
                                   init_params="",
                                   params=self.trainable_params,
                                   n_iter=1000,
                                   tol=1e-5)
                hmm.startprob_ = self.init_ws
                hmm.transmat_ = self.init_tps
                hmm.diffusion_coefficients_ = self.init_Ds
                # x_hmm = jumps, lengths = trajectory length in frames
                x_hmm, lengths = em.create_observation_sequence(jump_df)
                hmm.fit(x_hmm, lengths)
                self.model_stateprobs.append(hmm.startprob_)
                self.model_tps.append(hmm.transmat_)
                self.model_Ds.append(hmm.diffusion_coefficients_)
                self.model_scores.append(pd.DataFrame(hmm.evaluate(x_hmm, lengths)))

                state_sequence = hmm.predict(x_hmm, lengths)
                viterbi_df = pd.DataFrame({"idx": idx, "track_x":track_x, "jump": x_hmm[:, 0],
                                           "state": state_sequence+1,
                                           "D": np.take(hmm.diffusion_coefficients_, state_sequence)})
                self.state_sequences.append(viterbi_df)
                pbar.update(1)
        # get viterbi state occurrence
        frame = pd.concat(self.state_sequences)
        counts = frame["state"].value_counts(normalize=True)
        count_keys = list(counts.keys())
        count_keys.sort()
        zipped_lists = sorted(zip(list(counts), list(counts.keys())))
        counts_sorted = [percentage for percentage, count in zipped_lists]
        self.viterbi_counts = counts_sorted

    def get_average_results(self):
        """
        Get averaged hmm parameter values and print them.
        """
        self.mean_model_stateprobs = np.mean(self.model_stateprobs, axis=0)
        self.std_model_stateprobs = np.std(self.model_stateprobs, axis=0, ddof=1)
        self.sem_model_stateprobs = np.std(self.model_stateprobs, axis=0, ddof=1) / np.sqrt(len(self.names))
        self.mean_model_ws = np.mean(self.model_ws, axis=0)
        self.std_model_ws = np.std(self.model_ws, axis=0, ddof=1)
        self.sem_model_ws = np.std(self.model_ws, axis=0, ddof=1) / np.sqrt(len(self.names))
        self.mean_model_Ds = np.mean(self.model_Ds, axis=0)
        self.std_model_Ds = np.std(self.model_Ds, axis=0, ddof=1)
        self.sem_model_Ds = np.std(self.model_Ds, axis=0, ddof=1) / np.sqrt(len(self.names))
        self.mean_model_tps = np.mean(self.model_tps, axis=0)
        self.std_model_tps = np.std(self.model_tps, axis=0, ddof=1)
        self.sem_model_tps = np.std(self.model_tps, axis=0, ddof=1) / np.sqrt(len(self.names))
        print("Averaged HMM results:")
        for i in range(self.n_states):
            print("state probability" + str(i+1) + ":", str(self.mean_model_stateprobs[i]) + ";", "D" + str(i + 1) + " [um²/s]:", str(self.mean_model_Ds[i][0]/1000000))
        print("transition probabilities:")
        print(self.mean_model_tps)

    def violin_plot(self, dataframe, title, ylim, ylabel, xticks):
        if title == "transition probabilities":
            sns.set_palette(sns.color_palette("husl", 8))
        else:
            sns.set_palette(sns.color_palette(self.color_palette))
        fig = plt.figure()
        ax = sns.violinplot(data=dataframe, showmeans=True, inner="quartile",
                            meanprops={"marker": "s", "markerfacecolor": "white", "markeredgecolor": "0.25"})
        x_len = self.n_states ** 2 if title == "transition probabilities" else self.n_states
        ax = sns.scatterplot(x=[i for i in range(x_len)], y=dataframe.mean(), marker="X",
                             color=(19/255, 27/255, 30/255), edgecolor=(19/255, 27/255, 30/255), s=12)
        ax.set_title(title)
        ax.set(ylim=ylim)
        plt.ylabel(ylabel)
        if xticks:
            ax.set_xticklabels(xticks)
        plt.show()
        self.figures.append(fig)

    def plot_parameters(self):
        unicode_idx_subs = ["\u2081", "\u2082", "\u2083", "\u2084", "\u2085", "\u2086"]
        # diffusion coefficients
        Ds = []
        for i in range(self.n_states):
            Ds.append([d[i][0] for d in self.model_Ds])
        self.violin_plot(pd.DataFrame(Ds).T/1000000, "diffusion coefficients", (None, None), "diffusion coefficient [um\u00B2/s]",
                    ["D" + str(unicode_idx_subs[i]) for i in range(self.n_states)])
        # transition probabilities
        tps = [[] for i in range(self.n_states ** 2)]
        for tp in self.model_tps:
            c = 0
            for i in range(self.n_states):
                for j in range(self.n_states):
                    tps[c].append(tp[i][j])
                    c += 1
        x_labels = []
        for i in range(self.n_states):
            for j in range(self.n_states):
                x_labels.append([i + 1, j + 1])
        x_labels = [str(i[0]) + "\u279D" + str(i[1]) for i in x_labels]
        self.violin_plot(pd.DataFrame(tps).T, "transition probabilities", (None, None), "transition probabilitiy", x_labels)
        # state probabilities
        self.violin_plot(pd.DataFrame(self.model_stateprobs), "state probabilities", (0, 1), "state probability",
                         ["\u03C0" + str(unicode_idx_subs[i]) for i in range(self.n_states)])
        # weights
        self.violin_plot(pd.DataFrame(self.model_ws), "weights", (0, 1), "weights",
                         ["\u03C9" + str(unicode_idx_subs[i]) for i in range(self.n_states)])

    def get_jdd_weights(self):
        """
        Fit jdd with fixed diffusion coefficients to extract weights.
        """
        for c, jump_df in enumerate(self.jump_dfs):
            x_jmm = jump_df["jump_distance"].values
            jmm = em.JumpDistanceMixtureModel(n_components=self.n_states,
                                              degrees_of_freedom=4,
                                              tau=self.dt,
                                              init_params="",
                                              params="w")  # Notice, only the weights are optimized, diffusion is kept fix!
            # Since neither the weights nor the diffusion coefficients are initialized, we have to enter them manually.
            jmm.mu = em.postprocessing.calculate_expectation_value(diff_coeff=self.model_Ds[c],
                                                                   tau=self.dt,
                                                                   dof=4,
                                                                   sigma=0,  # We do not expect any errors in this tutorial
                                                                   epsilon=0)  # We do not expect any errors in this tutorial
            jmm.weights = self.model_stateprobs[c]
            jmm.fit(x_jmm)
            self.model_ws.append(jmm._weights)

    def plot_cell_fit(self, cell, x_range):
        """
        Plot the jump distribution fit (individual Ds and weights) and print hmm parameters of cell.
        """
        target_idx = self.names.index(cell)

        print("HMM results of cell:")
        for i in range(self.n_states):
            print("state probability" + str(i+1) + ":", str(self.model_stateprobs[target_idx][i]) + ";", "D" + str(i + 1) + " [um²/s]:", str(self.model_Ds[target_idx][i][0]/1000000))
        print("transition probabilities:")
        print(self.model_tps[target_idx])

        model_df = pd.DataFrame({"r": np.arange(0, x_range, 1), "superposition": np.zeros(x_range)})
        for i in range(self.n_states):
            unimodal_judi_model = em.JumpDistanceModel(
                diffusion_coefficient=self.model_Ds[target_idx][i],
                degrees_of_freedom=4,
                tau=self.dt)
            component_string = str(" $\\omega$ = %.2f, D = %.2e [um²/s]" % (
            self.model_ws[target_idx][i],
            self.model_Ds[target_idx][i]/1000000))
            model_df[component_string] = self.model_ws[target_idx][i] * unimodal_judi_model.pdf(
                distance=model_df["r"])
            model_df["superposition"] = model_df["superposition"] + model_df[component_string]

        sns.kdeplot(data=self.jump_dfs[target_idx], x="jump_distance", fill=True, bw_adjust=0.1, clip=[0, x_range])
        sns.lineplot(data=model_df.melt(id_vars=['r']), x="r", y="value", color="black", style="variable")
        plt.title(cell)
        plt.xlabel("jump distance [nm]")
        plt.ylabel("density")
        plt.show()

    def get_trainable_params(self, params):
        """
        based on boolean list, return trainable parameters (s=weights, t=transition prob, d=diffusion coeffs).
        """
        trainable_params = ""
        for c, i in enumerate(params):
            if c == 0 and i:
                trainable_params += "s"
            elif c == 1 and i:
                trainable_params += "t"
            elif c == 2 and i:
                trainable_params += "d"
        self.trainable_params = trainable_params

    @staticmethod
    def percentage_rounded(value):
        """
        Example: 0.4056357 -> 40.56.
        """
        exponent = 4
        while True:
            if value >= 0.5 * 10 ** (-exponent):
                value *= 10 ** exponent
                value = int(np.round(value))
                value = value / 10 ** (exponent - 2)
                return value
            exponent += 1

    @staticmethod
    def tp_px_mapping(tp, row, column, min_log=-5, min_px=0.5, max_px=5):
        """
        Get logarithmic value of tp, transform it to positive values in range 0-1 and clip it to reasonable arrow px sizes.
        """
        log_tp = np.log10(tp)
        if log_tp <= min_log:
            return min_px
        log_tp -= min_log  # positive range
        log_tp /= -min_log  # 0-1
        return max_px * log_tp + min_px * (1 - log_tp)  # sth between max & min px

    @staticmethod
    def gv_edge_color_gradient(c1, c2, res=50):
        """
        Interpolation between two colors c1 and c2 with res steps.
        """
        color_list = ""
        c1 = Color.from_hex(c1)
        c2 = Color.from_hex(c2)
        cDiff = c2 - c1
        weight = str(1.0 / (res + 1.0))
        color_list += c1.to_hex() + ";" + weight
        for i in range(res):
            color_list += ":"
            color_list += (c1 + cDiff.scalar_mult(i / res)).to_hex()
            color_list += ";" + weight
        return color_list

    def plot_state_transition_diagram(self, choose_weights, save_dir):
        if choose_weights == "jdd fit weights":
            weights = self.mean_model_ws
        elif choose_weights == "occurrence":
            weights = self.mean_viterbi_counts
        elif choose_weights == "state probabilities":
            weights = self.mean_model_stateprobs

        min_node_size = 1.5  # label are too large for node, the smallest % has to fit in the node
        mult_with_node_size = min_node_size / min(weights)
        edge_fontsize = "10"
        dot = Digraph(comment="state transition diagram", format="svg")  # access format!
        float_precision = "%.3f"
        var_width = True
        colored_edges = True

        # return the state % between 0-100 %
        mean_population_rounded = [str(self.percentage_rounded(x)) for x in weights]
        # state population is represented as area of circle, get its radius as input: A = pi * r^2 -> r = (A/pi)**0.5
        mean_node_size = list(
            map(lambda x: float_precision % (float(x) * mult_with_node_size / math.pi) ** (0.5), weights))
        diffusions = self.mean_model_Ds
        # represent D in um^2/s with certain float precision
        diffusions = list(map(lambda x: str(float_precision % x), diffusions))

        for i in range(len(diffusions)):
            dot.node(str(i + 1), mean_population_rounded[i] + "%",
                     color="black", fillcolor=self.color_palette[i], fontcolor="black",
                     style="filled", shape="circle", fixedsize="shape", width=mean_node_size[i], pos="0,0!")

        columns = [i for i in range(self.n_states)] * self.n_states
        rows = [[i] * len(self.mean_model_Ds) for i in range(len(self.mean_model_Ds))]
        rows = [item for sublist in rows for item in sublist]

        tps = []
        for i in range(np.shape(self.mean_model_tps)[0]):
            for j in range(np.shape(self.mean_model_tps)[1]):
                tps.append(self.mean_model_tps[i][j])

        for tp, column, row in zip(tps, rows, columns):
            label_name = " " + str(self.percentage_rounded(tp)) + "%"
            dot.edge(str(column + 1), str(row + 1), label=" " + label_name,
                     color=(self.gv_edge_color_gradient(self.color_palette[column], self.color_palette[row],
                                                   25) if colored_edges else "black"),
                     fontsize=edge_fontsize, style="filled",
                     penwidth=(str(self.tp_px_mapping(tp, row, column)) if var_width else "1"))
        dot.render(save_dir + r"\\State_transition_diagram", view=True)

    def D_to_MSD(self, diff):
        return 4*diff*self.dt

    def MSD_to_D(self, MSD):
        return MSD/4*self.dt

    def check_immobility(self, immobile_threshold):
        for c, diff in enumerate(self.mean_model_Ds):
            epsilon = em.postprocessing.static_error(apparent_msd_d0=self.D_to_MSD(diff), dof=4)[0]
            diff = diff[0]/1000000
            immobile_result = "immobile" if epsilon <= immobile_threshold else "not immobile"
            print("The localization precision of D{:.0f} = {:.5f} um²/s by NeNA is {:.2f} nm \u279D {:s}".format(c, diff, epsilon, immobile_result))

    def diffusion_correction(self, error):
        """
        Calculate the expected diffusion coefficient corrected for static and dynamic error.
        """
        self.model_Ds_corr = copy.deepcopy(self.model_Ds)
        for target_idx, diffs in enumerate(self.model_Ds):
            for d_idx, diff in enumerate(diffs):
                D_corrected = em.postprocessing.calculate_diffusion_coefficient(expected_value=self.D_to_MSD(diff),
                                                                        tau=self.dt,
                                                                        dof=4,
                                                                        sigma=self.dt,
                                                                        epsilon=error)
                self.model_Ds_corr[target_idx][d_idx] = D_corrected
        # calculate statistic
        self.mean_model_Ds_corr = np.mean(self.model_Ds_corr, axis=0)
        self.std_model_Ds_corr = np.std(self.model_Ds_corr, axis=0, ddof=1)
        self.sem_model_Ds_corr = np.std(self.model_Ds_corr, axis=0, ddof=1) / (self.n_states)**0.5

        for i in range(self.n_states):
            print("D_corr" + str(i+1) + " [\u00B5m²/s]:", str(self.mean_model_Ds_corr[i][0]/1000000))
        # display as plot
        unicode_idx_subs = ["\u2081", "\u2082", "\u2083", "\u2084", "\u2085", "\u2086"]
        Ds = []
        for i in range(self.n_states):
            Ds.append([d[i][0] for d in self.model_Ds_corr])
        self.violin_plot(pd.DataFrame(Ds).T/1000000, "corrected diffusion coefficients", (None, None), "corrected diffusion coefficient [um\u00B2/s]",
                    ["D" + str(unicode_idx_subs[i]) for i in range(self.n_states)])

    def run(self, params, choose_weights, input_dir):
        self.names, self.jump_dfs = self.read_files()
        self.get_trainable_params(params)
        self.fit_hmm()
        self.get_jdd_weights()
        self.get_average_results()
        self.plot_state_transition_diagram(choose_weights, input_dir)
        self.plot_parameters()

    def run_immob_check(self, immob_threshold):
        self.check_immobility(immob_threshold)

    def run_correction(self, epsilon):
        self.diffusion_correction(epsilon)
