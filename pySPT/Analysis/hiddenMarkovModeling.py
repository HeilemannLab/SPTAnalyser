"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Initial parameter estimation and analysis of hidden markov modeling.
"""

import os
import ermine as em
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm_notebook as tqdm

# TODO: Check D units ym² or nm²?

class InitHMM():
    """Get inital parameters for HMM by fitting jump distance distributions"""
    def __init__(self):
        self.dir_path, self.dt, self.n_start, self.n_end = "", 0.0, 0, 0
        self.names, self.jump_dfs = [], []
        self.model_dfs, self.information_criteria_plot = [], []
        self.n_states, self.n_BIC, self.n_AIC, self.n_AICc = [], [], [], []
        self.models_diffusion_coefficients, self.models_weights = [], []
        self.mean_diffs, self.mean_weights = [], []

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
            diff_coeff = self.models_diffusion_coefficients[i]
            weights = self.models_weights[i]
            self.mean_diffs.append(np.mean(self.models_diffusion_coefficients[i], axis=0))
            self.mean_weights.append(np.mean(self.models_weights[i], axis=0))
            print("n states =", n)
            for c, (diff, w) in enumerate(zip(self.mean_diffs[i], self.mean_weights[i]), 1):
                print("average w"+str(c)+":", str(w)+";", "average D"+str(c)+" [nm²/s]:", str(diff[0]))
            print("\n")

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
            component_string = str(" $\\omega$ = %.2f, D = %.2e [nm²/s]" % (
            self.models_weights[state_idx][target_idx][i],
            self.models_diffusion_coefficients[state_idx][target_idx][i]))
            model_df[component_string] = self.models_weights[state_idx][target_idx][i] * unimodal_judi_model.pdf(
                distance=model_df["r"])
            model_df["superposition"] = model_df["superposition"] + model_df[component_string]

        sns.kdeplot(data=self.jump_dfs[target_idx], x="jump_distance", fill=True, bw_adjust=0.1, clip=[0, x_range])
        sns.lineplot(data=model_df.melt(id_vars=['r']), x="r", y="value", color="black", style="variable")
        plt.title(cell)
        plt.xlabel("jump distance [nm]")
        plt.ylabel("density")
        plt.show()

