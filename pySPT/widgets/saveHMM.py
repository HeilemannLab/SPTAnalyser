"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Save hmm results.
"""

import os
import pandas as pd
import numpy as np
import ermine as em
import seaborn as sns
import matplotlib.pyplot as plt



class SaveInitHMM():
    def __init__(self, save_path, input_path, init_hmm, figure_format, x_axis):
        self.save_path = save_path
        self.input_path = input_path
        self.init_hmm = init_hmm
        self.x_axis = x_axis
        self.figure_format = figure_format

    def create_folder(self):
        os.mkdir(self.save_path)

    def save_init_params(self):
        for c, n_states in enumerate(range(self.init_hmm.n_start, self.init_hmm.n_end+1)):
            n_states_Ds = self.init_hmm.models_diffusion_coefficients[c]
            n_states_weights = self.init_hmm.models_weights[c]

            Ds = [[] for _ in range(len(n_states_Ds[0]))]
            ws = [[] for _ in range(len(n_states_Ds[0]))]
            for j in range(len(self.init_hmm.names)):
                for i in range(len(n_states_Ds[0])):
                    Ds[i].append(n_states_Ds[j][i][0])
                for i in range(len(n_states_Ds[0])):
                    ws[i].append(n_states_weights[j][i])

            weight_names = ["weight" + str(i) for i in range(1, n_states + 1)]
            D_names = ["D" + str(i) + " [nm2/s]" for i in range(1, n_states + 1)]

            df = pd.DataFrame()
            df["name"] = self.init_hmm.names
            for column_name, w in zip(weight_names, ws):
                df[column_name] = w
            for column_name, D in zip(D_names, Ds):
                df[column_name] = D

            row = ["average"]
            row.extend([i for i in self.init_hmm.mean_weights[c]])
            row.extend([i[0] for i in self.init_hmm.mean_diffs[c]])

            df.loc[len(df)] = row

            save_name = self.save_path + "\\initial_params_" + str(n_states) + "states.csv"
            df.to_csv(save_name, index=False)

    def save_information_criteria(self):
        save_name = self.save_path + "\\information_criteria.csv"
        self.init_hmm.model_dfs.to_csv(save_name, index=False)
        self.init_hmm.information_criteria_plot.savefig(self.save_path + "\\information_criteria." + self.figure_format,
                                                        format=self.figure_format, transparent=True, bbox_inches="tight")

    def save_jdd_fits(self):
        plt.ioff()
        for target_idx in range(len(self.init_hmm.names)):
            for n_states in range(self.init_hmm.n_start, self.init_hmm.n_end + 1):
                fig = plt.figure()
                state_idx = range(self.init_hmm.n_start, self.init_hmm.n_end + 1).index(n_states)
                model_df = pd.DataFrame({"r": np.arange(0, self.x_axis, 1), "superposition": np.zeros(self.x_axis)})
                for i in range(n_states):
                    unimodal_judi_model = em.JumpDistanceModel(
                        diffusion_coefficient=self.init_hmm.models_diffusion_coefficients[state_idx][target_idx][i],
                        degrees_of_freedom=4,
                        tau=self.init_hmm.dt)
                    component_string = str(" $\\omega$ = %.2f, D = %.2e [nm²/s]" % (
                        self.init_hmm.models_weights[state_idx][target_idx][i],
                        self.init_hmm.models_diffusion_coefficients[state_idx][target_idx][i]))
                    model_df[component_string] = self.init_hmm.models_weights[state_idx][target_idx][i] *\
                                                 unimodal_judi_model.pdf(distance=model_df["r"])
                    model_df["superposition"] = model_df["superposition"] + model_df[component_string]

                sns.kdeplot(data=self.init_hmm.jump_dfs[target_idx], x="jump_distance", fill=True, bw_adjust=0.1, clip=[0, self.x_axis])
                sns.lineplot(data=model_df.melt(id_vars=['r']), x="r", y="value", color="black", style="variable")
                plt.title(self.init_hmm.names[target_idx])
                plt.xlabel("jump distance [nm]")
                plt.ylabel("density")

                subfolder = self.input_path + "\\SPTAnalyser_" + os.path.splitext(os.path.splitext(self.init_hmm.names[target_idx])[0])[0]
                if not os.path.isdir(subfolder):
                    os.mkdir(subfolder)
                if not os.path.isdir(subfolder + "\\hmm"):
                    os.mkdir(subfolder + "\\hmm")
                fig.savefig(subfolder + "\\hmm" + "\\jdd_fit_" + str(n_states) + "." + self.figure_format,
                            format=self.figure_format, transparent=True, bbox_inches="tight")
                plt.close(fig)

    def save(self):
        self.create_folder()
        self.save_init_params()
        self.save_information_criteria()
        self.save_jdd_fits()
