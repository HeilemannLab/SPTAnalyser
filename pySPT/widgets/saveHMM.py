"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Save hmm results.
"""

import os
import math
import numpy as np
import pandas as pd
import ermine as em
import matplotlib.pyplot as plt
import seaborn as sns
from graphviz import Digraph


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
                    Ds[i].append(n_states_Ds[j][i][0] / 1000000)
                for i in range(len(n_states_Ds[0])):
                    ws[i].append(n_states_weights[j][i])

            weight_names = ["weight" + str(i) for i in range(1, n_states + 1)]
            D_names = ["D" + str(i) + " [um2/s]" for i in range(1, n_states + 1)]

            df = pd.DataFrame()
            df["name"] = self.init_hmm.names
            for column_name, w in zip(weight_names, ws):
                df[column_name] = w
            for column_name, D in zip(D_names, Ds):
                df[column_name] = D

            row = ["average"]
            row.extend([i for i in self.init_hmm.mean_weights[c]])
            row.extend([i[0]/1000000 for i in self.init_hmm.mean_diffs[c]])

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
                    component_string = str(" $\\omega$ = %.2f, D = %.2e [um²/s]" % (
                        self.init_hmm.models_weights[state_idx][target_idx][i],
                        self.init_hmm.models_diffusion_coefficients[state_idx][target_idx][i]/1000000))
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


class SaveHMM():
    def __init__(self, save_path, input_path, hmm, figure_format, x_axis):
        self.save_path = save_path
        self.input_path = input_path
        self.hmm = hmm
        self.x_axis = x_axis
        self.figure_format = figure_format

    def create_folder(self):
        if not os.path.isdir(self.save_path):
            os.mkdir(self.save_path)

    def save_jdd_fits(self):
        plt.ioff()
        for target_idx in range(len(self.hmm.names)):
            fig = plt.figure()
            model_df = pd.DataFrame({"r": np.arange(0, self.x_axis, 1), "superposition": np.zeros(self.x_axis)})
            for i in range(self.hmm.n_states):
                unimodal_judi_model = em.JumpDistanceModel(
                    diffusion_coefficient=self.hmm.model_Ds[target_idx][i],
                    degrees_of_freedom=4,
                    tau=self.hmm.dt)
                component_string = str(" $\\omega$ = %.2f, D = %.2e [um²/s]" % (
                    self.hmm.model_ws[target_idx][i],
                    self.hmm.model_Ds[target_idx][i]/1000000))
                model_df[component_string] = self.hmm.model_ws[target_idx][i] *\
                                             unimodal_judi_model.pdf(distance=model_df["r"])
                model_df["superposition"] = model_df["superposition"] + model_df[component_string]

            sns.kdeplot(data=self.hmm.jump_dfs[target_idx], x="jump_distance", fill=True, bw_adjust=0.1, clip=[0, self.x_axis])
            sns.lineplot(data=model_df.melt(id_vars=['r']), x="r", y="value", color="black", style="variable")
            plt.title(self.hmm.names[target_idx])
            plt.xlabel("jump distance [nm]")
            plt.ylabel("density")

            subfolder = self.input_path + "\\SPTAnalyser_" + os.path.splitext(os.path.splitext(self.hmm.names[target_idx])[0])[0]
            if not os.path.isdir(subfolder):
                os.mkdir(subfolder)
            if not os.path.isdir(subfolder + "\\hmm"):
                os.mkdir(subfolder + "\\hmm")
            fig.savefig(subfolder + "\\hmm" + "\\jdd_fit_hmm_diff." + self.figure_format,
                        format=self.figure_format, transparent=True, bbox_inches="tight")
            plt.close(fig)

    def save_hmm_param_plots(self):
        figure_names = ["diffusion_coefficients", "transition_probabilities", "state_probabilities", "weights",
                        "diffusion_coefficients_corr"]
        for figure, name in zip(self.hmm.figures, figure_names):
            figure.savefig(self.save_path + "\\" + name + "." + self.figure_format,
                           format=self.figure_format, transparent=True, bbox_inches="tight")

    def save_state_transition_diagram(self, choose_weights):
        if choose_weights == "jdd fit weights":
            weights = self.hmm.mean_model_ws
        elif choose_weights == "occurrence":
            weights = self.hmm.mean_viterbi_counts
        elif choose_weights == "state probabilities":
            weights = self.hmm.mean_model_stateprobs

        min_node_size = 1.5  # label are too large for node, the smallest % has to fit in the node
        mult_with_node_size = min_node_size / min(weights)
        edge_fontsize = "10"
        dot = Digraph(comment="state transition diagram", format=self.figure_format)  # access format!
        float_precision = "%.3f"
        var_width = True
        colored_edges = True

        # return the state % between 0-100 %
        mean_population_rounded = [str(self.hmm.percentage_rounded(x)) for x in weights]
        # state population is represented as area of circle, get its radius as input: A = pi * r^2 -> r = (A/pi)**0.5
        mean_node_size = list(
            map(lambda x: float_precision % (float(x) * mult_with_node_size / math.pi) ** (0.5), weights))
        diffusions = self.hmm.mean_model_Ds
        # represent D in um^2/s with certain float precision
        diffusions = list(map(lambda x: str(float_precision % x), diffusions))

        for i in range(len(diffusions)):
            dot.node(str(i + 1), mean_population_rounded[i] + "%",
                     color="black", fillcolor=self.hmm.color_palette[i], fontcolor="black",
                     style="filled", shape="circle", fixedsize="shape", width=mean_node_size[i], pos="0,0!")

        columns = [i for i in range(self.hmm.n_states)] * self.hmm.n_states
        rows = [[i] * len(self.hmm.mean_model_Ds) for i in range(len(self.hmm.mean_model_Ds))]
        rows = [item for sublist in rows for item in sublist]

        tps = []
        for i in range(np.shape(self.hmm.mean_model_tps)[0]):
            for j in range(np.shape(self.hmm.mean_model_tps)[1]):
                tps.append(self.hmm.mean_model_tps[i][j])

        for tp, column, row in zip(tps, rows, columns):
            label_name = " " + str(self.hmm.percentage_rounded(tp)) + "%"
            dot.edge(str(column + 1), str(row + 1), label=" " + label_name,
                     color=(self.hmm.gv_edge_color_gradient(self.hmm.color_palette[column], self.hmm.color_palette[row],
                                                   25) if colored_edges else "black"),
                     fontsize=edge_fontsize, style="filled",
                     penwidth=(str(self.hmm.tp_px_mapping(tp, row, column)) if var_width else "1"))
        dot.render(self.save_path + r"\\state_transition_diagram." + self.figure_format, view=False)

    def save_scores(self):
        x = pd.concat(self.hmm.model_scores)
        x.insert(0, "name", self.hmm.names)
        x.to_csv(self.save_path + "\\model_scores.csv", index=False)

    def save_viterbi(self):
        for target_idx, sequence in enumerate(self.hmm.state_sequences):
            subfolder = self.input_path + "\\SPTAnalyser_" + \
                        os.path.splitext(os.path.splitext(self.hmm.names[target_idx])[0])[0]
            if not os.path.isdir(subfolder):
                os.mkdir(subfolder)
            if not os.path.isdir(subfolder + "\\hmm"):
                os.mkdir(subfolder + "\\hmm")
            sequence.to_csv(subfolder + "\\hmm" + "\\state_sequence.csv", index=False)

    def save_hmm_params(self):
        Ds = [[] for _ in range(len(self.hmm.model_Ds))]
        ws = [[] for _ in range(len(self.hmm.model_ws))]
        stateprobs = [[] for _ in range(len(self.hmm.model_stateprobs))]
        stateocc = [[] for _ in range(len(self.hmm.viterbi_counts))]
        if self.hmm.model_Ds_corr:
            Ds_corr = [[] for _ in range(len(self.hmm.model_Ds_corr))]

        for j in range(len(self.hmm.names)):
            for i in range(self.hmm.n_states):
                Ds[i].append(self.hmm.model_Ds[j][i][0] / 1000000)
            for i in range(self.hmm.n_states):
                ws[i].append(self.hmm.model_ws[j][i])
            for i in range(self.hmm.n_states):
                stateprobs[i].append(self.hmm.model_stateprobs[j][i])
            for i in range(self.hmm.n_states):
                stateocc[i].append(self.hmm.viterbi_counts[j][i])
            if self.hmm.model_Ds_corr:
                for i in range(self.hmm.n_states):
                    Ds_corr[i].append(self.hmm.model_Ds_corr[j][i][0] / 1000000)
        tps = [[] for i in range(self.hmm.n_states ** 2)]
        for tp in self.hmm.model_tps:
            c = 0
            for i in range(self.hmm.n_states):
                for j in range(self.hmm.n_states):
                    tps[c].append(tp[i][j])
                    c += 1

        weight_names = ["weight" + str(i) for i in range(1, self.hmm.n_states + 1)]
        stateprobs_names = ["state prob" + str(i) for i in range(1, self.hmm.n_states + 1)]
        occurrence_names = ["occurrence" + str(i) for i in range(1, self.hmm.n_states + 1)]
        D_names = ["D" + str(i) + " [um2/s]" for i in range(1, self.hmm.n_states + 1)]
        D_corr_names = ["D_corr" + str(i) + " [um2/s]" for i in range(1, self.hmm.n_states + 1)]
        tp_names = []
        for i in range(self.hmm.n_states):
            for j in range(self.hmm.n_states):
                tp_names.append([i + 1, j + 1])
        tp_names = [str(i[0]) + "\u279D" + str(i[1]) for i in tp_names]

        df = pd.DataFrame()
        df["name"] = self.hmm.names
        for column_name, w in zip(weight_names, ws):
            df[column_name] = w
        for column_name, sp in zip(stateprobs_names, stateprobs):
            df[column_name] = sp
        for column_name, oc in zip(occurrence_names, stateocc):
            df[column_name] = oc
        for column_name, tp in zip(tp_names, tps):
            df[column_name] = tp
        for column_name, D in zip(D_names, Ds):
            df[column_name] = D
        if self.hmm.model_Ds_corr:
            for column_name, D_corr in zip(D_corr_names, Ds_corr):
                df[column_name] = D_corr

        row = ["average"]
        row.extend([i for i in self.hmm.mean_model_ws])
        row.extend([i for i in self.hmm.mean_model_stateprobs])
        row.extend(list(self.hmm.mean_viterbi_counts))
        mean_tps = []
        for i in range(self.hmm.n_states):
            for j in range(self.hmm.n_states):
                mean_tps.append(self.hmm.mean_model_tps[i][j])
        row.extend(mean_tps)
        row.extend([i[0]/1000000 for i in self.hmm.mean_model_Ds])
        if self.hmm.model_Ds_corr:
            row.extend([i[0] / 1000000 for i in self.hmm.mean_model_Ds_corr])
        df.loc[len(df)] = row
        save_name = self.save_path + "\\hmm_params.csv"
        df.to_csv(save_name, index=False)

    def save(self, choose_weights):
        self.create_folder()
        self.save_hmm_param_plots()
        self.save_state_transition_diagram(choose_weights)
        self.save_jdd_fits()
        self.save_scores()
        self.save_viterbi()
        self.save_hmm_params()
