"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Count the transitions of diffusion types of segments within trajectories.
"""


import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats


class TransitionCounts:
    """Transition counting between diffusion states applied to a sample."""
    def __init__(self, file_name, h5_file, tracked_file, n_diffusion_types, mask):
        """
        :param h5_file: h5 file object.
        :param tracked_file: tracked file from swift as pandas frame.
        :param n_diffusion_types: Number of diffusion type models (3=immobile, confined, free; 4=immobile, confined, free, "notype")
        :param mask: Segments < mask are ignored.
        :param cell_size: Size in um².
        :param min_length: Only trajectories >= min_length are considered.
        """
        self.h5_file = h5_file
        self.tracked_file = tracked_file
        self.n_diffusion_types = n_diffusion_types
        self.mask = mask
        self.cell_size = h5_file["settings"]["settings"][()][0][0][3]
        self.segment_ids, self.segment_booleans = self.extract_values("rossier", "rossierStatistics", 0, 5)
        self.segment_types = [self.booleans_to_type(str(i)) for i in self.segment_booleans]  # list of diff types
        self.segs_of_tracks, self.track_ids = self.get_segs_of_tracks()
        self.counts, self.counts_size_norm = self.count_transitions()
        self.trajectory_table = self.create_trajectory_table(os.path.splitext(file_name)[0])

    def extract_values(self, folder, frame, idx_min, idx_max):
        """
        Read dataframe of h5 file, get segment id and segment booleans describing diffusion type.
        """
        segment_ids = []
        segment_booleans = []
        rows = self.h5_file[folder][frame][()]
        for row in rows:
            segment_id = []
            segment_boolean = []
            for i in range(idx_min, idx_max):
                if i == idx_min:
                    segment_id.append(row[i])
                else:
                    segment_boolean.append(row[i])
            segment_ids.append(segment_id[0])
            segment_booleans.append(segment_boolean)
        return segment_ids, segment_booleans

    def booleans_to_type(self, boolean):
        """
        Convert boolean to str, regard immobile and notype separately or combined.
        """
        if self.n_diffusion_types == 4:
            dct = {str([1, 0, 0, 0]): "immobile",
                   str([0, 1, 0, 1]): "confined",
                   str([0, 0, 1, 1]): "free",
                   str([0, 0, 0, 0]): "notype"}
        elif self.n_diffusion_types == 3:
            dct = {str([1, 0, 0, 0]): "immobile",
                   str([0, 1, 0, 1]): "confined",
                   str([0, 0, 1, 1]): "free",
                   str([0, 0, 0, 0]): "immobile"}
        type = dct[boolean]
        return type

    def get_segs_of_tracks(self):
        """
        For each track, get its seg ids:
        - if at least one segment is classified.
        - filter segments shorter than mask value.
        """
        segs_of_tracks, track_ids = [], []
        track_file = self.tracked_file[["track.id", "seg.id", "seg.lifetime"]].sort_values(["track.id", "seg.id"])
        track_file = track_file[track_file["seg.lifetime"] > 0].to_numpy()
        for i in range(np.max(track_file, axis=0)[0] + 1):
            filtered_segs = []
            any_classified_seg = False
            track = track_file[track_file[:,0] == i]
            for seg_id in np.unique(track[:,1]):
                seg_lifetime = np.unique(track[track[:,1] == seg_id][:,2])[0]
                if seg_lifetime + 1 > self.mask:
                    filtered_segs.append(seg_id)
                    if seg_id in self.segment_ids:
                        any_classified_seg = True
            if any_classified_seg:
                segs_of_tracks.append(np.asarray(filtered_segs))
                track_ids.append(i)
        return segs_of_tracks, track_ids

    def count_transitions(self):
        """
        Get a list of counts, idx refers to transition type of chosen transition dct, last entry is transition to None.
        Only trajectories are regarded that contain at least one segment of filtered dataset.
        """
        if self.n_diffusion_types == 4:
            transition_dct = {"['immobile', 'immobile']": 0,
                              "['immobile', 'confined']": 1,
                              "['immobile', 'free']": 2,
                              "['immobile', 'notype']": 3,
                              "['confined', 'immobile']": 4,
                              "['confined', 'confined']": 5,
                              "['confined', 'free']": 6,
                              "['confined', 'notype']": 7,
                              "['free', 'immobile']": 8,
                              "['free', 'confined']": 9,
                              "['free', 'free']": 10,
                              "['free', 'notype']": 11,
                              "['notype', 'immobile']": 12,
                              "['notype', 'confined']": 13,
                              "['notype', 'free']": 14,
                              "['notype', 'notype']": 15}
        elif self.n_diffusion_types == 3:
            transition_dct = {"['immobile', 'immobile']": 0,
                              "['immobile', 'confined']": 1,
                              "['immobile', 'free']": 2,
                              "['confined', 'immobile']": 3,
                              "['confined', 'confined']": 4,
                              "['confined', 'free']": 5,
                              "['free', 'immobile']": 6,
                              "['free', 'confined']": 7,
                              "['free', 'free']": 8}
        counts = [0 for _ in range(len(transition_dct)+1)]

        for track in self.segs_of_tracks:
            if len(track) > 1:  # trajectories with at least two segments.
                seg_diffusion_types = []
                for i in track:
                    if i in self.segment_ids:
                        seg_diffusion_types.append(self.segment_types[self.segment_ids.index(i)])
                    else:
                        seg_diffusion_types.append(None)
                for i in range(0, len(seg_diffusion_types)-1):
                    transition = seg_diffusion_types[i:i+2]
                    if None in transition:
                        counts[-1] += 1
                    else:
                        counts[transition_dct[str(transition)]] += 1
        counts_size_norm = [i / self.cell_size for i in counts]
        return counts, counts_size_norm

    def create_trajectory_table(self, target_name):
        """
        Store transition counting information in table with columns:
        target.name, track.id, seg.id, seg.state, seg.lifetime, track.lifetime, x [nm], y [nm]
        """
        data = self.tracked_file[["track.id", "seg.id", "seg.lifetime", "track.lifetime", "frame", "x [nm]", "y [nm]"]]
        data = data[data["track.id"].isin(self.track_ids)]  # filter for chosen track ids
        data.insert(0, "target.name", [target_name]*len(data))  # add target.name column
        data = data.sort_values(["track.id", "seg.id", "frame"])

        seg_diff_states = []
        type_dict = {self.segment_ids[i]: self.segment_types[i] for i in range(len(self.segment_ids))}
        for _, row in data.iterrows():
            seg_id = row["seg.id"]
            t = type_dict[seg_id] if seg_id in type_dict else "none"
            seg_diff_states.append(t)
        data.insert(3, "seg.state", seg_diff_states)  # add seg.state column
        return data


def save_counts(transition_counts, file_names, save_dir):
    """
    Save transition counts per target as txt-file.
    """
    counts = [i.counts for i in transition_counts]
    sizes = [i.cell_size for i in transition_counts]
    if len(counts[0]) == 10:
        transitions = ["i-i", "i-c", "i-f", "c-i", "c-c", "c-f", "f-i", "f-c", "f-f", "x-none"]
    elif len(counts[0]) == 17:
        transitions = ["i-i", "i-c", "i-f", "i-n", "c-i", "c-c", "c-f", "c-n", "f-i", "f-c", "f-f", "f-n",
                       "n-i", "n-c", "n-f", "n-n", "x-none"]
    data = [list(e) for e in zip(*counts)]
    out_file_name = save_dir + "\\transition_counts.txt"
    header = "name\tsize[\u00B5m\u00B2]\t" + "\t".join(transitions)
    dtype = [("col1", "U"+str(max([len(i) for i in file_names]))), ("col2", float)]
    fmt = ["%"+str(max([len(i) for i in file_names]))+"s", "%.4e"]
    for i in range(len(transitions)):
        dtype.append(("col"+str(i+3), float))
        fmt.append("%.4e")
    save_data = np.zeros(np.array(file_names).shape, dtype=dtype)
    save_data["col1"] = np.array(file_names)
    save_data["col2"] = np.array(sizes)
    for i in range(len(transitions)):
        save_data["col"+str(i+3)] = np.array(data[i])
    np.savetxt(out_file_name, X=save_data, fmt=tuple(fmt), header=header, delimiter="\t")


def save_trajectory_tables(tables, save_dir):
    """
    Save trajectory table as txt-file.
    """
    data = pd.concat(tables)
    data.to_csv(save_dir + "\\trajectory_table.txt", sep="\t", index=False)


def save_mask_value(mask_val, save_dir):
    """
    Store mask value as txt-file.
    """
    file = open(save_dir + "\\mask_val.txt", "w+")
    file.write("Mask value: " + mask_val)


class Statistic:
    """Transition count information pulled from available samples."""
    def __init__(self):
        self.counts = []  # counts
        self.trajectory_tables = []
        self.sizes = []  # target sizes
        self.n_diffusion_types = 3  # number of diffusion models
        self.save_dir = ""

    def vis_counts(self, counts="absolute", norm="global", ylim=[None,None]):
        """
        Visualize distribution of transition counts.
        :param counts: Take absolute count values (=absolute) or normalize counts by cell size (=size_norm)
        :param norm: Normalize counts per cell (=global) or per diffusion type (=split) or not at all (=absolute)
        :param ylim: min and max y-axis-limit of visualization
        """
        if counts == "size_norm":
            data = []
            for count, size in zip(self.counts, self.sizes):
                data.append([i/size for i in count][:-1])  # without None
        elif counts == "absolute":
            data = [i[:-1] for i in self.counts]  # without None

        # choose color palette (for the 3 or 4 diffusion states)
        rgb_colors = [(64, 105, 224), (33, 140, 33), (225, 140, 0), (140, 0, 140)]
        rgb_colors = [[c / 255 for c in color] for color in rgb_colors]
        if len(data[0]) == 9:
            colors = []
            [colors.extend(3 * [i]) for i in rgb_colors[:-1]]
            colors = [tuple(i) for i in colors]
            transitions = ["i-i", "i-c", "i-f", "c-i", "c-c", "c-f", "f-i", "f-c", "f-f"]
        elif len(data[0]) == 16:
            colors = []
            [colors.extend(4 * [i]) for i in rgb_colors]
            colors = [tuple(i) for i in colors]
            transitions = ["i-i", "i-c", "i-f", "i-n", "c-i", "c-c", "c-f", "c-n", "f-i", "f-c", "f-f", "f-n",
                           "n-i", "n-c", "n-f", "n-n"]

        if norm == "global":
            data = [target / np.sum(target) for target in data]
            y_label = "normalized counts per cell" if counts == "absolute" else "normalized counts per area"
        elif norm == "split":
            data = []
            for target in self.counts:
                split = [target[i:i + self.n_diffusion_types] for i in range(0, len(target), self.n_diffusion_types)]
                split_norm_vals = [np.sum(i) for i in split]
                norm_vals = []
                for i in split_norm_vals:
                    norm_vals.extend(i for _ in range(self.n_diffusion_types))
                split_norms = [i / norm_val for i, norm_val in zip(target, norm_vals)]
                data.append(split_norms)
            y_label = "normalized counts per diffusion state"
        elif norm == "absolute":
            y_label = "counts per cell" if counts == "absolute" else "counts per area"

        data_dct = {}
        for c, transition in enumerate(transitions):
            data_dct[transition] = [i[c] for i in data]
        data_df = pd.DataFrame(data_dct)

        # visualizations
        save_name = "transition_counts_" + counts + "_" + "norm_" + norm
        self.plot_counts_violinplot(data_df, y_label, colors, ylim, save_name)
        self.significance_test(data, transitions, save_name)
        print("*"*48)
        # saving
        data_df.to_csv(self.save_dir + "\\" + save_name + ".txt", sep="\t", index=False)

    def plot_counts_violinplot(self, data_df,  y_label, colors, ylim, save_name):
        """Display counts as violinplots"""
        sns.violinplot(data=data_df, palette=colors, scale="width", cut=0, inner="quartile")
        sns.scatterplot(x=[i for i in range(9)], y=data_df.mean(), marker="X", color=(19/255,27/255,30/255),
                        edgecolor=(19/255,27/255,30/255), s=12)
        plt.xlabel("transition types"); plt.ylabel(y_label);
        plt.ylim([ylim[0], ylim[1]])
        plt.savefig(self.save_dir + "\\distributions_violin_qr_" + save_name + ".svg")
        plt.show()

    def plot_counts_matrix(self, data_df, title, save_name, average="mean", cmap="YlGnBu_r", vmin=0, vmax=1):
        """Display average counts as matrix"""
        x_axis_labels = ["i", "c", "f", "n"][:self.n_diffusion_types]
        y_axis_labels = ["i", "c", "f", "n"][:self.n_diffusion_types]

        if average == "mean":
            data = data_df.mean().to_numpy().reshape(self.n_diffusion_types, self.n_diffusion_types)
        else:
            data = data_df.median().to_numpy().reshape(self.n_diffusion_types, self.n_diffusion_types)

        ax = sns.heatmap(data=data, xticklabels=x_axis_labels, yticklabels=y_axis_labels, annot=True, vmin=vmin,
                         vmax=vmax, cmap=cmap)
        ax.xaxis.tick_top()
        plt.title(average + " " + title)
        plt.savefig(self.save_dir + "\\matrix_" + save_name + ".svg")
        plt.show()

    def segment_lengths_plot(self, xlim=[0,None], ylim=[None, None]):
        """Lengths of segments, visualized as histogram"""
        seg_lengths = []
        max_len = 0
        for cell in self.trajectory_tables:
            seg_lengths_cell = [k+1 for k, g in itertools.groupby(cell["seg.lifetime"])]
            if max(seg_lengths_cell) > max_len:
                max_len = max(seg_lengths_cell)
            seg_lengths.append(seg_lengths_cell)
        histograms = [np.histogram(i, bins=[x for x in range(max_len)]) for i in seg_lengths]
        mean_counts = np.mean([i[0] for i in histograms], axis=0)
        norm_mean_counts = mean_counts / np.sum(mean_counts)
        SEM_error = np.std([i[0] for i in histograms], axis=0, ddof=1)/(len(seg_lengths))**0.5 / np.sum(mean_counts)

        plt.bar(histograms[0][1][:-1], norm_mean_counts, yerr=SEM_error, width=1,
                edgecolor="black", capsize=0.5)
        plt.xlabel("segment length [frame]"); plt.ylabel("normalized counts"); plt.title("Segment lengths per cell")
        plt.xlim([xlim[0],xlim[1]])
        plt.ylim([ylim[0],ylim[1]])
        plt.savefig(self.save_dir + "\\segment_lengths.svg")
        plt.show()
        # saving
        header = [i.iloc[0,0] for i in self.trajectory_tables]
        pd.DataFrame(seg_lengths).transpose().to_csv(self.save_dir + "\\segment_lengths.txt", sep="\t", header=header, index=False)

    def segments_per_trajectory_plot(self, xlim=[0,None], ylim=[0, None]):
        """Number of segments per trajectory, visualize as histogram"""
        n_segs = []
        max_n = 0
        for target in self.trajectory_tables:
            n_segs_target = []
            for track_id in target["track.id"].unique():
                n = len(target[target["track.id"] == track_id]["seg.id"].unique())
                if n > max_n:
                    max_n = n
                n_segs_target.append(n)
            n_segs.append(n_segs_target)
        histograms = [np.histogram(i, bins=[x for x in range(max_n)]) for i in n_segs]
        mean_counts = np.mean([i[0] for i in histograms], axis=0)
        norm_mean_counts = mean_counts / np.sum(mean_counts)
        SEM_error = np.std([i[0] for i in histograms], axis=0, ddof=1)/(len(n_segs))**0.5 / np.sum(mean_counts)
        plt.bar(histograms[0][1][:-1], norm_mean_counts, yerr=SEM_error, width=1, edgecolor="black", capsize=0.5, color="tab:blue")
        plt.xlabel("number of segments per track"); plt.ylabel("normalized counts")
        plt.xticks(np.arange(0, max_n + 1, 1.0))
        plt.xlim([xlim[0],xlim[1]])
        plt.ylim(([ylim[0],ylim[1]]))
        plt.savefig(self.save_dir + "\\n_segments_per_track.svg")
        plt.show()
        # save
        header = [i.iloc[0,0] for i in self.trajectory_tables]
        pd.DataFrame(n_segs).transpose().to_csv(self.save_dir + "\\n_segments_per_track.txt", sep="\t", header=header, index=False)

    def transitions_wo_none_plot(self, ylim=[0,1]):
        """Average percentage of transitions involving not classified vs without classified segments, visualized as violinplot"""
        data = self.counts
        data = [target / np.sum(target) for target in data]
        data = [[np.sum(counts[:-1]), counts[-1]] for counts in data]
        data_dct = {}
        transitions = ["without none", "with none"]
        for c, transition in enumerate(transitions):
            data_dct[transition] = [i[c] for i in data]

        sns.violinplot(data=pd.DataFrame(data_dct), inner="quartile")
        plt.xlabel(""); plt.ylabel("normalized counts"); plt.title("Transitions without / with none")
        plt.ylim([ylim[0], ylim[1]])
        plt.savefig(self.save_dir + "\\transitions_wo_none_violin.svg")
        plt.show()
        # saving
        pd.DataFrame(data_dct).to_csv(self.save_dir + "\\transitions_wo_none.txt", sep="\t", index=False)

    def significance_test(self, data, transition_names, save_name):
        data = [list(e) for e in zip(*data)]
        # normality tests
        stats_shapiro, ps_shapiro, results_shapiro, stats_ks, ps_ks, results_ks = [], [], [], [], [], []
        for transition, transition_name in zip(data, transition_names):
            # Shapiro-Wilk, W=1, p>0.05 -> norm
            stat_shapiro, p_shapiro = scipy.stats.shapiro(transition)
            result_shapiro = "norm" if p_shapiro > 0.05 else "not norm"
            stats_shapiro.append(stat_shapiro); ps_shapiro.append(p_shapiro); results_shapiro.append(result_shapiro)
            # Kolmogorov-Smirnov, D=0, p>0.05 -> norm
            stat_ks, p_ks = scipy.stats.kstest(transition, "norm")
            result_ks = "norm" if p_ks > 0.05 else "not norm"
            stats_ks.append(stat_ks); ps_ks.append(p_ks); results_ks.append(result_ks)

        out_file_name = self.save_dir + "\\normality_test_" + save_name + ".txt"
        header = "transition type\tShapiro statistic\tShapiro p\tShapiro result\tKolmogorov-Smirnov statistic\t" \
                 "Kolmogorov-Smirnov p\tKolmogorov-Smirnov result\t"
        save_data = np.zeros(np.array(transition_names).size, dtype=[("col1", "U3"), ("col2", float), ("col3", float),
                                                                     ("col4", "U8"), ("col5", float), ("col6", float),
                                                                     ("col7", "U8")])
        save_data["col1"] = np.array(transition_names)
        save_data["col2"] = np.array(stats_shapiro)
        save_data["col3"] = np.array(ps_shapiro)
        save_data["col4"] = np.array(results_shapiro)
        save_data["col5"] = np.array(stats_ks)
        save_data["col6"] = np.array(ps_ks)
        save_data["col7"] = np.array(results_ks)
        np.savetxt(out_file_name, X=save_data, fmt=("%3s", "%.4e", "%.4e", "%.8s", "%.4e", "%.4e", "%.8s"),
                   header=header, delimiter="\t")

        # paired t-test, wilcoxon signed-rank test, p<0.05 -> *, p<0.01 -> **, p<0.001 -> ***
        unique_transition_names, stats_t, ps_t, results_t, stats_wilc, ps_wilc, results_wilc = [], [], [], [], [], [], []
        for i, i_name in zip(data, transition_names):
            for j, j_name in zip(data, transition_names):
                if not i == j:
                    if [i_name, j_name] not in unique_transition_names and [j_name, i_name] not in unique_transition_names:
                        unique_transition_names.append([i_name, j_name])
                        stat_t, p_t = scipy.stats.ttest_rel(i, j)
                        stats_t.append(stat_t); ps_t.append(p_t); results_t.append(self.p_to_sig_results(p_t))
                        stat_wilc, p_wilc = scipy.stats.wilcoxon(i, j)
                        stats_wilc.append(stat_wilc); ps_wilc.append(p_wilc); results_wilc.append(self.p_to_sig_results(p_wilc))
        unique_transition_names = [str(i[0]) + "-" + str(i[1]) for i in unique_transition_names]

        out_file_name = self.save_dir + "\\significance_test_" + save_name + ".txt"
        header = "targets\tWilcoxon statistic\tWilcoxon p\tWilcoxon result\tpaired t-test statistic\t" \
                 "paired t-test p\tpaired t-test result\t"
        save_data = np.zeros(np.array(unique_transition_names).size, dtype=[("col1", "U7"), ("col2", float),
                                                                            ("col3", float), ("col4", "U4"),
                                                                            ("col5", float), ("col6", float),
                                                                            ("col7", "U4")])
        save_data["col1"] = np.array(unique_transition_names)
        save_data["col2"] = np.array(stats_t)
        save_data["col3"] = np.array(ps_t)
        save_data["col4"] = np.array(results_t)
        save_data["col5"] = np.array(stats_wilc)
        save_data["col6"] = np.array(ps_wilc)
        save_data["col7"] = np.array(results_wilc)
        np.savetxt(out_file_name, X=save_data, fmt=("%7s", "%.4e", "%.4e", "%.4s", "%.4e", "%.4e", "%.4s"),
                   header=header, delimiter="\t")

    def p_to_sig_results(self, p_value, sig_levels=(0.05, 0.01, 0.001)):
        if p_value < sig_levels[2]:
            return "***"
        elif p_value < sig_levels[1]:
            return "**"
        elif p_value < sig_levels[0]:
            return "*"
        else:
            return "n.s."

    def visualize_cell_trajectories(self):
        pass
