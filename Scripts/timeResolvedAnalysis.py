"""
@author: Alexander Niedrig
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
Calculates mean values over timeframes, plots them, runs statistical tests, and rearranges input data in an output file
"""
import configparser
import math
import os
import shutil
import sys
import time
import warnings

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as scy


class IncorrectConfigException(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)


def get_matching_files(directory, target, exclusion_string):
    """
    Search in directory and its subdirectories for files with target in their name but missing an exclusion string
    :param directory: all files in this directory & subdirectories are checked
    :param target: substring of filename
    :param exclusion_string: excludes files that have any of the strings in this list
    :return: List of matching file paths
    """
    matching_files = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if target.lower() in name.lower():
                if not any([True for string in exclusion_string if string.lower() in name.lower()]):
                    matching_files.append(os.path.join(path, name))
    return matching_files


def parent_string(substring, parents, exclusionstr):
    """
    goes through a list of strings, searching for a specific substring
    :param substring: the string to search for
    :param parents: a list of potential parent strings
    :param exclusionstr: a list of strings that must not be in the parent string
    :return: the parent string
    """
    parent = [a for a in parents if
              (substring in a and exclusionstr not in a)]
    if len(parent) > 1:
        raise IncorrectConfigException("multiple options when looking for " + substring)
    elif len(parent) == 0:
        raise IncorrectConfigException("no file found when looking for " + substring)
    elif len(parent) == 1:
        parent = parent[0]
    return parent


def sort_cells(cells):
    """
    :param cells: a list of cell names
    :return: cells sorted in alphanumerical order taking multi digit numbers into account
    """
    digits = [[]]
    for cell in cells:
        number = cell.split("_")[-1].split(".")[0]
        while True:
            try:
                digits[len(number)].append(cell)
                break
            except IndexError:
                digits.append([])
                continue
    sorted = []
    for dig in digits:
        dig.sort()
        sorted += dig
    return sorted


def insert_error(meanframe, sdframe, semframe):
    """
    :param meanframe: a dataframe of calculated mean values
    :param sdframe: a dataframe of calculated standard deviations
    :param a dataframe of calculated standard errors of means
    :return: a dataframe with SD and SEM columns inserted after their mean
    """
    sdframe = sdframe.iloc[:, 2:]
    semframe = semframe.iloc[:, 2:]
    i = 0
    while i < sdframe.shape[1]:
        name = sdframe.columns[i]
        sd = sdframe.iloc[:, i]
        sem = semframe.iloc[:, i]
        sd.rename(name + '_SD', inplace=True)
        sem.rename(name + '_SEM', inplace=True)
        meanframe.insert(i * 3 + 3, name + '_SD', sd)
        meanframe.insert(i * 3 + 4, name + '_SEM', sem)
        i += 1
    return meanframe


def calc_mean_over_cS(binned_data, attribute):
    """
    :param binned_data: a list of binned dataframes
    :param attribute: the attribute to calculate the values over
    :return: numpy stack of the mean, standard deviation and standard error of mean of a certain attribute
    """
    attributes = pd.DataFrame()
    for i, df in enumerate(binned_data):
        columns = [col for col in df.columns if attribute in col and not col.endswith(("_SD", '_SEM'))]
        for col in columns:
            attributes = pd.concat([attributes, df[col]], axis=1)
    stats = np.column_stack([attributes.mean(axis=1), attributes.std(axis=1), attributes.sem(axis=1)])
    return stats


def plot_by_time(dataframes, attribute, t_lig, ligand_name, binned_data, error_type, clr):
    """
    :param dataframes: all datafames to be plotted
    :param attribute: the attribute to be plotted
    :param t_lig: time the ligand was added
    :param ligand_name: name of the ligand
    :param binned_data: list of binned data per coverslip
    :param error_type: whether to use standard deviation or SEM
    :param clr: color to plot the datapoints in
    :return: the created plot
    plots the given dataframes into a graph and creates bars of the mean values given over a timeframe
    """
    fix, ax = plt.subplots()
    # plots dots
    colors = []
    for frame in dataframes:
        color = ax.plot(frame.iloc[:, 1] / 60, frame[attribute], marker='o', ms=2, linestyle='None', color=clr)[
            0].get_color()
    # plots mean bars and errors
    stats = calc_mean_over_cS(binned_data, attribute)
    left = 0
    longestdf = ''
    for df in binned_data:
        if len(df) > len(longestdf):
            longestdf = df
    for i, row in enumerate(stats):
        # reads the error after it has been calculated in the calc_mean_over_cS function
        if error_type == 'SEM':
            error = row[2]
        else:
            error = row[1]
        if i == len(stats) - 1:
            right = float(longestdf.iloc[i, 1].split('-')[1]) / 60
        else:
            right = (((float(longestdf.iloc[i + 1, 1].split('-')[0]) - float(
                longestdf.iloc[i, 1].split('-')[1])) / 2) + float(longestdf.iloc[i, 1].split('-')[1])) / 60
        ax.barh(row[0],
                width=right - left,
                height=ax.get_ylim()[1] * 0.005, left=left, alpha=1,
                align='center', color='black')
        ax.barh(row[0],
                width=right - left,
                height=error * 2, left=left, alpha=0.5,
                align='center', color='grey')
        left = right
    # Y axis labeling
    atr = attribute.split('_')[0]
    if atr == 'P':
        ax.set_ylabel(attribute + ' [%]')
    if atr == 'D':
        ax.set_ylabel(attribute + r' [$\mu m^2 s^{-1}$]')
    if atr == 'L':
        ax.set_ylabel(attribute + ' [frames]')
    if atr == 'N':
        ax.set_ylabel(attribute)
    if atr == 'confinement_radius':
        ax.set_ylabel('confinement_radius' + r'$\mu m$')
    ax.set_xlim(left=-0.2)
    ax.set_xlabel('Time [s]')

    # add ligand bar
    if ligand_name == '':
        ligand_name = '+LIG'
    if len(t_lig) > 0 and t_lig[-1] == 's':
        ax.axvline(float(t_lig[:-1]) / 60, color='black')
        ax.text(float(t_lig[:-1]) / 60, ax.get_ylim()[1], ligand_name, ha='center', va='bottom', color='black')

    return plt


def plot_by_cells(dataframes, attribute, t_lig, ligand_name, binned_data, error_type, clr):
    """
    :param dataframes: all datafames to be plotted
    :param attribute: the attribute to be plotted
    :param t_lig: time the ligand was added
    :param ligand_name: name of the ligand
    :param binned_data: list of binned data per coverslip
    :param error_type: whether to use standard deviation or SEM
    :param clr: color to plot the datapoints in
    :return: the created plot
    plots the given dataframes into a graph and creates bars of the mean values given over a frame of cells
    """
    fig, ax = plt.subplots()

    # Plot dots and error bars
    for frame in dataframes:
        ax.plot([float(cell_name.split("_")[-1]) for cell_name in frame.iloc[:, 0]], frame[attribute], marker='o',
                linestyle='None', color=clr)[0].get_color()

    # Plot mean bars and errors
    stats = calc_mean_over_cS(binned_data, attribute)
    left = 0
    lens = [int(j.iloc[-1, 0].split('-')[1]) for j in binned_data]
    max_index = lens.index(max(lens))
    maxframe = binned_data[max_index]
    for i, row in enumerate(stats):
        if error_type == 'SEM':
            error = row[2]
        else:
            error = row[1]
        right = float(maxframe.iloc[i, 0].split('-')[1]) + 0.5
        ax.barh(row[0], width=right - left, height=ax.get_ylim()[1] * 0.005, left=left, alpha=1, align='center',
                color='black')
        ax.barh(row[0], width=right - left, height=error * 2, left=left, alpha=0.5, align='center', color='grey')
        left = right

    # Y axis labeling
    atr = attribute.split('_')[0]
    if atr == 'P':
        ax.set_ylabel(attribute + ' [%]')
    if atr == 'D':
        ax.set_ylabel(attribute + r' [$\mu m^2 s^{-1}$]')
    if atr == 'L':
        ax.set_ylabel(attribute + ' [frames]')
    if atr == 'N':
        ax.set_ylabel(attribute)
    if atr == 'confinement_radius':
        ax.set_ylabel('confinement_radius' + r' [$\mu m$]')
    ax.set_xlim(left=-0.2)
    ax.set_xlabel('Cells')

    if len(t_lig) > 0 and t_lig[-1] == 'c':
        ax.axvline((float(t_lig[:-1])), color='black')
        if ligand_name == '':
            ligand_name = '+LIG'
        ax.text((float(t_lig[:-1])), ax.get_ylim()[1], ligand_name, ha='center', va='bottom', color='black')

    return plt


def bin_data_cells(frame, bin_size):
    """
    :param frame: given dataframe
    :param bin_size: what interals to create the means over
    :return: frame with the means, SD and SEM for each time interval
    compiles data to ranges of means by their cell number, creates statistics over them.
    """
    meanframe = pd.DataFrame(columns=frame.columns)
    sdframe = pd.DataFrame(columns=frame.columns)
    semframe = pd.DataFrame(columns=frame.columns)
    cs_name = frame.iloc[0, 0].split('_')[0] + '_' + frame.iloc[0, 0].split('_')[-3] + '_'
    binnumber = 0
    # appends to the frames
    rowindex = 0
    complete = False
    while not complete:
        start = rowindex
        bin = pd.DataFrame(columns=frame.columns)
        while float(frame.iloc[rowindex, 0].split('_')[
                        -1]) <= binnumber * bin_size:  # runs until the time is greater than that of the current bin, appends data to a df
            bin.loc[len(bin)] = frame.iloc[rowindex, :]
            rowindex += 1
            if rowindex == len(frame):  # if the end of the df is reached, end iteration
                complete = True
                break
        if len(bin) == 0:
            binnumber += 1
            continue

        # adds mean, cell name interval, and time interval
        meanframe.loc[len(meanframe)] = bin.iloc[:, 2:].mean()
        meanframe.iloc[-1, 0] = cs_name + "cell_" + str((binnumber - 1) * bin_size) + "-" + str(binnumber * bin_size)
        meanframe.iloc[-1, 1] = str(int(bin.iloc[0, 1])) + "-" + str(int(bin.iloc[-1, 1]))
        # adds SD, cell name interval, and time interval
        sdframe.loc[len(sdframe)] = bin.iloc[:, 2:].std()
        sdframe.iloc[-1, 0] = cs_name + "cell_" + str((binnumber - 1) * bin_size) + "-" + str(binnumber * bin_size)
        sdframe.iloc[-1, 1] = str(int(bin.iloc[0, 1])) + "-" + str(int(bin.iloc[-1, 1]))
        # adds SEM, cell name interval, and time interval
        semframe.loc[len(semframe)] = bin.iloc[:, 2:].sem()
        semframe.iloc[-1, 0] = cs_name + "cell_" + str((binnumber - 1) * bin_size) + "-" + str(binnumber * bin_size)
        semframe.iloc[-1, 1] = str(int(bin.iloc[0, 1])) + "-" + str(int(bin.iloc[-1, 1]))

    binnumber += 1

    outframe = insert_error(meanframe, sdframe, semframe)
    return outframe


def bin_data_time(frame, bin_size):
    """
    :param frame: given dataframe
    :param bin_size: what intervals to create the means over
    :return: frame with the means, SD and SEM for each time interval
    compiles data to ranges of means by their time index, creates statistics over them.
    """
    meanframe = pd.DataFrame(columns=frame.columns)
    sdframe = pd.DataFrame(columns=frame.columns)
    semframe = pd.DataFrame(columns=frame.columns)
    cs_name = frame.iloc[0, 0].split('_')[0] + '_' + frame.iloc[0, 0].split('_')[-3] + '_'
    complete = False
    binnumber = 1
    rowindex = 0
    while not complete:
        start = rowindex
        bin = pd.DataFrame(columns=frame.columns)
        while frame.iloc[
            rowindex, 1] < binnumber * bin_size:  # runs until the time is greater than that of the current bin, appends data to a df
            bin.loc[len(bin)] = frame.iloc[rowindex, :]
            rowindex += 1
            if rowindex == len(frame):  # if the end of the df is reached, end iteration
                complete = True
                break
        if len(bin) == 0:
            binnumber += 1
            continue
        # adds mean and cell name interval
        meanframe.loc[len(meanframe)] = bin.iloc[:, 2:].mean()
        meanframe.iloc[-1, 0] = cs_name + "cell_" + str(start) + "-" + str(start + len(bin))
        # adds mean and cell name interval
        sdframe.loc[len(sdframe)] = bin.iloc[:, 2:].std()
        sdframe.iloc[-1, 0] = cs_name + "cell_" + str(start) + "-" + str(start + len(bin))
        # adds mean and cell name interval
        semframe.loc[len(semframe)] = bin.iloc[:, 2:].sem()
        semframe.iloc[-1, 0] = cs_name + "cell_" + str(start) + "-" + str(start + len(bin))
        # adds time intervals
        if len(bin) > 0:
            meanframe.iloc[-1, 1] = str((binnumber - 1) * bin_size) + "-" + str((binnumber) * bin_size)
            sdframe.iloc[-1, 1] = str((binnumber - 1) * bin_size) + "-" + str((binnumber) * bin_size)
            semframe.iloc[-1, 1] = str((binnumber - 1) * bin_size) + "-" + str((binnumber) * bin_size)

        binnumber += 1
    outframe = insert_error(meanframe, sdframe, semframe)
    return outframe


def test_by_cell(frames, bin_size, ligand, alpha, p1, p2, p3):
    """
    Runs statistic test over bins of a certain number of cells
    :param frames: dataframes to test
    :param bin_size: size of the bins to test (in cells)
    :param ligand: whether a ligand was used and thus whether significance tests are required
    :param alpha: alpha interval for normality test results
    :param p1: * interval for significance
    :param p2: ** interval for significance
    :param p3: *** interval for significance
    :return: dataframes of the test results
    """
    normality_frames = {}
    if ligand:
        significance_frames = {}
    cell_count=[]
    bin_number = int(len(frames[0]) / bin_size)
    if len(frames[0]) % bin_size != 0:
        bin_number += 1
    for j in range(bin_number):
        tempFrame = pd.DataFrame()
        for frame in frames:
            binframe=pd.DataFrame(columns=tempFrame.columns)
            for index,row in frame.iterrows():
                if float(row[0].split('_')[-1]) > j * bin_size and float(row[0].split('_')[-1]) <= (j+1) * bin_size:
                    binframe=binframe.append(row,ignore_index=True)
            tempFrame = pd.concat([tempFrame, binframe],
                                  axis=0)  # all data of a bin is collected here
        try:
            cell_count[j]+= len(tempFrame)
        except IndexError:
            cell_count.append(len(tempFrame))
        if len(tempFrame) < 3:
            for name in tempFrame.columns[2:]:
                normality_frames[name].loc[len(normality_frames[name])] = pd.Series(
                    [j + 1, str(j * bin_size + 1) + '-' + str((j + 1) * bin_size), 'NaN', 'NaN',
                     'Dataset too small (' + str(len(tempFrame)) + ')',
                     'NaN', 'NaN',
                     'Dataset too small (' + str(len(tempFrame)) + ')']).values
                significance_frames[name].loc[len(significance_frames[name])] = pd.Series(
                    ['1-' + str(j + 1), str(j * bin_size + 1) + '-' + str((j + 1) * bin_size), 'NaN', 'NaN',
                     'Dataset too small (' + str(len(tempFrame)) + ')',
                     'NaN', 'NaN',
                     'Dataset too small (' + str(len(tempFrame)) + ')']).values
            continue
        for name in tempFrame.columns[2:]:
            tempFrame2 = tempFrame.dropna(subset=[name])
            if name not in normality_frames.keys():  # creates new frame if necessary
                compare_frame = tempFrame2
                norm_frame = pd.DataFrame(
                    columns=['BinNumber', 'Cellrange', 'number of cells','Shapiro statistic', 'Shapiro p', 'Shapiro result',
                             'Kolmogorov-Smirnov statistic', 'Kolmogorov-Smirnov p',
                             'Kolmogorov-Smirnov result'])
                normality_frames[name] = norm_frame
            # runs the statistic tests
            shapstat, ps_value, kolstat, pk_value = normality_tests(tempFrame2, name)
            if ps_value < alpha:
                sr = 'not norm'
            else:
                sr = 'norm'
            if pk_value < alpha:
                kr = 'not norm'
            else:
                kr = 'norm'
            normality_frames[name].loc[len(normality_frames[name])] = pd.Series(
                [j + 1, str(j * bin_size + 1) + '-' + str((j + 1) * bin_size),str(len(tempFrame2)) + ' '+str(cell_count[j] - len(tempFrame2))+' were dropped due to NaN entries', shapstat, ps_value, sr, kolstat,
                 pk_value, kr]).values

            if ligand:
                if name not in significance_frames.keys():
                    sign_frame = pd.DataFrame(
                        columns=['Compared Bins', 'Cellrange', 'number of cells (control: '+str(cell_count)+' cells)', 'paired tTest statistic', 'paired tTest p',
                                 'paired tTest result',
                                 'Wilcoxon-signed-rank statistic', 'Wilcoxon p',
                                 'Wilcoxon result'])
                    significance_frames[name] = sign_frame
                else:
                    mstat, pm_value, wilstat, pw_value = significance_tests(tempFrame2,
                                                                            compare_frame,
                                                                            name)
                    if pm_value == 'NaN':
                        tr = 'test not applicable'
                    elif pm_value < p1:
                        if pm_value < p2:
                            if pm_value < p3:
                                tr = '***'
                            else:
                                tr = '**'
                        else:
                            tr = '*'
                    else:
                        tr = 'no significant difference'
                    if pw_value == 'NaN':
                        wr = 'test not applicable'
                    elif pw_value < p1:
                        if pw_value < p2:
                            if pw_value < p3:
                                wr = '***'
                            else:
                                wr = '**'
                        else:
                            wr = '*'
                    else:
                        wr = 'no significant difference'
                    significance_frames[name].loc[len(significance_frames[name])] = pd.Series(
                        ['1-' + str(j + 1), str(j * bin_size + 1) + '-' + str((j + 1) * bin_size),str(len(tempFrame2)) + ' '+str(cell_count[j] -len(tempFrame2))+' were dropped due to NaN entries', mstat, pm_value, tr,
                         wilstat, pw_value, wr]).values

    if ligand:
        return normality_frames, significance_frames
    else:
        return normality_frames


def test_by_time(frames, bin_size, ligand, alpha, p1, p2, p3):
    """
    Runs statistic test over bins of a certain number of cells
    :param frames: dataframes to test
    :param bin_size: size of the bins to test (in seconds)
    :param ligand: whether a ligand was used and thus whether significance tests are required
    :param alpha: alpha interval for normality test results
    :param p1: * interval for significance
    :param p2: ** interval for significance
    :param p3: *** interval for significance
    :return: dataframes of the test results
    """
    normality_frames = {}
    if ligand:
        significance_frames = {}

    max_times = []
    for frame in frames:
        max_times.append(max(frame.iloc[:, 1]))
    max_time = max(max_times)
    bin_number = int(max_time / bin_size + 1)
    if max_time % bin_size != 0:
        bin_number += 1
    for i in range(1, bin_number):
        tempFrame = pd.DataFrame()
        for frame in frames:
            rowindex = 0
            while frame.iloc[rowindex, 1] < i * bin_size:
                if frame.iloc[rowindex, 1] >= (i - 1) * bin_size:
                    tempFrame = tempFrame.append(frame.iloc[rowindex, :])
                rowindex += 1
                if rowindex == len(frame):
                    break
        if len(tempFrame) < 3:
            for name in tempFrame.columns[2:]:
                normality_frames[name].loc[len(normality_frames[name])] = pd.Series(
                    [i + 1, str(min(tempFrame['Time'])) + '-' + str(max(tempFrame['Time'])), 'NaN', 'NaN',
                     'Dataset too small', 'NaN', 'NaN', 'Dataset too small']).values
                if ligand:
                    significance_frames[name].loc[len(significance_frames[name])] = pd.Series(
                        ['1-' + str(i), str(min(tempFrame['Time'])) + '-' + str(max(tempFrame['Time'])), 'NaN', 'NaN',
                         'Dataset too small', 'NaN', 'NaN',
                         'Dataset too small']).values
            continue
        for name in tempFrame.columns[2:]:
            if name not in normality_frames.keys():
                compare_frame = tempFrame
                norm_frame = pd.DataFrame(
                    columns=['BinNumber', 'Timerange', 'Shapiro statistic', 'Shapiro p', 'Shapiro result',
                             'Kolmogorov-Smirnov statistic', 'Kolmogorov-Smirnov p',
                             'Kolmogorov-Smirnov result'])
                normality_frames[name] = norm_frame
            shapstat, ps_value, kolstat, pk_value = normality_tests(tempFrame.iloc[:, 2:], name)
            if ps_value < alpha:
                sr = 'not norm'
            else:
                sr = 'norm'
            if pk_value < alpha:
                kr = 'not norm'
            else:
                kr = 'norm'
            normality_frames[name].loc[len(normality_frames[name])] = pd.Series(
                [i + 1, str(min(tempFrame['Time'])) + '-' + str(max(tempFrame['Time'])), shapstat, ps_value, sr,
                 kolstat, pk_value, kr]).values
            if ligand:
                if name not in significance_frames.keys():
                    sign_frame = pd.DataFrame(
                        columns=['Compared Bins', 'Timerange', 'Mann-Whitney U statistic', 'Mann-Whitney U p',
                                 'Mann-Whitney U result',
                                 'Wilcoxon-signed-rank statistic', 'Wilcoxon p',
                                 'Wilcoxon result'])
                    significance_frames[name] = sign_frame
                else:
                    tstat, pt_value, wilstat, pw_value = significance_tests(tempFrame.iloc[:, 2:],
                                                                            compare_frame,
                                                                            name)
                    if pt_value == 'NaN':
                        tr = 'test not applicable'
                    elif pt_value < p1:
                        if pt_value < p2:
                            if pt_value < p3:
                                tr = '***'
                            else:
                                tr = '**'
                        else:
                            tr = '*'
                    else:
                        tr = 'no significant difference'
                    if pw_value == 'NaN':
                        wr = 'test not applicable'
                    elif pw_value < p1:
                        if pw_value < p2:
                            if pw_value < p3:
                                wr = '***'
                            else:
                                wr = '**'
                        else:
                            wr = '*'
                    else:
                        wr = 'no significant difference'
                    significance_frames[name].loc[len(significance_frames[name])] = pd.Series(
                        ['1-' + str(i), str(min(tempFrame['Time'])) + '-' + str(max(tempFrame['Time'])), tstat,
                         pt_value, tr, wilstat, pw_value, wr]).values

    if ligand:
        return normality_frames, significance_frames
    else:
        return normality_frames


def normality_tests(dataframe, attribute):
    """
    :param dataframe: a given dataframe to test the normality for
    :param attribute: what attribute to test for
    :return: test statistics and values
    """
    shapstat, ps_value = scy.shapiro(dataframe[attribute])
    kolstat, pk_value = scy.shapiro(dataframe[attribute])
    return shapstat, ps_value, kolstat, pk_value


def significance_tests(frame1, frame2, attribute):
    """
    :param frame1: distribution to compare to
    :param frame2: distribution to compare
    :param attribute: attribute to compare
    :return: test statistics and values
    """
    if len(frame1) == len(frame2):
        try:
            mstat, pm_value = scy.ttest_rel(frame1[attribute], frame2[attribute])
        except ValueError:
            mstat, pm_value = 'not applicable to data', 'NaN'
        try:
            wilstat, pw_value = scy.wilcoxon(frame1[attribute], frame2[attribute])
        except ValueError:
            wilstat, pw_value = 'not applicable to data', 'NaN'
    else:
        mstat, pm_value, wilstat, pw_value = 'uneven sample size: ' + str(len(frame1)) + " vs " + str(
            len(frame2)), 'NaN', 'uneven sample size: ' + str(len(frame1)) + " vs " + str(len(frame2)), 'NaN'
    return mstat, pm_value, wilstat, pw_value


def compile_columns(frames, columns, renameColumns):
    """
    compiles columns of the same name in a single dataframe
    :param frames: list of dataframes
    :param columns: list of column names
    :param renameColumns: whether to rename columns to prevent duplicate column names
    :return: a single dataframe consisting of all matching columns
    in the given list
    """
    compiled_frame = pd.DataFrame()
    for df in frames:
        date = df.iloc[0, 0].split('_')[0] + '_' + df.iloc[0, 0].split('_')[-3] + ': '
        common_columns = [column for column in df.columns if any(sub in column for sub in columns)]
        comp_df = df[common_columns]

        if renameColumns:
            comp_df = comp_df.rename(columns={col: date + col for col in comp_df.columns})
        compiled_frame = pd.concat([compiled_frame, comp_df], axis=1)
    if not renameColumns:
        return comp_df
    return compiled_frame


def generate_shortname(value):
    """
    :param value: value of a dataframe field
    :return: the value without a colon at the end
    """
    if type(value) != str:
        return str(value)
    elif value == '':
        return 'None'
    elif value[-1] == '.':
        return value[:-1]
    else:
        return value


def cut_cellnames(name):
    """
    :param name: name of a cell field
    :return: the cell number
    """

    return name.split("_")[-1]


def output_folder(file, folder, datasets):
    """
    writes data into a sheet in a group of a h5 file
    :param file: the h5 file to write to
    :param folder: name for the folder to write the data into
    :param datasets: the dataframe to write
    """
    fold = file.create_group(folder)
    for i, set in enumerate(
            datasets):  # set consists of a 2 entry list with the name of the set in the first and the dataframe in the second entry
        if folder != 'metadata':
            set[1].columns = rename_columns(set[1].columns.tolist())
        for col in set[1].columns:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                set[1][col] = set[1][col].apply(generate_shortname)
        compound_dtype = np.dtype([(a, h5py.special_dtype(vlen=str)) for a in set[1].columns])
        tab = fold.create_dataset(set[0], (len(set[1]),), dtype=compound_dtype)
        for i in range(set[1].shape[1]):
            data_array = np.array(set[1].iloc[:, i], dtype=compound_dtype)
            tab[set[1].columns[i]] = data_array


def rename_columns(old_col_names):
    """
    Adds units to column names
    """
    new_col_names = []
    for col in old_col_names:
        if col.split(' ')[-1].startswith('Time'):
            new_col_names.append(col + ' [s]')
        elif col.split(' ')[-1].startswith('P'):
            new_col_names.append(col + ' [%]')
        elif col.split(' ')[-1].startswith('D'):
            new_col_names.append(col + ' [um^2s^-1]')
        elif col.split(' ')[-1].startswith('L'):
            new_col_names.append(col + ' [frames]')
        elif col.split(' ')[-1].startswith('confinement'):
            new_col_names.append(col + ' [um]')
        else:
            new_col_names.append(col)
    return new_col_names


def load_cs(sorted, coverslip, file, tiffiles, coverslip_data):
    """
    loads the data of a coverslip into the coverslip_data variable
    :param sorted: a correctly sorted list of cell names
    :param coverslip: a correctly sorted list of cell names in a coverslip
    :param file: list of h5 files
    :param tiffiles: list of the tif files for the timestamps
    :param coverslip_data: list of already loaded data to add to
    :return:coverslip_data with the data added to it
    """
    for cell in coverslip:
        hdf5global = h5py.File(parent_string(cell, file, "metadata"), "r")
        try:
            float(coverslip_data[sorted.index(coverslip)].iloc[0, 1])
        except IndexError:
            csStartTime = os.path.getmtime(parent_string(cell, tiffiles, "metadata"))
        statglobal = hdf5global["statistics"]["statistics_3"][0, 0]
        P_immobile = float(statglobal[0])
        P_confined = float(statglobal[1])
        P_free = float(statglobal[2])
        D_immobile = float(statglobal[5])
        if math.isnan(D_immobile):
            DG_immobile = 0
        else:
            DG_immobile = D_immobile
        D_confined = float(statglobal[6])
        if math.isnan(D_confined):
            DG_confined = 0
        else:
            DG_confined = D_confined
        D_free = float(statglobal[7])
        if math.isnan(D_free):
            DG_free = 0
        else:
            DG_free = D_free
        D_global = (DG_immobile * P_immobile + DG_confined * P_confined + DG_free * P_free) * 0.01
        L_immobile = float(statglobal[11])
        L_confined = float(statglobal[12])
        L_free = float(statglobal[13])
        if math.isnan(L_immobile):
            LG_immobile = 0
        else:
            LG_immobile = L_immobile
        if math.isnan(D_confined):
            LG_confined = 0
        else:
            LG_confined = L_confined
        if math.isnan(D_free):
            LG_free = 0
        else:
            LG_free = L_free
        L_global = (LG_immobile * P_immobile + LG_confined * P_confined + LG_free * P_free) * 0.01
        N_confined = 0
        N_free = 0
        N_immobile = 0
        confinementRadiiSum = 0
        DE_immobile = float(statglobal[8])
        DE_confined = float(statglobal[9])
        DE_free = float(statglobal[10])
        LE_immobile = float(statglobal[14])
        LE_confined = float(statglobal[15])
        LE_free = float(statglobal[16])
        for i, j in enumerate(hdf5global["rossier"]["rossierStatistics"][()]):
            if j[2] == 1:
                N_confined += 1
                confinementRadiiSum += j[7]
            elif j[3] == 1:
                N_free += 1
            elif j[1] == 1:
                N_immobile += 1
            elif j[4] == 1:
                N_immobile += 1
        if N_confined > 0:
            confinementRadii = confinementRadiiSum / N_confined
        else:
            confinementRadii = 0
        N_global = N_immobile + N_confined + N_free
        coverslip_data[sorted.index(coverslip)].loc[
            len(coverslip_data[sorted.index(coverslip)])] = cell, os.path.getmtime(parent_string(cell,
                                                                                                 tiffiles,
                                                                                                 'metadata')) - csStartTime, P_immobile, P_confined, P_free, D_global, D_immobile, D_confined, D_free, L_global, L_immobile, L_confined, L_free, N_global, N_immobile, N_confined, N_free, confinementRadii, DE_immobile, DE_confined, DE_free, LE_immobile, LE_confined, LE_free
    if coverslip_data[sorted.index(coverslip)].empty:
        raise IncorrectConfigException('Coverslide number ' + str(
            sorted.index(coverslip) + 1) + ' created an empty dataframe. Please check your data')
    return coverslip_data


def fourset_output(save_dir, outputfile, dataset, coverslip_data, use_timestamps, plotcolor, t_lig, ligand_name,
                   error_type, binned):
    """
    :param save_dir: where to save
    :param outputfile: h5 file to output into
    :param dataset: the set inside the h5 file to output into
    :param coverslip_data: the data to output
    :param use_timestamps: whether timestamps were used
    :param plotcolor: what color to plot
    :param t_lig: when the ligand was added
    :param ligand_name: the name of the ligand
    :param error_type: whether to use SD or SEM
    :param binned: whether you're outputting binned data
    outputs data that has 4 diffusion types ('global', 'immobile', 'confined', 'free')
    """
    # determins attribute name
    for folder in ['D', 'L', 'N']:
        if folder == 'D':
            name = 'diffusion_coefficients'
        if folder == 'L':
            name = 'segment_lengths'
        if folder == 'N':
            name = 'number_of_segments'
        data = []
        if binned:
            os.mkdir(save_dir + '\\TRplots\\' + name)
            os.mkdir(save_dir + '\\TRplots\\' + name + '\\pngs')
            os.mkdir(save_dir + '\\TRplots\\' + name + '\\svgs')
            os.mkdir(save_dir + '\\TRplots\\' + name + '\\pdfs')
        # writes output by signal type
        for signal_type in ['global', 'immobile', 'confined', 'free']:
            if binned:
                data.append(compile_columns(dataset, ["Cell Name", "Time", folder + "_" + signal_type,
                                                      folder + '_' + signal_type + '_SD',
                                                      folder + '_' + signal_type + '_SEM'], True))
                if use_timestamps:
                    plot = plot_by_time(coverslip_data, folder + '_' + signal_type, t_lig, ligand_name, dataset,
                                        error_type, plotcolor)
                else:
                    plot = plot_by_cells(coverslip_data, folder + '_' + signal_type, t_lig, ligand_name, dataset,
                                         error_type,
                                         plotcolor)
                plot.savefig(save_dir + '\\TRplots\\' + name + '\\pngs\\' + signal_type + '.png', dpi=300)
                plot.savefig(save_dir + '\\TRplots\\' + name + '\\svgs\\' + signal_type + '.svg', dpi=300)
                plot.savefig(save_dir + '\\TRplots\\' + name + '\\pdfs\\' + signal_type + '.pdf', dpi=300)

            else:
                data.append(compile_columns(dataset, ["Cell Name", "Time", folder + "_" + signal_type,
                                                      folder + 'E_' + signal_type], True))
        output_folder(outputfile, name,
                      [['global', data[0]], ['immobile', data[1]], ['confined', data[2]], ['free', data[3]]])


def stack_data(dataframes):
    """
    :param dataframes: list of dataframes
    :return: a dataframe consisting of the frames stacked atop each other
    """
    columns = ["Cell Name", "Time", "P_immobile", "P_confined", "P_free", "D_global", "D_immobile", "D_confined",
               "D_free", "L_global", "L_immobile", "L_confined", "L_free", "N_global", "N_immobile", "N_confined",
               "N_free", 'confinement_radius', 'DE_immobile', 'DE_confined', 'DE_free', 'LE_immobile', 'LE_confined',
               'LE_free']
    for df in dataframes:
        new_col = {}
        for i, col in enumerate(df.columns):
            new_col[col] = columns[i]
        df.rename(columns=new_col)
    return [pd.concat(dataframes, axis=0)]


def main(config_path):
    start_time = time.time()
    csPaths = []
    paths = []
    files = []
    use_timestamps = False
    config = configparser.ConfigParser()
    config.sections()
    config.read(config_path)

    try:
        if config["USE_TIMESTAMPS"]["use_timestamps"].lower() == "true":
            use_timestamps = True
    except KeyError:
        raise IncorrectConfigException("Section USE_TIMESTAMPS missing in config.")
    try:
        bin_size_cells = int(config["BIN_SIZE"]["cells"])
        bin_size_time = float(config["BIN_SIZE"]["minutes"]) * 60 + float(config["BIN_SIZE"]["seconds"])
    except KeyError:
        raise IncorrectConfigException("Section BIN_SIZE missing in config.")
    try:
        if len([key for key in config["CS_DIRS"]]):
            for key in config["CS_DIRS"]:
                csPaths.append(config["CS_DIRS"][key])
        else:
            raise IncorrectConfigException("No coverslip directory defined in config.")
    except KeyError:
        raise IncorrectConfigException("Section CS_DIRS missing in config.")

    tiffiles = []
    for dir in csPaths:
        tiffiles += get_matching_files(dir + "\\cells\\tifs", "cell", ["_dl", 'metadata'])

    csNames = []
    for cs in csPaths:
        csNames.append('_'.join(os.listdir(cs + "\\cells\\tifs")[0].split("\\")[-1].split('_')[:-2]))

    try:
        if len([key for key in config["GLOBAL_DIR"]]):
            for key in config["GLOBAL_DIR"]:
                paths.append(config["GLOBAL_DIR"][key])
        else:
            raise IncorrectConfigException("No GLOBAL directory defined in config.")
    except KeyError:
        raise IncorrectConfigException("Section GLOBAL_DIR missing in config.")

    try:
        save_dir = config["SAVE_DIR"]["svdir"]
    except KeyError:
        raise IncorrectConfigException("Parameter svdir missing in config.")

    try:
        alpha = float(config["STAT_SETTINGS"]["alpha_norm"])
    except KeyError:
        raise IncorrectConfigException("Parameter alpha_norm missing in config.")

    try:
        if config["STAT_SETTINGS"]["run_stats"].lower() == "true":
            run_stats = True
        else:
            run_stats = False
    except KeyError:
        raise IncorrectConfigException("Section USE_TIMESTAMPS missing in config.")

    try:
        p1 = float(config["STAT_SETTINGS"]["alpha_sign_1"])
    except KeyError:
        raise IncorrectConfigException("Parameter alpha_sign_1 missing in config.")
    try:
        p2 = float(config["STAT_SETTINGS"]["alpha_sign_2"])
    except KeyError:
        raise IncorrectConfigException("Parameter alpha_sign_1 missing in config.")
    try:
        p3 = float(config["STAT_SETTINGS"]["alpha_sign_3"])
    except KeyError:
        raise IncorrectConfigException("Parameter alpha_sign_1 missing in config.")

    try:
        plotcolor = config["PLOT_SETTINGS"]["dot_color"]
        if plotcolor == '':
            plotcolor = '#FFA500'
    except KeyError:
        raise IncorrectConfigException("Parameter dot_color missing in config.")

    try:
        t_lig = config["PLOT_SETTINGS"]["ligand_index"]
        if t_lig == '':
            ligand = False
        else:
            ligand = True
            # adds identificatior to the end of the ligand index
            if use_timestamps:
                t_lig += 's'
            else:
                t_lig += 'c'
    except KeyError:
        raise IncorrectConfigException("Parameter t_lig missing in config.")
    try:
        ligand_name = config["PLOT_SETTINGS"]["ligand_name"]
    except KeyError:
        raise IncorrectConfigException("Parameter ligand_name missing in config.")
    try:
        error_type = config["PLOT_SETTINGS"]["error_type"]
    except KeyError:
        raise IncorrectConfigException("Parameter error_type missing in config.")

    # writes the location of the h5 files into columns
    for path in paths:
        if type(path) is float:
            raise IncorrectConfigException("mismatched number of files")
        else:
            files = get_matching_files(path, ".h5", ["statistics.h5"])
    # sets up output folder
    try:
        os.mkdir(save_dir)
    except FileExistsError:
        pass
    try:
        os.mkdir(save_dir + '\\timeResolvedAnalysis')
    except FileExistsError:
        shutil.rmtree(save_dir + '\\timeResolvedAnalysis')
        os.mkdir(save_dir + '\\timeResolvedAnalysis')
    save_dir += '\\timeResolvedAnalysis'

    # writes the .tif files into a dictionary with the key being their coverslip name
    input_files = {c: [] for c in csNames}
    for h5 in files:
        filename = h5.split("\\")[-1][:-2]
        tif = parent_string(filename, tiffiles,
                            "metadata")  # goes through the tif files and looks for the one with the right name
        coverslipname = '_'.join(tif.split("\\")[-1].split('_')[:-2])
        input_files[coverslipname].append(filename)
    # sorts the files in correct order, taking multi digit numbers into account e.g. 10 comes after 2
    sorted = []
    for cs in input_files.keys():
        sorted.append(sort_cells(input_files[cs]))
    # prepares dataframes for each coverslip
    coverslip_data = [pd.DataFrame(
        columns=["Cell Name", "Time", "P_immobile", "P_confined", "P_free", "D_global", "D_immobile", "D_confined",
                 "D_free", "L_global", "L_immobile", "L_confined", "L_free", "N_global", "N_immobile", "N_confined",
                 "N_free", "confinement_radius", "DE_immobile", "DE_confined", "DE_free", "LE_immobile", "LE_confined",
                 "LE_free"]) for cs in csNames]

    # writes all data into dataframes and creates their own binned dataframe, then appends them to the binned_data list
    binned_data = []
    for coverslip in sorted:  # coverslip is a sorted list of all cells in a coverslip
        coverslip_data = load_cs(sorted, coverslip, files, tiffiles, coverslip_data)
        if use_timestamps:
            binned_mean = bin_data_time(coverslip_data[sorted.index(coverslip)], bin_size_time)
        else:
            binned_mean = bin_data_cells(coverslip_data[sorted.index(coverslip)], bin_size_cells)
        binned_data += [binned_mean]
    stacked_data = stack_data(coverslip_data)

    largest_bindex = 0
    for i, bin in enumerate(binned_data):
        if len(bin) > len(binned_data[largest_bindex]):
            largest_bindex = i

    global_mean = pd.DataFrame(binned_data[largest_bindex].iloc[:, 0:2])
    for attribute in ["P_immobile", "P_confined", "P_free", "D_global", "D_immobile", "D_confined",
                      "D_free", "L_global", "L_immobile", "L_confined", "L_free", "N_global", "N_immobile",
                      "N_confined",
                      "N_free", "confinement_radius"]:
        global_mean = pd.concat([global_mean, pd.DataFrame(calc_mean_over_cS(binned_data, attribute),
                                                           columns=[attribute, attribute + "_SD", attribute + "_SEM"])],
                                axis=1)
    for i, row in enumerate(global_mean["Cell Name"]):
        global_mean.iloc[i, 0] = row.split("_")[-1]
    global_mean.rename(columns={'Cell Names': 'cell range'})

    # output
    outputFile = h5py.File(save_dir + '\\stats.h5', 'w')
    outputFile_bin = outputFile.create_group('bin')
    outputFile_raw = outputFile.create_group('raw')
    outputFile_stacked = outputFile.create_group('stacked')
    # sheet with global means output
    output_folder(outputFile, 'global means',
                  [['global means', global_mean]])

    # raw data output
    output_folder(outputFile_raw, 'fractions',
                  [['immobile', compile_columns(coverslip_data, ["Cell Name", "Time", "P_immobile"], True)],
                   ['confined', compile_columns(coverslip_data, ["Cell Name", "Time", "P_confined"], True)],
                   ['free', compile_columns(coverslip_data, ["Cell Name", "Time", "P_free"], True)]])

    # binned data output
    # writes fractions and confinement_radii
    output_folder(outputFile_bin, 'fractions',
                  [['immobile', compile_columns(binned_data, ["Cell Name", "Time", "P_immobile"], True)],
                   ['confined', compile_columns(binned_data, ["Cell Name", "Time", "P_confined"], True)],
                   ['free', compile_columns(binned_data, ["Cell Name", "Time", "P_free"], True)]])
    # stacked output
    output_folder(outputFile_stacked, 'fractions',
                  [['immobile', compile_columns(stacked_data, ["Cell Name", "Time", "P_immobile"], False)],
                   ['confined', compile_columns(stacked_data, ["Cell Name", "Time", "P_confined"], False)],
                   ['free', compile_columns(stacked_data, ["Cell Name", "Time", "P_free"], False)]])

    name = 'fractions'
    os.mkdir(save_dir + '\\TRplots')
    os.mkdir(save_dir + '\\TRplots\\' + name)
    os.mkdir(save_dir + '\\TRplots\\' + name + '\\pngs')
    os.mkdir(save_dir + '\\TRplots\\' + name + '\\svgs')
    os.mkdir(save_dir + '\\TRplots\\' + name + '\\pdfs')
    for signal_type in ['immobile', 'confined', 'free']:
        if use_timestamps:
            plot = plot_by_time(coverslip_data, 'P_' + signal_type, t_lig, ligand_name, binned_data, error_type,
                                plotcolor)
        else:
            plot = plot_by_cells(coverslip_data, 'P_' + signal_type, t_lig, ligand_name, binned_data, error_type,
                                 plotcolor)
        plot.savefig(save_dir + '\\TRplots\\' + name + '\\pngs\\' + signal_type + '.png', dpi=300)
        plot.savefig(save_dir + '\\TRplots\\' + name + '\\svgs\\' + signal_type + '.svg', dpi=300)
        plot.savefig(save_dir + '\\TRplots\\' + name + '\\pdfs\\' + signal_type + '.pdf', dpi=300)

    output_folder(outputFile_raw, 'confinement_radii',
                  [['confinement_radii',
                    compile_columns(coverslip_data, ["Cell Name", "Time", 'confinement_radius'], True)]])
    output_folder(outputFile_bin, 'confinement_radii',
                  [['confinement_radii',
                    compile_columns(binned_data, ["Cell Name", "Time", 'confinement_radius'], True)]])
    output_folder(outputFile_stacked, 'confinement_radii',
                  [['confinement_radii',
                    compile_columns(stacked_data, ["Cell Name", "Time", 'confinement_radius'], False)]])

    name = 'confinement_radii'
    os.mkdir(save_dir + '\\TRplots\\' + name)
    if use_timestamps:
        plot = plot_by_time(coverslip_data, 'confinement_radius', t_lig, ligand_name, binned_data, error_type,
                            plotcolor)
    else:
        plot = plot_by_cells(coverslip_data, 'confinement_radius', t_lig, ligand_name, binned_data, error_type,
                             plotcolor)
    plot.savefig(save_dir + '\\TRplots\\' + name + '\\confinement_radii.png', dpi=300)
    plot.savefig(save_dir + '\\TRplots\\' + name + '\\confinement_radii.svg', dpi=300)
    plot.savefig(save_dir + '\\TRplots\\' + name + '\\confinement_radii.pdf', dpi=300)

    # writes the output files for diff coef, seg lengths, and number of segments
    # raw
    fourset_output(save_dir, outputFile_raw, coverslip_data, coverslip_data, use_timestamps, plotcolor, t_lig,
                   ligand_name,
                   error_type, False)
    # binned
    fourset_output(save_dir, outputFile_bin, binned_data, coverslip_data, use_timestamps, plotcolor, t_lig, ligand_name,
                   error_type, True)

    # stacked
    fourset_output(save_dir, outputFile_stacked, stacked_data, coverslip_data, use_timestamps, plotcolor, t_lig,
                   ligand_name,
                   error_type, False)

    metadata = pd.DataFrame(
        columns=['use_timestamps', 'n_cells', 'Plot_Error_type', 'Plot_Ligand_time', 'stat_alpha', 'stat_p1',
                 'stat_p2', 'stat_p3'])
    metadata.loc[len(metadata)] = pd.Series(
        [use_timestamps, str(int(bin_size_time / 60)) + 'm' + str(bin_size_time % 60) + 's', error_type, t_lig,
         alpha, p1, p2, p3]).values

    if use_timestamps:
        metadata = pd.DataFrame(
            columns=['use_timestamps', 'n_cells', 'Plot_Error_type', 'Plot_Ligand_time', 'stat_alpha', 'stat_p1',
                     'stat_p2', 'stat_p3'])
        metadata.loc[len(metadata)] = pd.Series(
            [use_timestamps, str(int(bin_size_time / 60)) + 'm' + str(bin_size_time % 60) + 's', error_type, t_lig,
             alpha, p1, p2, p3]).values
    else:
        metadata = pd.DataFrame(
            columns=['use_timestamps', 'n_cells', 'Plot_Error_type', 'Plot_Ligand_index', 'stat_alpha', 'stat_p1',
                     'stat_p2', 'stat_p3'])
        metadata.loc[len(metadata)] = pd.Series(
            [use_timestamps, str(bin_size_cells), error_type, t_lig,
             alpha, p1, p2, p3]).values

    output_folder(outputFile, 'metadata',
                  [['metadata', metadata]])
    try:
        output_folder(outputFile, 'metadata',
                      [['metadata', metadata]])

    except ValueError:
        pass
    outputFile.close()
    # runs the statistic tests
    if use_timestamps:
        run_stats = False
    if run_stats:
        os.mkdir(save_dir + '\\tests')
        os.mkdir(save_dir + '\\tests\\normality')
        if use_timestamps:
            if ligand:
                normFrames, signFrames = test_by_time(coverslip_data, bin_size_time, ligand, alpha, p1, p2, p3)
            else:
                normFrames = test_by_time(coverslip_data, bin_size_time, ligand, alpha, p1, p2, p3)
        else:
            if ligand:
                normFrames, signFrames = test_by_cell(coverslip_data, bin_size_cells, ligand, alpha, p1, p2, p3)
            else:
                normFrames = test_by_cell(coverslip_data, bin_size_cells, ligand, alpha, p1, p2, p3)
        for key in normFrames.keys():
            key2 = key
            if key.split('_')[0] == 'P':
                name = 'fractions'
            elif key.split('_')[0] == 'D':
                name = 'diffusion_coefficients'
            elif key.split('_')[0] == 'L':
                name = 'segment_lengths'
            elif key.split('_')[0] == 'N':
                name = 'number_of_segments'
            elif key.split('_')[0] == 'confinement':
                name = 'confinement_radius'
                key2 = ""
            else:
                continue
            try:
                os.mkdir(save_dir + '\\tests\\normality\\' + name)
            except FileExistsError:
                pass
            normFrames[key].to_csv(
                save_dir + '\\tests\\normality\\' + name + '\\test_normality_' + name + "_" + key2.split('_')[
                    -1] + '.csv', index=False)
        if ligand:
            os.mkdir(save_dir + '\\tests\\significance')
            for key in signFrames.keys():
                key2 = key
                if key.split('_')[0] == 'P':
                    name = 'fractions'
                elif key.split('_')[0] == 'D':
                    name = 'diffusion_coefficients'
                elif key.split('_')[0] == 'L':
                    name = 'segment_lengths'
                elif key.split('_')[0] == 'N':
                    name = 'number_of_segments'
                elif key.split('_')[0] == 'confinement':
                    name = 'confinement_radius'
                    key2 = ""
                else:
                    continue
                try:
                    os.mkdir(save_dir + '\\tests\\significance\\' + name)
                except FileExistsError:
                    pass
                if key == key2:
                    signFrames[key].to_csv((
                                                       save_dir + '\\tests\\significance\\' + name + '\\test_significance_' + name + "_" +
                                                       key2.split('_')[
                                                           -1] + '.csv'), index=False)
                else:
                    signFrames[key].to_csv((
                            save_dir + '\\tests\\significance\\' + name + '\\test_significance_' + name + '.csv'), index=False)

    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except FileExistsError:
        print("Usage: python timeResolvedAnalysis.py your_config_file.ini")
