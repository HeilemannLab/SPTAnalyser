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
    Search in directory and its subdirectories for files with target in their name but missing the exclusion string

    :param directory: all files in this directory & subdirectories are checked
    :param target: substring of filename
    :param exclusion_string: excludes files with this string in the name
    :return: List of matching file paths
    """
    matching_files = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if target.lower() in name.lower():
                if not any([True for string in exclusion_string if string.lower() in name.lower()]):
                    matching_files.append(os.path.join(path, name))
    return matching_files


def parentString(substring, parents, exclusionstr):
    """
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
    :param a list of cells
    :return: cells sorted in alphanumerical order taking multi digit numbers into account"""
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


def insertError(meanframe, sdframe, semframe):
    '''
    :param meanframe: a dataframe of calculated mean values
    :param sdframe: a dataframe of calculated standard deviations
    :param a dataframe of calculated standard errors of means
    :return: a dataframe with SD and SEM columns inserted after their mean
    '''
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


def calcMeanOverCS(binned_data, attribute):
    '''
    :param binned_data: a list of binned dataframes
    :param attribute: the attribute to calculate the values over
    :return: numpy stac of the the mean, standard deviation and standard error of mean of a certain attribute'''
    attributes = pd.DataFrame()
    for i, df in enumerate(binned_data):
        columns = [col for col in df.columns if attribute in col and not col.endswith(("_SD", '_SEM'))]
        for col in columns:
            attributes[f"{i}_{col}"] = df[col]
    stats = np.column_stack([attributes.mean(axis=1), attributes.std(axis=1), attributes.sem(axis=1)])
    return stats


def plotByTime(dataframes, attribute, t_lig, ligand_name, binned_data, error_type, clr):
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
    stats = calcMeanOverCS(binned_data, attribute)
    left = 0
    for i, row in enumerate(stats):
        # reads the error after it has been inserted using insertError
        if error_type == 'SEM':
            error = row[2]
        else:
            error = row[1]
        if i == len(stats) - 1:
            right = float(binned_data[0].iloc[i, 1].split('-')[1]) / 60
        else:
            right = (((float(binned_data[0].iloc[i + 1, 1].split('-')[0]) - float(
                binned_data[0].iloc[i, 1].split('-')[1])) / 2) + float(binned_data[0].iloc[i, 1].split('-')[1])) / 60
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
    if atr == 'confinementR':
        ax.set_ylabel('Confinement Radius' + r'$\mu m$')
    ax.set_xlim(left=-0.2)
    ax.set_xlabel('Time [s]')

    # add ligand bar
    if ligand_name == '':
        ligand_name = '+LIG'
    if len(t_lig) > 0 and t_lig[-1] == 's':
        ax.axvline(float(t_lig[:-1]) / 60, color='black')
        ax.text(float(t_lig[:-1]) / 60, ax.get_ylim()[1], ligand_name, ha='center', va='bottom', color='black')


    return plt


def plotByCells(dataframes, attribute, t_lig, ligand_name, binned_data, error_type, clr):
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
        if len(t_lig) > 0 and t_lig[-1] == 'c':
            ax.axvline((float(t_lig[:-1]) +0.5), color='black')

    # Plot mean bars and errors
    stats = calcMeanOverCS(binned_data, attribute)
    left = 0
    lens = [int(j.iloc[-1, 0].split('-')[1]) for j in binned_data]
    max_index = lens.index(max(lens))
    maxframe = binned_data[max_index]
    for i, row in enumerate(stats):
        if error_type == 'SEM':
            error = row[2]
        else:
            error = row[1]
        right = float(maxframe.iloc[i, 0].split('-')[1])
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
    if atr == 'confinementR':
        ax.set_ylabel('Confinement Radius' + r' [$\mu m$]')
    ax.set_xlim(left=-0.2)
    ax.set_xlabel('Cells')

    return plt


def binDataCells(frame, bin_size):
    '''
    :param frame: given dataframe
    :param bin_size: what interals to create the means over
    :return: frame with the means, SD and SEM for each time interval
    compiles data to ranges of means by their cell number, creates statistics over them.
    '''
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

        #adds mean, cell name interval, and time interval
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

    outframe = insertError(meanframe, sdframe, semframe)
    return outframe


def binDataTime(frame, bin_size):
    '''
    :param frame: given dataframe
    :param bin_size: what interals to create the means over
    :return: frame with the means, SD and SEM for each time interval
    compiles data to ranges of means by their time index, creates statistics over them.
    '''
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
        #adds time intervals
        if len(bin) > 0:
            meanframe.iloc[-1, 1] = str((binnumber - 1) * bin_size) + "-" + str((binnumber) * bin_size)
            sdframe.iloc[-1, 1] = str((binnumber - 1) * bin_size) + "-" + str((binnumber) * bin_size)
            semframe.iloc[-1, 1] = str((binnumber - 1) * bin_size) + "-" + str((binnumber) * bin_size)

        binnumber += 1
    outframe = insertError(meanframe, sdframe, semframe)
    return outframe


def testByCell(frames, bin_size, ligand, alpha, p1, p2, p3):
    """Runs statistic test over bins of a certain number of cells"""
    normality_frames = {}
    if ligand:
        significance_frames = {}

    bin_number = int(len(frames[0]) / bin_size)
    if len(frames[0]) % bin_size != 0:
        bin_number += 1
    for j in range(bin_number):
        tempFrame = pd.DataFrame()
        for frame in frames:
            tempFrame = tempFrame.append(
                frame.iloc[j * bin_size:(j + 1) * bin_size, :])  # all data of a bin is collected here
        if len(tempFrame) < 3:
            for name in tempFrame.columns[2:]:
                normality_frames[name].loc[len(normality_frames[name])] = pd.Series(
                    [j + 1, str(j * bin_size + 1) + '-' + str((j + 1) * bin_size), 'NaN', 'NaN', 'Dataset too small',
                     'NaN', 'NaN',
                     'Dataset too small']).values
                significance_frames[name].loc[len(significance_frames[name])] = pd.Series(
                    [j + 1, str(j * bin_size + 1) + '-' + str((j + 1) * bin_size), 'NaN', 'NaN', 'Dataset too small',
                     'NaN', 'NaN',
                     'Dataset too small']).values
            continue
        for name in tempFrame.columns[2:]:

            if name not in normality_frames.keys():  # creates new frame if necessary
                compare_frame = tempFrame
                norm_frame = pd.DataFrame(
                    columns=['BinNumber', 'Cellrange', 'Shapiro statistic', 'Shapiro p', 'Shapiro result',
                             'Kolmogorov-Smirnov statistic', 'Kolmogorov-Smirnov p',
                             'Kolmogorov-Smirnov result'])
                normality_frames[name] = norm_frame
            # runs the statistic tests
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
                [j + 1, str(j * bin_size + 1) + '-' + str((j + 1) * bin_size), shapstat, ps_value, sr, kolstat,
                 pk_value, kr]).values

            if ligand:
                if name not in significance_frames.keys():
                    sign_frame = pd.DataFrame(
                        columns=['Compared Bins', 'Cellrange', 'Mann-Whitney U statistic', 'Mann-Whitney U p', 'Mann-Whitney U result',
                                 'Wilcoxon-signed-rank statistic', 'Wilcoxon p',
                                 'Wilcoxon result'])
                    significance_frames[name] = sign_frame
                else:
                    mstat, pm_value, wilstat, pw_value = significance_tests(tempFrame.iloc[:, 2:],
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
                        ['1-' + str(j + 1), str(j * bin_size + 1) + '-' + str((j + 1) * bin_size), mstat, pm_value, tr,
                         wilstat, pw_value, wr]).values

    if ligand:
        return normality_frames, significance_frames
    else:
        return normality_frames


def testByTime(frames, bin_size, ligand, alpha, p1, p2, p3):
    """runs statistical tests over a bin of a certain timeframe"""
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
                        [i + 1, str(min(tempFrame['Time'])) + '-' + str(max(tempFrame['Time'])), 'NaN', 'NaN',
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
                        columns=['Compared Bins', 'Timerange', 'Mann-Whitney U statistic', 'Mann-Whitney U p', 'Mann-Whitney U result',
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
                        ['1-' + str(i + 1), str(min(tempFrame['Time'])) + '-' + str(max(tempFrame['Time'])), tstat,
                         pt_value, tr, wilstat, pw_value, wr]).values

    if ligand:
        return normality_frames, significance_frames
    else:
        return normality_frames


def compileColumns(frames, columns,renameColumns):
    """
    :param frames: list of dataframes
    :param columns: list of column names
    :param renameColumns:
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
        return value
    elif value[-1] == '.':
        return value[:-1]
    else:
        return value

def outputFolder(file, folder, datasets):
    """
    :param file: the h5 file to write to
    :param folder: name for the folder to write the data into
    :param datasets: the dataframe to write
    """
    fold = file.create_group(folder)
    for i, set in enumerate(
            datasets):  # set consists of a 2 entry list with the name of the set in the first and the dataframe in the second entry
        for col in set[1].columns:
            set[1][col] = set[1][col].apply(generate_shortname)
        compound_dtype = np.dtype([(a, h5py.special_dtype(vlen=str)) for a in set[1].columns])
        tab = fold.create_dataset(set[0], (len(set[1]),), dtype=compound_dtype)
        for i in range(set[1].shape[1]):
            data_array = np.array(set[1].iloc[:, i], dtype=compound_dtype)
            tab[set[1].columns[i]] = data_array

def loadCS(sorted, coverslip, file, tiffiles, coverslipData):
    """
    :param sorted: a correctly sorted list of cell names
    :param coverslip: a correctly sorted list of cell names in a coverslip
    :param file: list of h5 files
    :param tiffiles: list of the tif files for the timestamps
    :param coverslipData: list of already loaded data to add to
    :return:
    """
    for cell in coverslip:
        hdf5global = h5py.File(parentString(cell, file, "metadata"), "r")
        try:
            float(coverslipData[sorted.index(coverslip)].iloc[0, 1])
        except IndexError:
            csStartTime = os.path.getmtime(parentString(cell, tiffiles, "metadata"))
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
        L_global = (L_immobile * P_immobile + L_confined * P_confined + L_free * P_free) * 0.01
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
        coverslipData[sorted.index(coverslip)].loc[
            len(coverslipData[sorted.index(coverslip)])] = cell, os.path.getmtime(parentString(cell,
                                                                                               tiffiles,
                                                                                               'metadata')) - csStartTime, P_immobile, P_confined, P_free, D_global, D_immobile, D_confined, D_free, L_global, L_immobile, L_confined, L_free, N_global, N_immobile, N_confined, N_free, confinementRadii, DE_immobile, DE_confined, DE_free, LE_immobile, LE_confined, LE_free
    if coverslipData[sorted.index(coverslip)].empty:
        raise IncorrectConfigException('Coverslide number ' + str(
            sorted.index(coverslip) + 1) + ' created an empty dataframe. Please check your data')
    return coverslipData


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
    try:
        mstat, pm_value = scy.mannwhitneyu(frame1[attribute], frame2[attribute])
    except ValueError:
        mstat, pm_value = 'not applicable to data', 'NaN'
    try:
        wilstat, pw_value = scy.wilcoxon(frame1[attribute], frame2[attribute])
    except ValueError:
        wilstat, pw_value = 'not applicable to data', 'NaN'

    return mstat, pm_value, wilstat, pw_value


def foursetOutput(save_dir, outputfile, dataset, coverslipData, use_timestamps, plotcolor, t_lig, ligand_name,
                  error_type, binned):
    """
    :param save_dir: where to save
    :param outputfile: h5 file to output into
    :param dataset: the set inside the h5 file to output into
    :param coverslipData: the data to output
    :param use_timestamps: whether timestamps were used
    :param plotcolor: what color to plot
    :param t_lig: when the ligand was added
    :param ligand_name: the name of the ligand
    :param error_type: whether to use SD or SEM
    :param binned: whether you're outputting binned data
    outputs data that has 4 diffusion types ('global', 'immobile', 'confined', 'free')
    """
    #determins attribute name
    for folder in ['D', 'L', 'N']:
        if folder == 'D':
            name = 'Diffusion Coefficients'
        if folder == 'L':
            name = 'Segment lengths'
        if folder == 'N':
            name = 'number of segments'
        data = []
        if binned:
            os.mkdir(save_dir + '\\TRplots\\' + name)
            os.mkdir(save_dir + '\\TRplots\\' + name + '\\pngs')
            os.mkdir(save_dir + '\\TRplots\\' + name + '\\svgs')
            os.mkdir(save_dir + '\\TRplots\\' + name + '\\pdfs')
        #writes output by signal type
        for signal_type in ['global', 'immobile', 'confined', 'free']:
            if binned:
                data.append(compileColumns(dataset, ["Cell Name", "Time", folder + "_" + signal_type,
                                                     folder + '_' + signal_type + '_SD',
                                                     folder + '_' + signal_type + '_SEM'],True))
                if use_timestamps:
                    plot = plotByTime(coverslipData, folder + '_' + signal_type, t_lig, ligand_name, dataset,
                                      error_type, plotcolor)
                else:
                    plot = plotByCells(coverslipData, folder + '_' + signal_type, t_lig, ligand_name, dataset,
                                       error_type,
                                       plotcolor)
                plot.savefig(save_dir + '\\TRplots\\' + name + '\\pngs\\' + signal_type + '.png', dpi=300)
                plot.savefig(save_dir + '\\TRplots\\' + name + '\\svgs\\' + signal_type + '.svg', dpi=300)
                plot.savefig(save_dir + '\\TRplots\\' + name + '\\pdfs\\' + signal_type + '.pdf', dpi=300)

            else:
                data.append(compileColumns(dataset, ["Cell Name", "Time", folder + "_" + signal_type,
                                                     folder + 'E_' + signal_type],True))
        outputFolder(outputfile, name,
                     [['global', data[0]], ['immobile', data[1]], ['confined', data[2]], ['free', data[3]]])


def stackData(dataframes):
    """
    :param dataframes: list of dataframes
    :return: a dataframe consisting of the frames stacked atop each other
    """
    columns = ["Cell Name", "Time", "P_immobile", "P_confined", "P_free", "D_global", "D_immobile", "D_confined",
               "D_free", "L_global", "L_immobile", "L_confined", "L_free", "N_global", "N_immobile", "N_confined",
               "N_free", "confinementR", 'DE_immobile', 'DE_confined', 'DE_free', 'LE_immobile', 'LE_confined',
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
        alpha = float(config["STAT_SETTINGS"]["alpha"])
    except KeyError:
        raise IncorrectConfigException("Parameter alpha missing in config.")

    try:
        if config["STAT_SETTINGS"]["run_stats"].lower() == "true":
            run_stats = True
        else:
            run_stats = False
    except KeyError:
        raise IncorrectConfigException("Section USE_TIMESTAMPS missing in config.")

    try:
        p1 = float(config["STAT_SETTINGS"]["p1"])
    except KeyError:
        raise IncorrectConfigException("Parameter p1 missing in config.")
    try:
        p2 = float(config["STAT_SETTINGS"]["p2"])
    except KeyError:
        raise IncorrectConfigException("Parameter p1 missing in config.")
    try:
        p3 = float(config["STAT_SETTINGS"]["p3"])
    except KeyError:
        raise IncorrectConfigException("Parameter p1 missing in config.")

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
            #adds identificatior to the end of the ligand index
            if use_timestamps:
                t_lig += 's'
            else:
                t_lig +='c'
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
        os.mkdir(save_dir + '\\TR_Analysis')
    except FileExistsError:
        shutil.rmtree(save_dir + '\\TR_Analysis')
        os.mkdir(save_dir + '\\TR_Analysis')
    save_dir += '\\TR_Analysis'

    # writes the .tif files into a dictionary with the key being their coverslip name
    input_files = {c: [] for c in csNames}
    for h5 in files:
        filename = h5.split("\\")[-1][:-2]
        tif = parentString(filename, tiffiles,
                           "metadata")  # goes through the tif files and looks for the one with the right name
        coverslipname = '_'.join(tif.split("\\")[-1].split('_')[:-2])
        input_files[coverslipname].append(filename)
    # sorts the files in correct order, taking multi digit numbers into account e.g. 10 comes after 2
    sorted = []
    for cs in input_files.keys():
        sorted.append(sort_cells(input_files[cs]))
    # prepares dataframes for each coverslip
    coverslipData = [pd.DataFrame(
        columns=["Cell Name", "Time", "P_immobile", "P_confined", "P_free", "D_global", "D_immobile", "D_confined",
                 "D_free", "L_global", "L_immobile", "L_confined", "L_free", "N_global", "N_immobile", "N_confined",
                 "N_free", "confinementR", 'DE_immobile', 'DE_confined', 'DE_free', 'LE_immobile', 'LE_confined',
                 'LE_free']) for cs in csNames]

    # writes all data into dataframes and creates their own binned dataframe, then appends them to the binnedData list
    binnedData = []
    for coverslip in sorted:  # coverslip is a sorted list of all cells in a coverslip
        coverslipData = loadCS(sorted, coverslip, files, tiffiles, coverslipData)
        if use_timestamps:
            binnedMean = binDataTime(coverslipData[sorted.index(coverslip)], bin_size_time)
        else:
            binnedMean = binDataCells(coverslipData[sorted.index(coverslip)], bin_size_cells)
        binnedData += [binnedMean]
    stacked_data = stackData(coverslipData)
    # output
    outputFile = h5py.File(save_dir + '\\stats.h5', 'w')
    outputFile_bin = outputFile.create_group('bin')
    outputFile_raw = outputFile.create_group('raw')
    outputFile_stacked = outputFile.create_group('stacked')
    # raw data output
    outputFolder(outputFile_raw, 'Fractions',
                 [['immobile', compileColumns(coverslipData, ["Cell Name", "Time", "P_immobile"],True)],
                  ['confined', compileColumns(coverslipData, ["Cell Name", "Time", "P_confined"],True)],
                  ['free', compileColumns(coverslipData, ["Cell Name", "Time", "P_free"],True)]])

    # binned data output
    # writes fractions and confinement radii
    outputFolder(outputFile_bin, 'Fractions',
                 [['immobile', compileColumns(binnedData, ["Cell Name", "Time", "P_immobile"],True)],
                  ['confined', compileColumns(binnedData, ["Cell Name", "Time", "P_confined"],True)],
                  ['free', compileColumns(binnedData, ["Cell Name", "Time", "P_free"],True)]])
    # stacked output
    outputFolder(outputFile_stacked, 'Fractions',
                 [['immobile', compileColumns(stacked_data, ["Cell Name", "Time", "P_immobile"],False)],
                  ['confined', compileColumns(stacked_data, ["Cell Name", "Time", "P_confined"],False)],
                  ['free', compileColumns(stacked_data, ["Cell Name", "Time", "P_free"],False)]])
    name = 'Fraction'
    os.mkdir(save_dir + '\\TRplots')
    os.mkdir(save_dir + '\\TRplots\\' + name)
    os.mkdir(save_dir + '\\TRplots\\' + name + '\\pngs')
    os.mkdir(save_dir + '\\TRplots\\' + name + '\\svgs')
    os.mkdir(save_dir + '\\TRplots\\' + name + '\\pdfs')
    for signal_type in ['immobile', 'confined', 'free']:
        if use_timestamps:
            plot = plotByTime(coverslipData, 'P_' + signal_type, t_lig, ligand_name, binnedData, error_type, plotcolor)
        else:
            plot = plotByCells(coverslipData, 'P_' + signal_type, t_lig, ligand_name, binnedData, error_type, plotcolor)
        plot.savefig(save_dir + '\\TRplots\\' + name + '\\pngs\\' + signal_type + '.png', dpi=300)
        plot.savefig(save_dir + '\\TRplots\\' + name + '\\svgs\\' + signal_type + '.svg', dpi=300)
        plot.savefig(save_dir + '\\TRplots\\' + name + '\\pdfs\\' + signal_type + '.pdf', dpi=300)

    outputFolder(outputFile_raw, 'Confinement Radii',
                 [['Confinement Radii', compileColumns(coverslipData, ["Cell Name", "Time", "confinementR"],True)]])
    outputFolder(outputFile_bin, 'Confinement Radii',
                 [['Confinement Radii', compileColumns(binnedData, ["Cell Name", "Time", "confinementR"],True)]])
    outputFolder(outputFile_stacked, 'Confinement Radii',
                 [['Confinement Radii', compileColumns(stacked_data, ["Cell Name", "Time", "confinementR"],False)]])

    name = 'Confinement Radii'
    os.mkdir(save_dir + '\\TRplots\\' + name)
    if use_timestamps:
        plot = plotByTime(coverslipData, 'confinementR', t_lig, ligand_name, binnedData, error_type, plotcolor)
    else:
        plot = plotByCells(coverslipData, 'confinementR', t_lig, ligand_name, binnedData, error_type, plotcolor)
    plot.savefig(save_dir + '\\TRplots\\' + name + '\\Confinement Radii.png', dpi=300)
    plot.savefig(save_dir + '\\TRplots\\' + name + '\\Confinement Radii.svg', dpi=300)
    plot.savefig(save_dir + '\\TRplots\\' + name + '\\Confinement Radii.pdf', dpi=300)

    # writes the output files for diff coef, seg lengths, and number of segments
    # raw
    foursetOutput(save_dir, outputFile_raw, coverslipData, coverslipData, use_timestamps, plotcolor, t_lig, ligand_name,
                  error_type, False)
    # binned
    foursetOutput(save_dir, outputFile_bin, binnedData, coverslipData, use_timestamps, plotcolor, t_lig, ligand_name,
                  error_type, True)

    # stacked
    foursetOutput(save_dir, outputFile_stacked, stacked_data, coverslipData, use_timestamps, plotcolor, t_lig,
                  ligand_name,
                  error_type, False)

    metadata = pd.DataFrame(
        columns=['use_timestamps', 'Bin_size', 'Plot_Error_type', 'Plot_Ligand_time', 'stat_alpha', 'stat_p1',
                 'stat_p2', 'stat_p3'])
    metadata.loc[len(metadata)] = pd.Series(
        [use_timestamps, str(int(bin_size_time / 60)) + 'm' + str(bin_size_time % 60) + 's', error_type, t_lig,
         alpha, p1, p2, p3]).values
    outputFolder(outputFile, 'metadata',
                 [['metadata', metadata]])
    if use_timestamps:
        metadata = pd.DataFrame(
            columns=['use_timestamps', 'Bin_size', 'Plot_Error_type', 'Plot_Ligand_time', 'stat_alpha', 'stat_p1',
                     'stat_p2', 'stat_p3'])
    else:
        metadata = pd.DataFrame(
            columns=['use_timestamps', 'Bin_size', 'Plot_Error_type', 'Plot_Ligand_index', 'stat_alpha', 'stat_p1',
                     'stat_p2', 'stat_p3'])
    metadata.loc[len(metadata)] = pd.Series(
        [use_timestamps, str(int(bin_size_time / 60)) + 'm' + str(bin_size_time % 60) + 's', error_type, t_lig,
         alpha, p1, p2, p3]).values
    try:
        outputFolder(outputFile, 'metadata',
                     [['metadata', metadata]])
    except ValueError:
        pass
    outputFile.close()
    # runs the statistic tests
    if run_stats:
        os.mkdir(save_dir + '\\tests')
        os.mkdir(save_dir + '\\tests\\normality')
        if use_timestamps:
            if ligand:
                normFrames, signFrames = testByTime(coverslipData, bin_size_time, ligand, alpha, p1, p2, p3)
            else:
                normFrames = testByTime(coverslipData, bin_size_time, ligand, alpha, p1, p2, p3)
        else:
            if ligand:
                normFrames, signFrames = testByCell(coverslipData, bin_size_cells, ligand, alpha, p1, p2, p3)
            else:
                normFrames = testByCell(coverslipData, bin_size_cells, ligand, alpha, p1, p2, p3)
        for key in normFrames.keys():
            if key.split('_')[0] == 'P':
                name = 'Fractions'
            if key.split('_')[0] == 'D':
                name = 'Diffusion Coefficients'
            if key.split('_')[0] == 'L':
                name = 'Segment lengths'
            if key.split('_')[0] == 'N':
                name = 'number of segments'
            if key.split('_')[0] == 'confinementR':
                name = 'confinement Radius'
            try:
                os.mkdir(save_dir + '\\tests\\normality\\' + name)
            except FileExistsError:
                pass
            normFrames[key].to_csv(
                save_dir + '\\tests\\normality\\' + name + '\\' + ''.join(key.split('_')[:2]) + '.csv', index=False)
        if ligand:
            os.mkdir(save_dir + '\\tests\\significance')
            for key in signFrames.keys():
                if key.split('_')[0] == 'P':
                    name = 'Fractions'
                if key.split('_')[0] == 'D':
                    name = 'Diffusion Coefficients'
                if key.split('_')[0] == 'L':
                    name = 'Segment lengths'
                if key.split('_')[0] == 'N':
                    name = 'number of segments'
                if key.split('_')[0] == 'confinementR':
                    name = 'confinement Radius'
                try:
                    os.mkdir(save_dir + '\\tests\\significance\\' + name)
                except FileExistsError:
                    pass
                signFrames[key].to_csv(
                    save_dir + '\\tests\\significance\\' + name + '\\' + ''.join(key.split('_')[:2]) + '.csv',
                    index=False)

    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    try:
        cfg_path = sys.argv[1]
        main(cfg_path)
    except FileExistsError:
        print("Usage: python TR_Analysis.py your_config_file.ini")
