'''
READ ME
    accessing_track(trial, track)
    accessing_trial(trial)
    accessing_extra_info(info_field)

DESCRIOPTION
    accessing_track(trial, track) gives the data of a single track of a given trial.
    accessing_trial(trial) gives the data of a whole trial
    accessing_extra_info(info_field'). Directly loaded from the matlab file.

PARAMETERS
    info_field in accessing_extra_info(): Fil in one of the extra information field names options:
        'Condition', 'Day', 'TimeStart' or 'Catches'.


OUTPUT
    accessing_track(trial, track) gives an output like this:
        [x_coordinates], [y_coordinates], [z_coordinates], [time]

    accessing_trail(trial) gives an output like this:
        [ [header, [x_coordinates], [y_coordinates], [z_coordinates], [time], exe... ]

    accessing_extra_info(info_field') put outs a list of the info asked in a list form. Every trial has one of each 'extra_info's'

'''

import os
import csv
import sys
import h5py
import numpy as np

from loading_matlab_file import path_matlab_file1, path_csv_folder1
basemap_csv_path  = os.path.join(path_csv_folder1,'database_csv')

def accessing_track(trial, track):
    track_file = os.path.join(basemap_csv_path, f'Trial_{trial}/Trial_{trial}_Track_{track}.csv')
    if os.path.isfile(track_file):
        with open(track_file, newline='') as csvfile:
            dataset = list(csv.reader(csvfile))

            x_coordinates, y_coordinates, z_coordinates, time_coordinates = [], [], [], []

            for i in range(1,5):
                string = dataset[i][1]
                string = string.strip('()')
                string = string.split(',')
                for value in string:
                    if i == 1:
                        time_coordinates.append(float(value))
                    if i == 2:
                        x_coordinates.append(float(value))
                    if i == 3:
                        y_coordinates.append(float(value))
                    if i == 4:
                        z_coordinates.append(float(value))
            msg = [x_coordinates, y_coordinates, z_coordinates, time_coordinates]
    else:
        msg = 'Error: This file does not exist.'
    return msg

def accessing_trial(trial):
    csv.field_size_limit(sys.maxsize)
    trial_map = os.path.join(basemap_csv_path,f'Trial_{trial}')
    total_trial_data = []
    for file_name in os.listdir(trial_map):
        track_file = os.path.join(trial_map, file_name) #looping through every file in a trial map
        if os.path.isfile(track_file):
            with open(track_file, newline='') as csvfile:
                dataset = list(csv.reader(csvfile))
                header = dataset[0][1]
                x_coordinates, y_coordinates, z_coordinates, time_coordinates = [], [], [], []

                for i in range(1,5):
                    string = dataset[i][1]
                    string = string.strip('()')
                    string = string.split(',')
                    for value in string:
                        if not value == '': # to solve an error (ValueError: could not convert string to float: '')
                            if i == 1:
                                time_coordinates.append(float(value))
                            if i == 2:
                                x_coordinates.append(float(value))
                            if i == 3:
                                y_coordinates.append(float(value))
                            if i == 4:
                                z_coordinates.append(float(value))
                total_trial_data.append([header, x_coordinates, y_coordinates, z_coordinates, time_coordinates])
    return total_trial_data


def accessing_extra_info(info_field):
    if info_field == 'Condition':
        with (h5py.File(path_matlab_file1, 'r') as mat_file): # open math file
            database_info = mat_file['Database']['Trial']['Condition']
            info = []
            for trial_info in range(database_info.shape[0]):
                ref_trial_info = database_info[trial_info, 0]
                trial_info_data = mat_file[ref_trial_info]
                non_string_data = np.array(trial_info_data).flatten()
                characters = []
                for x in non_string_data:
                    characters.append(chr(x))
                    string_data = ''.join(characters)
                info.append(string_data)
    else:
        with (h5py.File(path_matlab_file1, 'r') as mat_file): # open math file
            database_info = mat_file['Database']['Trial']['Info']
            info = []
            for trial_info in range(database_info.shape[0]):
                ref_trial_info = database_info[trial_info, 0]
                trial_info_data = mat_file[ref_trial_info]
                info.append(list(np.array(trial_info_data[info_field]).flatten())) #info_field as a string is put in as the name
            tmp = []
            for i in info:
                tmp.append(i[0])
            info = tmp
    return info

