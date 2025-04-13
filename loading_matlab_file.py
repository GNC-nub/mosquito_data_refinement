''' RUN THIS FILE FIRST + ADD OWN PATHS!
READ ME
    matlab_to_csv_map(filename)
    filtering_nan(lst)
DESCRIPTION
    matlab_to_csv_map(path_matlab_file, path_csv_folder)
        This is a function that takes the whole S3 matlab dataset and turns it into a csv files folder on the given path.
        The tracks are turned into a dictionary, then put into a pandas dataframe.
        It uses the references of the nested matlab dataset to collect the actual data.
        It creates a map with maps of every trial, in every trial map there are track csv files for easy access.
        Every track from every trail is put into a column with a header named like this: Trail_1_Track_2

        Trial 59 was different to the rest of the tracks, because it had only one track.
        Within the track 59 there were no groups only the actual data, so there is a different piece of code in the
        function for just trail 59.
    filtering_nan(lst)
        This function is used to filter out all the nan from the end and the beginning of each track.
        All the nan's in the middle of actual data are kept within the csv data.

PARAMETERS
    path_matlab_file is the path to the original matlab file (string).
        I used: '/Users/nubia/Desktop/Thesis/Datasets/S3_-_Database.mat', but put in your own path.
    path_csv_folder is the path to where the new datastructure should be stored (string)
        For example in my case I stored it on my desktop: '/Users/nubia/Desktop'

LIMITATIONS
    This only works for this specific dataset as it uses the names of the headers in the matlab dataset (S3_-_Database.mat)

STRUCTURE
    with, for-loops, if-else statement

OUTPUT
    When running this file csv files are made for every track and stored in this directory:
    '.../database_csv/Trial_{trial_number}/Trial_{trial_number}_Track_{track_number}.csv'

'''

import h5py
import numpy as np
import pandas as pd
import os





#Fill in your own paths here and run this file first! Do this before running main.py!

#path to the matlab file:
path_matlab_file1 = '/Users/nubia/Desktop/Thesis/Datasets/Database_time.mat'



#path to where you want to store the new data structure (for example a desktop):
path_csv_folder1 = '/Users/nubia/Desktop/final_thesis_testing'





def filtering_nan(lst):
    start_index = None
    end_index = None

    for i in range(len(lst)):
        if not np.isnan(lst[i]):
            end_index = i
    for i in range(len(lst) - 1, -1, -1):
        if not np.isnan(lst[i]):
            start_index = i

    if start_index == None or end_index == None:
        lst_filterd = []
    else:
        lst_filterd  = lst[start_index:end_index + 1]
    return lst_filterd

def matlab_to_csv_map(path_matlab_file, path_csv_folder):
    with (h5py.File(path_matlab_file, 'r') as mat_file): # open math file
        trial = mat_file['Database']['Trial']['Tracks']

        # making a map on the directory_csv_folder
        basemap_path  = os.path.join(path_csv_folder,'database_csv')
        os.makedirs(basemap_path, exist_ok=True)

        for trial_number in range(trial.shape[0]):
            print(trial_number+1)
            new_map = os.path.join(basemap_path, f'Trial_{trial_number+1}')
            os.makedirs(new_map, exist_ok=True) # make a map per trial

            ref_trial = trial[trial_number, 0] # Idk wherefore the 0 is (found it with trial and error working)
            trial_data = mat_file[ref_trial]
            x_data_group = trial_data['x']
            y_data_group = trial_data['y']
            z_data_group = trial_data['z']
            time_data_group = trial_data['time']

            if trial_number == 58:
                dictionary = {}
                x_vals, y_vals, z_vals, time_vals = [], [], [], []

                x_vals.append(np.array(x_data_group).flatten())
                y_vals.append(np.array(y_data_group).flatten())
                z_vals.append(np.array(z_data_group).flatten())
                time_vals.append(np.array(time_data_group).flatten())

                # from here i'm trying to filter out all the nan's
                x_list = list(x_vals[0])
                y_list = list(y_vals[0])
                z_list = list(z_vals[0])
                time_list = list(time_vals[0])

                filtered_x = filtering_nan(x_list)
                filtered_y = filtering_nan(y_list)
                filtered_z = filtering_nan(z_list)
                filtered_time = filtering_nan(time_list)

                if not filtered_x == []:  # again assuming that if one of the filterd lists is empty than all of them are
                    x_tuple = tuple(filtered_x)
                    y_tuple = tuple(filtered_y)
                    z_tuple = tuple(filtered_z)
                    time_tuple = tuple(filtered_time)

                    dictionary[f'Trial_{59}_Track_{1}'] = {
                        'x': x_tuple,
                        'y': y_tuple,
                        'z': z_tuple,
                        'time': time_tuple
                    }
                    df = pd.DataFrame(dictionary)
                    file_path = os.path.join(new_map, f'Trial_{59}_Track_{1}.csv')
                    df.to_csv(file_path)


            else:
                for track_number in range(x_data_group.shape[0]):
                    dictionary = {}
                    ref_x_data = x_data_group[track_number, 0] # Idk what for the 0 is
                    ref_y_data = y_data_group[track_number, 0]
                    ref_z_data = z_data_group[track_number, 0]
                    ref_time_data = time_data_group[track_number, 0]

                    x_vals, y_vals, z_vals, time_vals = [], [], [], []

                    x_data = mat_file[ref_x_data]
                    y_data = mat_file[ref_y_data]
                    z_data = mat_file[ref_z_data]
                    time_data = mat_file[ref_time_data]

                    x_vals.append(np.array(x_data).flatten())
                    y_vals.append(np.array(y_data).flatten())
                    z_vals.append(np.array(z_data).flatten())
                    time_vals.append(np.array(time_data).flatten())

                    # from here i'm trying to filter out all the nan's
                    x_list = list(x_vals[0])
                    y_list = list(y_vals[0])
                    z_list = list(z_vals[0])
                    time_list = list(time_vals[0])

                    filtered_x = filtering_nan(x_list)
                    filtered_y = filtering_nan(y_list)
                    filtered_z = filtering_nan(z_list)
                    filtered_time = filtering_nan(time_list)

                    if not filtered_x == []: #again assuming that if one of the filterd lists is empty than all of them are
                        x_tuple = tuple(filtered_x)
                        y_tuple = tuple(filtered_y)
                        z_tuple = tuple(filtered_z)
                        time_tuple = tuple(filtered_time)

                        dictionary[f'Trial_{trial_number + 1}_Track_{track_number + 1}'] = {
                            'x': x_tuple,
                            'y': y_tuple,
                            'z': z_tuple,
                            'time': time_tuple

                        }
                        df = pd.DataFrame(dictionary)
                        file_path = os.path.join(new_map, f'Trial_{trial_number + 1}_Track_{track_number + 1}.csv')
                        df.to_csv(file_path)


if __name__ =='__main__':
    matlab_to_csv_map(path_matlab_file1, path_csv_folder1)

