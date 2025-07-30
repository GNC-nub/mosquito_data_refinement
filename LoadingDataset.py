# This loading_dictionary loads the data from matlab into a dictionary with every header like this:
'''
[f'Trial_{trial_number + a}_Track_{track_number + 1}'] = {
                            'x': x_tuple,
                            'y': y_tuple,
                            'z': z_tuple,
                            'time': time_tuple
                        }
'''
# Also loading it into a dataframe after that
# function loading_dataframe loads the data as a dataframe

path_matlab_file1 = '/Users/nubia/Desktop/Thesis_2.0/dataset/Database.mat' #replace this directory with your own!

import h5py
import numpy as np
import pandas as pd
import math



# this function deletes all the edge cases of Nan's
# so from [Nan, Nan, 2, 3, NaN, 4, 5, 6, Nan] --> [2, 3, Nan, 4, 5, 6]
def filtering_nan(lst):
    start_index = None
    end_index = None

    for i in range(len(lst)):
        if not math.isnan(lst[i]):
            end_index = i
    for i in range(len(lst) - 1, -1, -1):
        if not math.isnan(lst[i]):
            start_index = i

    if start_index == None or end_index == None:
        lst_filterd = []
    else:
        lst_filterd = lst[start_index:end_index + 1]
    return lst_filterd

def loading(path_matlab_file):
    with (h5py.File(path_matlab_file, 'r') as mat_file): # open math file
        trial = mat_file['Database']['Trial']['Tracks']
        print("--Start loading dataset--")
        dictionary = {}
        for trial_number in range(trial.shape[0]):
            if trial_number < 58: # to exclude trial 59 trial 60 needs to become trial 59
                a = 1
            else:
                a = 0
            if trial_number != 58: # Deleting trial 59!
                
                if ((trial_number+a) % 20) == 0:
                    print(f'Trial: {trial_number+a}')
                ref_trial = trial[trial_number, 0] # Idk wherefore the 0 is (found it with trial and error working)
                trial_data = mat_file[ref_trial]
                x_data_group = trial_data['x']
                y_data_group = trial_data['y']
                z_data_group = trial_data['z']
                time_data_group = trial_data['time']
                trial_key = f'Trial_{trial_number + a}'
                if trial_key not in dictionary:
                    dictionary[trial_key] = pd.DataFrame(columns=['x', 'y', 'z', 'time'])
                for track_number in range(x_data_group.shape[0]):
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
                    if len(filtered_x) != (len(filtered_y) or len(filtered_y) or len(filtered_time)):
                        print('ERROR') #Checks if the lengths of all the coodinates are the same.

                    if not filtered_x == []: # Assuming that if one of the filtered lists is empty, then all of them are
                        x_tuple = tuple(filtered_x)
                        y_tuple = tuple(filtered_y)
                        z_tuple = tuple(filtered_z)
                        time_tuple = tuple(filtered_time)
                        
                        dictionary[trial_key].at[f'Track_{track_number + 1}', 'x'] = x_tuple
                        dictionary[trial_key].at[f'Track_{track_number + 1}', 'y'] = y_tuple
                        dictionary[trial_key].at[f'Track_{track_number + 1}', 'z'] = z_tuple
                        dictionary[trial_key].at[f'Track_{track_number + 1}', 'time'] = time_tuple
        print("--Finished loading dataset--")
    return dictionary

def loading_dataframe():
    df = loading(path_matlab_file1)
    return df

