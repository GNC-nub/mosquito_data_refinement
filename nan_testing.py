'''
READ ME
    First load load_matlab_file or put it the right path to the matlab file!!
DESCRIPTION
    A script to find out what percentage the dataset is made up of nan's.
    The dataset is made up of 11,78 % of nan's.
    26,18 % of the tracks only contain nan's
OUTPUT
    The percentage of nan's in the whole dataset.
    The amount of tracks with only nan's
    The number of real datapoints and
    The number of total datapoints including nan's
    The percentage of tracks that only contain nan's


    Output, so I don't have to run it again:
        The nan count is 1266603.0.
        The datapoints count containing actual data is 9480653.0.
        The total is 42989024.
        The percentage of nan on the whole dataset is 2.9463404426208886 %.
        The number of tracks containing only nan's is 9873.0.
        The total number of tracks is 37705.0.
        The percentage of tracks containing only nan's is 26.184856119878003 %
'''


import h5py
import numpy as np
import pandas as pd

from loading_matlab_file import path_matlab_file1

def matlab_to_dataframe(filename):
    dict = {}
    with (h5py.File(filename, 'r') as mat_file):
        trial = mat_file['Database']['Trial']['Tracks']

        for trial_number in range(trial.shape[0]):
            print(trial_number)
            ref_trial = trial[trial_number, 0] # Idk wherefore the 0 is (found it with trial and error working)
            trial_data = mat_file[ref_trial]

            x_data_group = trial_data['x']
            y_data_group = trial_data['y']
            z_data_group = trial_data['z']
            time_data_group = trial_data['z']

            if trial_number == 58:
                x_vals, y_vals, z_vals, time_vals = [], [], [], []

                x_vals.append(np.array(x_data_group).flatten())
                y_vals.append(np.array(y_data_group).flatten())
                z_vals.append(np.array(z_data_group).flatten())
                time_vals.append(np.array(time_data_group).flatten())

                x_tuple = tuple(x_vals[0])
                y_tuple = tuple(y_vals[0])
                z_tuple = tuple(z_vals[0])
                time_tuple = tuple(time_vals[0])

                dict[f'Trial_{trial_number + 1}_Track_{track_number + 1}'] = {
                    'x': x_tuple,
                    'y': y_tuple,
                    'z': z_tuple,
                    'time': time_tuple
                }
                df_show = pd.DataFrame(dict[f'Trial_{trial_number + 1}_Track_{track_number + 1}'])
                #print(df_show.head())
            else:
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

                    x_tuple = tuple(x_vals[0])
                    y_tuple = tuple(y_vals[0])
                    z_tuple = tuple(z_vals[0])
                    time_tuple = tuple(time_vals[0])

                    dict[f'Trial_{trial_number+1}_Track_{track_number + 1}'] = {
                        'x': x_tuple,
                        'y': y_tuple,
                        'z': z_tuple,
                        'time': time_tuple
                    }

        df = pd.DataFrame(dict)
        print(df.head())
    return df

nan_count = 0
datapoint_count = 0
total_points = 0
nantrack_count = 0
total_tracks = 0

df_dataset = matlab_to_dataframe(path_matlab_file1)

#illiteriate through every track, through every coordinate (x,y,z):
for row in range(df_dataset.shape[0]):
    for column in range(df_dataset.shape[1]):
        track = df_dataset.iloc[row, column]
        total_tracks += 1
        for point in track:
            total_points += 1
            x = True
            if pd.isna(point):
                nan_count += 1
            else:
                datapoint_count += 1
                x = False
        if x == True:
            nantrack_count += 1


#Every x,y,z, time tuples is still seperate tacks, so to acccount for this i devide by 3.
nantrack_count /= 4
total_tracks /= 4
nan_count /= 4
datapoint_count /= 4

print(f'The nan count is {nan_count}.\nThe datapoints count containing actual data is {datapoint_count}. \nThe total is {total_points}.')
print(f'The percentage of nan on the whole dataset is {(nan_count/total_points)*100} %.')
print(f"The number of tracks containing only nan's is {nantrack_count}.\nThe total number of tracks is {total_tracks}.")
print(f"The percentage of tracks containing only nan's is {(nantrack_count/total_tracks)*100} %")
