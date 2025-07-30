# Testing file

from LoadingDataset import *

# A function to test if all the datapoints are correctly loaded.
# Testing if the number of datapoints are the same in the matlab file and in the CSV file on my desktop

from loading_matlab_file import *
from nan_testing import *
def nr_datapoints_matlab(path_matlab_file):
    df_dataset = matlab_to_dataframe(path_matlab_file)
    total_tracks, total_points, nan_count, datapoint_count, nantrack_count = 0, 0, 0, 0, 0
    # Illiteriate through every track, through every coordinate (x,y,z):
    for row in range(df_dataset.shape[0]):
        for column in range(df_dataset.shape[1]):
            track = df_dataset.iloc[row, column]
            total_tracks += 1
            nan_bool = True
            for point in track:
                total_points += 1
                if pd.isna(point):
                    nan_count += 1
                else:
                    datapoint_count += 1
                    nan_bool = False
            if nan_bool == True:
                nantrack_count += 1
    return total_tracks, total_points, nan_count, datapoint_count, nantrack_count

def nr_datapoints_csv():
    print("--Loading CSV file--")
    total_tracks, total_points, nan_count, datapoint_count, nantrack_count = 0, 0, 0, 0, 0
    for trial_num in range(1, 65):
        if (trial_num % 10) == 0:
            print(f'Trial: {trial_num}')
        elif trial_num == 64:
            print('--Trial loading done--')
        trial = accessing_trial(trial_num)
        # Illiteriate through every track, through every coordinate (x,y,z):
        for track in trial:
            nan_bool = True
            total_tracks += 1
            for coordinates in track[1:]:
                for point in coordinates:
                    total_points += 1
                    if math.isnan(point) or np.isnan(point):
                        nan_count += 1
                    else:
                        datapoint_count += 1
                        nan_bool = False
            if nan_bool == True:
                nantrack_count += 1
    return total_tracks, total_points, nan_count, datapoint_count, nantrack_count


#Testing datapoints numbers of LoadingDataset.py

def nr_datapoints_dataframe():
    total_tracks, total_points, nan_count, datapoint_count, nantrack_count = 0, 0, 0, 0, 0
    df_dataset = loading_dataframe()
    for row in range(df_dataset.shape[0]):
        for column in range(df_dataset.shape[1]):
            track = df_dataset.iloc[row, column]
            total_tracks += 1
            nan_bool = True
            for point in track:
                total_points += 1
                if pd.isna(point):
                    nan_count += 1
                else:
                    datapoint_count += 1
                    nan_bool = False
            if nan_bool == True:
                nantrack_count += 1
    return total_tracks, total_points, nan_count, datapoint_count, nantrack_count


def nr_datapoints_test(path_matlab_file):
    tot_tracks_ml, tot_points_ml, nans_ml, datapoint_ml, nantracks_ml = nr_datapoints_matlab(path_matlab_file)
    tot_tracks_csv, tot_points_csv, nan_csv, datapoints_csv, nantrack_csv =nr_datapoints_csv()
    tot_tracks_df, tot_points_df, nan_df, datapoints_df, nantrack_df = nr_datapoints_dataframe()
    print(f'\nTOTAL TRACKS\nTotal tracks matlab = {tot_tracks_ml}\nTotal tracks csv ={tot_tracks_csv}\nTotal tracks df = {tot_tracks_df}')
    dif_tot_tracks = tot_tracks_ml - tot_tracks_csv
    dif_tot_tracks_df = tot_tracks_ml - tot_tracks_df
    print(f'Difference mat - CSV = {dif_tot_tracks}\nDifference mat - df = {dif_tot_tracks_df}')

    print(f'\nTOTAL POINTS\nTotal points matlab = {tot_points_ml}\nTotal points csv = {tot_points_csv}\nTotal points in df = {tot_points_df}')
    dif_tot_points = tot_points_ml - tot_points_csv
    dif_tot_points_df = tot_points_ml - tot_points_df
    print(f'Difference mat - CSV = {dif_tot_points}\nDifference mat - df = {dif_tot_points_df}')

    print(f'\nDATA VS NaN points\nDatapoints matlab = {datapoint_ml}\nNaNpoints matlab = {nans_ml}')
    print(f'Datapoints CSV = {datapoints_csv}\nNanpoints CSV = {nan_csv}\nDatapoints in df = {datapoints_df}\nNaNpoints in df ={nan_df}')
    dif_datapoints = datapoint_ml - datapoints_csv
    dif_nan = nans_ml - nan_csv
    dif_datapoints_df = datapoint_ml - datapoints_df
    dif_nan_df = nans_ml - nan_df
    print(f'Differance in datapoints (mat - CSV) = {dif_datapoints}\nDifferance in NaNs (mat - CSV) = {dif_nan}')
    print(f'Differance in datapoints (mat - df) = {dif_datapoints_df}\nDifferance in NaNs (mat - df) = {dif_nan_df}')

    print(f'\nNaNtrack count\nMatlab = {nantracks_ml}\nCSV = {nantrack_csv}')
    dif_nantrack = nantracks_ml - nantrack_csv
    print(f'Differance = {dif_nantrack}')
    dif_nantrack_df = nantracks_ml - nantrack_df
    print(f'Differance = {dif_nantrack_df}')

nr_datapoints_test(path_matlab_file1)


''' The result of this is that there are 2169 less total points in the df than in the matlab file. 
I assume most f this is due to filtering out nan's. 
Total tracks remains the same. 
There are more NaN's in the df and less actual datapoints :( IDK WHY 
'''

