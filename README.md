## 1. Loading the correct packages 
Make sure the following packages are pip installed in your working environment / interpreter: 
- numpy
- h5py
- sys
- os
- csv
- matplotlib 
- scripy 
- statsmodels 
- pandas


## 2. Loading the dataset
The first step is to load the dataset from the matlab file to a csv folder data structure for easy, direct and/or partial access.

1. Download the correct matlab file from the source and coppy your directory. 
2. Open the 'loading_matlab_file.py' file.
3. Add your own directory of the matlab file on the correct line in the 'loading_matlab_file.py' file.


    DESCRIPTION 'loading_matlab_file.py' 
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

## 3. Running main.py and nan_testing.py 
Make sure you completed step one, before running any other files! 

For generating the information about the nan data in the original matlab file code run: nan_testing.py 

For generating all the figures run main.py. There are short discriptions above every set of figure's. 
You can togle on and off the figures by adding a # in fromt of the line. 

Main.py calls mostly on the ClassMosquito.py file where three classes are written. 

    DESCRIPTION 'main.py' 
    This script puts out all the data needed for the BSc Thesis of Nubia Middelkoop.
    This code uses 3 classes from the file ClassMosquito.py:
        Track (to load single tracks), Trials (to load single trials) and Dataset (to load the whole dataset).
    PARAMETERS 'ClassMosquito.py' 
        The class Dataset does not take any parameters
        The class Trial takes only one integer as an input: The trial nuber
        The class Track takes two integers parameters: the trial number and the track number
            Some tracks don't exist in this dataset,then an error message arises.

## 4. The other files: accessing_data.py and landing_capturing_area.py 
Other files are called on by the ClassMosquito.py file. 

In accessing_data.py file the following functions are present: 
    accessing_track(trial, track); 
    accessing_trial(trial); 
    accessing_extra_info(info_field). 

These functions are called on in ClassMosquito.py to access the data from the csv folder structure in the created path. 

In landing_capturing_area.py the follwing functions are present and called on in ClassMosquito.py: 
    landing_area(x, y, z)
    capturing_area(x, y, z)

They both take 3 float numbers as a single coordinate point to see if this point is within the landing or capturing region. 
