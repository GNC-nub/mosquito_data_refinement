'''
READ ME
    In this file are all the functions that are used within all the classes of ClassMosquito.py.
    These are general functions like opening a specific csv file, calculating the landing area,
    or generating the visualisation of the trap in 2 or 3D.

    The following functions are within this file:
        landing_area(x, y, z)
        capturing_area(x, y, z)
        accessing_track(trial, track)
        accessing_trial(trial)
        accessing_extra_info(info_field)
        getTrap()
        getTrap2D()



READ ME for:
    landing_area(x, y, z)
    capturing_area(x, y, z)
DESCRIPTION
    These functions determine if a single coordinate is within a capturing or a landing area

PARAMETERS
    Both functions take three float numbers (x, y and z) form a  single coordinate.

OUTPUT
    If the coordinate is within the capturing or respectively landing area the output is True, otherwise it is False.

RESULTS
    In the whole dataset 3197 mosquitos are landing. (9.38 %)
    In the whole dataset 1194 mosquitos are captured. (3.5 %)



READ ME for:
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



READ ME for:
    getTrap(body_lower_z=-0.38,body_upper_z=-0.083,inlet_upper_z=0,body_radius=0.15,inlet_radius=0.055)
    getTrap2D(body_lower_z=-0.38, body_upper_z=-0.083, body_radius=0.15, inlet_radius=0.055)
DISCRIPTION:
    Generates the correct values to create a trap in a figure either in 2D or 3D.


'''





import numpy as np
import os
import csv
import sys
import h5py
from loading_matlab_file import path_matlab_file1, path_csv_folder1
basemap_csv_path  = os.path.join(path_csv_folder1,'database_csv')


def landing_area_side(x, y, z, boundary=0.03, trap_height=0.388, trap_radius=0.15, inlet_height=0.083,
                      inlet_radius=0.055):
    r = np.sqrt((x ** 2) + (y ** 2))
    landing = False

    # landing_area of the inlet
    if -boundary < z < 0:
        if inlet_radius < r < inlet_radius + boundary:
            landing = True
    elif -(inlet_height - boundary) < z < -boundary:
        if inlet_radius - boundary < r < inlet_radius + boundary:
            landing = True
    # landing_area of the body side
    elif -trap_height < z < -(inlet_height + boundary):
        if trap_radius - boundary < r < trap_radius + boundary:
            landing = True
    return landing


def landing_area_top(x, y, z, boundary=0.03, trap_radius=0.15, inlet_height=0.083,
                     inlet_radius=0.055):
    r = np.sqrt((x ** 2) + (y ** 2))
    landing = False
    if -(inlet_height + boundary) < z < -(inlet_height - boundary):
        if inlet_radius - boundary < r < trap_radius + boundary:
            landing = True
    return landing

def landing_area(x, y, z, specific_area = 'whole', boundary=0.03):
    boolean = False
    if specific_area == 'whole':
        if landing_area_top(x, y, z, boundary=boundary) or landing_area_side(x, y, z, boundary=boundary):
            boolean = True
    elif specific_area == 'top':
        if landing_area_top(x, y, z, boundary=boundary):
            boolean = True
    elif specific_area == 'side':
        if landing_area_side(x, y, z, boundary=boundary):
            boolean = True
    return boolean

def capturing_area(x, y, z, boundary = 0.03, inlet_radius = 0.055):
    distance = np.sqrt((x**2)+(y**2))
    capture = False
    if 0 < z < boundary:
        if distance < (inlet_radius + boundary):
            capture = True
    elif -boundary < z < 0:
        if distance < inlet_radius:
            capture = True
    return capture



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




def getTrap(body_lower_z=-0.38,body_upper_z=-0.083,inlet_upper_z=0,body_radius=0.15,inlet_radius=0.055):
    theta = np.linspace(0, 2 * np.pi, 50)

    z_body = np.linspace(body_lower_z, body_upper_z, 50)
    theta_grid_body, z_grid_body = np.meshgrid(theta, z_body)
    x_grid_body = body_radius * np.cos(theta_grid_body)
    y_grid_body = body_radius * np.sin(theta_grid_body)

    z_inlet = np.linspace(body_upper_z, inlet_upper_z, 50)
    theta_grid_inlet, z_grid_inlet = np.meshgrid(theta, z_inlet)
    x_grid_inlet = inlet_radius * np.cos(theta_grid_inlet)
    y_grid_inlet = inlet_radius * np.sin(theta_grid_inlet)

    return x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet


def getTrap2D(body_lower_z=-0.38, body_upper_z=-0.083, body_radius=0.15, inlet_radius=0.055):
    inlet_r = [0, inlet_radius, inlet_radius, 0]
    inlet_z = [0, 0, body_upper_z, body_upper_z]

    body_r = [0, body_radius, body_radius, 0]
    body_z = [body_upper_z, body_upper_z, body_lower_z, body_lower_z]
    return inlet_r, inlet_z, body_r, body_z