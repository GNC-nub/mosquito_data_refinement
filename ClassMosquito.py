'''
READ ME
DISCRIPTION
    This code uses 3 classes from the file ClassMosquito.py:
        Track (to load single tracks), Trials (to load single trials) and Dataset (to load the whole dataset).
PARAMETERS
    The class Dataset does not take any parameters
    The class Trial takes only one integer as an input: The trial nuber
    The class Track takes two integers parameters: the trial number and the track number
        Some tracks don't exist in this dataset,then an error message arises.
'''
import numpy as np

from accessing_data import *
from landing_capturing_area import *

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import f_oneway, sem, norm
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import pandas as pd
from statistics import median

class Track:
    def __init__(self, trial_num, track_num): # if you want to load a track, fill in both trial and track, otherwise only whole trial
        self.total_mosquitos_per_trial = 50
        self.trial_num = trial_num
        self.track_num = track_num
        coordinate_list_track = accessing_track(trial_num, track_num)
        self.x = coordinate_list_track[0]
        self.y = coordinate_list_track[1]
        self.z = coordinate_list_track[2]
        self.time = coordinate_list_track[3]
        condition = accessing_extra_info('Condition')[trial_num]
        self.condition = condition
        self.header = f'Trial_{trial_num+1}_Track_{track_num}'

    def getTrack(self):
        return [self.x, self.y, self.z, self.time]

    def getTrap(self, body_lower_z=-0.38, body_upper_z=-0.083, inlet_upper_z=0, body_radius=0.15, inlet_radius=0.055):
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

    def getTrap2D(self, body_lower_z=-0.38, body_upper_z=-0.083, body_radius=0.15, inlet_radius=0.055):
        inlet_r = [0, inlet_radius, inlet_radius, 0]
        inlet_z = [0, 0, body_upper_z, body_upper_z]

        body_r = [0, body_radius, body_radius, 0]
        body_z = [body_upper_z, body_upper_z, body_lower_z, body_lower_z]
        return inlet_r, inlet_z, body_r, body_z



    def functionVelocity(self, x, delta_t):
        v = []
        v1 = ( -3 * x[0] + 4 * x[1] - x[2] )  / 2 * delta_t
        vn = (3 * x[-1] - 4 * x[-2] + x[-3] ) / 2 * delta_t
        v.append(v1)
        for i in range(1, len(x)):
            v_mid = ( x[i] - x[i-1] ) / delta_t
            v.append(v_mid)
        v.append(vn)
        return v

    def getVelocityTrack(self):
        x, y, z, t = self.getTrack()
        delta_t = t[1] = t[0]
        velocity_x = self.functionVelocity(x, delta_t)
        velocity_y = self.functionVelocity(y, delta_t)
        velocity_z = self.functionVelocity(z, delta_t)
        return velocity_x, velocity_y, velocity_z, delta_t

    def getAcelerationTrack(self):
        vx, vy, vz, delta_t = self.getVelocityTrack()
        acceleration_x = self.functionVelocity(vx, delta_t)
        acceleration_y = self.functionVelocity(vy, delta_t)
        acceleration_z = self.functionVelocity(vz, delta_t)
        return acceleration_x, acceleration_y, acceleration_z, delta_t


    def lastCoordinatesTrack(self):
        last_x =self.x[-1]
        last_y = self.y[-1]
        last_z = self.z[-1]
        last_time = self.time[-1]
        return [last_x, last_y, last_z, last_time]

    def plotTrack(self):
        last_x = self.x[-1]
        last_y = self.y[-1]
        last_z = self.z[-1]
        last_time = self.time[-1]

        ax = plt.figure().add_subplot(projection='3d')
        ax.plot(self.x, self.y, self.z)

        ax.scatter(last_x, last_y, last_z, color='r', marker='o')
        #plt.figtext(0.1, 0.9,
        #            f'The last coordinate is ({round(last_x, 3)}, {round(last_y, 3)}, {round(last_z, 3)}) at {round(last_time)} s',
        #            size=10, zorder=1)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('3D (x, y, z) plot of a single track')

        x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet = self.getTrap()
        ax.plot_surface(x_grid_body, y_grid_body, z_grid_body, alpha=0.5, color='b')
        ax.plot_surface(x_grid_inlet, y_grid_inlet, z_grid_inlet, alpha=0.5, color='b')
        ax.set_aspect('equal', adjustable='box')

        plt.show()

    def plotTheta2DTrack(self):
        r = []
        last_r, last_z = np.sqrt(self.x[-1]**2 + self.y[-1]**2), self.z[-1]
        for i in range(len(self.x)):
            r.append(np.sqrt(self.x[i]**2 + self.y[i]**2))
        plt.plot(r, self.z)
        plt.scatter(last_r, last_z, color='r', marker='o')
        inlet_r, inlet_z, body_r, body_z = self.getTrap2D()
        # Plot trap
        plt.fill(inlet_r, inlet_z, color='purple')
        plt.fill(body_r, body_z, color='purple')
        plt.xlim(0, 0.8)
        plt.ylim(-0.5, 0.5)
        plt.xlabel('r')
        plt.ylabel('z')
        plt.title('2D (r, z) plot of a single track')
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

class Trial:
    def __init__(self, trial_num):  # if you want to load a track, fill in both trial and track, otherwise only whole trial
        self.total_mosquitos_per_trial = 50
        self.trial_num = trial_num
        self.coordinate_list_trial = None
        self.condition = accessing_extra_info('Condition')[trial_num -1]
        self.header = f'Trial_{trial_num + 1}'
        self.landingpoints = None
        self.take_off_points = None
        self.num_tracks = None


    def initiateCoordinateList(self):
        self.coordinate_list_trial = accessing_trial(self.trial_num)  # [ [header, [x_coordinates], [y_coordinates], [z_coordinates], [time] ], exe... ]

    def initializeLandingPoints(self):
        landingpoints_list = []
        last_coordinates = self.lastCoordinatesTrial()
        for coordinate in last_coordinates:
            header, [x, y, z, time] = coordinate
            if landing_area(x, y, z) == True:
                landingpoints_list.append([x, y, z, time])
        self.landingpoints = landingpoints_list

    def initializeTakeOffPoints(self):
        take_off_points = []
        last_coordinates = self.firstCoordinatesTrial()
        for coordinate in last_coordinates:
            header, [x, y, z, time] = coordinate
            if landing_area(x, y, z) == True:
                take_off_points.append([x, y, z, time])
        self.take_off_points = take_off_points

    def getTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        return self.coordinate_list_trial  # [ [header, [x_coordinates], [y_coordinates], [z_coordinates], [time] ] , exe... ]

        # first/last coordinates
    def lastCoordinatesTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        storage = []
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            last_coordinates = [x[-1], y[-1], z[-1], time[-1]]
            storage.append([header, last_coordinates])
        return storage

    def firstCoordinatesTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        storage = []
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            first_coordinates = [x[0], y[0], z[0], time[0]]
            storage.append([header, first_coordinates])
        return storage

    def getNumTracks(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        num_tracks = len(self.coordinate_list_trial)
        return num_tracks

    def getEndTimeTrial(self):
        last_time_point = 0
        for track in self.getTrial():
            time = track[4]
            end_point = time[-1]
            if end_point > last_time_point:
                last_time_point = end_point
        return last_time_point

    def getStartTimeTrial(self):
        first_time_point = 1500
        for track in self.getTrial():
            time = track[4]
            first_point = time[0]
            if first_point < first_time_point:
                first_time_point = first_point
        return first_time_point

    def getDurationTrial(self):
        return self.getEndTimeTrial() - self.getStartTimeTrial()

    def getTrap(self, body_lower_z=-0.38, body_upper_z=-0.083, inlet_upper_z=0, body_radius=0.15, inlet_radius=0.055):
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

    def getTrap2D(self, body_lower_z=-0.38, body_upper_z=-0.083, body_radius=0.15, inlet_radius=0.055):
        inlet_r = [0, inlet_radius, inlet_radius, 0]
        inlet_z = [0, 0, body_upper_z, body_upper_z]

        body_r = [0, body_radius, body_radius, 0]
        body_z = [body_upper_z, body_upper_z, body_lower_z, body_lower_z]
        return inlet_r, inlet_z, body_r, body_z



#resting time
    def getRestingPairsTimesPoints(self, threshold = 0.02):
        self.initializeLandingPoints()
        self.initializeTakeOffPoints()
        potential_pairs = []

        for i_land, landing_point in enumerate(self.landingpoints):
            x_land, y_land, z_land, time_land = landing_point

            for i_takeoff, take_off_point in enumerate(self.take_off_points):
                    x_take, y_take, z_take, time_take = take_off_point
                    dx, dy, dz, dtime = (x_take - x_land), (y_take - y_land), (z_take - z_land), (time_take - time_land)

                    resting_time = dtime
                    distance = np.sqrt((dx ** 2) + (dy ** 2) + (dz ** 2))

                    if distance < threshold and resting_time > 0:
                        mean_x = (x_land + x_take) / 2
                        mean_y = (y_land + y_take) / 2
                        mean_z = (z_land + z_take) / 2
                        resting_point = [mean_x, mean_y, mean_z]
                        potential_pairs.append((distance, resting_time, i_land, i_takeoff, resting_point))

        potential_pairs.sort()

        pairs = {}
        used_takeoff_points = set()

        resting_times = []
        resting_points = []

        for distance, resting_time, i_land, i_takeoff, resting_point in potential_pairs:
            if i_takeoff not in used_takeoff_points:
                pairs[i_land] = i_takeoff
                used_takeoff_points.add(i_takeoff)
                resting_times.append(resting_time)
                resting_points.append(resting_point)
        return pairs, resting_times, resting_points

    def getRestingTimeTrial(self, radius = 0.02):
        pairs, resting_times, resting_points = self.getRestingPairsTimesPoints(radius)
        return resting_times
    def getRestingPointsTrial(self, radius = 0.02):
        pairs, resting_times, resting_points = self.getRestingPairsTimesPoints(radius)
        return resting_points
    def getRestingPairsTrial(self, radius =0.02):
        pairs, resting_times, resting_points = self.getRestingPairsTimesPoints(radius)
        return pairs
    def countRestingPairsTrial(self, radius = 0.02):
        return len(self.getRestingTimeTrial(radius))

# not really used, but keeping it for now: checks if thre are any rsting times that are more than 1320.
    def getFalseRestingTimeTrial(self, radius = 0.02, resting_time_limit = 1320): #based on the fact that most trials are from 180 s to 1500 s
        count = 0
        for resting_time in self.getRestingTimeTrial(radius):
            if resting_time > resting_time_limit:
                count += 1
        return count


#landing/capture
    def getMeasuredCatchesTrial(self):
        return int(accessing_extra_info("Catches")[self.trial_num - 1])

    def countSimulatedCatchesTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        count = 0
        for track in self.coordinate_list_trial:
            x = track[1][-1]
            y = track[2][-1]
            z = track[3][-1]
            if capturing_area(x, y, z) == True:
                count += 1
        return count

    def countLandingsTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        count = 0
        for track in self.coordinate_list_trial:
            x = track[1][-1]
            y = track[2][-1]
            z = track[3][-1]
            if landing_area(x, y, z) == True:
                count += 1
        return count

    def countTakeOffsTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        count = 0
        for track in self.coordinate_list_trial:
            x = track[1][0]
            y = track[2][0]
            z = track[3][0]
            if landing_area(x, y, z) == True:
                count += 1
        return count


#Take_off analysis
    def countLandingToCaptureTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        count = 0
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            if landing_area(x[0], y[0], z[0]) == True and capturing_area(x[-1], y[-1], z[-1]) == True:
                    count += 1
        return count

    def countLandingAgainTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        count = 0
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            if landing_area(x[0], y[0], z[0]) == True and landing_area(x[-1], y[-1], z[-1]) == True:
                count += 1
        return count

    def getCoordinatesLandingToCapture(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        take_off_coordinates = []
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            if landing_area(x[0], y[0], z[0]) == True and capturing_area(x[-1], y[-1], z[-1]) == True:
                take_off_coordinates.append([x[0], y[0], z[0]])
        return take_off_coordinates

    def getCoordinatesLandingAgain(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        take_off_coordinates = []
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            if landing_area(x[0], y[0], z[0]) == True and landing_area(x[-1], y[-1], z[-1]) == True:
                take_off_coordinates.append([x[0], y[0], z[0]])
        return take_off_coordinates

    def getTracksLandingAgain(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        landing_again_tracks = []
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            if landing_area(x[0], y[0], z[0]) == True and landing_area(x[-1], y[-1], z[-1]) == True:
                landing_again_tracks.append([x, y, z])
        return landing_again_tracks

    def getTracksLandingtoCapture(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        landing_to_capture_tracks = []
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            if landing_area(x[0], y[0], z[0]) == True and capturing_area(x[-1], y[-1], z[-1]) == True:
                landing_to_capture_tracks.append([x, y, z])
        return landing_to_capture_tracks



#plots
    def plotTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()

        ax = plt.figure().add_subplot(projection='3d')

        for track in self.coordinate_list_trial:
            x = track[1]
            y = track[2]
            z = track[3]

            last_x = x[-1]
            last_y = y[-1]
            last_z = z[-1]

            ax.plot(x, y, z)
            ax.scatter(last_x, last_y, last_z, color= 'r', marker='o', s = 1)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('a) 3d (x, y, z) plot of all the tracks in one trial')

        #plot trap
        x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet = self.getTrap()
        ax.plot_surface(x_grid_body, y_grid_body, z_grid_body, alpha=0.5, color='b')
        ax.plot_surface(x_grid_inlet, y_grid_inlet, z_grid_inlet, alpha=0.5, color='b')
        ax.set_aspect('equal', adjustable='box')
        plt.show()

    def plotLandingPointsTrial(self):
        ax = plt.figure().add_subplot(projection='3d')
        for track in self.lastCoordinatesTrial():
            last_x, last_y, last_z, last_time = track[1]
            if landing_area(last_x, last_y, last_z) == True:
                ax.scatter(last_x, last_y, last_z, color='r', marker='o')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet = self.getTrap()
        ax.plot_surface(x_grid_body, y_grid_body, z_grid_body, alpha=0.5, color='b')
        ax.plot_surface(x_grid_inlet, y_grid_inlet, z_grid_inlet, alpha=0.5, color='b')

        plt.show()

    def plotTheta2DTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            r = []
            last_r, last_z = np.sqrt(x[-1]**2 + y[-1]**2),z[-1]
            for i in range(len(x)):
                r.append(np.sqrt(x[i]**2 + y[i]**2))
            plt.plot(r, z)
            plt.scatter(last_r, last_z, color='r', marker='o', s = 1)
        # Plot trap
        inlet_r, inlet_z, body_r, body_z = self.getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple')
        plt.fill(body_r, body_z, color='purple')
        plt.xlim(0, 0.8)
        plt.ylim(-0.5, 0.5)
        plt.xlabel('r')
        plt.ylabel('z')
        plt.title('b) 2D (r, z) plot of all the tracks in one trial')
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

    def plotThetaLandingPointsTrial(self):
        r_list =[]
        z_list =[]
        self.initializeLandingPoints()
        for coordinate in self.landingpoints:
            x, y, z, time = coordinate
            r = np.sqrt(x**2 + y**2)
            r_list.append(r)
            z_list.append(z)
        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)
        plt.hist2d(r_list, z_list, bins=[5, 30], cmap=custom_cmap)
        plt.colorbar(label='Density')
        plt.xlabel('r coordinate')
        plt.ylabel('z coordinate')
        plt.title('Density Heatmap of Landing Coordinates')
        #plot trap
        inlet_r, inlet_z, body_r, body_z = self.getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple')
        plt.fill(body_r, body_z, color='purple')
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

    def plotTrackLandingAgainTrial(self, nr = 8):
        ax = plt.figure().add_subplot(projection='3d')
        tracks = self.getTracksLandingAgain()
        x, y, z = tracks[nr]
        last_x, last_y, last_z = x[-1], y[-1], z[-1]
        ax.plot(x, y, z, color = 'green')
        ax.scatter(last_x, last_y, last_z, color='r', marker='o', s=10)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('a) A take-off track that lands again')

        #plot trap
        x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet = self.getTrap()
        ax.plot_surface(x_grid_body, y_grid_body, z_grid_body, alpha=0.5, color='b')
        ax.plot_surface(x_grid_inlet, y_grid_inlet, z_grid_inlet, alpha=0.5, color='b')
        ax.set_aspect('equal', adjustable='box')
        plt.show()

    def plotTrackLandingToCaptureTrial(self, nr = 1):
        ax = plt.figure().add_subplot(projection='3d')
        tracks = self.getTracksLandingtoCapture()
        x, y, z = tracks[nr]
        last_x, last_y, last_z = x[-1], y[-1], z[-1]
        ax.plot(x, y, z, color = 'green')
        ax.scatter(last_x, last_y, last_z, color='r', marker='o', s=10)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('b) A take-off track that leads to capture')

        #plot trap
        x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet = self.getTrap()
        ax.plot_surface(x_grid_body, y_grid_body, z_grid_body, alpha=0.5, color='b')
        ax.plot_surface(x_grid_inlet, y_grid_inlet, z_grid_inlet, alpha=0.5, color='b')
        ax.set_aspect('equal', adjustable='box')
        plt.show()

class Dataset:

    def __init__(self):  # if you want to load a track, fill in both trial and track, otherwise only whole trial
        self.total_mosquitos_per_trial = 50
        self.total_number_trials = 65

        self.coordinate_list_dataset = None
        self.trialobjects = None

        self.catches_without = None
        self.catches_with_heat = None
        self.catches_with_heat_water = None

        self.landing_without = None
        self.landing_with_heat = None
        self.landing_with_heat_water = None

        self.landing_all_trials = None
        self.catches_all_trials = None

        self.num_r_cells_matrix = 20
        self.num_z_cells_matrix = 50
        self.upper_z_coord_matrix = 0
        self.lower_z_coord_matrix = -0.5
        self.first_r_coord_matrix = 0
        self.last_r_coord_matrix = 0.2
        self.area_cell_matrix = 0.0001
        self.r_edges_matrix = np.linspace(self.first_r_coord_matrix, self.last_r_coord_matrix,
                                          self.num_r_cells_matrix + 1)
        self.z_edges_matrix = np.linspace(self.lower_z_coord_matrix, self.upper_z_coord_matrix,
                                          self.num_z_cells_matrix + 1)

#intitializing data
    def initializeCoordinateList(self):
        self.coordinate_list_dataset = []
        for i in range(1, self.total_number_trials + 1):
            coordinate_trial = accessing_trial(i)
            self.coordinate_list_dataset.append(coordinate_trial)

    def initializeMeasuredCatchesConditions(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        data_without = []
        data_with_heat = []
        data_with_heat_water = []
        for trial_object in self.getTrialObjects():
            if trial_object.condition == "without":
                data_without.append(trial_object.getMeasuredCatchesTrial())
            elif trial_object.condition == "with_heat":
                data_with_heat.append(trial_object.getMeasuredCatchesTrial())
            else:
                data_with_heat_water.append(trial_object.getMeasuredCatchesTrial())
        self.measured_catches_without = data_without
        self.measured_catches_with_heat = data_with_heat
        self.measured_catches_with_heat_water = data_with_heat_water

    def initializeSimulatedCatchesConditions(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        data_without = []
        data_with_heat = []
        data_with_heat_water = []
        for trial_object in self.getTrialObjects():
            if trial_object.condition == "without":
                data_without.append(trial_object.countSimulatedCatchesTrial())
            elif trial_object.condition == "with_heat":
                data_with_heat.append(trial_object.countSimulatedCatchesTrial())
            else:
                data_with_heat_water.append(trial_object.countSimulatedCatchesTrial())
        self.catches_without = data_without
        self.catches_with_heat = data_with_heat
        self.catches_with_heat_water = data_with_heat_water

    def initializeLandingConditions(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        data_without = []
        data_with_heat = []
        data_with_heat_water = []
        for trial_object in self.getTrialObjects():
            if trial_object.condition == "without":
                data_without.append(trial_object.countLandingsTrial())
            elif trial_object.condition == "with_heat":
                data_with_heat.append(trial_object.countLandingsTrial())
            else:
                data_with_heat_water.append(trial_object.countLandingsTrial())
        self.landing_without = data_without
        self.landing_with_heat = data_with_heat
        self.landing_with_heat_water = data_with_heat_water

    def initializeLandingAllTrials(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        landing_count_per_trial = []
        for trial_object in self.getTrialObjects():
            landing_count_per_trial.append(trial_object.countLandingsTrial())
        self.landing_all_trials = landing_count_per_trial

    def initializeMeasuredCatchesAllTrials(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        catches_count_per_trial = []
        for trial_object in self.getTrialObjects():
            catches_count_per_trial.append(trial_object.getMeasuredCatchesTrial())
        self.catches_all_trials = catches_count_per_trial

    def initializeSimulatedCatchesAllTrials(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        catches_count_per_trial = []
        for trial_object in self.getTrialObjects():
            catches_count_per_trial.append(trial_object.countSimulatedCatchesTrial())
        self.catches_all_trials = catches_count_per_trial

#getting data

    def getTrialObjects(self):
        object_array = []
        for i in range(1, self.total_number_trials + 1):
            obj = Trial(trial_num=i)
            object_array.append(obj)
        return object_array

    def getTrap(self, body_lower_z=-0.38, body_upper_z=-0.083, inlet_upper_z=0, body_radius=0.15, inlet_radius=0.055):
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

    def getTrap2D(self, body_lower_z=-0.38, body_upper_z=-0.083, body_radius=0.15, inlet_radius=0.055):
        inlet_r = [0, inlet_radius, inlet_radius, 0]
        inlet_z = [0, 0, body_upper_z, body_upper_z]

        body_r = [0, body_radius, body_radius, 0]
        body_z = [body_upper_z, body_upper_z, body_lower_z, body_lower_z]
        return inlet_r, inlet_z, body_r, body_z

    def getStartTimes(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        beginning_points = []
        for trial_object in self.trialobjects:
            if trial_object.trial_num != 59:
                beginning_points.append(trial_object.getStartTimeTrial())
        return beginning_points
    def getEndTimes(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        end_points = []
        for trial_object in self.trialobjects:
            if trial_object.trial_num != 59:
                end_points.append(trial_object.getEndTimeTrial())
        return end_points

    def getNumTracksPerTrial(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        num_tracks_per_trial = []
        for trial_object in self.trialobjects:
            num_tracks_per_trial.append(trial_object.getNumTracks())
        return num_tracks_per_trial

    def getAvarageLengthTracks(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        times = []
        for trial_object in self.trialobjects:
            trial_object.initiateCoordinateList()
            for track in trial_object.coordinate_list_trial:
                header, x, y, z, time = track
                length = time[-1] - time[0]
                times.append(length)
        return np.average(times)


# landing/take_off points
    def getlandingPointsTheta(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list = []
        z_list = []
        for trial_object in self.trialobjects:
            trial_object.initializeLandingPoints()
            for coordinate in trial_object.landingpoints:
                x, y, z, time = coordinate
                r = np.sqrt(x ** 2 + y ** 2)
                r_list.append(r)
                z_list.append(z)
        return r_list, z_list
    def getTakeOffPointsTheta(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list = []
        z_list = []
        for trial_object in self.trialobjects:
            trial_object.initializeTakeOffPoints()
            for coordinate in trial_object.take_off_points:
                x, y, z, time = coordinate
                r = np.sqrt(x ** 2 + y ** 2)
                r_list.append(r)
                z_list.append(z)
        return r_list, z_list
    def getlandingPointsThetaPerCondition(self, condition):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list = []
        z_list = []
        for trial_object in self.trialobjects:
            trial_object.initializeLandingPoints()
            if trial_object.condition == condition:
                for coordinate in trial_object.landingpoints:
                    x, y, z, time = coordinate
                    r = np.sqrt(x ** 2 + y ** 2)
                    r_list.append(r)
                    z_list.append(z)
        return r_list, z_list


# resting time
    def getAvarageRestingTime(self, lower_boundary = 0):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_list = []
        for trial_object in self.trialobjects:
            resting_time_trial = trial_object.getRestingTimeTrial()
            for point in resting_time_trial:
                if point > lower_boundary:
                    resting_time_list.append(point)
        mean = sum(resting_time_list)/len(resting_time_list)
        std = np.std(resting_time_list)
        return f'{mean}' + u"\u00B1" + f'{std}'

    def getMedianRestingTime(self, lower_boundary = 0):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_list = []
        for trial_object in self.trialobjects:
            resting_time_trial = trial_object.getRestingTimeTrial()
            for point in resting_time_trial:
                if point > lower_boundary:
                    resting_time_list.append(point)
        resting_time_median = median(resting_time_list)
        first_quartile = np.quantile(resting_time_list, 0.25)
        third_quartile = np.quantile(resting_time_list, 0.75)
        return first_quartile, resting_time_median, third_quartile

    def getLongestRestingTime(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_list = []
        for trial_object in self.trialobjects:
            resting_time_trial = trial_object.getRestingTimeTrial()
            for point in resting_time_trial:
                resting_time_list.append(point)
        longest_resting_time = max(resting_time_list)
        return longest_resting_time

    def countSmallRestingTimes(self, seconds):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_list = []
        for trial_object in self.trialobjects:
            resting_time_trial = trial_object.getRestingTimeTrial()
            for point in resting_time_trial:
                if point < seconds:
                    resting_time_list.append(point)
        return len(resting_time_list)

    def countTotalRestingTimes(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_list = []
        for trial_object in self.trialobjects:
            resting_time_trial = trial_object.getRestingTimeTrial()
            for point in resting_time_trial:
                    resting_time_list.append(point)
        return len(resting_time_list)

    def countLandingPoints(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        landing_points = 0
        for trial_object in self.trialobjects:
            landing_points += trial_object.countLandingsTrial()
        return landing_points

    def countCapturingPoints(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        capturing_points = 0
        for trial_object in self.trialobjects:
            capturing_points += trial_object.countSimulatedCatchesTrial()
        return capturing_points
    def getLastCoordinatesDataset(self):
        if self.coordinate_list_dataset == None:
            self.initializeCoordinateList()
        dataset_coordinates = []
        for i in range(1, self.total_number_trials+1):
            data = accessing_trial(i)
            trial_data = []
            for track in data:
                header, x, y, z, time = track
                last_coordinates = [x[-1], y[-1], z[-1], time[-1]]
                trial_data.append([header, last_coordinates])
            dataset_coordinates.append(trial_data)
        return dataset_coordinates # OUTPUT [ [ [ [header], [x_coordinate,  y_coordinate, z_coordinate, time_point] ], exe...] ]


#analysis
    def testAnovaConditions(self, data):
        without, with_heat, with_heat_water = data
        f_stat, p_value = f_oneway(without, with_heat, with_heat_water)
        return f_stat, p_value
    def testTukeysHSDConditions(self, data):
        without, with_heat, with_heat_water = data
        df = pd.DataFrame({
            'Values': np.concatenate([without, with_heat, with_heat_water]),
            'Condition': ['without'] * len(without) + ['with_heat'] * len(with_heat) + ['with_heat_water'] * len(with_heat_water)})
        tukey = pairwise_tukeyhsd(endog=df['Values'],groups=df['Condition'], alpha=0.05)
        tukey_results = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])
        significant_results = tukey_results[tukey_results['p-adj'] < 0.05]
        significant_pairs = []
        for i in range(len(significant_results)):
            group1 = significant_results.iloc[i]['group1']
            group2 = significant_results.iloc[i]['group2']
            significant_pairs.append([group1, group2])
        return significant_pairs

    def testNormalDistribution(self, data):
        plt.subplot(1,2,1)
        plt.hist(data, bins = 10)
        plt.xlabel('Data values')
        plt.ylabel('Density')
        plt.title('Histogram')
        plt.subplot(1,2,2)
        stats.probplot(data, dist="norm", plot=plt)
        plt.title("Q-Q Plot")
        plt.show()

    def testConfidenceInterval(self, data):
        mean = np.mean(data)
        se = sem(data)
        ci = 1.96 * se # for 95% confidence interval
        return f'{mean}'+u"\u00B1"+f'{ci}'


#calculating percentages
    def calculatingPercentagesLandingAgain(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        total_take_offs = 0
        land_again = 0
        for trial_object in self.trialobjects:
            land_again += trial_object.countLandingAgainTrial()
            total_take_offs += trial_object.countTakeOffsTrial()
        percentage_land_again = land_again / total_take_offs * 100
        return percentage_land_again

    def calculatingPercentagesLandingToCapture(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        total_take_offs = 0
        land_to_capture = 0
        for trial_object in self.trialobjects:
            land_to_capture += trial_object.countLandingToCaptureTrial()
            total_take_offs += trial_object.countTakeOffsTrial()
        percentage_land_to_capture = land_to_capture / total_take_offs * 100
        return percentage_land_to_capture


#plotting boxplot/histogram

    def plotDuration(self):
        fig, ax = plt.subplots()
        for i, (start, end) in enumerate(zip(self.getStartTimes(), self.getEndTimes())):
            ax.plot([start, end], [i, i], marker='o', color='teal')
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Trial number")
        ax.set_title("Trial Durations")
        plt.show()

    def plotViolinStartEndTimes(self):
        plt.subplot(1, 2, 1)
        data_start = self.getStartTimes()
        plt.violinplot(data_start)
        plt.title(f"Starting times")
        plt.xlabel('Start')
        plt.ylabel('Time (s)')

        plt.subplot(1, 2, 2)
        data_end = self.getEndTimes()
        plt.violinplot(data_end)
        plt.title(f"Ending times")
        plt.xlabel('End')
        plt.ylabel('Time (s)')
        plt.show()
        plt.show()


    def plotBoxplotCatchesConditions(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        if self.catches_without == None and self.catches_with_heat == None and self.catches_with_heat_water == None:
            self.initializeSimulatedCatchesConditions()
        data = [self.catches_without, self.catches_with_heat, self.catches_with_heat_water]

        #do the anova
        f_stat, p_value = self.testAnovaConditions(data)
        if p_value > 0.05:
            print('plotBoxplotCatchesConditions: No significant differences between any group.\n')
        plt.boxplot(data)
        plt.xticks([1, 2, 3], ['Without', 'With heat', 'With heat and water'])
        plt.title('b) Captures on M-tego trap per condition')
        plt.ylabel('Captures per trial')
        plt.show()

    def plotBoxplotLandingCondition(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        if self.landing_without == None and self.landing_with_heat == None and self.landing_with_heat_water == None:
            self.initializeLandingConditions()
        #anova
        data = [self.landing_without, self.landing_with_heat, self.landing_with_heat_water]
        f_stat, p_value = self.testAnovaConditions(data)
        if p_value < 0.05:
            list_significant_pairs = self.testTukeysHSDConditions(data)
            group1 = list_significant_pairs[0][0]
            group2 = list_significant_pairs[0][1]
            print(f'plotBoxplotLandingCondition: Significant differences between: {group1} and {group2}, with {p_value}.\n')
            if len(list_significant_pairs) > 1:
                print('More significant pairs')
        else:
            print('plotBoxplotLandingCondition: No significant differences between any group.\n')
        plt.boxplot(data)
        plt.xticks([1, 2, 3], ['Without', 'With heat', 'With heat and water'])
        plt.title('c) Landings on M-tego trap per condition')
        plt.xlabel('Condition')
        plt.ylabel('Landings per trial')
        plt.show()

    def plotBoxplotCatchesVsLandings(self):
        if self.landing_all_trials == None:
            self.initializeLandingAllTrials()
        if self.catches_all_trials == None:
            self.initializeSimulatedCatchesAllTrials()
        data = [self.catches_all_trials, self.landing_all_trials]
        #paired t test
        t_stat, p_value = stats.ttest_rel(self.catches_all_trials, self.landing_all_trials)
        if p_value < 0.05:
            print(f'plotBoxplotCatchesVsLandings: Significant, p-value = {p_value}.\n')
        else:
            print(f'plotBoxplotCatchesVsLandings: Not significant.\n')

        plt.boxplot(data)
        plt.xticks([1, 2], ['Captures', 'Landings'])
        plt.title('a) Captures and Landings per trial')
        plt.ylabel('Number of captures or landings per trial')
        plt.show()


    def plotHistogramRestingTime(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_list = []
        for trial_object in self.trialobjects:
            resting_time_trial = trial_object.getRestingTimeTrial()
            for point in resting_time_trial:
                resting_time_list.append(point)
        plt.hist(resting_time_list, bins=20)
        plt.title('Resting time histogram')
        plt.xlabel('Resting time (s)')
        plt.ylabel('Number of mosquitos resting')
        plt.show()

    def plotHistogramRestingTimeZoomedIn(self, low, high):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_list = []
        for trial_object in self.trialobjects:
            resting_time_trial = trial_object.getRestingTimeTrial()
            for point in resting_time_trial:
                if low < point < high:
                    resting_time_list.append(point)
        plt.hist(resting_time_list, bins=30)
        plt.title(f'A histogram of the resting times above {low} below {high} seconds')
        plt.xlabel('Resting time (s)')
        plt.ylabel("Number of mosquito's resting")
        plt.show()

    def plotHistogramRestingTimeCondition(self, condition):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_list = []
        for trial_object in self.trialobjects:
            condition_trial = trial_object.condition
            if condition_trial == condition:
                resting_time_trial = trial_object.getRestingTimeTrial()
                for point in resting_time_trial:
                    resting_time_list.append(point)
        plt.hist(resting_time_list, bins=20)
        plt.title(f"Resting time histogram of the trap '{condition}' short-range cues")
        plt.xlabel('Resting time (s)')
        plt.ylabel('Number of mosquitos resting')
        plt.show()

    def plotBoxplotRestingTimeCondition(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_without = []
        resting_time_with_heat = []
        resting_time_with_heat_water = []
        for trial_object in self.trialobjects:
            condition_trial = trial_object.condition
            resting_time_trial = trial_object.getRestingTimeTrial()
            if condition_trial == 'without':
                for point in resting_time_trial:
                    resting_time_without.append(point)
            elif condition_trial == 'with_heat':
                for point in resting_time_trial:
                    resting_time_with_heat.append(point)
            else:
                for point in resting_time_trial:
                    resting_time_with_heat_water.append(point)
        data = [resting_time_without, resting_time_with_heat, resting_time_with_heat_water]
        #anova
        f_stat, p_value = self.testAnovaConditions(data)
        if p_value < 0.05:
            list_significant_pairs = self.testTukeysHSDConditions(data)
            pair1 = f'{list_significant_pairs[0][0]} and {list_significant_pairs[0][1]}'
            pair2 = f'{list_significant_pairs[1][0]} and {list_significant_pairs[1][1]}'
            plt.figtext(0.01, 0.05, f'Significant differences between:\n{pair1}\n{pair2}', size=5)
        else:
            plt.figtext(0.01, 0.09, 'No significant differences\nbetween any group', size=5)

        plt.boxplot(data)
        plt.title(f"Boxplot of the resting time with different trap short-range cues")
        plt.xticks([1, 2, 3], ['Without', 'With heat', 'With heat and water'])
        plt.xlabel('Condition')
        plt.ylabel('Resting time (s)')
        plt.show()

    def plotViolinplotRestingTimeCondition(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_without = []
        resting_time_with_heat = []
        resting_time_with_heat_water = []
        for trial_object in self.trialobjects:
            condition_trial = trial_object.condition
            resting_time_trial = trial_object.getRestingTimeTrial()
            if condition_trial == 'without':
                for point in resting_time_trial:
                    resting_time_without.append(point)
            elif condition_trial == 'with_heat':
                for point in resting_time_trial:
                    resting_time_with_heat.append(point)
            else:
                for point in resting_time_trial:
                    resting_time_with_heat_water.append(point)
        data = [resting_time_without, resting_time_with_heat, resting_time_with_heat_water]
        # anova
        f_stat, p_value = self.testAnovaConditions(data)
        if p_value < 0.05:
            list_significant_pairs = self.testTukeysHSDConditions(data)
            pair1 = f'{list_significant_pairs[0][0]} and {list_significant_pairs[0][1]}'
            pair2 = f'{list_significant_pairs[1][0]} and {list_significant_pairs[1][1]}'
            plt.figtext(0.01, 0.05, f'Significant differences between:\n{pair1}\n{pair2}', size=5)
        else:
            plt.figtext(0.01, 0.09, 'No significant differences\nbetween any group', size=5)

        plt.violinplot(data)
        plt.title(f"The resting time with different trap short-range cues")
        plt.xticks([1, 2, 3], ['Without', 'With heat', 'With heat and water'])
        plt.xlabel('Condition')
        plt.ylabel('Resting time (s)')
        plt.show()


    def plotBoxplotRadiusAssociations(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        pairs1 = []
        pairs2 = []
        pairs3 = []
        pairs4 = []
        for trial_object in self.trialobjects:
            pairs1.append(trial_object.countRestingPairsTrial(0.01))
            pairs2.append(trial_object.countRestingPairsTrial(0.02))
            pairs3.append(trial_object.countRestingPairsTrial(0.03))
            pairs4.append(trial_object.countRestingPairsTrial(0.04))
        data = [pairs1, pairs2, pairs3, pairs4]
        plt.boxplot(data)
        plt.title('Nr of of associated landings with take offs per radius threshold')
        plt.xticks([1, 2, 3, 4], ['r = 0.01', 'r = 0.02', 'r = 0.03', 'r = 0.04'])
        plt.ylabel('nr of associated landings with take offs')
        plt.show()

    def plotBoxplotCountLandingTakeOffAccosiation(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        landings, takeoffs, associated_pairs = [], [], []
        for trial_object in self.trialobjects:
            landings.append(trial_object.countLandingsTrial())
            takeoffs.append(trial_object.countTakeOffsTrial())
            associated_pairs.append(trial_object.countRestingPairsTrial())
        data = [landings, takeoffs, associated_pairs]
        plt.boxplot(data)
        plt.title('Boxplot of the number of landings, take-offs and landing / take_off pairs')
        plt.xticks([1, 2, 3], ['Landings', 'Take-offs', 'associated pairs'])
        plt.ylabel('nr of associated landings with take offs')
        plt.show()




# creating matrixes
    def getMatrixNormilizingVolume(self):
        distance_row = np.linspace(0.005, self.last_r_coord_matrix - 0.005, self.num_r_cells_matrix)
        distance_matrix = np.tile(distance_row, (self.num_z_cells_matrix, 1))
        volume_matrix = distance_matrix * 2 * np.pi * self.area_cell_matrix
        volume_matrix_np = np.array(volume_matrix)
        return volume_matrix_np

    def getMatrixLandingPoints(self):
        r, z = self.getlandingPointsTheta()
        landingpoint_count_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r, z, bins=(
            self.r_edges_matrix, self.z_edges_matrix))
        landingpoint_count_matrix = landingpoint_count_matrix.T
        landingpoint_count_matrix_np = np.array(landingpoint_count_matrix)
        return landingpoint_count_matrix_np

    def getMatrixTakeOffPoints(self):
        r, z = self.getTakeOffPointsTheta()
        takeoffpoints_count_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r, z, bins=(
            self.r_edges_matrix, self.z_edges_matrix))
        takeoffpoints_count_matrix = takeoffpoints_count_matrix.T
        takeoffpoints_count_matrix_np = np.array(takeoffpoints_count_matrix)
        return takeoffpoints_count_matrix_np

    def getMatrixLandingPointsPerCondition(self, condition):
        r, z = self.getlandingPointsThetaPerCondition(condition)
        landingpoint_count_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r, z, bins=(
            self.r_edges_matrix, self.z_edges_matrix))
        landingpoint_count_matrix = landingpoint_count_matrix.T
        landingpoint_count_matrix_np = np.array(landingpoint_count_matrix)
        return landingpoint_count_matrix_np

    def getMatrixRestingTimes(self, lower_time_boundary = 0, upper_time_boundary = 1500):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list, z_list, w_list = [], [], []
        for trial_object in self.trialobjects:
            resting_times = trial_object.getRestingTimeTrial()
            points = trial_object.getRestingPointsTrial()
            for i, resting_time in enumerate(resting_times):
                if lower_time_boundary < resting_time < upper_time_boundary:
                    w_list.append(resting_time)
                    x, y, z = points[i]
                    r = np.sqrt(x ** 2 + y ** 2)
                    r_list.append(r)
                    z_list.append(z)
        resting_time_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r_list, z_list, bins=(
            self.r_edges_matrix, self.z_edges_matrix), weights=w_list)
        resting_time_matrix = resting_time_matrix.T
        resting_time_matrix_np = np.array(resting_time_matrix)

        resting_time_count_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r_list, z_list, bins=(
            self.r_edges_matrix, self.z_edges_matrix))
        resting_time_count_matrix = resting_time_count_matrix.T
        resting_time_count_matrix_np = np.array(resting_time_count_matrix)
        return resting_time_matrix_np, resting_time_count_matrix_np

    def getMatrixRestingTimePerCondition(self, condition):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list, z_list, w_list = [], [], []
        for trial_object in self.trialobjects:
            if trial_object.condition == condition:
                resting_times = trial_object.getRestingTimeTrial()
                points = trial_object.getRestingPointsTrial()

                for resting_time in resting_times:
                    w_list.append(resting_time)

                for coordinate in points:
                    x, y, z = coordinate
                    r = np.sqrt(x ** 2 + y ** 2)
                    r_list.append(r)
                    z_list.append(z)

        resting_time_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r_list, z_list, bins=(
            self.r_edges_matrix, self.z_edges_matrix), weights=w_list)
        resting_time_matrix = resting_time_matrix.T
        resting_time_matrix_np = np.array(resting_time_matrix)

        resting_time_count_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r_list, z_list, bins=(
            self.r_edges_matrix, self.z_edges_matrix))
        resting_time_count_matrix = resting_time_count_matrix.T
        resting_time_count_matrix_np = np.array(resting_time_count_matrix)
        return resting_time_matrix_np, resting_time_count_matrix_np

    def getMatrixCaptureProbability(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list, z_list, w_list = [], [], []
        for trial_object in self.trialobjects:
            coordinates = trial_object.getCoordinatesLandingToCapture()
            for coordinate in coordinates:
                x, y, z = coordinate
                r = np.sqrt(x ** 2 + y ** 2)
                r_list.append(r)
                z_list.append(z)
        landing_to_capture_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r_list, z_list, bins=(
            self.r_edges_matrix, self.z_edges_matrix))
        landing_to_capture_matrix = landing_to_capture_matrix.T
        landing_to_capture_matrix_np = np.array(landing_to_capture_matrix)
        return landing_to_capture_matrix_np

    def getMatrixLandingAgainProbability(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list, z_list, w_list = [], [], []
        for trial_object in self.trialobjects:
            coordinates = trial_object.getCoordinatesLandingAgain()
            for coordinate in coordinates:
                x, y, z = coordinate
                r =  np.sqrt(x ** 2 + y ** 2)
                r_list.append(r)
                z_list.append(z)
        landing_again_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r_list, z_list, bins=(
            self.r_edges_matrix, self.z_edges_matrix))
        landing_again_matrix = landing_again_matrix.T
        landing_again_matrix_np = np.array(landing_again_matrix)
        return landing_again_matrix_np

    def getMatrixVelocity(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()




#plotting 2D heatmaps

    def plotHeatmapLandingPoints(self):
        r_list, z_list = self.getlandingPointsTheta()
        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)
        plt.hist2d(r_list, z_list, bins=[10, 40], cmap=custom_cmap)
        plt.colorbar(label='Density')
        plt.xlabel('r coordinate')
        plt.ylabel('z coordinate')
        plt.title('Density Heatmap of Landing Coordinates Dataset')
        # plot trap
        inlet_r, inlet_z, body_r, body_z = self.getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple')
        plt.fill(body_r, body_z, color='purple')
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

    def plotHeatmapLandingPointNormilized(self):
        volume_matrix = self.getMatrixNormilizingVolume()
        landingpoint_matrix = self.getMatrixLandingPoints()

        density_matrix = landingpoint_matrix / volume_matrix

        #plot matrix
        fig = plt.figure()
        ax = fig.add_subplot(132, title='Density Heatmap: Landing points per volume', aspect='equal')
        X, Y = np.meshgrid(self.r_edges_matrix, self.z_edges_matrix)

        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)
        c = ax.pcolormesh(X, Y, density_matrix, cmap=custom_cmap)
        ax.set_facecolor('white')

        plt.colorbar(c, ax=ax, fraction=0.065, pad = 0.13, label='Density (points/m$^3$)')
        plt.xlabel('r coordinates')
        plt.ylabel('z coordinates')

        # plot trap
        inlet_r, inlet_z, body_r, body_z = self.getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple')
        plt.fill(body_r, body_z, color='purple')
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

    def plotHeatmapRestingPoints(self, lower_time_boundary = 0, upper_time_boundary = 1500, title = 'Density Heatmap: Resting points per volume'):
        volume_matrix = self.getMatrixNormilizingVolume()
        resting_time_matrix, resting_time_count_matrix = self.getMatrixRestingTimes(lower_time_boundary,
                                                                                    upper_time_boundary)
        density_matrix = resting_time_count_matrix / volume_matrix

        # plot matrix
        fig = plt.figure()
        ax = fig.add_subplot(132, title=title, aspect='equal')
        X, Y = np.meshgrid(self.r_edges_matrix, self.z_edges_matrix)

        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)
        c = ax.pcolormesh(X, Y, density_matrix, cmap=custom_cmap)
        ax.set_facecolor('white')

        plt.colorbar(c, ax=ax, fraction=0.065, pad=0.13, label='Density (points/m$^3$)')
        plt.xlabel('r coordinates')
        plt.ylabel('z coordinates')

        # plot trap
        inlet_r, inlet_z, body_r, body_z = self.getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple')
        plt.fill(body_r, body_z, color='purple')
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()  # nor #ormil

    def plotHeatmapRestingTimes(self, lower_time_boundary = 0, upper_time_boundary = 1500, title = 'Heatmap Resting Time'):
        volume_matrix = self.getMatrixNormilizingVolume()
        resting_time_matrix, resting_time_count_matrix = self.getMatrixRestingTimes(lower_time_boundary, upper_time_boundary)

        resting_time_count_matrix = np.nan_to_num(resting_time_count_matrix, nan=0.0)

        resting_time_matrix_norm = resting_time_matrix / volume_matrix
        resting_time_matrix_average = np.divide(resting_time_matrix_norm, resting_time_count_matrix,
                                                where=resting_time_count_matrix != 0)

        # plot matrix
        fig = plt.figure()
        ax = fig.add_subplot(132, title=title, aspect='equal')
        X, Y = np.meshgrid(self.r_edges_matrix, self.z_edges_matrix)

        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)
        c = ax.pcolormesh(X, Y, resting_time_matrix_average, cmap=custom_cmap)
        ax.set_facecolor('white')

        plt.colorbar(c, ax=ax, fraction=0.065, pad=0.13, label='Density resting time (s/m$^3$)')
        plt.xlabel('r coordinates')
        plt.ylabel('z coordinates')

        # plot trap
        inlet_r, inlet_z, body_r, body_z = self.getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple')
        plt.fill(body_r, body_z, color='purple')
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()


    def plotHeatmapRestingTimePerCondition(self):
        fig, ax = plt.subplots(1, 3, figsize=(10, 5))
        fig.suptitle('Density Heatmap: Average resting time per volume', fontsize = 16)
        X, Y = np.meshgrid(self.r_edges_matrix, self.z_edges_matrix)
        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)

        volume_matrix = self.getMatrixNormilizingVolume()
        conditions = ['without', 'with_heat', 'with_heat_water']
        for index, condition in enumerate(conditions):
            resting_time_matrix, resting_time_count_matrix = self.getMatrixRestingTimePerCondition(condition)
            resting_time_count_matrix = np.nan_to_num(resting_time_count_matrix, nan=0.0)
            resting_time_matrix_norm = resting_time_matrix / volume_matrix
            resting_time_matrix_average = np.divide(resting_time_matrix_norm, resting_time_count_matrix,
                                                where=resting_time_count_matrix != 0)
            c = ax[index].pcolormesh(X, Y, resting_time_matrix_average, cmap=custom_cmap)
            ax[index].set_xlabel('r coordinates')
            ax[0].set_ylabel('z coordinates')
            ax[index].set_facecolor('white')
            ax[index].set_title(f'{condition}')

            inlet_r, inlet_z, body_r, body_z = self.getTrap2D()
            ax[index].fill(inlet_r, inlet_z, color='purple')
            ax[index].fill(body_r, body_z, color='purple')
            ax[index].set_xlim(0, 0.3)
            ax[index].set_ylim(-0.45, 0.1)
            ax[index].set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        fig.colorbar(c, ax=ax, orientation = 'vertical', label='Density resting time (s/m$^3$)')
        plt.show()



    def plotHeatmapLandingToCaptureProbability(self): #klopt geen ene kut van nu
        volume_matrix = self.getMatrixNormilizingVolume()
        captured_coord_matrix = self.getMatrixCaptureProbability()
        takeoffpoints_count_matrix = self.getMatrixTakeOffPoints()

        matrix_norm = captured_coord_matrix / volume_matrix
        matrix_average_norm = np.divide(matrix_norm, takeoffpoints_count_matrix,
                                                where=takeoffpoints_count_matrix != 0)

        matrix_average_norm /= matrix_average_norm.max() # the max value stays 1. (so it keeps being a probaility)

        fig = plt.figure()
        ax = fig.add_subplot(132, title='Density Heatmap: Positional probability to get captured after landing', aspect='equal')
        X, Y = np.meshgrid(self.r_edges_matrix, self.z_edges_matrix)

        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)
        c = ax.pcolormesh(X, Y, matrix_average_norm, cmap=custom_cmap)
        ax.set_facecolor('white')

        plt.colorbar(c, ax=ax, fraction=0.065, pad=0.13, label='Probability')
        plt.xlabel('r coordinates')
        plt.ylabel('z coordinates')

        # plot trap
        inlet_r, inlet_z, body_r, body_z = self.getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple')
        plt.fill(body_r, body_z, color='purple')
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

    def plotHeatmapLandingAgainProbability(self): #klopt geen ene kut van nu
        volume_matrix = self.getMatrixNormilizingVolume()
        landing_again_matrix = self.getMatrixLandingAgainProbability()
        takeoffpoints_count_matrix = self.getMatrixTakeOffPoints()

        matrix_norm = landing_again_matrix / volume_matrix
        matrix_average_norm = np.divide(matrix_norm, takeoffpoints_count_matrix,
                                                where=takeoffpoints_count_matrix != 0)

        matrix_average_norm /= matrix_average_norm.max() # the max value stays 1. (so it keeps being a probaility)

        fig = plt.figure()
        ax = fig.add_subplot(132, title='Density Heatmap: Positional probability to land again after landing', aspect='equal')
        X, Y = np.meshgrid(self.r_edges_matrix, self.z_edges_matrix)

        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)
        c = ax.pcolormesh(X, Y, matrix_average_norm, cmap=custom_cmap)
        ax.set_facecolor('white')

        plt.colorbar(c, ax=ax, fraction=0.065, pad=0.13, label='Probability')
        plt.xlabel('r coordinates')
        plt.ylabel('z coordinates')

        # plot trap
        inlet_r, inlet_z, body_r, body_z = self.getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple')
        plt.fill(body_r, body_z, color='purple')
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

