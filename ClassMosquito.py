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



from supportive_functions import *

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import f_oneway, sem, norm
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import pandas as pd
import numpy as np
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

    def getTrack(self): # [x, y, z, time]
        return [self.x, self.y, self.z, self.time]

# Returns hoppings, landing_track, take_off_track
    def splitCoordinatesTrack(self, boundary=0.02):
        x, y, z, t = self.getTrack()
        in_run = False
        x_hop, y_hop, z_hop, t_hop = [], [], [], []
        hoppings = []
        landing_track =[]
        take_off_track = []
        walking_track = []
        for i in range(len(x)):
            if landing_area(x[i], y[i], z[i], boundary=boundary):
                if not in_run:
                    # Start time of a landing
                    in_run = True
                    x_hop = [x[i]]
                    y_hop = [y[i]]
                    z_hop = [z[i]]
                    t_hop = [t[i]]
                else:
                    x_hop.append(x[i])
                    y_hop.append(y[i])
                    z_hop.append(z[i])
                    t_hop.append(t[i])
            else:
                if in_run:
                    # End time of a landing
                    in_run = False
                    hoppings.append([x_hop, y_hop, z_hop, t_hop])
                    x_hop, y_hop, z_hop, t_hop = [], [], [], []
        if in_run:
            # When a track ends in landing, the landing still gets added
            hoppings.append([x_hop, y_hop, z_hop, t_hop])

        if hoppings:
            if len(hoppings) == 1 and landing_area(x[-1], y[-1], z[-1], boundary = boundary) and landing_area(x[0], y[0], z[0], boundary=boundary):
                walking_track = hoppings[0]
                hoppings.pop()
            else:
                if landing_area(x[-1], y[-1], z[-1], boundary = boundary):
                    landing_track = hoppings[-1]
                    hoppings.pop()
                if landing_area(x[0], y[0], z[0], boundary=boundary):
                    take_off_track = hoppings[0]
                    hoppings.pop(0)
        return hoppings, landing_track, take_off_track, walking_track

    def getAllTracksInlandingArea(self, boundary = 0.02):
        all_track = []
        hoppings, landing, take_off, walking_track = self.splitCoordinatesTrack(boundary=boundary)

        if walking_track:
            all_track.append(walking_track)

        if hoppings:
            all_track += hoppings

        if take_off:
            all_track.append(take_off)

        if landing:
            all_track.append(landing)
        return all_track


    # Returns a list of lists of coordinates. A list of all te hoppings.
    # --> [ [ [x], [y], [z], [t] ], [more hoppings], exe... ]
    def getHoppingCoordinatesTrack(self, boundary = 0.02):
        hoppings, landing_track, take_off_track, walking_track = self.splitCoordinatesTrack(boundary=boundary)
        return hoppings

    def getLandingTrack(self, boundary = 0.02):
        hoppings, landing_track, take_off_track, walking_track = self.splitCoordinatesTrack(boundary=boundary)
        return landing_track
    def getWalkingTrack(self, boundary = 0.02):
        walk = []
        hoppings, landing_track, take_off_track, walking_track = self.splitCoordinatesTrack(boundary=boundary)
        if walking_track:
            walk = walking_track
        return walk
    def getTakeOffTrack(self, boundary = 0.02):
        take_off = []
        hoppings, landing_track, take_off_track, walking_track = self.splitCoordinatesTrack(boundary=boundary)
        if take_off_track:
            take_off = take_off_track
        return take_off

# Returns a list of the landing points of this track
    # --> [ [x,y,z,t], exe.. ]
    def getHoppingLandingPointsTrack(self, boundary = 0.02):
        hoppings = self.getHoppingCoordinatesTrack(boundary=boundary)
        landing_points = []
        for hop in hoppings:
            x, y, z, t = hop
            landing_point = nearest_neighbor_to_trap_surface(x, y, z, t)
            landing_points.append(landing_point)
        return landing_points

    def getDiplacementHoppingsTrack(self, boundary = 0.02):
        hoppings = self.getHoppingCoordinatesTrack(boundary=boundary)
        distances = []
        if hoppings:
            for hop in hoppings:
                x, y, z, t = hop
                distance = np.sqrt( (x[0] - x[-1])**2 +
                                    (y[0] - y[-1])**2 +
                                    (z[0] - z[-1])**2
                                    )
                distances.append(distance)
        return distances

    def boolCaptureTrack(self):
        bool = False
        if capturing_area(self.x[-1], self.y[-1], self.z[-1]):
            bool = True
        return bool
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

    def lastCoordinateTrack(self):
        last_x =self.x[-1]
        last_y = self.y[-1]
        last_z = self.z[-1]
        last_time = self.time[-1]
        return [last_x, last_y, last_z, last_time]

    def firstCoorinatesTrack(self):
        return [self.x[0], self.y[0], self.z[0], self.time[0]]

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

        x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet = getTrap()
        ax.plot_surface(x_grid_body, y_grid_body, z_grid_body, alpha=0.5, color='b')
        ax.plot_surface(x_grid_inlet, y_grid_inlet, z_grid_inlet, alpha=0.5, color='b')
        ax.set_aspect('equal', adjustable='box')

        plt.show()

    def plot2DTrack(self):
        r = []
        last_r, last_z = np.sqrt(self.x[-1]**2 + self.y[-1]**2), self.z[-1]
        for i in range(len(self.x)):
            r.append(np.sqrt(self.x[i]**2 + self.y[i]**2))
        plt.plot(r, self.z)
        plt.scatter(last_r, last_z, color='r', marker='o')
        inlet_r, inlet_z, body_r, body_z = getTrap2D()
        # Plot trap
        plt.fill(inlet_r, inlet_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.fill(body_r, body_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.xlim(0, 0.8)
        plt.ylim(-0.5, 0.5)
        plt.xlabel('r')
        plt.ylabel('z')
        plt.title('2D (r, z) plot of a single track')
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

    def plotLandingTrack(self, boundary = 0.02):
        x, y, z, t = self.getTrack()
        ax = plt.axes(projection='3d')
        for i in range(1, len(x)):
            if landing_area(x[i], y[i], z[i], boundary = boundary):
                ax.plot([x[i-1], x[i]], [y[i-1], y[i]], [z[i-1], z[i]], color='r')
            else:
                ax.plot([x[i-1], x[i]], [y[i-1], y[i]], [z[i-1], z[i]], color='g')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('3D (x, y, z) plot of a single track')
        tracks = self.getAllTracksInlandingArea(boundary=boundary)

        for landing_point in tracks:
            x, y, z, t = landing_point
            ax.scatter(x, y, z, color = 'r', marker='o')


        x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet = getTrap()
        ax.plot_surface(x_grid_body, y_grid_body, z_grid_body, alpha=0.5, color='b')
        ax.plot_surface(x_grid_inlet, y_grid_inlet, z_grid_inlet, alpha=0.5, color='b')
        ax.set_aspect('equal', adjustable='box')
        plt.show()

    def plotLanding2DTrack(self, boundary = 0.02):
        x, y, z, t = self.getTrack()
        x, y = np.array(x), np.array(y)
        r = np.sqrt(x ** 2 + y ** 2)
        for i in range(1, len(x)):
            if landing_area(x[i], y[i], z[i], boundary = boundary):
                plt.plot([r[i-1], r[i]], [z[i-1], z[i]], color='r')
            else:
                plt.plot([r[i-1], r[i]], [z[i-1], z[i]], color='g')
        inlet_r, inlet_z, body_r, body_z = getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple', linewidth=0, alpha=0.5)
        plt.fill(body_r, body_z, color='purple', linewidth=0, alpha=0.5)
        plt.xlim(0, 0.8)
        plt.ylim(-0.5, 0.5)
        plt.xlabel('r')
        plt.ylabel('z')
        plt.title('2D (r, z) plot of a single track')
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

    def getRestingTimeTrack(self, boundary = 0.02):
        tracks = self.getAllTracksInlandingArea(boundary=boundary)
        print(tracks)
        resting_times = []
        for i, hop in enumerate(tracks):
            x, y, z, t = hop
            resting_time = t[-1] - t[0]
            resting_times.append(resting_time)
        return resting_times


class Trial:
    def __init__(self, trial_num):  # if you want to load a track, fill in both trial and track, otherwise only whole trial
        self.total_mosquitos_per_trial = 50
        self.trial_num = trial_num
        self.coordinate_list_trial = None
        if trial_num < 59:  # to exclude tria 59
            self.condition = accessing_extra_info('Condition')[trial_num -1]
        else:
            self.condition = accessing_extra_info('Condition')[trial_num]
        self.header = f'Trial_{trial_num + 1}'
        self.hopping_points = None
        self.landing_points = None
        self.take_off_points = None
        self.paired_points = None
        self.num_tracks = None
        self.track_objects = None

        self.hopping_tracks = None
        self.landing_tracks = None
        self.take_off_tracks = None
        self.walking_tracks = None

# All the coordinates of one trial N
    # --> [ [header, [x_coordinates], [y_coordinates], [z_coordinates], [time] ] , exe... ]
    def getTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        return self.coordinate_list_trial

# Creating objects for all the given tracks. N
    def getTrackObjects(self):
        object_array = []
        trial = self.getTrial()
        for track in trial:
            header, x, y, z, t = track
            track_num = int(header.split('_')[-1])
            obj = Track(trial_num=self.trial_num, track_num=track_num)
            object_array.append(obj)
        return object_array

# Initializes self.coordinate_list_trial, aka a list of all the point in a trial N
    # --> [ [header, [x_coordinates], [y_coordinates], [z_coordinates], [time] ], exe... ]
    def initiateCoordinateList(self):
        self.coordinate_list_trial = accessing_trial(self.trial_num)

# First/last coordinates #

# List of last coordinates of a trial N
    # --> [ [ header, [x, y, z, time] ], exe... ]
    def lastCoordinateTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        storage = []
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            last_coordinates = [x[-1], y[-1], z[-1], time[-1]]
            storage.append([header, last_coordinates])
        return storage

# List of first coordinates of a trial N
    # --> [ [ header, [x, y, z, time] ], exe... ]
    def firstCoordinateTrial(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        storage = []
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            first_coordinates = [x[0], y[0], z[0], time[0]]
            storage.append([header, first_coordinates])
        return storage

# Total numer of tracks of a trial N
    # --> amount
    def getNumTracks(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        num_tracks = len(self.coordinate_list_trial)
        return num_tracks

# End time of a trial N
    # --> amount
    def getEndTimeTrial(self):
        last_time_point = 0
        for track in self.getTrial():
            time = track[4]
            end_point = time[-1]
            if end_point > last_time_point:
                last_time_point = end_point
        return last_time_point

# Start time of a trial N
    # --> amount
    def getStartTimeTrial(self):
        first_time_point = 1500
        for track in self.getTrial():
            time = track[4]
            first_point = time[0]
            if first_point < first_time_point:
                first_time_point = first_point
        return first_time_point

# Duration of a trial N
    # --> amount
    def getDurationTrial(self):
        return self.getEndTimeTrial() - self.getStartTimeTrial()



#  Generates the most likely resting pairs in a dictionary, with lists of the resting time


# and the associated resting spots in a trial.
        # --> [paired_tracks], [resting_times], [resting_points]
    def generatePairs(self, radius = 0.02, boundary = 0.02):
        if self.track_objects == None:
            self.track_objects = self.getTrackObjects()

        take_off_tracks = []
        landing_tracks = []
        walking_tracks = []
        paired_tracks = []
        for track_object in self.track_objects:
            landing_tracks.append([track_object.track_num, track_object.getLandingTrack(boundary=boundary)])
            take_off_tracks.append([track_object.track_num, track_object.getTakeOffTrack(boundary=boundary)])
            walking_tracks.append([track_object.track_num, track_object.getWalkingTrack(boundary=boundary)])

        potential_new_landing = []
        potential_new_take_off =[]

        for i_walk, walkings in enumerate(walking_tracks):
            walking_num, walking = walkings
            if walking:
                x_walk, y_walk, z_walk, t_walk = walking
                x_walk_beg, y_walk_beg, z_walk_beg, time_walk_beg = x_walk[0], y_walk[0], z_walk[0], t_walk[0]
                x_walk_end, y_walk_end, z_walk_end, time_walk_end = x_walk[-1], y_walk[-1], z_walk[-1], t_walk[-1]
                for i_land, landings in enumerate(landing_tracks):
                    landing_num, landing = landings
                    if landing:
                        x_land, y_land, z_land, t_land = landing
                        x_land_end, y_land_end, z_land_end, time_land_end = x_land[-1], y_land[-1], z_land[-1], t_land[
                            -1]

                        dx, dy, dz, dtime = (x_walk_beg - x_land_end), (y_walk_beg - y_land_end), (
                                z_walk_beg - z_land_end), (time_walk_beg - time_land_end)
                        distance = np.sqrt((dx ** 2) + (dy ** 2) + (dz ** 2))
                        resting_at_merge = dtime
                        if distance < radius and resting_at_merge > 0:
                            potential_new_landing.append((distance,  i_land, i_walk, landing_num))
                for i_takeoff, take_offs in enumerate(take_off_tracks):
                    take_off_num, take_off = take_offs
                    if take_off:
                        x_take, y_take, z_take, t_take = take_off
                        x_take_beg, y_take_beg, z_take_beg, time_take_beg = x_take[0], y_take[0], z_take[0], t_take[0]
                        dx, dy, dz, dtime = (x_take_beg - x_walk_end), (y_take_beg - y_walk_end), (
                                z_take_beg - z_walk_end), (time_take_beg - time_walk_end)
                        distance = np.sqrt((dx ** 2) + (dy ** 2) + (dz ** 2))
                        resting_at_merge = dtime
                        if distance < radius and resting_at_merge > 0:
                            potential_new_take_off.append((distance, i_takeoff, i_walk, take_off_num))

        potential_new_landing.sort()
        potential_new_take_off.sort()
        used_walking_take_points = set()
        used_walking_landing_points = set()
        used_walking_tracks = set()
        new_walking_tracks = []

        for distance_land, i_land, i_walk_land, landing_num in potential_new_landing:
            for distance_take, i_takeoff, i_walk_take, take_off_num in potential_new_take_off:
                if i_walk_land == i_walk_take:
                    if i_land not in used_walking_landing_points and i_takeoff not in used_walking_take_points:
                        used_walking_landing_points.add(i_land)
                        used_walking_take_points.add(i_walk_take)
                        merged_track = [sum(axes, []) for axes in zip(landing_tracks[i_land][1], walking_tracks[i_walk_land][1], take_off_tracks[i_takeoff][1])]
                        landing_tracks.pop(i_land)
                        take_off_tracks.pop(i_takeoff)
                        walking_tracks.pop(i_walk_land)
                        paired_tracks.append(merged_track)

        for distance, i_land, i_walk, landing_num in potential_new_landing:
            if i_land not in used_walking_landing_points:
                used_walking_landing_points.add(i_land)
                used_walking_tracks.add(i_walk)
                merged_lan_track = [a + b for a, b in zip(landing_tracks[i_land][1], walking_tracks[i_walk][1])]
                landing_tracks.pop(i_land)
                walking_tracks.pop(i_walk)
                landing_tracks.append([landing_num, merged_lan_track])

        for distance, i_takeoff, i_walk, take_off_num in potential_new_take_off:
            if i_takeoff not in used_walking_take_points:
                used_walking_take_points.add(i_takeoff)
                used_walking_tracks.add(i_walk)
                merged_take_track = [a + b for a, b in zip(walking_tracks[i_walk][1], take_off_tracks[i_takeoff][1])]
                take_off_tracks.pop(i_takeoff)
                walking_tracks.pop(i_walk)
                take_off_tracks.append([take_off_num, merged_take_track])

        for track_num, walking in walking_tracks:
            if walking:
                new_walking_tracks.append(walking)

        potential_pairs = []
        for i_land, landings in enumerate(landing_tracks):
            landing_num, landing = landings
            if landing:
                x_land, y_land, z_land, t_land = landing
                x_land_beg, y_land_beg, z_land_beg, time_land_beg = x_land[0], y_land[0], z_land[0], t_land[0]
                x_land_end, y_land_end, z_land_end, time_land_end = x_land[-1], y_land[-1], z_land[-1], t_land[-1]
                for i_takeoff, take_offs in enumerate(take_off_tracks):
                    take_off_num, take_off = take_offs
                    if take_off:
                        x_take, y_take, z_take, t_take = take_off
                        x_take_beg, y_take_beg, z_take_beg, time_take_beg = x_take[0], y_take[0], z_take[0], t_take[0]
                        x_take_end, y_take_end, z_take_end, time_take_end = x_take[-1], y_take[-1], z_take[-1], t_take[-1]
                        dx, dy, dz, dtime = (x_take_beg - x_land_end), (y_take_beg - y_land_end), (
                                    z_take_beg - z_land_end), (time_take_beg - time_land_end)
                        resting_at_merge = dtime
                        total_resting_time = time_take_end - time_land_beg

                        distance = np.sqrt((dx ** 2) + (dy ** 2) + (dz ** 2))
                        if distance < radius and resting_at_merge > 0:
                            mean_x = (x_land_end + x_take_beg) / 2
                            mean_y = (y_land_end + y_take_beg) / 2
                            mean_z = (z_land_end + z_take_beg) / 2
                            resting_point = [mean_x, mean_y, mean_z]
                            potential_pairs.append((distance, total_resting_time, i_land, i_takeoff, resting_point))

        potential_pairs.sort()
        used_takeoff_points = set()
        used_landing_points = set()

        paired_resting_times = []
        paired_resting_points = []
        new_take_off_tracks = []
        new_landings_tracks = []

        for distance, total_resting_time, i_land, i_takeoff, resting_point in potential_pairs:
            if i_takeoff not in used_takeoff_points:
                used_takeoff_points.add(i_takeoff)
                used_landing_points.add(i_land)
                paired_resting_times.append(total_resting_time)
                paired_resting_points.append(resting_point)
                merge_track = [a + b for a, b in zip(landing_tracks[i_land][1], take_off_tracks[i_takeoff][1])]
                paired_tracks.append(merge_track)
        for i, item in enumerate(landing_tracks):
            if i not in used_landing_points and item[1]:
                new_landings_tracks.append(item[1])
        for i, item in enumerate(take_off_tracks):
            if i not in used_takeoff_points and item[1]:
                new_take_off_tracks.append(item[1])
        return paired_tracks, paired_resting_times, paired_resting_points, new_landings_tracks, new_take_off_tracks, new_walking_tracks

# Initializes a list in self.hoppings of all the hoppings in trial N
    # --> [ [ [[x], [y], [z], [t]], [[z], [y], [z], [t]] ... exe ]
    def initializeHoppingCoordinatesTrial(self, boundary=0.02):
        if self.track_objects == None:
            self.track_objects = self.getTrackObjects()
        hopping_trial = []
        for track_object in self.track_objects:
            if track_object.getHoppingCoordinatesTrack(boundary=boundary):
                hopping_trial += track_object.getHoppingCoordinatesTrack(boundary=boundary)
        self.hopping_tracks = hopping_trial

    def getHoppingsTrackTrial(self, boundary = 0.02):
        if not self.hopping_tracks:
            self.initializeHoppingCoordinatesTrial(boundary=boundary)
        return self.hopping_tracks

    def getLandingTracksTrial(self, radius = 0.02, boundary = 0.02):
        pairs, resting_times, resting_points, new_landings_tracks, new_take_off_tracks, new_walking_tracks = self.generatePairs(radius = radius,
                                                                                                            boundary=boundary)
        return new_landings_tracks

    def getTakeOffTracksTrial(self, radius = 0.02, boundary = 0.02):
        pairs, resting_times, resting_points, new_landings_tracks, new_take_off_tracks, new_walking_tracks = self.generatePairs(radius = radius,
                                                                                                            boundary=boundary)
        return new_take_off_tracks

    def getPairedTracksTrial(self, radius = 0.02, boundary = 0.02):
        pairs, resting_times, resting_points, new_landings_tracks, new_take_off_tracks, new_walking_tracks = self.generatePairs(radius,
                                                                                                            boundary=boundary)
        return pairs

    def getWalkingTracksTrial(self, radius = 0.02, boundary = 0.02):
        pairs, resting_times, resting_points, new_landings_tracks, new_take_off_tracks, new_walking_tracks = self.generatePairs(
            radius,
            boundary=boundary)
        return new_walking_tracks

    def initializeLandingTracksTrial(self, radius = 0.02, boundary = 0.02 ):
        self.landing_tracks = self.getLandingTracksTrial(radius=radius, boundary=boundary)

    def initializeTakeOffTracksTrial(self, radius = 0.02, boundary = 0.02):
        self.take_off_tracks = self.getTakeOffTracksTrial(radius=radius, boundary=boundary)

    def initializeWalkingTrackTrial(self, radius = 0.02, boundary = 0.02):
        self.walking_tracks = self.getWalkingTracksTrial(radius=radius, boundary=boundary)
    def initializeHoppingPoints(self, boundary = 0.02):
        if not self.hopping_tracks:
            self.initializeHoppingCoordinatesTrial(boundary = boundary)
        points = []
        for hop in self.hopping_tracks:
            x, y, z, t = hop
            points.append(nearest_neighbor_to_trap_surface(x, y, z, t))
        self.hopping_points = points

    def initializeLandingPoints(self, radius = 0.02, boundary = 0.02):
        if not self.landing_points:
            self.initializeLandingTracksTrial(radius= radius, boundary=boundary)
        landingpoints_list = []
        for track in self.landing_tracks:
            x, y, z, time = track
            landingpoints_list.append(nearest_neighbor_to_trap_surface(x, y, z, time))
        self.landing_points = landingpoints_list

    def initializeTakeOffPoints(self, radius= 0.02, boundary = 0.02):
        if not self.take_off_points:
            self.initializeTakeOffTracksTrial(radius= radius,boundary=boundary)
        landingpoints_list = []
        for track in self.take_off_points:
            x, y, z, time = track
            landingpoints_list.append(nearest_neighbor_to_trap_surface(x, y, z, time))
        self.take_off_points = landingpoints_list

    def initializePairedPoints(self, radius = 0.02, boundary= 0.02):
        self.paired_points = self.getRestingPointsPairsTrial(radius= radius,boundary=boundary)

# Only get the resting time list of a trial
    # --> [resting_times]
    def getRestingTimePairsTrial(self, radius = 0.02, boundary = 0.02):
        pairs, resting_times, resting_points, new_landings_tracks, new_take_off_tracks, new_walking_tracks = self.generatePairs(radius, boundary=boundary)
        return resting_times

    def getRestingTimeTakeOffsTrial(self, radius = 0.02, boundary = 0.02):
        if not self.take_off_tracks:
            self.initializeTakeOffTracksTrial(radius= radius,boundary=boundary)
        resting_times = []
        for take_off in self.take_off_tracks:
            x, y, z, t = take_off
            resting_time = (t[-1] - t[0])
            resting_times.append(resting_time)
        return resting_times

    def getRestingTimeLandingsTrial(self, radius = 0.02, boundary = 0.02):
        if not self.landing_tracks:
            self.initializeLandingTracksTrial(radius= radius,boundary=boundary)
        resting_times = []
        for landing in self.landing_tracks:
            x, y, z, t = landing
            resting_time = (t[-1] - t[0])
            resting_times.append(resting_time)
        return resting_times

    def getRestingTimeHoppingsTrial(self, boundary = 0.02):
        if not self.hopping_tracks:
            self.initializeHoppingCoordinatesTrial(boundary=boundary)
        resting_times = []

        for hop in self.hopping_tracks:
            x, y, z, t = hop
            duration = t[-1] - t[0]
            resting_times.append(duration)
        return resting_times

    def getRestingTimeWalkingsTrial(self, radius = 0.02, boundary = 0.02):
        if not self.walking_tracks:
            self.initializeWalkingTrackTrial(radius=radius,boundary=boundary)
        resting_times = []

        for hop in self.walking_tracks:
            x, y, z, t = hop
            duration = t[-1] - t[0]
            resting_times.append(duration)
        return resting_times

    def getRestingTimeTrial(self, radius=0.02, boundary=0.02):
        hoppings = self.getRestingTimeHoppingsTrial(boundary = boundary)
        landings = self.getRestingTimeLandingsTrial(radius=radius, boundary=boundary)
        take_offs = self.getRestingTimeTakeOffsTrial(radius=radius, boundary=boundary)
        pairs = self.getRestingTimePairsTrial(radius = radius, boundary=boundary)
        walking = self.getRestingTimeWalkingsTrial(radius = radius, boundary=boundary)
        return hoppings + landings + take_offs + pairs + walking

    # Only get the resting points list of a trial
    # --> [resting_points]
    def getRestingPointsPairsTrial(self, radius = 0.02, boundary = 0.02):
        pairs, resting_times, resting_points, new_landings_tracks, new_take_off_tracks, new_walking_tracks = self.generatePairs(radius, boundary=boundary)
        return resting_points


# Get the total number of associated (landing -- take-off) pairs of a trial
    # --> amount
    def countPairsTrial(self, radius = 0.02, boundary = 0.02):
        return len(self.getPairedTracksTrial(radius=radius, boundary=boundary))

    def countWalkingTracksTrial(self, boundary = 0.02):
        if self.track_objects == None:
            self.track_objects = self.getTrackObjects()
        count = 0
        for track_object in self.track_objects:
             if track_object.getWalkingTrack(boundary=boundary):
                 count += 1
        return count


# PLotting resting times #

    def plotRestingTimesViolinTrial(self, radius = 0.02, boundary = 0.02):
        hoppings = self.getRestingTimeHoppingsTrial(boundary=boundary)
        landings = self.getRestingTimeLandingsTrial(radius=radius, boundary=boundary)
        take_offs = self.getRestingTimeTakeOffsTrial(radius=radius, boundary=boundary)
        pairs = self.getRestingTimePairsTrial(radius=radius, boundary=boundary)
        walkings = self.getRestingTimeWalkingsTrial(radius=radius, boundary=boundary)
        all_resting_times = self.getRestingTimeTrial(radius=radius, boundary=boundary)
        data = [hoppings, landings, take_offs, walkings, pairs, all_resting_times]
        titles = ['Hoppings', 'Landings', 'Take-offs', 'walkings', 'Resting pairs', 'All resting times']
        fig, axs = plt.subplots(nrows=1, ncols=len(data))
        for i, ax in enumerate(axs):
            vp = ax.violinplot([data[i]], showmeans=True)
            ax.set_title(titles[i])
            ax.set_ylabel('Resting Times in Seconds (s)')
            ax.set_xticks([])
            count = len(data[i])
            y_min, y_max = ax.get_ylim()
            y_pos = y_min - 0.05 * (y_max - y_min)
            ax.text(0.95, y_pos, f'n={count}', ha='center', va='top', fontsize=9)

        plt.suptitle(f'Resting Times for Trial {self.trial_num}', fontsize = 14)
        plt.tight_layout()
        plt.show()


    def plotDisplacementViolin(self, radius = 0.02, boundary=0.02):

        return








    # Landing and Capture rates #

# Get total of the catches for a trial measured by Cribellier et al. (2020)
    # --> amount
    def getMeasuredCatchesTrial(self):
        return int(accessing_extra_info("Catches")[self.trial_num - 1])

# Get total of the catches of a trial measured by this program. If the track ends up in the defined capture area it is cached.
    # --> amount
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

# Get total of the landings of a trial measured by this program. If the track ends up in the defined landing area it is a landing.
    # --> amount
    def countLandingsTrial(self, boundary = 0.02):
        if self.hoppings == None:
            self.initializeHoppingCoordinatesTrial()

        return len(self.hoppings_track(boundary=boundary))

# Get total of the take-offs of a trial measured by this program. If the track begins up in the defined landing area it is a take-off.
    # --> amount
    def countTakeOffsTrial(self, boundary =0.02):
        return len(self.getTakeOffTracksTrial(boundary=boundary))

# What happens after take-off analysis #

# Get total amount of tracks in a trial that begin in take-off and end in capture
    # --> amount
    def countLandingToCaptureTrial(self, boundary = 0.02):
        if self.track_objects == None:
            self.track_objects = self.getTrackObjects()
        count = 0
        for track_object in self.track_objects:
            num_hop = len(track_object.getAllTracksInlandingArea(boundary=boundary))
            if num_hop > 1:
                if track_object.boolCaptureTrack():
                    count += 1
        return count

# Get total amount of tracks in a trial that begin and end in landing
# Landing again
    # --> amount
    def countLandingAgainTrial(self, boundary =0.02):
        if self.track_objects == None:
            self.track_objects = self.getTrackObjects()
        count = 0
        for track_object in self.track_objects:
            num_hop = len(track_object.getAllTracksInlandingArea(boundary=boundary))
            if num_hop > 1:
                land_again = num_hop - 1
                count += land_again
        return count

# Get list of coordinates that begin in take-ff and end in capture of a trial
    # --> [[x, y, z], exe... ]
    def getCoordinatesLandingToCapture(self, boundary = 0.02):
        if self.track_objects == None:
            self.track_objects = self.getTrackObjects()
        take_off_coordinates = []
        for track_object in self.track_objects:
            if capturing_area(track_object.x[-1], track_object.y[-1], track_object.z[-1], boundary=boundary):
                x, y, z, t = track_object.getHoppingLandingPointsTrack(boundary=boundary)[-1]
                take_off_coordinates.append([x, y, z])
        return take_off_coordinates

# Get list of coordinates that begin and end in landing of a trial ! NOT UPDATED ±
# Landing again
    # --> [[x, y, z], exe... ]
    def getCoordinatesLandingAgain(self, boundary = 0.03):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        take_off_coordinates = []
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            if landing_area(x[0], y[0], z[0], boundary = boundary) == True and landing_area(x[-1], y[-1], z[-1], boundary = boundary) == True:
                take_off_coordinates.append([x[0], y[0], z[0]])
        return take_off_coordinates

# Get a list of the whole tracks that land again
    # --> [[ [x], [y], [z] ], exe... ]
    def getTracksLandingAgain(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        landing_again_tracks = []
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            if landing_area(x[0], y[0], z[0]) == True and landing_area(x[-1], y[-1], z[-1]) == True:
                landing_again_tracks.append([x, y, z])
        return landing_again_tracks

# Get a list of the whole tracks that begin in take-off and end in capture
    # --> [[ [x], [y], [z] ], exe... ]
    def getTracksLandingtoCapture(self):
        if self.coordinate_list_trial == None:
            self.initiateCoordinateList()
        landing_to_capture_tracks = []
        for track in self.coordinate_list_trial:
            header, x, y, z, time = track
            if landing_area(x[0], y[0], z[0]) == True and capturing_area(x[-1], y[-1], z[-1]) == True:
                landing_to_capture_tracks.append([x, y, z])
        return landing_to_capture_tracks


# Plotting  #

# Plot all the tracks of a trial in one 3D figure
    # --> 3D plot
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

        # plot trap
        x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet = getTrap()
        ax.plot_surface(x_grid_body, y_grid_body, z_grid_body, alpha=0.5, color='b')
        ax.plot_surface(x_grid_inlet, y_grid_inlet, z_grid_inlet, alpha=0.5, color='b')
        ax.set_aspect('equal', adjustable='box')
        plt.show()

# Plot all the landing points of a trial in one 3D figure
    # --> 3D plot
    def plotLandingPointsTrial(self):
        ax = plt.figure().add_subplot(projection='3d')
        for track in self.lastCoordinateTrial():
            last_x, last_y, last_z, last_time = track[1]
            if landing_area(last_x, last_y, last_z) == True:
                ax.scatter(last_x, last_y, last_z, color='r', marker='o')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet = getTrap()
        ax.plot_surface(x_grid_body, y_grid_body, z_grid_body, alpha=0.5, color='b')
        ax.plot_surface(x_grid_inlet, y_grid_inlet, z_grid_inlet, alpha=0.5, color='b')

        plt.show()

# Plot all the tracks of a trial in one 2D figure
    # --> 2D plot
    def plot2DTrial(self):
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
        inlet_r, inlet_z, body_r, body_z = getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.fill(body_r, body_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.xlim(0, 0.8)
        plt.ylim(-0.5, 0.5)
        plt.xlabel('r')
        plt.ylabel('z')
        plt.title('b) 2D (r, z) plot of all the tracks in one trial')
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

# Plot all the landing points of a trial in one 2D figure
    # --> 2D plot
    def plot2DLandingPointsTrial(self):
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
        inlet_r, inlet_z, body_r, body_z = getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.fill(body_r, body_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

#PLot a single track that begins and ends in landing
    # --> 3D plot
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
        x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet = getTrap()
        ax.plot_surface(x_grid_body, y_grid_body, z_grid_body, alpha=0.5, color='b')
        ax.plot_surface(x_grid_inlet, y_grid_inlet, z_grid_inlet, alpha=0.5, color='b')
        ax.set_aspect('equal', adjustable='box')
        plt.show()

# PLot a single track that begins in take-off and ends in capture
    # --> 3D plot
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
        x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet = getTrap()
        ax.plot_surface(x_grid_body, y_grid_body, z_grid_body, alpha=0.5, color='b')
        ax.plot_surface(x_grid_inlet, y_grid_inlet, z_grid_inlet, alpha=0.5, color='b')
        ax.set_aspect('equal', adjustable='box')
        plt.show()

class Dataset:

    def __init__(self):  # if you want to load a track, fill in both trial and track, otherwise only whole trial
        self.total_mosquitos_per_trial = 50
        self.total_number_trials = 64

        self.coordinate_list_dataset = None
        self.trialobjects = None

        self.catches_without = None
        self.catches_with_heat = None
        self.catches_with_heat_water = None

        self.measured_catches_without = None
        self.measured_catches_with_heat = None
        self.measured_catches_with_heat_water = None

        self.landing_without = None
        self.landing_with_heat = None
        self.landing_with_heat_water = None

        self.landing_all_trials = None
        self.catches_all_trials = None
        self.measured_catches_all_trials = None


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

# Initializing #

# Initialize self.coordinate_list_dataset
    # --> [ [ [header, x_coordinates, y_coordinates, z_coordinates, time_coordinates] , exe ... ] ]
    def initializeCoordinateList(self):
        self.coordinate_list_dataset = []
        for i in range(1, self.total_number_trials + 1):
            coordinate_trial = accessing_trial(i)
            self.coordinate_list_dataset.append(coordinate_trial)

# Initialize all the catches measured by Cribellier et al. (2020) in a lists per short-range condition.   !NOT IN USE!
    # --> 3 lists: self.measured_catches_without, self.measured_catches_with_heat, self.measured_catches_with_heat_water
    # ---> [ num, num, num ]
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

# !NOT IN USE!
# Initialize all the catches measured by Cribellier et al. (2020) in a lists: self.measured_catches_all_trials
    # ---> [ num, num, num ]
    def initializeMeasuredCatchesAllTrials(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        catches_count_per_trial = []
        for trial_object in self.getTrialObjects():
            catches_count_per_trial.append(trial_object.getMeasuredCatchesTrial())
        self.measured_catches_all_trials = catches_count_per_trial


# Initialize all the simulated (by this program) catches per short-range cue condition
    # --> 3 lists: self.catches_without, self.catches_with_heat, self.catches_with_heat_water
    # ---> [ num, num, num ]
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

# Initialize all the simulated (by this program) catches in a list: self.catches_all_trials
    # ---> [ num, num, num ]
    def initializeSimulatedCatchesAllTrials(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        catches_count_per_trial = []
        for trial_object in self.getTrialObjects():
            catches_count_per_trial.append(trial_object.countSimulatedCatchesTrial())
        self.catches_all_trials = catches_count_per_trial

# Initialize all the simulated (by this program) landings per short-range cue condition
    # --> 3 lists: self.landing_without, self.landing_with_heat, self.landing_with_heat_water
    # ---> [ num, num, num ]
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

# Initialize all the simulated (by this program) landings in a list: self.landings_all_trials
    # ---> [ num, num, num ]
    def initializeLandingAllTrials(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        landing_count_per_trial = []
        for trial_object in self.getTrialObjects():
            landing_count_per_trial.append(trial_object.countLandingsTrial())
        self.landing_all_trials = landing_count_per_trial

# Get the trial information as objects in a list
    # --> list of all the trials (inc their information)
    def getTrialObjects(self):
        object_array = []
        for i in range(1, self.total_number_trials + 1):
            obj = Trial(trial_num=i)
            object_array.append(obj)
        return object_array


# Duration #


# Get start times of all the trials
    # --> list
    def getStartTimes(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        beginning_points = []
        for trial_object in self.trialobjects:
            if trial_object.trial_num != 59:
                beginning_points.append(trial_object.getStartTimeTrial())
        return beginning_points

# Get end times of all the trials
    # ---> list
    def getEndTimes(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        end_points = []
        for trial_object in self.trialobjects:
            if trial_object.trial_num != 59:
                end_points.append(trial_object.getEndTimeTrial())
        return end_points

# Get number of tracks per trial
    # --> list
    def getNumTracksPerTrial(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        num_tracks_per_trial = []
        for trial_object in self.trialobjects:
            num_tracks_per_trial.append(trial_object.getNumTracks())
        return num_tracks_per_trial

# Get average duration of all the tracks per trial
    # --> list (length 64)
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

    def countWalkingTracks(self, boundary=0.02):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        count = 0
        for trial_object in self.trialobjects:
            count += trial_object.countWalkingTracksTrial(boundary=boundary)
        return count
# Landing -- take_off #

# Get landing coordinates in 2D in two lists: r and z
    # ---> list r , list z
    def getlandingPointsTheta(self, boundary = 0.03):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list = []
        z_list = []
        for trial_object in self.trialobjects:
            trial_object.initializeLandingPoints(boundary=boundary)
            for coordinate in trial_object.landingpoints:
                x, y, z, time = coordinate
                r = np.sqrt(x ** 2 + y ** 2)
                r_list.append(r)
                z_list.append(z)
        return r_list, z_list

# Get take off coordinates in 2D in two lists: r and z
    # ---> list r , list z
    def getTakeOffPointsTheta(self, boundary = 0.03):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list = []
        z_list = []
        for trial_object in self.trialobjects:
            trial_object.initializeTakeOffPoints(boundary= boundary)
            for coordinate in trial_object.take_off_points:
                x, y, z, time = coordinate
                r = np.sqrt(x ** 2 + y ** 2)
                r_list.append(r)
                z_list.append(z)
        return r_list, z_list

# Get landing coordinates in 2D, specify the condition
    # --> list r, list z
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


# Resting time #

# Get average resting time in a string
    # --> mean +- standard deviation
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

# Get the median resting time
    # --> 1st quartile, median, 3rd quartile
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

# Get the longest resting time in the whole dataset
    # --> num / amount
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

# Count the num of resting times below the boundary given
    # --> num / amount
    def countSmallRestingTimes(self, boundary_seconds):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_list = []
        for trial_object in self.trialobjects:
            resting_time_trial = trial_object.getRestingTimeTrial()
            for point in resting_time_trial:
                if point < boundary_seconds:
                    resting_time_list.append(point)
        return len(resting_time_list)

# Count the total resting times in the whole dataset
    # --> num / amount
    def countTotalRestingTimes(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        resting_time_list = []
        for trial_object in self.trialobjects:
            resting_time_trial = trial_object.getRestingTimeTrial()
            for point in resting_time_trial:
                    resting_time_list.append(point)
        return len(resting_time_list)

# Count the total landing point in the whole dataset
    # --> num / amount
    def countLandingPoints(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        landing_points = 0
        for trial_object in self.trialobjects:
            landing_points += trial_object.countLandingsTrial()
        return landing_points

# Count the total captured mosquitoes in the whole dataset
    # --> num / amount
    def countCapturingPoints(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        capturing_points = 0
        for trial_object in self.trialobjects:
            capturing_points += trial_object.countSimulatedCatchesTrial()
        return capturing_points

# Get a list of all the last coordinates in the dataset  !NOT IN USE!
    # -->  [ [ [ [header], [x_coordinate,  y_coordinate, z_coordinate, time_point] ], exe...] ]
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
        return dataset_coordinates


# Statistics #

# Does the anova test between the three short-range cue conditions
    # Input:  [ [without], [with_heat], [with_heat_water] ]
    # Output: F_stat, p_value
    def testAnovaConditions(self, data):
        without, with_heat, with_heat_water = data
        f_stat, p_value = f_oneway(without, with_heat, with_heat_water)
        return f_stat, p_value

# Does the TurkeysHSD test between the three short-range cue conditions
    # Input:  [ [without], [with_heat], [with_heat_water] ]
    # Output: Significant pairs
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

# Plots the normal distribution and a Q-Q plot of the given data (list)
    # Input: list
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

# Gives the confidence interval of the given data (list)
    # Input: list
    # Output: mean +- confidence interval (string)
    def testConfidenceInterval(self, data):
        mean = np.mean(data)
        se = sem(data)
        ci = 1.96 * se # for 95% confidence interval
        return f'{mean}'+u"\u00B1"+f'{ci}'


# Landing again / landing --> capture #

# Gives the percentage of landings after take-off (from the total take-offs)
    # --> num / amount
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

# Gives the percentage of captures after take-off (from the total take-offs)
    # --> num / amount
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


# Plotting durations #

# Plots the durations of all the trials in a plot to visualize te dataset
    # --> plot with vertical lines, x = time, y = trials
    def plotDuration(self):
        fig, ax = plt.subplots()
        for i, (start, end) in enumerate(zip(self.getStartTimes(), self.getEndTimes())):
            ax.plot([start, end], [i, i], marker='o', color='teal')
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Trial number")
        ax.set_title("Trial Durations")
        plt.show()

# Plots start and end times in a violin plot
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

# PLots a boxplot of the amount the catches per condition
    # Prints significance in the terminal
    def plotBoxplotCatchesConditions(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        if self.catches_without == None:
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

# Plots boxplot of the amount of landings per condition
    # Prints significance in the terminal
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

# Plots boxplot of the amount of landings per condition
    # Prints significance in the terminal
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


# Plotting resting times #

# PLots histogram of resting times of the whole dataset
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

# Plots a histogram of the resting times between a given time boundary (in seconds)
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

# PLot a histogram of the resting time per given condition
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

# PLot a boxplot of the resting times differentiating between the three short-range cue conditions
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
        # anova
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

# PLot a violin plot of the resting times differentiating between the three short-range cue conditions
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



# creating matrices #

# Get the volume matrix, that is used to normalize a matrix by volume.
    # --> matrix
    def getMatrixNormilizingVolume(self):
        distance_row = np.linspace(0.005, self.last_r_coord_matrix - 0.005, self.num_r_cells_matrix)
        distance_matrix = np.tile(distance_row, (self.num_z_cells_matrix, 1))
        volume_matrix = distance_matrix * 2 * np.pi * self.area_cell_matrix
        volume_matrix_np = np.array(volume_matrix)
        return volume_matrix_np

# Get the matrix of the landing points
    # --> matrix
    def getMatrixLandingPoints(self, boundary = 0.03):
        r, z = self.getlandingPointsTheta(boundary = boundary)
        landingpoint_count_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r, z, bins=(
            self.r_edges_matrix, self.z_edges_matrix))
        landingpoint_count_matrix = landingpoint_count_matrix.T
        landingpoint_count_matrix_np = np.array(landingpoint_count_matrix)
        return landingpoint_count_matrix_np

# Get the matrix of the take-off points
    # --> matrix
    def getMatrixTakeOffPoints(self, boundary = 0.03):
        r, z = self.getTakeOffPointsTheta(boundary = boundary)
        takeoffpoints_count_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r, z, bins=(
            self.r_edges_matrix, self.z_edges_matrix))
        takeoffpoints_count_matrix = takeoffpoints_count_matrix.T
        takeoffpoints_count_matrix_np = np.array(takeoffpoints_count_matrix)
        return takeoffpoints_count_matrix_np

# Get the matrix of the landing points of a given short-range cue condition
    # --> matrix
    def getMatrixLandingPointsPerCondition(self, condition):
        r, z = self.getlandingPointsThetaPerCondition(condition)
        landingpoint_count_matrix, r_edges_hist, z_edges_hist = np.histogram2d(r, z, bins=(
            self.r_edges_matrix, self.z_edges_matrix))
        landingpoint_count_matrix = landingpoint_count_matrix.T
        landingpoint_count_matrix_np = np.array(landingpoint_count_matrix)
        return landingpoint_count_matrix_np

# Get the matrix of the resting times, with a possibility to change the upper/lower time boundary
    # --> 2 matrices: resting_time_matrix_np, resting_time_count_matrix_np
    def getMatrixRestingTimes(self, lower_time_boundary = 0, upper_time_boundary = 1500, boundary = 0.03):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list, z_list, w_list = [], [], []
        for trial_object in self.trialobjects:
            resting_times = trial_object.getRestingTimeTrial(boundary=boundary)
            points = trial_object.getRestingPointsTrial(boundary=boundary)
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

# get matrix of the resting time of a given condition
    # --> 2 matrices: resting_time_matrix_np, resting_time_count_matrix_np
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

# Get the matrix of the probability that a mosquito that takes-off in a certain space ends in capture.
    # --> matrix
    def getMatrixCaptureProbability(self, boundary = 0.03):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list, z_list, w_list = [], [], []
        for trial_object in self.trialobjects:
            coordinates = trial_object.getCoordinatesLandingToCapture(boundary = boundary)
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

# Get the matrix of the probability that a mosquito that takes-off in a certain space lands again.
    # --> matrix
    def getMatrixLandingAgainProbability(self, boundary = 0.03):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        r_list, z_list, w_list = [], [], []
        for trial_object in self.trialobjects:
            coordinates = trial_object.getCoordinatesLandingAgain(boundary=boundary)
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


# 2D heatmaps #

# Plot heatmap of all the landing points (not normalized by volume)
    # --> heatmap
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
        inlet_r, inlet_z, body_r, body_z = getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple', edgecolor = 'none', alpha = 0.5)
        plt.fill(body_r, body_z, color='purple', edge_color = 'none', alpha = 0.5)
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

# Plot heatmap of all the landing points, normalized by volume!
    # --> heatmap
    def plotHeatmapLandingPointNormilized(self, boundary = 0.03):
        volume_matrix = self.getMatrixNormilizingVolume()
        landingpoint_matrix = self.getMatrixLandingPoints(boundary = boundary)

        density_matrix = landingpoint_matrix / volume_matrix

        # plot matrix
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
        inlet_r, inlet_z, body_r, body_z = getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.fill(body_r, body_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

# Plots heatmap of all the resting points
    # --> heatmap
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
        inlet_r, inlet_z, body_r, body_z = getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.fill(body_r, body_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()  # nor #ormil

# Plots heatmap of all the resting times
    # --> heatmap
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
        inlet_r, inlet_z, body_r, body_z = getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.fill(body_r, body_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

# Plots heatmap of the resting times per short-range cue condition
    # --> heatmap
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

            inlet_r, inlet_z, body_r, body_z = getTrap2D()
            ax[index].fill(inlet_r, inlet_z, color='purple', linewidth = 0, alpha = 0.5)
            ax[index].fill(body_r, body_z, color='purple', linewidth = 0, alpha = 0.5)
            ax[index].set_xlim(0, 0.3)
            ax[index].set_ylim(-0.45, 0.1)
            ax[index].set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        fig.colorbar(c, ax=ax, orientation = 'vertical', label='Density resting time (s/m$^3$)')
        plt.show()


#  PLot heatmap of the probability that if a mosquito that takes-off in a spot it ends in capture
    # --> heatmap
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
        inlet_r, inlet_z, body_r, body_z = getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.fill(body_r, body_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()

#  PLot heatmap of the probability that if a mosquito that takes-off in a spot it lands again
    # --> heatmap
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
        inlet_r, inlet_z, body_r, body_z = getTrap2D()
        plt.fill(inlet_r, inlet_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.fill(body_r, body_z, color='purple', linewidth = 0, alpha = 0.5)
        plt.xlim(0, 0.3)
        plt.ylim(-0.45, 0.1)
        plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
        plt.show()



# Sensitivity analysis / association analysis #

# PLot a boxplot differentiating between different possible radius and the number of landing--take-off associations
    def plotBoxplotRadiusAssociations(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        pairs1 = []
        pairs2 = []
        pairs3 = []
        pairs4 = []
        for trial_object in self.trialobjects:
            pairs1.append(trial_object.countPairsTrial(0.01))
            pairs2.append(trial_object.countPairsTrial(0.02))
            pairs3.append(trial_object.countPairsTrial(0.03))
            pairs4.append(trial_object.countPairsTrial(0.04))
        data = [pairs1, pairs2, pairs3, pairs4]
        plt.boxplot(data)
        plt.title('Nr of of associated landings with take offs per radius threshold')
        plt.xticks([1, 2, 3, 4], ['r = 0.01', 'r = 0.02', 'r = 0.03', 'r = 0.04'])
        plt.ylabel('nr of associated landings with take offs')
        plt.show()

# PLot boxplot showing the differences in amount between the landings, take-offs and landing--take-off associated pairs.
    def plotBoxplotCountLandingTakeOffAccosiation(self):
        if self.trialobjects == None:
            self.trialobjects = self.getTrialObjects()
        landings, takeoffs, associated_pairs = [], [], []
        for trial_object in self.trialobjects:
            landings.append(trial_object.countLandingsTrial())
            takeoffs.append(trial_object.countTakeOffsTrial())
            associated_pairs.append(trial_object.countPairsTrial())
        data = [landings, takeoffs, associated_pairs]
        plt.boxplot(data)
        plt.title('Boxplot of the number of landings, take-offs and landing / take_off pairs')
        plt.xticks([1, 2, 3], ['Landings', 'Take-offs', 'associated pairs'])
        plt.ylabel('nr of associated landings with take offs')
        plt.show()

# PLot sub heatmaps with different boundary options to see if the landing point results change
    # --> 4 subplots from 0.01 to 0.04 m width
    def plotHeatmapLandingPointsBoundaryAssociationTest(self):
        volume_matrix = self.getMatrixNormilizingVolume()

        boundaries = [0.01, 0.02, 0.03, 0.04]
        titles = ['1 cm', '2 cm', '3 cm', '4 cm']

        # plot matrices
        fig, axs = plt.subplots(2, 2)
        axs = axs.flatten()
        X, Y = np.meshgrid(self.r_edges_matrix, self.z_edges_matrix)
        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)

        for i, ax in enumerate(axs):
            density_matrix = self.getMatrixLandingPoints(boundary=boundaries[i]) / volume_matrix
            heatmap = ax.pcolormesh(X, Y, density_matrix, cmap=custom_cmap)
            ax.set_title(titles[i])
            ax.set_xlabel('r')
            ax.set_ylabel('z')
            ax.set_facecolor('white')
            inlet_r, inlet_z, body_r, body_z = getTrap2D()
            fig.colorbar(heatmap, ax=ax, fraction=0.065, pad=0.13, label='Density (points/m$^3$)')
            ax.fill(inlet_r, inlet_z, color='purple', linewidth=0, alpha=0.5)
            ax.fill(body_r, body_z, color='purple', linewidth=0, alpha=0.5)
            ax.set_xlim(0, 0.3)
            ax.set_ylim(-0.45, 0.1)
            ax.set_aspect('equal', adjustable='box')

        fig.suptitle('Heatmaps with different width boundary around the trap \nLanding points per volume')
        fig.tight_layout()
        plt.show()


# PLot sub heatmaps with different boundary options to see if the resting points result change
    # --> 4 subplots from 0.01 to 0.04 m width
    def plotHeatmapRestingPointsBoundaryAssociationTest(self, lower_time_boundary = 0, upper_time_boundary = 1500):
        volume_matrix = self.getMatrixNormilizingVolume()

        boundaries = [0.01, 0.02, 0.03, 0.04]
        titles = ['1 cm', '2 cm', '3 cm', '4 cm']

        # plot matrices
        fig, axs = plt.subplots(2, 2)
        axs = axs.flatten()
        X, Y = np.meshgrid(self.r_edges_matrix, self.z_edges_matrix)
        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)

        for i, ax in enumerate(axs):
            resting_time_matrix_np, resting_time_count_matrix_np = self.getMatrixRestingTimes(lower_time_boundary,upper_time_boundary, boundary=boundaries[i])
            density_matrix =  resting_time_count_matrix_np/ volume_matrix
            heatmap = ax.pcolormesh(X, Y, density_matrix, cmap=custom_cmap)
            ax.set_title(titles[i])
            ax.set_xlabel('r')
            ax.set_ylabel('z')
            ax.set_facecolor('white')
            inlet_r, inlet_z, body_r, body_z = getTrap2D()
            fig.colorbar(heatmap, ax=ax, fraction=0.065, pad=0.13, label='Density (points/m$^3$)')
            ax.fill(inlet_r, inlet_z, color='purple', linewidth=0, alpha=0.5)
            ax.fill(body_r, body_z, color='purple', linewidth=0, alpha=0.5)
            ax.set_xlim(0, 0.3)
            ax.set_ylim(-0.45, 0.1)
            ax.set_aspect('equal', adjustable='box')

        fig.suptitle('Heatmaps with different width boundary around the trap \nResting points per volume')
        fig.tight_layout()
        plt.show()

# PLot sub heatmaps with different area_boundary options to see if the resting time result change

    def plotHeatmapRestingTimesBoundaryAssociationTest(self, lower_time_boundary = 0, upper_time_boundary = 1500):
        volume_matrix = self.getMatrixNormilizingVolume()

        boundaries = [0.01, 0.02, 0.03, 0.04]
        titles = ['1 cm', '2 cm', '3 cm', '4 cm']

        # plot matrices
        fig, axs = plt.subplots(2, 2)
        axs = axs.flatten()
        X, Y = np.meshgrid(self.r_edges_matrix, self.z_edges_matrix)
        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)

        for i, ax in enumerate(axs):
            resting_time_matrix_np, resting_time_count_matrix_np = self.getMatrixRestingTimes(lower_time_boundary,upper_time_boundary,boundary=boundaries[i])
            resting_time_count_matrix = np.nan_to_num(resting_time_count_matrix_np, nan=0.0)
            resting_time_matrix_norm = resting_time_matrix_np / volume_matrix
            resting_time_matrix_average = np.divide(resting_time_matrix_norm, resting_time_count_matrix,
                                                    where=resting_time_count_matrix != 0)
            heatmap = ax.pcolormesh(X, Y, resting_time_matrix_average, cmap=custom_cmap)
            ax.set_title(titles[i])
            ax.set_xlabel('r')
            ax.set_ylabel('z')
            ax.set_facecolor('white')
            inlet_r, inlet_z, body_r, body_z = getTrap2D()
            fig.colorbar(heatmap, ax=ax, fraction=0.065, pad=0.13, label='Density (points/m$^3$)')
            ax.fill(inlet_r, inlet_z, color='purple', linewidth=0, alpha=0.5)
            ax.fill(body_r, body_z, color='purple', linewidth=0, alpha=0.5)
            ax.set_xlim(0, 0.3)
            ax.set_ylim(-0.45, 0.1)
            ax.set_aspect('equal', adjustable='box')

        fig.suptitle('Heatmaps with different width boundary around the trap \nResting times per volume')
        fig.tight_layout()
        plt.show()

# PLot sub heatmaps with different area_boundary options to see if the probability to land again results change
    def plotHeatmapLandingToCaptureProbabilityAssociationTest(self):
        volume_matrix = self.getMatrixNormilizingVolume()
        boundaries = [0.01, 0.02, 0.03, 0.04]
        titles = ['1 cm', '2 cm', '3 cm', '4 cm']

        # plot matrices
        fig, axs = plt.subplots(2, 2)
        axs = axs.flatten()
        X, Y = np.meshgrid(self.r_edges_matrix, self.z_edges_matrix)
        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)

        for i, ax in enumerate(axs):
            captured_coord_matrix = self.getMatrixCaptureProbability(boundary=boundaries[i])
            takeoffpoints_count_matrix = self.getMatrixTakeOffPoints(boundary=boundaries[i])

            matrix_norm = captured_coord_matrix / volume_matrix
            matrix_average_norm = np.divide(matrix_norm, takeoffpoints_count_matrix,
                                            where=takeoffpoints_count_matrix != 0)
            matrix_average_norm /= matrix_average_norm.max()  # the max value stays 1. (so it keeps being a probability)
            heatmap = ax.pcolormesh(X, Y, matrix_average_norm, cmap=custom_cmap)
            ax.set_title(titles[i])
            ax.set_xlabel('r')
            ax.set_ylabel('z')
            ax.set_facecolor('white')
            inlet_r, inlet_z, body_r, body_z = getTrap2D()
            fig.colorbar(heatmap, ax=ax, fraction=0.065, pad=0.13, label='Probability')
            ax.fill(inlet_r, inlet_z, color='purple', linewidth=0, alpha=0.5)
            ax.fill(body_r, body_z, color='purple', linewidth=0, alpha=0.5)
            ax.set_xlim(0, 0.3)
            ax.set_ylim(-0.45, 0.1)
            ax.set_aspect('equal', adjustable='box')

        fig.suptitle('Heatmaps with different width boundary around the trap \nProbability of a take-off to end in capture')
        fig.tight_layout()
        plt.show()


# PLot sub heatmaps with different area_boundary options to see if the probability to get captured after take-off change
    def plotHeatmapLandingAgainProbabilityAssociationTest(self):
        volume_matrix = self.getMatrixNormilizingVolume()
        boundaries = [0.01, 0.02, 0.03, 0.04]
        titles = ['1 cm', '2 cm', '3 cm', '4 cm']

        # plot matrices
        fig, axs = plt.subplots(2, 2)
        axs = axs.flatten()
        X, Y = np.meshgrid(self.r_edges_matrix, self.z_edges_matrix)
        colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]  # White -> Yellow -> Red -> Dark Red
        custom_cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)
        for i, ax in enumerate(axs):
            landing_again_matrix = self.getMatrixLandingAgainProbability(boundary=boundaries[i])
            takeoffpoints_count_matrix = self.getMatrixTakeOffPoints(boundary=boundaries[i])
            matrix_norm = landing_again_matrix / volume_matrix
            matrix_average_norm = np.divide(matrix_norm, takeoffpoints_count_matrix,
                                            where=takeoffpoints_count_matrix != 0)
            matrix_average_norm /= matrix_average_norm.max()  # the max value stays 1. (so it keeps being a probaility)
            heatmap = ax.pcolormesh(X, Y, matrix_average_norm, cmap=custom_cmap)
            ax.set_title(titles[i])
            ax.set_xlabel('r')
            ax.set_ylabel('z')
            ax.set_facecolor('white')
            inlet_r, inlet_z, body_r, body_z = getTrap2D()
            fig.colorbar(heatmap, ax=ax, fraction=0.065, pad=0.13, label='Probability')
            ax.fill(inlet_r, inlet_z, color='purple', linewidth=0, alpha=0.5)
            ax.fill(body_r, body_z, color='purple', linewidth=0, alpha=0.5)
            ax.set_xlim(0, 0.3)
            ax.set_ylim(-0.45, 0.1)
            ax.set_aspect('equal', adjustable='box')

        fig.suptitle(
            'Heatmaps with different width boundary around the trap \nProbability of a take-off to land again')
        fig.tight_layout()
        plt.show()