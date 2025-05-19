# In this file I'm going to generate a new csv datastructure to generate the pairs and link the tracks together.

trial_num = 1

import ClassMosquito
from loading_matlab_file import *
from supportive_functions import *



def makePairedCSVDataset(path_csv_folder, radius = 0.02, boundary = 0.02):
    basemap_paired_path = os.path.join(path_csv_folder, 'paired_database_csv')
    os.makedirs(basemap_paired_path, exist_ok=True)
    basemap_boundary_tracks_path = os.path.join(path_csv_folder, 'boundary_tracks_csv')
    os.makedirs(basemap_boundary_tracks_path, exist_ok=True)

    for trial_num in range(1, 64):
        trial = ClassMosquito.Trial(trial_num)
        paired_tracks, new_landings_tracks, new_take_off_tracks, new_walking_tracks, stitch_num_land_take, stitch_num_land_walk_take, stitch_num_land_walk, stitch_num_walk_take, altered_track_nums = trial.generatePairsForCSV(radius=radius, boundary=boundary)

        new_trial_map = os.path.join(basemap_paired_path, f'Trial_{trial_num}')
        os.makedirs(new_trial_map, exist_ok=True)
        dictionary = {}
        for track_object in trial.getTrackObjects():
            if track_object.track_num not in altered_track_nums:
                x, y, z, t = track_object.getTrack()
                dictionary[f'Trial_{trial_num}_Track_{track_object.track_num}'] = {
                    'x': x,
                    'y': y,
                    'z': z,
                    'time': t
                    }
                df = pd.DataFrame(dictionary)
                file_path = os.path.join(new_trial_map, f'Trial_{trial_num}_Track_{track_object.track_num}.csv')
                df.to_csv(file_path)
            else:
                merged = []
                found = False
                i = 0
                while not found:
                    if i < len(stitch_num_land_take) and track_object.track_num == stitch_num_land_take[i][0]:
                        for track_object2 in trial.getTrackObjects():
                            if track_object2.track_num == stitch_num_land_take[i][1]:
                                merged = [a + b for a, b in zip(track_object.getTrack(), track_object2.getTrack())]
                                found = True
                    elif i < len(stitch_num_land_walk) and track_object.track_num == stitch_num_land_walk[i][0]:
                        for track_object2 in trial.getTrackObjects():
                            if track_object2.track_num == stitch_num_land_walk[i][1]:
                                merged = [a + b for a, b in zip(track_object.getTrack(), track_object2.getTrack())]
                                found = True
                    elif i < len(stitch_num_walk_take) and track_object.track_num == stitch_num_walk_take[i][0]:
                        for track_object2 in trial.getTrackObjects():
                            if track_object2.track_num == stitch_num_walk_take[i][1]:
                                merged = [a + b for a, b in zip(track_object.getTrack(), track_object2.getTrack())]
                                found = True
                    elif i < len(stitch_num_land_walk_take) and track_object.track_num == stitch_num_land_walk_take[i][0]:
                        for track_object2 in trial.getTrackObjects():
                            if track_object2.track_num == stitch_num_land_walk_take[i][1]:
                                for track_object3 in trial.getTrackObjects():
                                    if track_object3.track_num == stitch_num_land_walk_take[i][2]:
                                        merged_intermediate = [a + b for a, b in zip(track_object.getTrack(), track_object2.getTrack())]
                                        merged = [a + b for a, b in zip(merged_intermediate, track_object3.getTrack())]
                                        found = True
                    i += 1
                if merged:
                    x, y, z, t = merged
                    dictionary[f'Trial_{trial_num}_Track_{track_object.track_num}'] = {
                        'x': x,
                        'y': y,
                        'z': z,
                        'time': t
                    }
                    df = pd.DataFrame(dictionary)
                    file_path = os.path.join(new_trial_map, f'Trial_{trial_num}_Track_{track_object.track_num}.csv')
                    df.to_csv(file_path)
                else:
                    raise Exception("Something went wrong")

        new_boundary_trial_map = os.path.join(basemap_boundary_tracks_path, f'Trial_{trial_num}')
        os.makedirs(new_boundary_trial_map, exist_ok=True)
        for track_object in trial.getTrackObjects():
            if track_object.track_num not in altered_track_nums:
                x, y, z, t = track_object.getTrack()
                in_run = False
                x_hop, y_hop, z_hop, t_hop = [], [], [], []
                num_hops = 0
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
                            num_hops += 1
                            dictionary[f'Trial_{trial_num}_Track_{track_object.track_num}_Hop_{num_hops}'] = {
                                'x': x_hop,
                                'y': y_hop,
                                'z': z_hop,
                                'time': t_hop
                            }
                            x_hop, y_hop, z_hop, t_hop = [], [], [], []
                    if in_run:
                        # When a track ends in landing, the landing still gets added
                        num_hops += 1
                        dictionary[f'Trial_{trial_num}_Track_{track_object.track_num}_Hop_{num_hops}'] = {
                            'x': x_hop,
                            'y': y_hop,
                            'z': z_hop,
                            'time': t_hop
                        }

            for i, track in enumerate(paired_tracks):
                x, y, z, t = track
                dictionary[f'Trial_{trial_num}_Paired_Track_{i}'] = {
                    'x': x,
                    'y': y,
                    'z': z,
                    'time': t
                }
            for i, track in enumerate(new_landings_tracks):
                x, y, z, t = track
                dictionary[f'Trial_{trial_num}_Landing_Track_{i}'] = {
                    'x': x,
                    'y': y,
                    'z': z,
                    'time': t
                }
            for i, track in enumerate(new_take_off_tracks):
                x, y, z, t = track
                dictionary[f'Trial_{trial_num}_Takeoff_Track_{i}'] = {
                    'x': x,
                    'y': y,
                    'z': z,
                    'time': t
                }
            for i, track in enumerate(new_walking_tracks):
                x, y, z, t = track
                dictionary[f'Trial_{trial_num}_Walking_Track_{i}'] = {
                    'x': x,
                    'y': y,
                    'z': z,
                    'time': t
                }

            df = pd.DataFrame(dictionary)
            file_path = os.path.join(new_boundary_trial_map, f'Trial_{trial_num}_Tracks.csv')
            df.to_csv(file_path)


if __name__ =='__main__':
    makePairedCSVDataset(path_csv_folder=path_csv_folder1)

