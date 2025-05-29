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

    for trial_num in range(1, 65):
        if trial_num == 65:
            trial = ClassMosquito.Trial(trial_num)
            paired_tracks, new_landings_tracks, new_take_off_tracks, new_walking_tracks, stitch_num_land_take, stitch_num_land_walk_take, stitch_num_land_walk, stitch_num_walk_take, altered_track_nums = trial.generatePairsForCSV(radius=radius, boundary=boundary)

            new_trial_map = os.path.join(basemap_paired_path, f'Trial_{trial_num}')
            os.makedirs(new_trial_map, exist_ok=True)
            track_objects1 = trial.getTrackObjects()
            for track_object in track_objects1:
                track_num1 = track_object.track_num
                if track_num1 not in altered_track_nums:
                    print(f'stage0: Trial {trial_num}, track {track_num1} is normal')
                    x, y, z, t = track_object.getTrack()
                    dictionary = {}
                    dictionary[f'Trial_{trial_num}_Track_{track_num1}'] = {
                        'x': x,
                        'y': y,
                        'z': z,
                        'time': t
                        }
                    df = pd.DataFrame(dictionary)
                    file_path = os.path.join(new_trial_map, f'Trial_{trial_num}_Track_{track_num1}.csv')
                    df.to_csv(file_path)
            print('stage1')
            for i in range(len(stitch_num_land_take)):
                print(f'Trial {trial_num}, stage1.1')
                track1 = ClassMosquito.Track(trial_num, stitch_num_land_take[i][0])
                track2 = ClassMosquito.Track(trial_num, stitch_num_land_take[i][1])
                merged = [a + b for a, b in zip(track1.getTrack(), track2.getTrack())]
                x, y, z, t = merged
                dictionary = {}
                title = f'Trial_{trial_num}_land_Track_{stitch_num_land_take[i][0]}_take_Track_{stitch_num_land_take[i][1]}'
                dictionary[title] = {
                    'x': x,
                    'y': y,
                    'z': z,
                    'time': t
                }
                df = pd.DataFrame(dictionary)
                file_path = os.path.join(new_trial_map, f'{title}.csv')
                df.to_csv(file_path)
            for i in range(len(stitch_num_land_walk)):
                print(f'Trial {trial_num}, stage1.2')
                track1 = ClassMosquito.Track(trial_num, stitch_num_land_walk[i][0])
                track2 = ClassMosquito.Track(trial_num, stitch_num_land_walk[i][1])
                merged = [a + b for a, b in zip(track1.getTrack(), track2.getTrack())]
                x, y, z, t = merged
                dictionary = {}
                title = f'Trial_{trial_num}_land_Track{stitch_num_land_walk[i][0]}_walk_Track_{stitch_num_land_walk[i][1]}'
                dictionary[title] = {
                    'x': x,
                    'y': y,
                    'z': z,
                    'time': t
                }
                df = pd.DataFrame(dictionary)
                file_path = os.path.join(new_trial_map, f'{title}.csv')
                df.to_csv(file_path)
            for i in range(len(stitch_num_walk_take)):
                print(f'Trial {trial_num}, stage1.3')
                track1 = ClassMosquito.Track(trial_num, stitch_num_walk_take[i][0])
                track2 = ClassMosquito.Track(trial_num, stitch_num_walk_take[i][1])
                merged = [a + b for a, b in zip(track1.getTrack(), track2.getTrack())]
                x, y, z, t = merged
                dictionary = {}
                title = f'Trial_{trial_num}_walk_Track{stitch_num_walk_take[i][0]}_take_Track_{stitch_num_walk_take[i][1]}'
                dictionary[title] = {
                    'x': x,
                    'y': y,
                    'z': z,
                    'time': t
                }
                df = pd.DataFrame(dictionary)
                file_path = os.path.join(new_trial_map, f'{title}.csv')
                df.to_csv(file_path)
            for i in range(len(stitch_num_land_walk_take)):
                print(f'Trial {trial_num}, stage1.4')
                track1 = ClassMosquito.Track(trial_num, stitch_num_land_walk_take[i][0])
                track2 = ClassMosquito.Track(trial_num, stitch_num_land_walk_take[i][1])
                track3 = ClassMosquito.Track(trial_num, stitch_num_land_walk_take[i][2])
                merged1 = [a + b for a, b in zip(track1.getTrack(), track2.getTrack())]
                merged2 = [a + b for a, b in zip(merged1, track3.getTrack())]
                x, y, z, t = merged2
                dictionary = {}
                title = f'Trial_{trial_num}_land_Track{stitch_num_land_walk_take[i][0]}_walk_Track{stitch_num_land_walk_take[i][1]}_take_Track_{stitch_num_land_walk_take[i][2]}'
                dictionary[title] = {
                    'x': x,
                    'y': y,
                    'z': z,
                    'time': t
                }
                df = pd.DataFrame(dictionary)
                file_path = os.path.join(new_trial_map, f'{title}.csv')
                df.to_csv(file_path)
            new_boundary_trial_map = os.path.join(basemap_boundary_tracks_path, f'Trial_{trial_num}')
            os.makedirs(new_boundary_trial_map, exist_ok=True)
            print(f'Trial {trial_num},stage3')
            for track_object in trial.getTrackObjects():
                if track_object.track_num not in altered_track_nums:
                    x, y, z, t = track_object.getTrack()
                    in_run = False
                    x_hop, y_hop, z_hop, t_hop = [], [], [], []
                    num_hops = 0
                    for i in range(len(x)):
                        if landing_area(x[i], y[i], z[i], boundary=boundary):
                            print(f'Trial {trial_num}, stage4')
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
                                title = f'Trial_{trial_num}_Track_{track_object.track_num}_Hop_{num_hops}'
                                dictionary = {}
                                dictionary[title] = {
                                    'x': x_hop,
                                    'y': y_hop,
                                    'z': z_hop,
                                    'time': t_hop
                                }
                                df = pd.DataFrame(dictionary)
                                file_path = os.path.join(new_boundary_trial_map, f'{title}.csv')
                                df.to_csv(file_path)
                                x_hop, y_hop, z_hop, t_hop = [], [], [], []
                        if in_run:
                            # When a track ends in landing, the landing still gets added
                            num_hops += 1
                            title = f'Trial_{trial_num}_Track_{track_object.track_num}_Hop_{num_hops}'
                            dictionary = {}
                            dictionary[title] = {
                                'x': x_hop,
                                'y': y_hop,
                                'z': z_hop,
                                'time': t_hop
                            }
                            df = pd.DataFrame(dictionary)
                            file_path = os.path.join(new_boundary_trial_map, f'{title}.csv')
                            df.to_csv(file_path)

                for i, track in enumerate(paired_tracks):
                    print(f'Trial {trial_num}, stage5.1')
                    x, y, z, t = track
                    title = f'Trial_{trial_num}_Paired_Track_{i}'
                    dictionary = {}
                    dictionary[title] = {
                        'x': x,
                        'y': y,
                        'z': z,
                        'time': t
                    }
                    df = pd.DataFrame(dictionary)
                    file_path = os.path.join(new_boundary_trial_map, f'{title}.csv')
                    df.to_csv(file_path)
                for i, track in enumerate(new_landings_tracks):
                    print(f'Trial {trial_num}, stage5.2')
                    x, y, z, t = track
                    title = f'Trial_{trial_num}_Landing_Track_{i}'
                    dictionary = {}
                    dictionary[title] = {
                        'x': x,
                        'y': y,
                        'z': z,
                        'time': t
                    }
                    df = pd.DataFrame(dictionary)
                    file_path = os.path.join(new_boundary_trial_map, f'{title}.csv')
                    df.to_csv(file_path)
                for i, track in enumerate(new_take_off_tracks):
                    print(f'Trial {trial_num}, stage5.3')
                    x, y, z, t = track
                    title = f'Trial_{trial_num}_Takeoff_Track_{i}'
                    dictionary = {}
                    dictionary[title] = {
                        'x': x,
                        'y': y,
                        'z': z,
                        'time': t
                    }
                    df = pd.DataFrame(dictionary)
                    file_path = os.path.join(new_boundary_trial_map, f'{title}.csv')
                    df.to_csv(file_path)
                for i, track in enumerate(new_walking_tracks):
                    print(f'Trial {trial_num}, stage5.4')
                    x, y, z, t = track
                    title = f'Trial_{trial_num}_Walking_Track_{i}'
                    dictionary = {}
                    dictionary[title] = {
                        'x': x,
                        'y': y,
                        'z': z,
                        'time': t
                    }
                    df = pd.DataFrame(dictionary)
                    file_path = os.path.join(new_boundary_trial_map, f'{title}.csv')
                    df.to_csv(file_path)
            print(f'Trial {trial_num}, stage6')

if __name__ =='__main__':
    makePairedCSVDataset(path_csv_folder=path_csv_folder1)

