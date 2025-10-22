import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap



## ---------------    Landing Functions ---------------
# Now this function can handle whole arrays
def landing_area_side(x, y, z, boundary=0.02, trap_height=0.388, trap_radius=0.15, inlet_height=0.083,
                      inlet_radius=0.055):
    r = np.sqrt(x**2 + y**2)

    landing = (
        ((-boundary < z) & (z < 0) &
         (inlet_radius < r) & (r < inlet_radius + boundary))
        |
        ((-(inlet_height - boundary) < z) & (z < -boundary) &
         (inlet_radius - boundary < r) & (r < inlet_radius + boundary))
        |
        ((-trap_height < z) & (z < -(inlet_height + boundary)) &
         (trap_radius - boundary < r) & (r < trap_radius + boundary))
    )

    return landing  # boolean array


def landing_area_top(x, y, z, boundary=0.02, trap_radius=0.15, inlet_height=0.083,
                     inlet_radius=0.055):
    r = np.sqrt(x**2 + y**2)
    landing = (
        (-(inlet_height + boundary) < z) & (z < -(inlet_height - boundary)) &
        (inlet_radius - boundary < r) & (r < trap_radius + boundary)
    )
    return landing  # boolean array


def landing_area_array(x, y, z, boundary=0.02):
    return landing_area_top(x, y, z, boundary=boundary) | landing_area_side(x, y, z, boundary=boundary)

# ---- Evaluates if the x, y, z POINT input is within the landing area at the side of the trap ----
    
def landing_area_point_side(x, y, z, boundary=0.03, trap_height=0.388, trap_radius=0.15, inlet_height=0.083,
                      inlet_radius=0.055):
    r = np.sqrt(x**2 + y**2)
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


def landing_area_point_top(x, y, z, boundary=0.03, trap_radius=0.15, inlet_height=0.083,
                     inlet_radius=0.055):
    r = np.sqrt(x**2 + y**2)
    landing = False
    if -(inlet_height + boundary) < z < -(inlet_height - boundary):
        if inlet_radius - boundary < r < trap_radius + boundary:
            landing = True
    return landing


def landing_area_point(x, y, z, specific_area = 'whole', boundary=0.03):
    boolean = False
    if specific_area == 'whole':
        if landing_area_point_top(x, y, z, boundary=boundary) or landing_area_point_side(x, y, z, boundary=boundary):
            boolean = True
    elif specific_area == 'top':
        if landing_area_point_top(x, y, z, boundary=boundary):
            boolean = True
    elif specific_area == 'side':
        if landing_area_point_side(x, y, z, boundary=boundary):
            boolean = True
    return boolean



# ----- transform list of 3D x, y, z (or x, y, z points) to r, z 2D ----
def transformation_2D(x, y, z):
    """
    Convert 3D Cartesian coordinates (x, y, z) into 2D cylindrical coordinates (r, z).
    Works for scalars, lists, or NumPy arrays.
    """
    # Convert to NumPy arrays but preserve scalars
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    z = np.atleast_1d(z)

    r = np.sqrt(x**2 + y**2)

    # If the input was scalar, return scalars
    if r.size == 1:
        return float(r), float(z)
    return r, z


#  ------------- traps ----------------

def getTrap(body_lower_z=-0.38, body_upper_z=-0.083, inlet_upper_z=0,
            body_radius=0.15, inlet_radius=0.055, num_points=50):
    
    theta = np.linspace(0, 2 * np.pi, num_points)

    # Body cylinder
    z_body = np.linspace(body_lower_z, body_upper_z, num_points)
    theta_grid_body, z_grid_body = np.meshgrid(theta, z_body)
    x_grid_body = body_radius * np.cos(theta_grid_body)
    y_grid_body = body_radius * np.sin(theta_grid_body)

    # Inlet cylinder/cone
    z_inlet = np.linspace(body_upper_z, inlet_upper_z, num_points)
    theta_grid_inlet, z_grid_inlet = np.meshgrid(theta, z_inlet)
    x_grid_inlet = inlet_radius * np.cos(theta_grid_inlet)
    y_grid_inlet = inlet_radius * np.sin(theta_grid_inlet)

    return x_grid_body, y_grid_body, z_grid_body, x_grid_inlet, y_grid_inlet, z_grid_inlet

def getTrap2D(body_lower_z=-0.38, body_upper_z=-0.083,
              body_radius=0.15, inlet_radius=0.055):
    
    # Inlet polygon (r,z)
    inlet_r = np.array([0, inlet_radius, inlet_radius, 0])
    inlet_z = np.array([0, 0, body_upper_z, body_upper_z])

    # Body polygon (r,z)
    body_r = np.array([0, body_radius, body_radius, 0])
    body_z = np.array([body_upper_z, body_upper_z, body_lower_z, body_lower_z])

    return inlet_r, inlet_z, body_r, body_z

def plotTrap(ax):
    x_body, y_body, z_body, x_inlet, y_inlet, z_inlet = getTrap()

    # Body surface
    ax.plot_surface(x_body, y_body, z_body, alpha=0.3, color="purple", edgecolor="none")
    # Inlet surface
    ax.plot_surface(x_inlet, y_inlet, z_inlet, alpha=0.3, color="purple", edgecolor="none")


## ----------------- Resting points ------------------

def distance_to_cylinder_surface(point, cylinder):
    x, y, z = point
    radius, z_min, z_max = cylinder
    radial_distance = np.abs(np.sqrt(x**2 + y**2) - radius)

    if z_min <= z <= z_max:
        combined_distance = radial_distance
    else:
        vertical_distance = min(np.abs(z - z_min), np.abs(z - z_max))
        if radial_distance < 0:
            combined_distance = vertical_distance
        else:
            combined_distance = np.sqrt(radial_distance**2 + vertical_distance**2)
    return combined_distance

def find_resting_points_in_track(x, y, z, boundary = 0.02):
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    mask = landing_area_array(x, y, z, boundary=boundary)  # vectorized mask
    resting_points = []

    if not np.any(mask):
        return resting_points

    # Identify contiguous True segments
    inside = np.diff(np.concatenate(([0], mask.astype(int), [0])))
    start_indices = np.where(inside == 1)[0]
    end_indices = np.where(inside == -1)[0]

    # Define trap cylinders
    cylinders = [(0.15, -0.38, -0.083), (0.055, -0.083, 0)]

    for start, end in zip(start_indices, end_indices):
        segment_x = x[start:end]
        segment_y = y[start:end]
        segment_z = z[start:end]

        # Compute distances to trap surface for each point in the segment
        distances = np.array([
            min(distance_to_cylinder_surface((segment_x[i], segment_y[i], segment_z[i]), cyl) for cyl in cylinders)
            for i in range(len(segment_x))
        ])

        min_idx = np.argmin(distances)
        resting_points.append((segment_x[min_idx], segment_y[min_idx], segment_z[min_idx]))

    return resting_points


def resting_points(df, start_trial_num=1, end_trial_num=65, boundary = 0.2):
    all_x, all_y, all_z = [], [], []

    for trial_num in range(start_trial_num, end_trial_num):
        for _, track_data in df[f'Trial_{trial_num}'].iterrows():
            points = find_resting_points_in_track(track_data['x'], track_data['y'], track_data['z'], boundary=boundary)
            for px, py, pz in points:
                all_x.append(px)
                all_y.append(py)
                all_z.append(pz)

    return np.array(all_x), np.array(all_y), np.array(all_z)
    

## ------------ Resting times ---------------


def find_resting_times_in_track(x, y, z, time, boundary=0.02):
    """
    Returns: list of [x, y, z, duration]. Aka resting points with durations.
    """
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    time = np.asarray(time)
    mask = landing_area_array(x, y, z, boundary=boundary)  # vectorized mask
    resting_times = []
    
    if not np.any(mask):
        return resting_times

    inside = np.diff(np.concatenate(([0], mask.astype(int), [0])))
    start_indices = np.where(inside == 1)[0]
    end_indices = np.where(inside == -1)[0]

    # Trap cylinders (radius, z_min, z_max)
    cylinders = [(0.15, -0.38, -0.083), (0.055, -0.083, 0)]

    for start, end in zip(start_indices, end_indices):
        # Require at least two samples to get positive duration
        if end - start <= 1:
            continue
        dur = float(time[end - 1] - time[start])
        if dur <= 0:
            continue
        
        segment_x = x[start:end]
        segment_y = y[start:end]
        segment_z = z[start:end]

        # Compute distances to trap surface for each point in the segment
        distances = np.array([
            min(distance_to_cylinder_surface((segment_x[i], segment_y[i], segment_z[i]), cyl) for cyl in cylinders)
            for i in range(len(segment_x))
        ])

        min_idx = np.argmin(distances)
        resting_times.append((segment_x[min_idx], segment_y[min_idx], segment_z[min_idx], dur))
    return resting_times

def resting_times(df, start_trial_num=1, end_trial_num=65, boundary=0.2):
    all_x, all_y, all_z, all_durations = [], [], [], []

    for trial_num in range(start_trial_num, end_trial_num):
        for _, track_data in df[f'Trial_{trial_num}'].iterrows():
            points = find_resting_times_in_track(track_data['x'], track_data['y'], track_data['z'], track_data['time'], boundary=boundary)
            for px, py, pz, dur in points:
                all_x.append(px)
                all_y.append(py)
                all_z.append(pz)
                all_durations.append(dur)

    return np.array(all_x), np.array(all_y), np.array(all_z), np.array(all_durations)


## ------------- Density Matrices 2D ----------------


def get_volume_matrix(num_r_cells=20, num_z_cells=50, first_r_coord=0, last_r_coord=0.2,
                      lower_z_coord=-0.5, upper_z_coord=0, area_cell=0.0001):
    distance_row = np.linspace(first_r_coord + 0.005, last_r_coord - 0.005, num_r_cells)
    distance_matrix = np.tile(distance_row, (num_z_cells, 1))
    return distance_matrix * 2 * np.pi * area_cell


def density_matrix_2d_normalized(x, y, z, num_r_cells=20, num_z_cells=50,
                            first_r_coord=0, last_r_coord=0.2,
                            lower_z_coord=-0.5, upper_z_coord=0):
    r = np.sqrt(x**2 + y**2)
    r_edges = np.linspace(first_r_coord, last_r_coord, num_r_cells + 1)
    z_edges = np.linspace(lower_z_coord, upper_z_coord, num_z_cells + 1)
    hist, _, _ = np.histogram2d(r, z, bins=(r_edges, z_edges))
    hist = hist.T  # match z-r orientation

    volume_matrix = get_volume_matrix(num_r_cells, num_z_cells, first_r_coord, last_r_coord,
                                      lower_z_coord, upper_z_coord)
    density_matrix = hist / volume_matrix
    return density_matrix, r_edges, z_edges


def density_matrix_2d_weighted_normalized(x, y, z, weights, num_r_cells=20, num_z_cells=50,
                            first_r_coord=0, last_r_coord=0.2,
                            lower_z_coord=-0.5, upper_z_coord=0):
    r = np.sqrt(x**2 + y**2)
    r_edges = np.linspace(first_r_coord, last_r_coord, num_r_cells + 1)
    z_edges = np.linspace(lower_z_coord, upper_z_coord, num_z_cells + 1)
    hist, _, _ = np.histogram2d(r, z, bins=(r_edges, z_edges), weights=weights)
    hist = hist.T  # match z-r orientation

    volume_matrix = get_volume_matrix(num_r_cells, num_z_cells, first_r_coord, last_r_coord,
                                      lower_z_coord, upper_z_coord)
    density_matrix = hist / volume_matrix
    return density_matrix, r_edges, z_edges


def plot_density_heatpmap(density_matrix, r_edges, z_edges, dataset_name = 'df', boundary = 0.02):
    colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]
    cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)

    plt.figure(figsize=(6, 8))
    plt.imshow(density_matrix, origin='lower', aspect='auto',
               extent=[r_edges[0], r_edges[-1], z_edges[0], z_edges[-1]],
               cmap=cmap)

    plt.colorbar(label='Density per unit volume')
    inlet_r, inlet_z, body_r, body_z = getTrap2D()

    # Fill the trap body
    plt.fill(body_r, body_z, color='purple', linewidth = 0, alpha = 0.5)

    # Fill the inlet
    plt.fill(inlet_r, inlet_z, color='purple', linewidth = 0, alpha = 0.5)

    plt.xlim(0, 0.3)
    plt.ylim(-0.4, 0.1)
    plt.xlabel('r coordinate')
    plt.ylabel('z coordinate')
    plt.title(f'Resting Points Density Heatmap ({dataset_name}, boundary = {boundary})')
    plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
    plt.show()


def plot_density_heatpmap_restingtime(density_matrix, r_edges, z_edges, dataset_name = 'df', boundary = 0.02):
    colors = [(1, 1, 1), (1, 0.8, 0), (1, 0, 0), (0.5, 0, 0)]
    cmap = LinearSegmentedColormap.from_list("custom_red_hot", colors)

    plt.figure(figsize=(6, 8))
    plt.imshow(density_matrix, origin='lower', aspect='auto',
               extent=[r_edges[0], r_edges[-1], z_edges[0], z_edges[-1]],
               cmap=cmap)

    plt.colorbar(label='Resting time per unit volume')
    inlet_r, inlet_z, body_r, body_z = getTrap2D()

    plt.fill(body_r, body_z, color='purple', linewidth = 0, alpha = 0.5)
    plt.fill(inlet_r, inlet_z, color='purple', linewidth = 0, alpha = 0.5)

    plt.xlim(0, 0.3)
    plt.ylim(-0.4, 0.1)
    plt.xlabel('r coordinate')
    plt.ylabel('z coordinate')
    plt.title(f'Resting Times Density Heatmap ({dataset_name}, boundary = {boundary})')
    plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is square
    plt.show()
