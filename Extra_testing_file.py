# testing file extra: 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec


# ==== Your functions pasted here ====
from CommonFunctions import * 
from matplotlib.colors import LinearSegmentedColormap

df_path = '/Users/nubia/Desktop/Thesis_2.0/dataset/df.joblib'
merged_df_path = '/Users/nubia/Desktop/Thesis_2.0/dataset/merged_df.joblib'

import joblib 

df = joblib.load(merged_df_path)
import numpy as np
import matplotlib.pyplot as plt

# ---- params (these must match the ones used in your landing functions) ----
BOUNDARY = 0.02
TRAP_RADIUS = 0.15
INLET_RADIUS = 0.055
INLET_HEIGHT = 0.083
TRAP_HEIGHT = 0.388   # body height (positive number)
ZOOM_SIZE = 0.10      # total half-width/height for zoom window in r and z

# --- region masks that match your geometry definitions ---

def mask_inlet_side(x, y, z, boundary=BOUNDARY, inlet_height=INLET_HEIGHT, inlet_radius=INLET_RADIUS):
    r = np.sqrt(x**2 + y**2)
    return (
        ((-boundary < z) & (z < 0) &
         (inlet_radius < r) & (r < inlet_radius + boundary))
        |
        ((-(inlet_height - boundary) < z) & (z < -boundary) &
         (inlet_radius - boundary < r) & (r < inlet_radius + boundary))
    )

def mask_body_side(x, y, z, boundary=BOUNDARY,
                   trap_height=TRAP_HEIGHT, trap_radius=TRAP_RADIUS, inlet_height=INLET_HEIGHT):
    r = np.sqrt(x**2 + y**2)
    return (
        (-trap_height < z) & (z < -(inlet_height + boundary)) &
        (trap_radius - boundary < r) & (r < trap_radius + boundary)
    )

def mask_top_ring(x, y, z, boundary=BOUNDARY, trap_radius=TRAP_RADIUS,
                  inlet_height=INLET_HEIGHT, inlet_radius=INLET_RADIUS):
    r = np.sqrt(x**2 + y**2)
    return (
        (-(inlet_height + boundary) < z) & (z < -(inlet_height - boundary)) &
        (inlet_radius - boundary < r) & (r < trap_radius + boundary)
    )

# ---- helper: compute closest-to-wall point in a region (one point) ----

def closest_wall_point_in_region(x, y, z, region, boundary=BOUNDARY):
    """
    Returns (idx, r0, z0) of the closest-to-wall point within the chosen region.
    If no points lie in the region, returns (None, None, None).
    """
    x = np.asarray(x); y = np.asarray(y); z = np.asarray(z)
    r = np.sqrt(x**2 + y**2)

    if region == 'inlet':
        region_mask = mask_inlet_side(x, y, z, boundary=boundary)
        # distance to inlet cylinder surface (radial distance only; region already limits z)
        distances = np.abs(r - INLET_RADIUS)
    elif region == 'body':
        region_mask = mask_body_side(x, y, z, boundary=boundary)
        distances = np.abs(r - TRAP_RADIUS)
    elif region == 'top':
        region_mask = mask_top_ring(x, y, z, boundary=boundary)
        # closest “wall” here is either inlet radius or body radius (whichever is closer)
        distances = np.minimum(np.abs(r - INLET_RADIUS), np.abs(r - TRAP_RADIUS))
    else:
        raise ValueError("region must be one of {'inlet','body','top'}")

    if not np.any(region_mask):
        return None, None, None

    # choose the single point with min distance among masked points
    candidate_idxs = np.where(region_mask)[0]
    best_local = np.argmin(distances[region_mask])
    idx = candidate_idxs[best_local]
    r0 = r[idx]
    z0 = z[idx]
    return idx, r0, z0

# ---- main: plot one track with three zoomed panels ----

def plot_zoomed_track(df, trial_num=1, track_label='Track_1',
                      boundary=BOUNDARY, zoom_size=ZOOM_SIZE):
    """
    Makes a 1x3 subplot with zoomed views around:
      - inlet cylinder side
      - body cylinder side
      - top ring (just above the body, to the right / outside of the body)
    The track is shown in r-z space. One resting point (closest to wall) is shown in red per panel.
    """
    # get the track
    row = df[f'Trial_{trial_num}'].loc[track_label]
    x = np.array(row['x']); y = np.array(row['y']); z = np.array(row['z'])
    r = np.sqrt(x**2 + y**2)

    # trap polygons (r-z)
    inlet_r, inlet_z, body_r, body_z = getTrap2D()

    # figure
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    regions = [('inlet', 'Inlet side'),
               ('body',  'Body side'),
               ('top',   'Top ring (just above body)')]

    for ax, (region_key, region_title) in zip(axes, regions):
        # background trap
        ax.fill(body_r, body_z, color='purple', linewidth=0, alpha=0.5)
        ax.fill(inlet_r, inlet_z, color='purple', linewidth=0, alpha=0.5)

        # full track in r-z
        ax.plot(r, z, lw=1.5, alpha=0.7)

        # region mask (for visual cue)
        if region_key == 'inlet':
            region_mask = mask_inlet_side(x, y, z, boundary=boundary)
        elif region_key == 'body':
            region_mask = mask_body_side(x, y, z, boundary=boundary)
        else:
            region_mask = mask_top_ring(x, y, z, boundary=boundary)

        # highlight region points (optional)
        if np.any(region_mask):
            ax.scatter(r[region_mask], z[region_mask], s=8, alpha=0.6)

        # resting point in this region
        idx, r0, z0 = closest_wall_point_in_region(x, y, z, region_key, boundary=boundary)
        if idx is not None:
            ax.scatter([r0], [z0], s=50, color='red', zorder=5, label='resting point')
            ax.set_xlim(r0 - zoom_size, r0 + zoom_size)
            ax.set_ylim(z0 - zoom_size, z0 + zoom_size)
        else:
            # fallback zoom: around the region’s typical location
            if region_key == 'inlet':
                cx, cz = INLET_RADIUS, -INLET_HEIGHT/2
            elif region_key == 'body':
                cx, cz = TRAP_RADIUS, -(INLET_HEIGHT + (TRAP_HEIGHT - INLET_HEIGHT)/2)
            else:  # top ring
                cx, cz = (INLET_RADIUS + TRAP_RADIUS)/2, -INLET_HEIGHT
            ax.set_xlim(cx - zoom_size, cx + zoom_size)
            ax.set_ylim(cz - zoom_size, cz + zoom_size)
            ax.text(cx, cz, "No points in region", ha='center', va='center', fontsize=9)

        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('r')
        ax.set_ylabel('z')
        ax.set_title(region_title)

    plt.tight_layout()
    plt.show()
