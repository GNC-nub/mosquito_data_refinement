import numpy as np

def distance_to_cylinder_surface(point, cylinder):
    """
    Shortest Euclidean distance from a 3D point to the *surface* of a finite cylinder
    aligned with the z-axis, centered at the origin, with radius R and caps at z_min/z_max.
    """
    x, y, z = point
    radius, z_min, z_max = cylinder

    rho = np.sqrt(x**2 + y**2)
    radial_distance = np.abs(rho - radius)

    if z_min <= z <= z_max:
        return radial_distance
    else:
        vertical_distance = min(np.abs(z - z_min), np.abs(z - z_max))
        return np.sqrt(radial_distance**2 + vertical_distance**2)


def find_resting_points_in_track(x, y, z, t=None, boundary=0.02):
    """
    If `t` is not given:
        Returns [(x, y, z), ...]

    If `t` is given:
        Returns [(x, y, z, duration), ...] where duration is how long the
        track stayed inside the boundary segment containing that point.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    n = len(x)
    if not (len(y) == n and len(z) == n):
        raise ValueError("x, y, z must have the same length")

    if t is not None:
        t = np.asarray(t, dtype=float)
        if len(t) != n:
            raise ValueError("t must have same length as x/y/z")

    # Vectorized mask of being "inside" the boundary
    mask = landing_area_array(x, y, z, boundary=boundary)

    if not np.any(mask):
        return []

    # Find contiguous True segments
    edges = np.diff(np.concatenate(([0], mask.astype(int), [0])))
    start_indices = np.where(edges == 1)[0]
    end_indices = np.where(edges == -1)[0]

    # Trap cylinders
    cylinders = [(0.15, -0.38, -0.083), (0.055, -0.083, 0.0)]

    resting_points = []

    for start, end in zip(start_indices, end_indices):
        seg_slice = slice(start, end)
        segment_x = x[seg_slice]
        segment_y = y[seg_slice]
        segment_z = z[seg_slice]

        distances = np.array([
            min(
                distance_to_cylinder_surface((segment_x[i], segment_y[i], segment_z[i]), cyl)
                for cyl in cylinders
            )
            for i in range(len(segment_x))
        ])

        min_idx_local = int(np.argmin(distances))
        min_idx_global = start + min_idx_local

        if t is None:
            resting_points.append(
                (x[min_idx_global], y[min_idx_global], z[min_idx_global])
            )
        else:
            segment_t = t[seg_slice]
            duration = float(segment_t[-1] - segment_t[0])
            resting_points.append(
                (x[min_idx_global], y[min_idx_global], z[min_idx_global], duration)
            )

    return resting_points
