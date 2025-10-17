import numpy as np
import joblib
import matplotlib.pyplot as plt
from IPython.display import display
import ipywidgets as widgets

from CommonFunctions import getTrap, landing_area_array

DF_PATH = '/Users/nubia/Desktop/Thesis_2.0/dataset/df.joblib'
MERGED_DF_PATH = '/Users/nubia/Desktop/Thesis_2.0/dataset/merged_df.joblib'

def _load_data(df=None, merged_df=None):
    if df is None:
        df = joblib.load(DF_PATH)
    if merged_df is None:
        merged_df = joblib.load(MERGED_DF_PATH)
    return df, merged_df

def _get_track(df, trial_key='Trial_1', track_key='Track_3'):
    row = df[trial_key].loc[track_key]
    x = np.asarray(row['x'])
    y = np.asarray(row['y'])
    z = np.asarray(row['z'])
    t = np.asarray(row['time']) if 'time' in row and row['time'] is not None else np.arange(x.size, dtype=float)
    return x, y, z, t

def show_interactive_track(trial_key='Trial_1', track_key='Track_3', boundary=0.02, df=None, merged_df=None):
    """3D interactive visualization of a single track with a time slider."""
    df, merged_df = _load_data(df, merged_df)
    x, y, z, t = _get_track(df, trial_key, track_key)
    if x.size == 0:
        raise ValueError('Empty track data.')

    mask = landing_area_array(x, y, z, boundary=boundary)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    xgb, ygb, zgb, xgi, ygi, zgi = getTrap()
    ax.plot_surface(xgb, ygb, zgb, alpha=0.35, color='b')
    ax.plot_surface(xgi, ygi, zgi, alpha=0.35, color='b')

    ax.plot(x, y, z, color='gray', linewidth=1.2, alpha=0.9)

    point = ax.scatter([x[0]], [y[0]], [z[0]], s=40, c=('red' if mask[0] else 'blue'))
    path, = ax.plot(x[:1], y[:1], z[:1], color='orange', linewidth=2, alpha=0.9)

    ax.set_title(f'Interactive Track: {track_key} in {trial_key} (boundary={boundary})')
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    xr, yr, zr = np.ptp(x), np.ptp(y), np.ptp(z)
    ax.set_box_aspect((xr if xr else 1, yr if yr else 1, zr if zr else 1))
    ax.view_init(elev=20, azim=200)

    time_text = ax.text2D(0.02, 0.95, f't = {t[0]:.3f} s', transform=ax.transAxes)

    slider = widgets.IntSlider(value=0, min=0, max=len(t)-1, step=1, description='Index', continuous_update=True)
    play = widgets.Play(interval=50, value=0, min=0, max=len(t)-1, step=1, description='')  # <- fix
    play.style = {'description_width': '0px'}  # optional: hide label space
    widgets.jslink((play, 'value'), (slider, 'value'))
    time_label = widgets.Label(value=f't = {t[0]:.3f} s')


    def update(i):
        color = 'red' if mask[i] else 'blue'
        point._offsets3d = (np.array([x[i]]), np.array([y[i]]), np.array([z[i]]))
        point.set_color(color)
        path.set_data(x[:i+1], y[:i+1])
        path.set_3d_properties(z[:i+1])
        time_text.set_text(f't = {float(t[i]):.3f} s')
        time_label.value = f't = {float(t[i]):.3f} s'
        fig.canvas.draw_idle()

    def _on_change(change):
        if change['name'] == 'value':
            update(change['new'])

    slider.observe(_on_change, names='value')

    controls = widgets.HBox([play, slider, time_label])
    display(widgets.VBox([controls, fig.canvas]))

    update(0)
    return controls

def show_track3_trial1(boundary=0.02, df=None, merged_df=None):
    return show_interactive_track('Trial_1', 'Track_3', boundary=boundary, df=df, merged_df=merged_df)


