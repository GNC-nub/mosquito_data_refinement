'''
READ ME
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
'''

import numpy as np

def landing_area(x, y, z, landing_width=0.03, trap_height = 0.388, trap_radius = 0.15, inlet_height = 0.083, inlet_radius = 0.055):
    distance = np.sqrt((x**2)+(y**2))
    landing = False

    if -landing_width < z < 0:
        if inlet_radius < distance < inlet_radius + landing_width:
            landing = True
    elif -(inlet_height - landing_width) < z < -landing_width:
        if inlet_radius - landing_width < distance < inlet_radius + landing_width:
            landing = True
    elif -(inlet_height + landing_width) < z < -(inlet_height - landing_width):
        if inlet_radius - landing_width < distance < trap_radius + landing_width:
            landing = True
    elif -trap_height < z < -(inlet_height + landing_width):
        if trap_radius - landing_width < distance < trap_radius + landing_width:
            landing = True
    return landing

def capturing_area(x, y, z, capturing_width = 0.03, inlet_radius = 0.055):
    distance = np.sqrt((x**2)+(y**2))
    capture = False
    if 0 < z < capturing_width:
        if distance < (inlet_radius + capturing_width):
            capture = True
    elif -capturing_width < z < 0:
        if distance < inlet_radius:
            capture = True
    return capture


