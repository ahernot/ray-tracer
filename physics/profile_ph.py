# profile_ph v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# NO MODEL

import numpy as np
from scipy import interpolate


pH_profile = np.array([
    [0, 7.98],
    [-100, 7.88],
    [-200, 7.84],
    [-250, 7.87],
    [-300, 7.85],
    [-400, 7.79],
    [-450, 7.7],
    [-500, 7.7],
    [-600, 7.71],
    [-800, 7.73],
    [-1000, 7.72],
    [-2000, 7.73],
    [-3000, 7.74],
    [-10000, 7.74]
])

calc_pH = interpolate.interp1d(
    pH_profile.T[0][::-1],
    pH_profile.T[1][::-1],
    kind='linear',
    bounds_error=False,
    fill_value='extrapolate'
)
