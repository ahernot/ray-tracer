# profile_temperature v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# NO MODEL  # TODO: Quadratic or linear?

import numpy as np
from scipy import interpolate


temperature_profile = np.array([
    [0, 22],
    [-100, 21],
    [-200, 20],
    [-300, 16],
    [-400, 10],
    [-500, 4],
    [-600, 3.5],
    [-700, 3],
    [-800, 3],
    [-900, 3],
    [-1000, 3],
    [-1100, 3],
    [-1200, 3],
    [-1300, 3],
    [-1400, 3],
    [-1500, 3],
    [-10000, 3]
])

calc_T = interpolate.interp1d(
    temperature_profile.T[0][::-1],
    temperature_profile.T[1][::-1],
    kind='quadratic',
    bounds_error=False,
    fill_value='extrapolate'
)
