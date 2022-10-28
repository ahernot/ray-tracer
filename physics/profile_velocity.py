# profile_velocity v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np
from scipy import interpolate

from physics.profile_salinity import calc_S
from physics.profile_temperature import calc_T
from physics.model_velocity import sound_velocity_medwin


calc_c = lambda z: sound_velocity_medwin(z, calc_S(z), calc_T(z))  # vectorised

z = np.linspace(-10000, 0, 10000)
calc_dz_c = interpolate.interp1d(z, np.gradient(calc_c(z)), kind='quadratic')
