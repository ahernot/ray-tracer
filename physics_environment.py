# 2D Physics Environment (v1.1)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# TODO: add a lowres setting into Environment2D
# TODO: xx, yy = np.meshgrid(x, y) ?
# TODO: quadratic or linear interpolation?
# TODO: rename size to range


import numpy as np
from scipy import interpolate

from preferences import *

import physics
from physics.profile_salinity import calc_S
from physics.profile_temperature import calc_T
from physics.profile_pressure import calc_P
from physics.profile_ph import calc_pH

from physics.model_velocity import sound_velocity_medwin
from physics.model_rho import calc_rho
from physics.model_absorption import calc_dl_dG_coefs, calc_dl_dG
from physics.model_impedance import calc_Z


X_DEFAULT = 0.
Z_DEFAULT = 0.



class PhysicsEnvironment2D:

    def __init__ (self, range_min, range_max, **kwargs):
        """
        Frequency-invariant parameters
        Makes everything into a spacial field (time-invariant)
        :param range_min:
        :param range_max:

        kwargs
        :param res_x: horizontal resolution (in meters)
        :param res_z: vertical resolution (in meters)
        """
        # ! UNIDIRECTIONIAL FOR NOW (functions of z as only spatial dimension)

        self.res_x = kwargs.get('res_x', RES_X_DEFAULT)  # Default is x-invariant  # TODO: set default values
        self.res_z = kwargs.get('res_z', RES_Z_DEFAULT)  # TODO: set default values
        self.range_min = range_min
        self.range_max = range_max
        self.size = self.range_max - self.range_min

        n_x = int(self.size[0] // self.res_x) + 1 if self.res_x != 0 else 1  # Add one extra point to supersede resolution
        n_z = int(self.size[1] // self.res_z) + 1 if self.res_z != 0 else 1  # Add one extra point to supersede resolution
        self.x = np.linspace(self.range_min[0], self.range_max[0], n_x) if self.res_x != 0 else np.array([X_DEFAULT])  # Check if x-invariant
        self.z = np.linspace(self.range_min[1], self.range_max[1], n_z) if self.res_z != 0 else np.array([Z_DEFAULT])  # Check if z-invariant

        # Generate
        self.__generate()

    def __generate (self):

        # Already interpolated (and interpolations not needed)
        self.S = physics.profile_salinity.calc_S (self.z)
        self.T = physics.profile_temperature.calc_T (self.z)
        self.P = physics.profile_pressure.calc_P (self.z)
        self.pH = physics.profile_ph.calc_pH (self.z)

        # Velocity
        self.c = physics.model_velocity.sound_velocity_medwin (self.z, self.S, self.T)
        self.calc_c = interpolate.interp1d (self.z, self.c, kind='quadratic', bounds_error=False, fill_value='extrapolate')
        self.dz_c = np.gradient(self.c) / self.res_z
        self.calc_dz_c = interpolate.interp1d (self.z, self.dz_c, kind='quadratic', bounds_error=False, fill_value='extrapolate')

        # Rho
        self.rho = physics.model_rho.calc_rho (self.z, self.S, self.T, self.P)
        self.calc_rho = interpolate.interp1d (self.z, self.rho, kind='quadratic', bounds_error=False, fill_value='extrapolate')

        # Absorption (keep full frequency resolution, interpolate spatially)
        self.__dl_dG_coefs = physics.model_absorption.calc_dl_dG_coefs (self.z, self.T, self.S, self.pH)
        self.__calc_dl_dG_coefs_interp = interpolate.interp1d (self.z, self.__dl_dG_coefs, axis=1, kind='quadratic', bounds_error=False, fill_value='extrapolate')
        self.calc_dl_dG = lambda f, z: physics.model_absorption.calc_dl_dG (f, z, self.__calc_dl_dG_coefs_interp)

        # Impedance
        self.Z = physics.model_impedance.calc_Z(self.rho, self.c)
        self.calc_Z = interpolate.interp1d(self.z, self.Z, kind='quadratic', bounds_error=False, fill_value='extrapolate')
