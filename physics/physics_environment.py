# 2D Physics Environment (v1.0)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# TODO: store cartesian product with coordinates for each sampled point (useful?)


import numpy as np
from scipy import interpolate

import physics
from physics.profile_salinity import calc_S
from physics.profile_temperature import calc_T
from physics.profile_pressure import calc_P
from physics.profile_ph import calc_pH

from physics.model_velocity import sound_velocity_medwin
from physics.model_rho import calc_rho

from physics.model_absorption import calc_dz_dG_coefs


X_DEFAULT = 0.
Z_DEFAULT = 0.


class PhysicsEnvironment2D:

    def __init__ (self, range_min, range_max, **kwargs):
        """
        Frequency-invariant parameters
        Makes everything into a spacial field (time-invariant)
        :param range_min:
        :param range_max:

        kwargs:
        :param res_x: horizontal resolution (in meters)
        :param res_z: vertical resolution (in meters)
        """

        # ! UNIDIRECTIONIAL FOR NOW (functions of z as only spatial dimension)

        self.__calc_S = physics.profile_salinity.calc_S  # Salinity as a function of z  # TODO: x, z
        self.__calc_T = physics.profile_temperature.calc_T  # Temperature as a function of z  # TODO: x, z
        self.__calc_P = physics.profile_pressure.calc_P  # Pressure as a function of z  # TODO: x, z
        self.__calc_pH = physics.profile_ph.calc_pH  # pH as a function of z  # TODO: x, z
        self.__calc_c = physics.model_velocity.sound_velocity_medwin
        self.__calc_rho = physics.model_rho.calc_rho
        self.__calc_dz_dG_coefs = physics.model_absorption.calc_dz_dG_coefs

        # self.res_x = kwargs.get('res_x', 0)  # Default is x-invariant  # TODO: set default values
        self.res_z = kwargs.get('res_z', 10)  # TODO: set default values
        self.range_min = range_min
        self.range_max = range_max

        # self.x = np.arange(self.range_min[0], self.range_max[0], self.res_x) if self.res_x != 0 else np.array([X_DEFAULT])  # Check if x-invariant
        self.z = np.arange(self.range_min[1], self.range_max[1], self.res_z) if self.res_z != 0 else np.array([Z_DEFAULT])  # Check if z-invariant



    def generate_params (self):

        # Already interpolated (and interpolations not needed)
        self.S = self.__calc_S(self.z)
        self.T = self.__calc_T(self.z)
        self.P = self.__calc_P(self.z)
        self.pH = self.__calc_pH(self.z)

        # Velocity
        self.c = self.__calc_c (self.z, self.S, self.T)
        self.calc_c = interpolate.interp1d(self.z, self.c, kind='quadratic')
        self.dz_c = np.gradient(self.c)
        self.calc_dz_c = interpolate.interp1d(self.z, self.dz_c, kind='quadratic')

        # Rho
        self.rho = self.__calc_rho(self.z, self.S, self.T, self.P)
        self.calc_rho = interpolate.interp1d(self.z, self.rho, kind='quadratic')

        # Absorption (keep full frequency resolution, interpolate spatially)
        self.__dz_dG_coefs = self.__calc_dz_dG_coefs(self.z, self.T, self.S, self.pH)
        self.__calc_dz_dG_coefs_interp = interpolate.interp1d(self.z, self.__dz_dG_coefs, kind='quadratic')
        self.calc_dz_dG = lambda f, z: self.__calc_dz_dG(f, z, self.__calc_dz_dG_coefs_interp)
        
        


    def extend (self, range_min, range_max):
        # Extend environment (shouldn't need to be used)

        # extend params

        # extend v, rho, absorption

        # re-interpolate
        pass
