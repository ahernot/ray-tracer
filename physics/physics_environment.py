# 2D Physics Environment (v1.0)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# PARAMS: Salinity, Pressure, Temperature, pH
#       Velocity, Rho(density), Absorption


import numpy as np

import physics
from physics.model_absorption import calc_absorption_dB


X_DEFAULT = 0.
Z_DEFAULT = 0.


class PhysicsEnvironment2D:

    def __init__ (self, range_min, range_max, **kwargs):
        """
        :param range_min:
        :param range_max:

        kwargs:
        :param res_x: horizontal resolution (in meters)
        :param res_z: vertical resolution (in meters)
        """

        self.res_x = kwargs.get('res_x', 0)  # Default is x-invariant  # TODO: set default values
        self.res_z = kwargs.get('res_z', 100)  # TODO: set default values
        self.range_min = range_min
        self.range_max = range_max

        self.x = np.arange(self.range_min[0], self.range_max[0], self.res_x) if self.res_x != 0 else np.array([X_DEFAULT])  # Check if x-invariant
        self.z = np.arange(self.range_min[1], self.range_max[1], self.res_z) if self.res_z != 0 else np.array([Z_DEFAULT])  # Check if z-invariant

        # cartesian product with coordinates for each sampled point

        # 2D interpolations


    def generate_params (self):
        pass


    def extend (self, range_min, range_max):
        # Extend environment
        pass

