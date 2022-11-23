# 2D Environment (v2.0)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# TODO: self.__source = source
# TODO: function to model ray speed as a function of x and z
# TODO: current as additional displacement (vector field)


import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate

from preferences import *
from physics_environment import PhysicsEnvironment2D



class Environment2D:

    def __init__ (self, floor, ceiling, range_min, range_max, **kwargs):
        """
        2D Environment
        :param floor:
        :param ceiling:
        
        kwargs
        :param res_x:
        :param res_z:
        """

        # Initialise range and resolution
        self.res_x = kwargs.get('res_x', RES_X_DEFAULT)
        self.res_z = kwargs.get('res_z', RES_Z_DEFAULT)
        self.range_min = range_min
        self.range_max = range_max
        self.size = self.range_max - self.range_min

        # Generate physics environment
        self.penv = PhysicsEnvironment2D(self.range_min, self.range_max, res_x=self.res_x, res_z=self.res_z)
        self.penv.generate()

        # Process environment bounds
        self.floor = floor
        self.ceil = ceiling
        self.__floor_sampled = self.floor(self.penv.x)
        self.__ceil_sampled = self.ceil(self.penv.x)
        self.__dx_floor_sampled = np.gradient(self.__floor_sampled) / self.res_x
        self.__dx_ceil_sampled = np.gradient(self.__ceil_sampled) / self.res_x
        self.dx_floor = interpolate.interp1d(self.penv.x, self.__dx_floor_sampled, kind='quadratic', bounds_error=True)  # floor is static
        self.dx_ceil = interpolate.interp1d(self.penv.x, self.__dx_ceil_sampled, kind='quadratic', bounds_error=True)  # ceiling is static
        self.floor_avg = np.mean(self.__floor_sampled)
        self.ceil_avg = np.mean(self.__ceil_sampled)
        

    def __repr__ (self):
        return self.__class__.__name__ + ' object'

    def plot (self, fig, **kwargs):
        plt.plot(self.penv.x, self.__ceil_sampled, figure=fig, **kwargs)
        plt.plot(self.penv.x, self.__floor_sampled, figure=fig, **kwargs)



# class EnvironmentStatic (Environment):
#     def __init__ (self, source: Point2D, floor, ceiling):
#         pass
