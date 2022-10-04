import numpy as np
from scipy import interpolate

from geometry import Point2D, Vector2D


class Environment2D:

    def __init__ (self, floor=None, ceiling=None):
        # self.__source = source
        # function to model ray speed as a function of x and z
        # current as additional displacement (vector field)

        # TODO: Use discrete matrices instead (set resolution) for x,t derivation and time-moving limits (floor, …)

        #self.range = {'x': (0, 100000), 'z': (-10000, 10000)}
        self.range_min = Vector2D (0, -5500)
        self.range_max = Vector2D (100000, 0)

        self.floor = floor
        self.ceil = ceiling

        res = 10000
        x = np.linspace(self.range_min.x, self.range_max.x, res)
        self.dx_floor = interpolate.interp1d(x, np.gradient(self.floor(x)), kind='quadratic')  # floor is static
        self.dx_ceil = interpolate.interp1d(x, np.gradient(self.ceil(x)), kind='quadratic')  # ceiling is static

        # reflection power coefficient as a function of distance


# class EnvironmentStatic (Environment):
#     def __init__ (self, source: Point2D, floor, ceiling):
#         pass
