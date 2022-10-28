import numpy as np
from scipy import interpolate


class Environment2D:

    def __init__ (self, floor=None, ceiling=None):
        # self.__source = source
        # function to model ray speed as a function of x and z
        # current as additional displacement (vector field)

        # TODO: Use discrete matrices instead (set resolution) for x,t derivation and time-moving limits (floor, …)

        #self.range = {'x': (0, 100000), 'z': (-10000, 10000)}
        self.range_min = np.array([0, -5500])
        self.range_max = np.array([100000, 0])
        self.size = self.range_max - self.range_min

        self.floor = floor
        self.ceil = ceiling

        res = 10000
        mult = res / self.range_max[0]
        x = np.linspace(self.range_min[0], self.range_max[0], res)
        self.dx_floor = interpolate.interp1d(x, np.gradient(self.floor(x)) * mult, kind='quadratic')  # floor is static
        self.dx_ceil = interpolate.interp1d(x, np.gradient(self.ceil(x)) * mult, kind='quadratic')  # ceiling is static

        # reflection power coefficient as a function of distance


# class EnvironmentStatic (Environment):
#     def __init__ (self, source: Point2D, floor, ceiling):
#         pass
