import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate


class Environment2D:

    def __init__ (self, floor=None, ceiling=None, **kwargs):
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

        self.__res = kwargs.get('res', 10000)
        mult = self.__res / self.range_max[0]

        self.__x = np.linspace(self.range_min[0], self.range_max[0], self.__res)  # x resolution across environment  # TODO: use PhysicsEnvironment2D.x
        self.__floor_sampled = self.floor(self.__x)
        self.__ceil_sampled = self.ceil(self.__x)

        self.dx_floor = interpolate.interp1d(self.__x, np.gradient(self.__floor_sampled) * mult, kind='quadratic')  # floor is static
        self.dx_ceil = interpolate.interp1d(self.__x, np.gradient(self.__ceil_sampled) * mult, kind='quadratic')  # ceiling is static

        # reflection power coefficient as a function of distance


    def plot (self, fig, **kwargs):
        plt.plot(self.__x, self.__ceil_sampled, figure=fig, **kwargs)
        plt.plot(self.__x, self.__floor_sampled, figure=fig, **kwargs)



# class EnvironmentStatic (Environment):
#     def __init__ (self, source: Point2D, floor, ceiling):
#         pass
