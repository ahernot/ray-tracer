# 2D Ray Propagation Model (v3.0)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

from numpy import np

from geometry import Point2D, Vector2D
from environment import Environment2D


DX_MAX_DEFAULT = 0.1
DZ_MAX_DEFAULT = 0.1
N_STEPS_MAX = 100000
N_REBOUNDS_MAX = -1  # -1=infinity


# TODO New version: each ray is subdivided into segments and with each rebound new subclasses are created as a tree (can have multiple branches per rebound)

class Ray2D:

    def __init__ (self, env: Environment2D, source: Point2D, angle, **kwargs):
        """
        :param env:
        :param source:
        :param angle: Casting angle, in radians
        """

        self.__source: Point2D = source
        self.__env: Environment2D = env
        self.__angle_init = angle
        self.__is_propagated = False

        self.n_steps_max = kwargs.get('n_steps_max', N_STEPS_MAX)
        self.n_rebounds_max = kwargs.get('n_rebounds_max', N_REBOUNDS_MAX)

        # Minimum resolutions
        self.dx_max = kwargs.get('dx_max', DX_MAX_DEFAULT)
        self.dz_max = kwargs.get('dz_max', DZ_MAX_DEFAULT)

        # Initialise
        self.X = np.array([self.__source.x, ])
        self.Z = np.array([self.__source.z, ])


    def __repr__ (self):
        return f'Ray object'  # TODO: improve repr

    def propagate (self):  # TODO: add propagation params
        if self.__is_propagated: return

        # Initialise
        angle = -1 * self.__angle_init + (np.pi / 2)  # convert angle notation
        x, z = self.__source.coordinates()
        dx_z = 1 / np.tan(angle)
        dxdx_z = 0  # no initial curvature
        # Initialise solver
        c0 = self.calc_c(0)
        mult = -1 * np.power(c0 / np.sin(angle), 2)  # differential equation multiplier



        for i in range(self.n_steps_max):

            # Calculate new integration segment with angle from previous step
            if np.abs(dx_z) > self.dx_z_max:  # steeper slope => limit dz=dz_max
                dz = np.sign(dx_z) * self.dz_max
                dx = dz / dx_z
            else:  # smaller slope => limit dx=dx_max
                dx = self.dx_max
                dz = dx_z * dx
            # Calculate new point (before reflections)
            x_new = x + dx
            z_new = z + dz



            # Check backwards propagation (if not allowed): TODO

            # Check reflections
            # floor
            # ceiling

            # Check simulation bounds (!!!! USE NEW X, Z AFTER REFLECTION EVENT)
            if x not in range (self.__env.range_min.x, self.__env.range_max.x) or z not in range (self.__env.range_min.z, self.__env.range_max.z):
                print('out of bounds. stopping')
                # TODO: plot last point at x = xmax or z = zmax
                break

        

            # Add new point
            self.X = np.concatenate((self.X, np.array([x_new, ])), axis=0)
            self.Z = np.concatenate((self.Z, np.array([z_new, ])), axis=0)
            x = x_new
            z = z_new
            # Calculate new point's properties
            c = self.calc_c (z_new)
            dz_c = self.calc_dz_c (z_new)
            # Update derivatives
            dxdx_z = mult * dz_c / np.power(c, 3)
            dx_z += dxdx_z * dx

        # # Generate interpolated path function
        # self.Z_func = interpolate.interp1d(self.X, self.Z, kind='linear')

        self.__is_propagated = True
        