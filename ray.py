# 2D Ray Propagation Model (v3.0)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np

from scipy import interpolate  # TODO
from scipy.optimize import fsolve
from scipy.misc import derivative  # TODO

from geometry import Point2D, Vector2D, intersect
from physics.absorption_model import calc_absorption_dB
from physics.profile_velocity import calc_c, calc_dz_c
from physics.profile_salinity import calc_S
from physics.profile_temperature import calc_T
from physics.profile_ph import calc_pH
from physics.surface_reflection_model import calc_refcoef_surface
from physics.reflection import reflection_coefficient
from environment import Environment2D


DX_MAX_DEFAULT = 0.1
DZ_MAX_DEFAULT = 0.1
N_STEPS_MAX = 100000
N_REBOUNDS_MAX = -1  # -1=infinity


# TODO New version: each ray is subdivided into segments and with each rebound new subclasses are created as a tree (can have multiple branches per rebound)
# TODO: generalise as calc_c(x, z) and derive using numpy => use matrices for calculations (set resolution) / use 2D functions for calculations

class Ray2D:
    # Only one child ray per reflection: highest intensity specular reflection only

    def __init__ (self, env: Environment2D, source: Point2D, angle, **kwargs):
        """
        :param env: Simulation environment
        :param source: Source point
        :param angle: Casting angle, in radians
        """

        # Solver functions
        self.calc_der = kwargs.get('calc_der', derivative)
        self.func_solve = kwargs.get('func_solve', fsolve)

        self.__source: Point2D = source
        self.__env: Environment2D = env
        self.__angle_init = angle
        self.__is_propagated = False

        self.n_steps_max = kwargs.get('n_steps_max', N_STEPS_MAX)
        self.n_rebounds_max = kwargs.get('n_rebounds_max', N_REBOUNDS_MAX)

        # Minimum resolutions
        self.dx_max = kwargs.get('dx_max', DX_MAX_DEFAULT)
        self.dz_max = kwargs.get('dz_max', DZ_MAX_DEFAULT)
        self.dx_z_max = np.abs(self.dz_max / self.dx_max)

        # Initialise
        self.X = np.array([self.__source.x, ])
        self.Z = np.array([self.__source.z, ])


    def __repr__ (self):
        return f'Ray object'  # TODO: improve repr

    def propagate (self, calc_c=calc_c, calc_dz_c=calc_dz_c):
        if self.__is_propagated: return

        # Initialise
        angle = -1 * self.__angle_init + (np.pi / 2)  # convert angle notation
        x, z = self.__source.coordinates()
        dx_z = 1 / np.tan(angle)
        dxdx_z = 0  # no initial curvature
        # Initialise solver
        c0 = calc_c(0)
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

            # print(f'DEBUG: {x}, {z} => {x_new}, {z_new}')



            # Check backwards propagation (if not allowed): TODO


            # Check reflections > NOTE: intersections are calculated with a fineness of dx (no complex segment intersection algorithm followed by a solver)
            
            if self.__env.floor:
                fx, fxn = self.__env.floor(x), self.__env.floor(x_new)
                if z_new < fxn or (abs(z - fx) < 2 * self.dz_max and intersect(Point2D(x, z), Point2D(x_new, z_new), Point2D(x, fx), Point2D(x_new, fxn))):
                    x_new = float( self.func_solve( lambda x1: fx - dx_z * (x1 - x) - z, x0=x ))
                    z_new = self.__env.floor(x_new)
                    k = np.array([1., dx_z])  # Direction of incident ray
                    u = np.array([1., self.__env.dx_floor(x_new)])  # Direction of floor
                    n = np.array([-1*u[1], u[0]])  # Normal of floor, going up
                    l = np.dot(k, u)*u - np.dot(k, n)*n  # Direction of reflected ray
                    dx_z = np.dot(l, np.array([0, 1])) / np.dot(l, np.array([1, 0]))
        
            if self.__env.ceil:
                cx, cxn = self.__env.ceil(x), self.__env.ceil(x_new)
                if z_new > cxn or (abs(z - cx) < 2 * self.dz_max and intersect(Point2D(x, z), Point2D(x_new, z_new), Point2D(x, cx), Point2D(x_new, cxn))):
                    x_new = float( self.func_solve( lambda x1: cx - dx_z * (x1 - x) - z, x0=x ))
                    z_new = self.__env.ceil(x_new)
                    k = np.array([1., dx_z])  # Direction of incident ray
                    u = np.array([1., self.__env.dx_ceil(x_new)])  # Direction of floor
                    n = np.array([u[1], -1*u[0]])  # Normal of floor, going down
                    l = np.dot(k, u)*u - np.dot(k, n)*n  # Direction of reflected ray
                    dx_z = np.dot(l, np.array([0, 1])) / np.dot(l, np.array([1, 0]))

            # Check simulation bounds (!!!! USE NEW X, Z AFTER REFLECTION EVENT)
            if x_new < self.__env.range_min.x or x_new > self.__env.range_max.x:
                # TODO: plot last point at x = xmax or z = zmax
                print('DEBUG: out of bounds (x-axis)')
                break
            if z_new < self.__env.range_min.z or z_new > self.__env.range_max.z:
                # TODO: plot last point at x = xmax or z = zmax
                print('DEBUG: out of bounds (z-axis)')
                break

        

            # Add new point
            self.X = np.concatenate((self.X, np.array([x_new, ])), axis=0)
            self.Z = np.concatenate((self.Z, np.array([z_new, ])), axis=0)
            x = x_new
            z = z_new
            # Calculate new point's properties
            c = calc_c (z_new)
            dz_c = calc_dz_c (z_new)
            # Update derivatives for next integration segment
            dxdx_z = mult * dz_c / np.power(c, 3)
            dx_z += dxdx_z * dx

        # # Generate interpolated path function
        # self.Z_func = interpolate.interp1d(self.X, self.Z, kind='linear')

        self.__is_propagated = True