# 2D Ray Propagation Model (v3.0)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np

from scipy import interpolate  # TODO
from scipy.optimize import fsolve
from scipy.misc import derivative  # TODO

from physics.absorption_model import calc_absorption_dB
from physics.profile_velocity import calc_c, calc_dz_c
from physics.profile_salinity import calc_S
from physics.profile_temperature import calc_T
from physics.profile_ph import calc_pH
from physics.surface_reflection_model import calc_refcoef_surface
from physics.reflection import reflection_coefficient
from environment import Environment2D


DX_MAX_DEFAULT = 10.
DZ_MAX_DEFAULT = 0.1
N_STEPS_MAX = 100000
N_REBOUNDS_MAX = -1  # -1=infinity


# TODO New version: each ray is subdivided into segments and with each rebound new subclasses are created as a tree (can have multiple branches per rebound)
# TODO: generalise as calc_c(x, z) and derive using numpy => use matrices for calculations (set resolution) / use 2D functions for calculations

class Ray2D:

    def __init__ (self, env: Environment2D, source: np.ndarray, angle, **kwargs):
        """
        :param env: Simulation environment
        :param source: Source point
        :param angle: Casting angle (from horizontal), in radians
        """

        # Solver functions
        self.calc_der = kwargs.get('calc_der', derivative)
        self.func_solve = kwargs.get('func_solve', fsolve)

        self.__source: np.ndarray = source
        self.__env: Environment2D = env
        self.__angle_init = angle
        self.__is_propagated = False

        # Initialise
        self.XZ = np.expand_dims(self.__source.copy(), axis=0)


    def __repr__ (self):
        return f'Ray object'  # TODO: improve repr

    def propagate (self, calc_c=calc_c, calc_dz_c=calc_dz_c, **kwargs):
        verbose = kwargs.get('verbose', False)
        if self.__is_propagated:
            if verbose: print('ERROR: Ray already propagated')
            return
        
        self.backprop = kwargs.get('backprop', True)
        self.n_steps_max = kwargs.get('n_steps_max', N_STEPS_MAX)
        self.n_rebounds_max = kwargs.get('n_rebounds_max', N_REBOUNDS_MAX)
        self.n_rebounds = 0

        # Minimum resolutions
        self.dx_max = kwargs.get('dx_max', DX_MAX_DEFAULT)
        self.dz_max = kwargs.get('dz_max', DZ_MAX_DEFAULT)
        dx_z_max = np.abs(self.dz_max / self.dx_max)  # Decision slope between both min resolutions



        # Initialise
        P = self.__source
        angle = -1 * self.__angle_init + (np.pi / 2)  # convert angle notation
        dx_z = 1 / np.tan(angle)
        dxdx_z = 0  # no initial curvature
        x_dir = 1.  # forwards propagation
        k = np.array([x_dir, dx_z])
        # Initialise solver
        c0 = calc_c(0)
        mult = -1 * np.power(c0 / np.sin(angle), 2)  # differential equation multiplier



        for i in range(self.n_steps_max):

            # Enforce minimum resolution (k becomes this step's [dx, dz])
            if abs(dx_z) > dx_z_max: k *= self.dz_max / abs(dx_z)
            else: k *= self.dx_max

            # Unpack coordinates
            x, z = P
            x_new, z_new = P + k

            # print(f'DEBUG: {x}, {z} => {x_new}, {z_new}')
            
            if self.__env.floor and z_new < self.__env.floor(x_new):  # Calculate intersection point and new direction vector
                x_new = float( self.func_solve( lambda x1: self.__env.floor(x) - dx_z * (x1 - x) - z, x0=x ))
                z_new = self.__env.floor(x_new)
                P = np.array([x_new, z_new])
                u = np.array([1., self.__env.dx_floor(x_new)])  # Direction of floor
                n = np.array([-1*u[1], u[0]])  # Normal of floor, going up
                k = np.dot(k, u)*u - np.dot(k, n)*n  # Direction of reflected ray
                self.n_rebounds += 1
                if self.n_rebounds_max > -1 and self.n_rebounds > self.n_rebounds_max:
                    if verbose: print(f'DEBUG: Max number of rebounds reached ({self.n_rebounds_max})')
                    self.XZ = np.insert(self.XZ, i+1, P, axis=0)  # Add final point
                    break
                if verbose: print(f'DEBUG: #{self.n_rebounds} - Ground rebound. New dir: {k}')
        
            elif self.__env.ceil and z_new > self.__env.ceil(x_new):
                x_new = float( self.func_solve( lambda x1: self.__env.ceil(x) - dx_z * (x1 - x) - z, x0=x ))
                z_new = self.__env.ceil(x_new)
                P = np.array([x_new, z_new])
                u = np.array([1., self.__env.dx_ceil(x_new)])  # Direction of floor
                n = np.array([u[1], -1*u[0]])  # Normal of floor, going down
                k = np.dot(k, u)*u - np.dot(k, n)*n  # Direction of reflected ray
                self.n_rebounds += 1
                if self.n_rebounds_max > -1 and self.n_rebounds > self.n_rebounds_max:
                    if verbose: print(f'DEBUG: Max number of rebounds reached ({self.n_rebounds_max})')
                    self.XZ = np.insert(self.XZ, i+1, P, axis=0)  # Add final point
                    break
                if verbose: print(f'DEBUG: #{self.n_rebounds} - Surface rebound. New dir: {k}')
            
            else:
                P = np.array([x_new, z_new])

            # Check simulation bounds
            if x_new < self.__env.range_min.x:
                x_new = self.__env.range_min.x
                z_new = -1 * dx_z * (x_new - x) + z  # Only hit when going left (x_dir = -1)
                self.XZ = np.insert(self.XZ, i+1, np.array([x_new, z_new]), axis=0)  # Add final point
                if verbose: print('DEBUG: Out of bounds (x-axis min)')
                break 
            elif x_new > self.__env.range_max.x:
                x_new = self.__env.range_min.x
                z_new = dx_z * (x_new - x) + z  # Only hit when going right (x_dir = 1)
                self.XZ = np.insert(self.XZ, i+1, np.array([x_new, z_new]), axis=0)  # Add final point
                if verbose: print('DEBUG: Out of bounds (x-axis max)')
                break
            elif z_new < self.__env.range_min.z:
                z_new = self.__env.range_min.z
                x_new = x_dir * (z_new - z) / dx_z + x
                self.XZ = np.insert(self.XZ, i+1, np.array([x_new, z_new]), axis=0)  # Add final point
                if verbose: print('DEBUG: Out of bounds (z-axis min)')
                break
            elif z_new > self.__env.range_max.z:
                z_new = self.__env.range_max.z
                x_new = x_dir * (z_new - z) / dx_z + x
                self.XZ = np.insert(self.XZ, i+1, np.array([x_new, z_new]), axis=0)  # Add final point
                if verbose: print('DEBUG: Out of bounds (z-axis max)')
                break
            
            
            # Add new point
            self.XZ = np.insert(self.XZ, i+1, P, axis=0)
            # Calculate new point's properties
            c = calc_c (z_new)
            dz_c = calc_dz_c (z_new)
            # Unpack k
            x_dir = np.sign(k[0])
            dx_z = k[1] / k[0]
            # Update k for next integration segment
            dxdx_z = mult * dz_c / np.power(c, 3)
            dx_z += dxdx_z * k[0]
            k = np.array([x_dir, dx_z])

            if not self.backprop and x_dir < 0:
                if verbose: print('DEBUG: Backpropagation')
                break

            if i == self.n_steps_max - 1:
                if verbose: print(f'DEBUG: Maximum iterations reached ({self.n_steps_max})')

        # # Generate interpolated path function
        # self.Z_func = interpolate.interp1d(self.X, self.Z, kind='linear')

        self.__is_propagated = True
