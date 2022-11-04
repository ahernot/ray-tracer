# 2D Ray Propagation Model (v3.1)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# TODO: generalise as calc_c(x, z) and derive using numpy => use matrices for calculations (set resolution) / use 2D functions for calculations
# TODO: Add min and max values into self. for when using interpolated function, generated at the same time as the latter
# TODO: increase pace when far enough away from borders
# TODO: add a verbose_indent kwargs field to indent verbose when called from Simulation2D with verbose enabled
# TODO: pregenerate lower-res functions such as calc_absorption in fixed-size arrays to approximate the values if res=low selected
# TODO: add an absorption_max criteria (either in mult or in dB)

# http://www.sengpielaudio.com/calculator-FactorRatioLevelDecibel.htm



import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate
from scipy.optimize import fsolve
from scipy.misc import derivative

from preferences import *
import physics
from physics.model_impedance import calc_Z
from physics.model_reflection import calc_refcoef_surface, calc_refcoef_sediment
from physics.profile_absorption import calc_absorption_dB
from environment import Environment2D




class Ray2D:

    def __init__ (self, env: Environment2D, source: np.ndarray, freq, angle, **kwargs):
        """
        Initialise ray
        
        :param env: Simulation environment
        :param source: Source point
        :param angle: Casting angle (from horizontal), in radians

        kwargs
        :param calc_der: Derivative calculation function
        :param func_solve: Equation solver
        """

        # Solver functions
        self.calc_der = kwargs.get('calc_der', derivative)
        self.func_solve = kwargs.get('func_solve', fsolve)

        self.__is_propagated = False

        self.env: Environment2D = env
        self.source: np.ndarray = source
        self.freq = freq
        self.angle_init = angle
        self.stop_reason = None

        # Initialise
        self.XZ = np.expand_dims(self.source.copy(), axis=0)
        self.T = np.array([0.,])
        self.steps = 1

        self.range_min = np.zeros(2)
        self.range_max = np.zeros(2)



    def __repr__ (self):
        return f'Ray of frequency {self.freq}'  # TODO: improve repr

    def plot (self, fig):
        if self.__is_propagated:
            plt.plot(self.XZ[:,0], self.XZ[:,1], figure=fig)

    def propagate (self, **kwargs):
        """
        Propagate ray

        kwargs
        :param backprop: Allow backpropagation (default=True)
        :param calc_c:
        :param calc_dz_c:
        :param dx_max:
        :param dz_max:
        :param n_steps_max:
        :param n_rebounds_max: Maximum number of rebounds (default=infinite)
        :param verbose: Verbose (default=False)
        """

        verbose = kwargs.get('verbose', False)
        if self.__is_propagated:
            if verbose: print('ERROR: Ray already propagated')
            return
        
        self.backprop = kwargs.get('backprop', True)
        self.n_steps_max = kwargs.get('n_steps_max', N_STEPS_MAX)
        self.n_rebounds_max = kwargs.get('n_rebounds_max', N_REBOUNDS_MAX)
        self.n_rebounds = 0
        self.__rebounds = list()

        # Minimum resolutions
        self.dx_max = kwargs.get('dx_max', DX_MAX_DEFAULT)
        self.dz_max = kwargs.get('dz_max', DZ_MAX_DEFAULT)
        dx_z_max = np.abs(self.dz_max / self.dx_max)  # Decision slope between both min resolutions


        # Initialise
        P = self.source
        angle = -1 * self.angle_init + (np.pi / 2)  # convert angle notation
        dx_z = 1 / np.tan(angle)
        dxdx_z = 0  # no initial curvature
        x_dir = 1.  # forwards propagation
        k = np.array([x_dir, dx_z])
        # Initialise solver
        c = self.env.penv.calc_c(0)
        mult = -1 * np.power(c / np.sin(angle), 2)  # differential equation multiplier


        for i in range(self.n_steps_max):

            # Enforce minimum resolution (k becomes this step's [dx, dz])
            if abs(dx_z) > dx_z_max: k *= self.dz_max / abs(dx_z)
            else: k *= self.dx_max

            # Unpack coordinates
            x, z = P
            x_new, z_new = P_new = P + k

            # Check floor rebounds
            if self.env.floor and z_new < self.env.floor(x_new):  # Calculate intersection point and new direction vector
                x_new = float( self.func_solve( lambda x1: self.env.floor(x) - dx_z * (x1 - x) - z, x0=x ))
                z_new = self.env.floor(x_new)
                P_new = np.array([x_new, z_new])
                u = np.array([1., self.env.dx_floor(x_new)])  # Direction of floor
                n = np.array([-1*u[1], u[0]])  # Normal of floor, going up
                k = np.dot(k, u)*u - np.dot(k, n)*n  # Direction of reflected ray

                # Calculate reflection coefficient
                wavelength = c / self.freq
                angle = 1 / ((1 + (k[1]/k[0]) ** 2) ** 0.5)
                refcoef = calc_refcoef_surface(wavelength=wavelength, angle=angle)
                self.__rebounds.append({'step': i+1, 'gain_dB': 10 * np.log10(refcoef), 'surface': 'ground'})
                self.n_rebounds += 1
                
                if self.n_rebounds_max > -1 and self.n_rebounds > self.n_rebounds_max:
                    self.XZ = np.insert(self.XZ, i+1, P_new, axis=0)  # Add final point
                    if verbose: print(f'DEBUG: Max number of rebounds reached ({self.n_rebounds_max})')
                    self.stop_reason = 'max-rebounds'
                    break
                if verbose: print(f'DEBUG: #{self.n_rebounds} - Floor rebound. New dir: {k}')

            # Check ceiling rebounds
            elif self.env.ceil and z_new > self.env.ceil(x_new):
                x_new = float( self.func_solve( lambda x1: self.env.ceil(x) - dx_z * (x1 - x) - z, x0=x ))
                z_new = self.env.ceil(x_new)
                P_new = np.array([x_new, z_new])
                u = np.array([1., self.env.dx_ceil(x_new)])  # Direction of ceiling
                n = np.array([u[1], -1*u[0]])  # Normal of ceiling, going down
                k = np.dot(k, u)*u - np.dot(k, n)*n  # Direction of reflected ray

                # Calculate reflection coefficient
                wavelength = c / self.freq
                angle = 1 / ((1 + (k[1]/k[0]) ** 2) ** 0.5)
                refcoef = calc_refcoef_sediment(wavelength=wavelength, Zp0=calc_Z(z))  # Uses c from previous iteration
                self.__rebounds.append({'step': i+1, 'gain_dB': 10 * np.log10(refcoef), 'surface': 'water-surface'})
                self.n_rebounds += 1

                if self.n_rebounds_max > -1 and self.n_rebounds > self.n_rebounds_max:
                    self.XZ = np.insert(self.XZ, i+1, P_new, axis=0)  # Add final point
                    if verbose: print(f'DEBUG: Max number of rebounds reached ({self.n_rebounds_max})')
                    self.stop_reason = 'max-rebounds'
                    break
                if verbose: print(f'DEBUG: #{self.n_rebounds} - Ceiling rebound. New dir: {k}')
            
            else:
                P_new = np.array([x_new, z_new])

            # Check simulation bounds
            if x_new < self.env.range_min[0]:
                x_new = self.env.range_min[0]
                z_new = -1 * dx_z * (x_new - x) + z  # Only hit when going left (x_dir = -1)
                self.XZ = np.insert(self.XZ, i+1, np.array([x_new, z_new]), axis=0)  # Add final point
                if verbose: print('DEBUG: Out of bounds (x-axis min)')
                self.stop_reason = 'exit-xmin'
                break 
            elif x_new > self.env.range_max[0]:
                x_new = self.env.range_max[0]
                z_new = dx_z * (x_new - x) + z  # Only hit when going right (x_dir = 1)
                self.XZ = np.insert(self.XZ, i+1, np.array([x_new, z_new]), axis=0)  # Add final point
                if verbose: print('DEBUG: Out of bounds (x-axis max)')
                self.stop_reason = 'exit-xmax'
                break
            elif z_new < self.env.range_min[1]:
                z_new = self.env.range_min[1]
                x_new = x_dir * (z_new - z) / dx_z + x
                self.XZ = np.insert(self.XZ, i+1, np.array([x_new, z_new]), axis=0)  # Add final point
                if verbose: print('DEBUG: Out of bounds (z-axis min)')
                self.stop_reason = 'exit-zmin'
                break
            elif z_new > self.env.range_max[1]:
                z_new = self.env.range_max[1]
                x_new = x_dir * (z_new - z) / dx_z + x
                self.XZ = np.insert(self.XZ, i+1, np.array([x_new, z_new]), axis=0)  # Add final point
                if verbose: print('DEBUG: Out of bounds (z-axis max)')
                self.stop_reason = 'exit-zmax'
                break


            # Add new point
            self.XZ = np.insert(self.XZ, i+1, P_new, axis=0)
            P = P_new.copy()

            # Calculate new point's properties
            c = self.env.penv.calc_c (z_new)
            dz_c = self.env.penv.calc_dz_c (z_new)
            # Unpack k
            x_dir = np.sign(k[0])
            dx_z = k[1] / k[0]
            # Update k for next integration segment
            dxdx_z = mult * dz_c / np.power(c, 3)
            dx_z += dxdx_z * k[0]
            k = np.array([x_dir, dx_z])


            # Check backpropagation
            if not self.backprop and x_dir < 0:
                if verbose: print('DEBUG: Backpropagation')
                self.stop_reason = 'backprop'
                break
            
            # Check number of steps (verbose only)
            if i == self.n_steps_max - 1:
                if verbose: print(f'DEBUG: Maximum iterations reached ({self.n_steps_max})')
                self.stop_reason = 'max-iter'

        # Count simulation steps
        self.steps = self.XZ.shape[0]

        # Update range
        self.range_min = np.array([np.min(self.XZ[:, 0]), np.min(self.XZ[:, 1])])
        self.range_max = np.array([np.max(self.XZ[:, 0]), np.max(self.XZ[:, 1])])

        # Calculate integration segments
        self.dL = np.linalg.norm(np.diff(self.XZ, axis=0), axis=1)  # dL at each arrival point (excluding initial point)
        self.L = np.cumsum(np.insert(self.dL, 0, 0.))
        self.C = self.env.penv.calc_c(self.XZ[:, 1])  # velocity at each point
        self.dT = self.dL / self.C[:-1]  # dT at each arrival point (excluding initial point)
        self.T = np.cumsum(np.insert(self.dT, 0, 0.))

        # Calculate gain
        dG_dB = -1 * calc_absorption_dB(self.freq, self.XZ[:-1, 1]) / 1000 * self.dL  # dG_dB for each dL (decibels)
        self.G_dB = np.cumsum(np.insert(dG_dB, 0, 0.))  # Cumulative gain for each point (decibels)
        for rebound in self.__rebounds:
            G_dB_add = np.zeros(self.steps, dtype=float)
            G_dB_add[rebound['step']+1:] = rebound['gain_dB']
            self.G_dB += G_dB_add
        self.Tmult = np.power(10, self.G_dB / 10)  # Cumulative transmittance multiplier for each point

        # Generate interpolated path function
        self.calc_z = interpolate.interp1d(self.XZ[:, 0], self.XZ[:, 1], kind='linear')
        
        self.__is_propagated = True
