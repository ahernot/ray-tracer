# 2D Ray Propagation Model (v3.1)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# TODO: increase pace when far enough away from borders
# TODO: add a verbose_indent kwargs field to indent verbose when called from Simulation2D with verbose enabled
# TODO: add an absorption_max criteria (either in mult or in dB)
# TODO: decouple ray path and frequency: calc_absorption(*freqs) method

# http://www.sengpielaudio.com/calculator-FactorRatioLevelDecibel.htm



import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate
from scipy.optimize import fsolve

from preferences import *
from physics.model_reflection import calc_refcoef_surface, calc_refcoef_sediment
from environment import Environment2D




class Ray2D:

    def __init__ (self, env: Environment2D, source: np.ndarray, freq, angle, **kwargs):
        """
        Initialise ray
        
        :param env: Simulation environment
        :param source: Source point
        :param angle: Casting angle (from horizontal), in radians
        :param kwargs/target: Target point
        """

        self.__is_propagated = False
        self.stop_reason = None
        self.dist_to_target = None

        self.env: Environment2D = env
        self.source: np.ndarray = source
        self.target: np.ndarray = kwargs.get('target', None)
        self.freq = freq
        self.angle = angle

        # Initialise
        self.XZ = np.expand_dims(self.source.copy(), axis=0)
        self.T = np.array([0.,])
        self.steps = 1

        self.range_min = np.zeros(2)
        self.range_max = np.zeros(2)



    def __repr__ (self):
        return f'Ray of frequency {self.freq}'  # TODO: improve repr

    def plot (self, fig, **kwargs):
        if self.__is_propagated:
            plt.plot(self.XZ[:,0], self.XZ[:,1], figure=fig, **kwargs)

    def propagate (self, **kwargs):
        """
        Propagate ray
        :param kwargs/backprop: Allow backpropagation (default=True)
        :param kwargs/dx_max:
        :param kwargs/dz_max:
        :param kwargs/n_steps_max:
        :param kwargs/n_rebounds_max: Maximum number of rebounds (default=infinite)
        :param kwargs/verbose: Verbose (default=False)
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
        angle = -1 * self.angle + (np.pi / 2)  # convert angle notation
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


            # Check simulation horizontal bounds
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


            # Check floor rebounds
            if self.env.floor and z_new < self.env.floor(x_new):  # Calculate intersection point and new direction vector
                x_new = float( fsolve( lambda x1: self.env.floor(x) - dx_z * (x1 - x) - z, x0=x ))
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
                x_new = float( fsolve( lambda x1: self.env.ceil(x) - dx_z * (x1 - x) - z, x0=x ))
                z_new = self.env.ceil(x_new)
                P_new = np.array([x_new, z_new])
                u = np.array([1., self.env.dx_ceil(x_new)])  # Direction of ceiling
                n = np.array([u[1], -1*u[0]])  # Normal of ceiling, going down
                k = np.dot(k, u)*u - np.dot(k, n)*n  # Direction of reflected ray

                # Calculate reflection coefficient
                wavelength = c / self.freq
                angle = 1 / ((1 + (k[1]/k[0]) ** 2) ** 0.5)
                refcoef = calc_refcoef_sediment(wavelength=wavelength, Zp0=self.env.penv.calc_Z(z))  # Uses c from previous iteration
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


            # Check simulation vertical bounds
            if z_new < self.env.range_min[1]:
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
        dG = -1 * self.env.penv.calc_dz_dG(self.freq, self.XZ[:-1, 1]) / 1000 * self.dL  # dG_dB for each dL (decibels)
        self.G = np.cumsum(np.insert(dG, 0, 0.))  # Cumulative gain for each point (decibels)
        for rebound in self.__rebounds:
            G_add = np.zeros(self.steps, dtype=float)
            G_add[rebound['step']+1:] = rebound['gain_dB']
            self.G += G_add
        self.Tmult = np.power(10, self.G / 10)  # Cumulative transmittance multiplier for each point

        # Generate interpolated path function
        self.calc_z = interpolate.interp1d(self.XZ[:, 0], self.XZ[:, 1], kind='linear')

        # Calculate distance to target
        if self.target is not None:
            try: self.dist_to_target = abs(self.target[1] - self.calc_z(self.target[0]))
            except: self.dist_to_target = np.nan
        
        self.__is_propagated = True
    
    def reverse (self):
        pass  # TODO: Ray2D.reverse




# TODO: add __repr__

class RayPack2D:

    def __init__ (self):
        self.rays = list()  # Raw list of rays (stored as references to Ray2D objects)
        self.angles = dict()  # Dictionary of Ray2D objects indexed by angle
        self.freqs = dict()  # Dictionary of Ray2D objects indexed by frequency
        self.stop_reasons = dict()  # Dictionary of Ray2D objects indexed by stop reason
        self.dist = dict()
        self.angles_sorted = None
        self.dist_sorted = None
        self.freqs_sorted = None

        self.range_min_plot = np.zeros(2)
        self.range_max_plot = np.zeros(2)
        self.size_plot = np.zeros(2)

        self.n_rays = 0
        self.n_angles = 0
        self.n_freqs = 0

        self.spectrum_vals = dict()  # Samples of self.spectrum non-normalised distribution function
        self.spectrum_total = None
        self.spectrum_distrib = dict()  # xf

        self.energy_total = 1.
        self.__energy_norm = None
        self.ray_energy = None  # Dict of ray normalised energy unit per frequency

    def add (self, *rays: Ray2D):

        for ray in rays:
            # Update counters
            self.n_rays += 1
            self.n_angles += 0 if ray.angle in self.angles else 1
            self.n_freqs += 0 if ray.freq in self.freqs else 1

            # Add ray to database
            self.rays .append(ray)
            self.angles[ray.angle] = self.angles[ray.angle] + [ray] if ray.angle in self.angles else [ray]
            self.freqs[ray.freq] = self.freqs[ray.freq] + [ray] if ray.freq in self.freqs else [ray]
            self.stop_reasons[ray.stop_reason] = self.stop_reasons[ray.stop_reason] + [ray] if ray.stop_reason in self.stop_reasons else [ray]
            self.dist[ray.dist_to_target] = self.dist[ray.dist_to_target] + [ray] if ray.dist_to_target in self.dist else [ray]

            # Update plot range
            if ray.range_min[0] < self.range_min_plot[0]: self.range_min_plot[0] = ray.range_min[0]
            if ray.range_min[1] < self.range_min_plot[1]: self.range_min_plot[1] = ray.range_min[1]
            if ray.range_max[0] > self.range_max_plot[0]: self.range_max_plot[0] = ray.range_max[0]
            if ray.range_max[1] > self.range_max_plot[1]: self.range_max_plot[1] = ray.range_max[1]

        # Update plot size
        self.size_plot = self.range_max_plot - self.range_min_plot

        # Update sorted arrays
        self.dist_sorted = np.sort(list(self.dist.keys()))
        self.angles_sorted = np.sort(list(self.angles.keys()))
        self.freqs_sorted = np.sort(list(self.freqs.keys()))

    def regen_energy (self):
        # Regenerate normalised ray energy unit
        self.__energy_norm = self.energy_total / np.sum([len(self.freqs[freq]) * self.spectrum_distrib[freq] for freq in self.freqs])
        self.ray_energy = {freq: self.spectrum_distrib[freq] * self.__energy_norm for freq in self.freqs}
