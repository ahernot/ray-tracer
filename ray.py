# 2D Ray Propagation Model (v3.3)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# TODO: Increase pace when far enough away from borders
# TODO: Ray2D.reverse()
# TODO: RayPack2D.__repr__()
# IMPORTANT: Always have a margin of at least dx or dz between the simulation bound and the features (floor / ceiling)


import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate
from scipy.optimize import fsolve

from preferences import *
from physics.model_reflection import calc_refcoef_surface, calc_refcoef_sediment
from environment import Environment2D


# TODO: Remove target from Ray2D

class Ray2D:

    def __init__ (self, env: Environment2D, source: np.ndarray, angle, **kwargs):
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
        self.target: np.ndarray = kwargs.get('target', None)  # TODO: deprecate self.target globally
        self.angle = angle

        # Initialise
        self.XZ = np.expand_dims(self.source.copy(), axis=0)
        self.T = np.array([0.,])
        self.steps = 1  # TODO: rename to self.n_steps

        self.freqs = list()
        self.G = dict()
        self.Tmult = dict()  # TODO: deprecate
        # self.G_target = np.empty((0, 2))
        self.Tmult_target = np.empty((0, 2))  # TODO: deprecate

        self.range_min = np.zeros(2)
        self.range_max = np.zeros(2)

    def __repr__ (self):
        return f'Ray2D ({self.source}, {self.angle})'
        # repr_list = [
        #     f'{"Propagated " if self.__is_propagated else ""}{self.__class__.__name__} object cast from {self.source} {f"to target {self.target} " if self.target is not None else ""}at angle {self.angle} rad',
        #     f'\tStop reason: {self.stop_reason}' if self.stop_reason else '',
        #     f'\tDistance to target: {self.dist_to_target} m' if self.dist_to_target else ''
        # ]
        # return '\n'.join(repr_list)

    def plot (self, fig, **kwargs):  # TODO: Plot gain: self.L, self.G[freq]
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
        :param kwargs/verbose_depth: Verbose depth
        :param kwargs/verbose_depth_max: Verbose max depth
        """

        out_vals = ['exit-xmin', 'exit-xmax', 'exit-zmin', 'exit-zmax']
        def fit_to_bounds (x_prev, x_new, z_prev, z_new, dx_z):
            _out = 0

            # Check simulation horizontal bounds
            if x < self.env.range_min[0]:
                x_new = self.env.range_min[0]
                z_new = -1 * dx_z * (x_new - x_prev) + z_prev  # Only hit when going left (x_dir = -1)
                _out = 1
            elif x_new > self.env.range_max[0]:
                x_new = self.env.range_max[0]
                z_new = dx_z * (x_new - x_prev) + z_prev  # Only hit when going right (x_dir = 1)
                _out = 2
            
            # Check simulation vertical bounds
            if z_new < self.env.range_min[1]:
                z_new = self.env.range_min[1]
                x_new = x_dir * (z_new - z_prev) / dx_z + x_prev
                _out = 3
            elif z_new > self.env.range_max[1]:
                z_new = self.env.range_max[1]
                x_new = x_dir * (z_new - z_prev) / dx_z + x_prev
                _out = 4
            
            return x_new, z_new, _out
        

        # Verbose param
        self.verbose = kwargs.get('verbose', False)
        self.verbose_depth_max = kwargs.get('verbose_depth_max', 0)
        self.verbose_depth =  kwargs.get('verbose_depth', 0)
        if self.verbose: self.verbose_next = (self.verbose_depth < self.verbose_depth_max) or (self.verbose_depth_max == -1)
        self.__vi = '\t' * self.verbose_depth + self.__class__.__name__ + ': '

        # Exit if already completed
        if self.__is_propagated:
            if self.verbose: print(f'{self.__vi}ERROR: Ray already propagated')  # TODO: Warning
            return
        
        # Initialize variables (kwargs)
        self.backprop = kwargs.get('backprop', True)  # TODO: Rename to allow_backprop
        self.n_steps_max = kwargs.get('n_steps_max', N_STEPS_MAX)
        self.n_rebounds_max = kwargs.get('n_rebounds_max', N_REBOUNDS_MAX)
        self.n_rebounds = 0
        self.__rebounds = list()

        # Minimum resolutions (maximum steps)
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

        # Initialise solver parameters
        c = self.env.penv.calc_c(0)  # TODO: 0? or z0?
        mult = -1 * np.power(c / np.sin(angle), 2)  # differential equation multiplier



        # Main loop
        for i in range(self.n_steps_max):
            # Enforce minimum resolution (k becomes this step's [dx, dz])
            if abs(dx_z) > dx_z_max: k *= self.dz_max / abs(dx_z)
            else: k *= self.dx_max
            # Unpack coordinates
            x, z = P
            x_new, z_new = P_new = P + k



            # Check simulation horizontal bounds (for horizotnally-interpolated functions)
            x_new_temp, z_new_temp, out = fit_to_bounds (x, x_new, z, z_new, dx_z)
            if out in (1, 2):
                x_new, z_new = x_new_temp, z_new_temp
                self.stop_reason = out_val = out_vals[out-1]
                if self.verbose: print(f'{self.__vi}Stopping: Out of bounds ({out_val})')
                break



            # Project onto floor & ceiling (for rebound check)
            z_floor_projection = self.env.floor(x_new)
            z_ceil_projection = self.env.ceil(x_new)

            # Check floor rebounds
            if self.env.floor and z_new < z_floor_projection:
                x_new = float( fsolve( lambda x1: self.env.floor(x) - dx_z * (x1 - x) - z, x0=x ))
                z_new = z_floor_projection
                u = np.array([1., self.env.dx_floor(x_new)])  # Direction of floor
                n = np.array([-1*u[1], u[0]])  # Normal of floor, going up
                k = np.dot(k, u)*u - np.dot(k, n)*n  # Direction of reflected ray

                # Log rebound
                angle = 1 / ((1 + (k[1]/k[0]) ** 2) ** 0.5)
                self.__rebounds.append({'step': i+1, 'surface': 'sediment', 'velocity': c, 'angle': angle, 'position': P})
                self.n_rebounds += 1
                
                # Stop if last allowed rebound
                if self.n_rebounds_max > -1 and self.n_rebounds > self.n_rebounds_max:
                    x_new, z_new, out = fit_to_bounds (x, x_new, z, z_new, dx_z)  # Check bounds
                    P_new = np.array([x_new, z_new])
                    self.XZ = np.insert(self.XZ, i+1, P_new, axis=0)  # Add final point
                    if self.verbose: print(f'{self.__vi}Stopping: Max number of rebounds reached ({self.n_rebounds_max})')
                    self.stop_reason = 'max-rebounds'
                    break
                if self.verbose: print(f'{self.__vi}Rebound #{self.n_rebounds} - floor. New dir: {k}')

            # Check ceiling rebounds
            elif self.env.ceil and z_new > z_ceil_projection:
                x_new = float( fsolve( lambda x1: self.env.ceil(x) - dx_z * (x1 - x) - z, x0=x ))
                z_new = z_ceil_projection
                u = np.array([1., self.env.dx_ceil(x_new)])  # Direction of ceiling
                n = np.array([u[1], -1*u[0]])  # Normal of ceiling, going down
                k = np.dot(k, u)*u - np.dot(k, n)*n  # Direction of reflected ray

                # Log rebound
                angle = 1 / ((1 + (k[1]/k[0]) ** 2) ** 0.5)  # TODO: is this the right angle?
                self.__rebounds.append({'step': i+1, 'surface': 'water-surface', 'velocity': c, 'angle': angle, 'position': P})
                # TODO: Velocity stored as a 0-d array??
                self.n_rebounds += 1

                # Stop if last allowed rebound
                if self.n_rebounds_max > -1 and self.n_rebounds > self.n_rebounds_max:
                    x_new, z_new, out = fit_to_bounds (x, x_new, z, z_new, dx_z)  # Check bounds
                    P_new = np.array([x_new, z_new])
                    self.XZ = np.insert(self.XZ, i+1, P_new, axis=0)  # Add final point
                    if self.verbose: print(f'{self.__vi}Stopping: Max number of rebounds reached ({self.n_rebounds_max})')
                    self.stop_reason = 'max-rebounds'
                    break
                if self.verbose: print(f'{self.__vi}Rebound #{self.n_rebounds} - ceiling. New dir: {k}')



            # Check simulation bounds (all)
            x_new_temp, z_new_temp, out = fit_to_bounds (x, x_new, z, z_new, dx_z)
            if out:
                x_new, z_new = x_new_temp, z_new_temp
                self.stop_reason = out_val = out_vals[out-1]
                if self.verbose: print(f'{self.__vi}Stopping: Out of bounds ({out_val})')
                break



            # Add new point
            P_new = np.array([x_new, z_new])
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
                self.stop_reason = 'backprop'
                if self.verbose: print(f'{self.__vi}Stopping: Backpropagation')
                break
            # Check number of steps
            if i == self.n_steps_max - 1:
                self.stop_reason = 'max-iter'
                if self.verbose: print(f'{self.__vi}Stopping: Maximum iterations reached ({self.n_steps_max})')
                


        # Count simulation steps
        self.steps = self.XZ.shape[0]

        # Populate ray's effective range
        self.range_min = np.array([np.min(self.XZ[:, 0]), np.min(self.XZ[:, 1])])
        self.range_max = np.array([np.max(self.XZ[:, 0]), np.max(self.XZ[:, 1])])

        # Calculate integration segments
        self.dL = np.linalg.norm(np.diff(self.XZ, axis=0), axis=1)  # dL at each arrival point (excluding initial point)
        self.L = np.cumsum(np.insert(self.dL, 0, 0.))
        self.C = self.env.penv.calc_c(self.XZ[:, 1])  # Velocity at each point
        self.dT = self.dL / self.C[:-1]  # dT at each arrival point (excluding initial point)
        self.T = np.cumsum(np.insert(self.dT, 0, 0.))  # Time from ray launch at each point

        # Calculate time at target  # calculate step closest to target (on the x-axis)  # TODO: DISTANCE ON BOTH AXES (linalg.norm)
        self.step_target = np.argmin(np.abs(self.XZ[:, 0] - self.target[0])) if self.target is not None else None  # Step closest to target on the x-axis
        self.T_target = self.T[self.step_target] if self.target is not None else None

        # Generate interpolated path function
        self.calc_z = interpolate.interp1d(self.XZ[:, 0], self.XZ[:, 1], kind='linear', bounds_error=True)

        # Calculate distance to target  # TODO: better check interpolation range (if range_x < target.x) and find target step
        if self.target is not None:
            try: self.dist_to_target = abs(self.target[1] - self.calc_z(self.target[0]))
            except: self.dist_to_target = np.nan
        
        self.__is_propagated = True

    def populate (self, *freqs):

        for freq in freqs:
            # Check if freq already generated
            if freq in self.freqs: continue
            
            # Calculate base gain
            dG = -1 * self.env.penv.calc_dl_dG(freq, self.XZ[:-1, 1]) / 1000 * self.dL  # dG for each dL (decibels)
            G = np.cumsum(np.insert(dG, 0, 0.))  # Cumulative gain for each point (decibels)

            # Process rebounds
            for rebound in self.__rebounds:

                wavelength = rebound['velocity'] / freq
                if rebound['surface'] == 'sediment':
                    z = rebound['position'][1]
                    refcoef = calc_refcoef_sediment(wavelength=wavelength, Zp0=self.env.penv.calc_Z(z))
                elif rebound['surface'] == 'water-surface':
                    angle = rebound['angle']
                    refcoef = calc_refcoef_surface(wavelength=wavelength, angle=angle)

                G_add = np.zeros(self.steps, dtype=float)
                G_add[rebound['step']+1:] = 10 * np.log10(refcoef)
                G += G_add
            
            # Cumulative transmittance multiplier for each point
            Tmult = np.power(10, G / 10)

            # Save values
            self.G[freq] = G
            self.Tmult[freq] = Tmult
            self.freqs.append(freq)

            # Save target values  # TODO: move this to simulation? or raypack?
            if (self.target is not None) and (len(G) > self.step_target):
                # self.G_target = np.insert(self.G_target, 0, np.array([freq, G[self.step_target]]), axis=0)
                self.Tmult_target = np.insert(self.Tmult_target, 0, np.array([freq, Tmult[self.step_target]]), axis=0)
        
        # Sort target arrays
        # self.G_target = np.sort(self.G_target, axis=0)
        self.Tmult_target = np.sort(self.Tmult_target, axis=0)        
        self.calc_Tmult_target = interpolate.interp1d(
            np.concatenate( (-1 * self.Tmult_target[::-1, 0], self.Tmult_target[:, 0]) ),
            np.concatenate( (self.Tmult_target[::-1, 1], self.Tmult_target[:, 1]) ),
            kind='linear',
            bounds_error=True
        )  # TODO: currently only on x-coordinate of the target (dz vertical difference not accounted for)

    def populate_NEW (self, *freqs):

        for freq in freqs:
            # Check if freq already generated
            if freq in self.freqs: continue
            
            # Calculate base gain for each point
            dG = -1 * self.env.penv.calc_dl_dG(freq, self.XZ[:-1, 1]) / 1000 * self.dL  # dG for each dL (decibels)
            G = np.cumsum(np.insert(dG, 0, 0.))  # Cumulative gain for each point (decibels)

            # Process rebounds
            for rebound in self.__rebounds:

                wavelength = rebound['velocity'] / freq
                if rebound['surface'] == 'sediment':
                    z = rebound['position'][1]
                    refcoef = calc_refcoef_sediment(wavelength=wavelength, Zp0=self.env.penv.calc_Z(z))
                elif rebound['surface'] == 'water-surface':
                    angle = rebound['angle']
                    refcoef = calc_refcoef_surface(wavelength=wavelength, angle=angle)

                G_add = np.zeros(self.steps, dtype=float)
                G_add[rebound['step']+1:] = 10 * np.log10(refcoef)
                G += G_add

            # Save values
            self.G[freq] = G
            self.freqs.append(freq)
        
        # Sort frequencies
        self.freqs.sort()


    def closest_step_to_point (self, target: np.ndarray):  # Calculate closest approach
        d_target = np.sqrt(np.power(target[0]-self.XZ[:, 0], 2) + np.power(target[1]-self.XZ[:, 1], 2))
        return np.argmin(d_target)
    
    def point (self, step: int) -> np.ndarray:  # Point coordinates at step
        return self.XZ[step]

    def gain (self, step: int) -> np.ndarray:  # Gain array (for each populated freq) at step
        return np.array([ self.G[freq][step] for freq in self.freqs ])  # NOTE: self.freqs is sorted


    # def dist_to_point (self, target: np.ndarray): 
    #     # Calculate shortest distance to a point on the ray's path

    #     # TODO: needs to calculate for each phase of the path (between 2 rebounds), then choose point with largest amplitude at 1Hz
    #     # steps_closest_list = list()
    #     # step_range_start = 0
    #     # for rebound in self.__rebounds:
    #     #     step_range_stop = rebound['step']
    #     #     step_closest = np.argmin(d_target[step_range_start:step_range_stop]) + step_range_start
    #     #     steps_closest_list.append(step_closest)
    #     #     step_range_start = step_range_stop
    #     # # Add from last rebound to last step
    #     # steps_closest_list.append (np.argmin(d_target[step_range_start:]) + step_range_start)
    #     # print(steps_closest_list)
    #     # for s in steps_closest_list:
    #     #     print(self.XZ[s])
        
    #     step = self.closest_step_to_point(target)
    #     return np.linalg.norm(target - self.XZ[step])




class RayPack2D:

    def __init__ (self):
        self.rays = list()  # Raw list of rays (stored as references to Ray2D objects)
        self.stop_reasons = dict()  # Dictionary of Ray2D objects indexed by stop reason
        self.angles = dict()  # Dictionary of Ray2D objects indexed by angle
        self.dist = dict()  # Dictionary of Ray2D objects indexed by distance to target
        self.dist_sorted = None
        self.freqs = list()
        self.n_rays, self.n_angles, self.n_freqs = 0, 0, 0

        self.range_min_plot = np.zeros(2)
        self.range_max_plot = np.zeros(2)
        self.size_plot = np.zeros(2)

    def __repr__ (self):
        repr_list = [
            f'Raypack containing {self.n_rays} rays cast with {self.n_angles} angles'
        ]
        return '\n'.join(repr_list)

    def add (self, *rays: Ray2D):

        for ray in rays:
            if ray.angle in self.angles: continue  # Skip already generated rays
            # Update counters
            self.n_rays += 1
            self.n_angles += 0 if ray.angle in self.angles else 1

            # Add ray to database
            self.rays .append(ray)
            self.angles[ray.angle] = ray  # Unicity of ray per angle
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

    def populate (self, *freqs):
        
        for freq in freqs:

            # Check if freq already generated
            if freq in self.freqs: continue

            # Populate rays
            for ray in self.rays:
                ray.populate(*freqs)

            self.freqs .append(freq)
        self.n_freqs = len(self.freqs)
