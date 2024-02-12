import numpy as np
import matplotlib.pyplot as plt
from typing import Union

from scipy import interpolate
from scipy.optimize import fsolve

from preferences import *
from physics.model_reflection import calc_refcoef_surface, calc_refcoef_sediment
from environment import Environment2D



class Ray2D:

    def __init__ (self, env: Environment2D, source: np.ndarray, angle: Union[float, int], **kwargs):  # TODO: kwargs?
        self.env: Environment2D = env  # Reference pass
        self.source: np.ndarray = source  # TODO: cast to ndarray if type mismatch (needed for P in self.propagate)
        self.angle: Union[float, int] = angle  # Angle between the leading vector k and x (angle to the horizontal, signed)

        self.data = dict()  # TODO: store G, T, etc in a big dictionary
        self.XZ = np.expand_dims(self.source.copy(), axis=0)  # TODO: rename to self.xz?


    def propagate (self, **kwargs):  # TODO: propagate with conditions (no backpropagation, …)

        # Parse kwargs
        n_steps_max    = kwargs.get('n_steps_max', N_STEPS_MAX)
        n_rebounds_max = kwargs.get('n_rebounds_max', N_REBOUNDS_MAX)
        allow_backprop = kwargs.get('allow_backprop', False)
        dx_max         = kwargs.get('dx_max', DX_MAX_DEFAULT)
        dz_max         = kwargs.get('dz_max', DZ_MAX_DEFAULT)
        verbose        = kwargs.get('verbose', False)
        vi = ''  # TODO

        # Calculate decision slope between min resolutions on both axes
        dx_z_max = np.abs(dz_max / dx_max)  # Decision slope between both min resolutions




        # # Calculate initial parameters
        # P = self.source
        # angle = -1 * self.angle + (np.pi / 2)  # convert angle notation  # TODO: don't? cf equation
        # dx_z = 1 / np.tan(angle)  # Initial slope
        # dxdx_z = 0  # Initial curvature
        # x_dir = 1.  # forwards propagation  # TODO: needed?
        # k = np.array([x_dir, dx_z])


        # # Initialize propagation equation parameters
        # c = self.env.penv.calc_c(0)  # TODO: 0? or z0? => z0 right?
        # mult = -1 * np.power(c / np.sin(angle), 2)  # Differential equation multiplier


        
        # TODO: self.stop_reason?


        # TODO: preset unit vectors
        ex, ey = np.identity(2) .astype(np.float64)

        # Compute initial leading vector
        dx_z = np.tan(self.angle)  # Initial slope
        k = np.array([1., dx_z])
        k /= np.linalg.norm(k)  # TODO: better way of doing this? with scalar product or something
        # TODO: x_dir and z_dir are the sign of kx and kz

        # Compute propagation constant
        c = self.env.penv.calc_c(self.source[1])
        const = np.cos(self.angle) / c

        xz = self.source .astype(np.float64)

        # Initialize counters
        n_rebounds = 0
        for i in range(n_steps_max):

            if verbose: print(xz)

            # Compute new leading vector
            c = self.env.penv.calc_c(xz[1])
            kx = const * c  # Equation solving of np.dot(k, ex) / c
            kz = -1 * np.sqrt(1 - np.power(kx, 2))  # TODO: sign (rebound => flip the sign of kz or kx based on angle of bounce surface)
            k = np.array([kx, kz])  # TODO: where is the sign gone?

            # TODO: dx_max, dz_max => apply dynamic resolution (try to get largest dl every time)

            # TODO: compute new point
            xz += k
            self.XZ = np.insert(self.XZ, i+1, xz, axis=0)






        # # Propagation loop
        # for i in range(n_steps_max):

        #     # Enforce minimum resolution (k becomes this step's [dx, dz])
        #     if abs(dx_z) > dx_z_max: k *= dz_max / abs(dx_z)
        #     else: k *= dx_max

        #     # Unpack coordinates
        #     x, z = P
        #     x_new, z_new = P_new = P + k
            
        #     # Check simulation horizontal bounds (for horizontally-interpolated functions??????????what does that mean)  # TODO: define exit function with verbose and break elsewhere?
        #     x_new_temp, z_new_temp, out = self.env.fit_to_bounds (x, x_new, z, z_new, dx_z)
        #     if out in ('exit-xmin', 'exit-xmax'):
        #         x_new, z_new = x_new_temp, z_new_temp
        #         if verbose: print(f'{vi}Stopping: Out of bounds ({out}) at 0')
        #         break

        #     # Check for rebounds
        #     # TODO

        #     # Check simulation bounds (all)
        #     x_new_temp, z_new_temp, out = self.env.fit_to_bounds (x, x_new, z, z_new, dx_z)
        #     if out:
        #         x_new, z_new = x_new_temp, z_new_temp  # TODO: add to self.XZ?
        #         if verbose: print(f'{vi}Stopping: Out of bounds ({out}) at 1')
        #         break

        #     # Add new point
        #     P_new = np.array([x_new, z_new])
        #     self.XZ = np.insert(self.XZ, i+1, P_new, axis=0)
        #     P = P_new.copy()

        #     # Calculate new point's properties
        #     c = self.env.penv.calc_c (z_new)
        #     dz_c = self.env.penv.calc_dz_c (z_new)
        
        #     # Unpack k
        #     x_dir = np.sign(k[0])
        #     dx_z = k[1] / k[0]
        
        #     # Update k for next integration segment
        #     dxdx_z = mult * dz_c / np.power(c, 3)
        #     dx_z += dxdx_z * k[0]
        #     k = np.array([x_dir, dx_z])  # Non-normalized yet (min resolution enforced at start of loop)

        #     # Check backpropagation
        #     if not allow_backprop and x_dir < 0:  # TODO: x_dir needed? or just do x_new - x
        #         stop_reason = 'backprop'  # TODO
        #         if verbose: print(f'{vi}Stopping: Backpropagation')
        #         break
            
        #     # Check number of steps
        #     if i == n_steps_max - 1:
        #         self.stop_reason = 'max-iter'
        #         if verbose: print(f'{vi}Stopping: Maximum iterations reached ({self.n_steps_max})')


        

        # #! Data aggregation
        # self.n_steps = self.XZ.shape[0]
        # self.range_min = np.array([np.min(self.XZ[:, 0]), np.min(self.XZ[:, 1])])  # Populate ray's effective range
        # self.range_max = np.array([np.max(self.XZ[:, 0]), np.max(self.XZ[:, 1])])  # Populate ray's effective range

        # # Calculate aggregated data
        # self.data['c'] = self.env.penv.calc_c(self.XZ[:, 1])  # Velocity at each point
        # self.data['dl'] = np.linalg.norm(np.diff(self.XZ, axis=0), axis=1)  # dL at each arrival point (excluding initial point => shape of n_steps-1)
        # self.data['l'] = np.cumsum(np.insert(self.data['dl'], 0, 0.))
        # self.data['dt'] = self.data['dl'] / self.data['c'][:-1]  # dT at each arrival point (excluding initial point)
        # self.data['t'] = np.cumsum(np.insert(self.data['dt'], 0, 0.))  # Time from ray launch at each point

        # # Generate interpolated path function
        # # self.calc_z = interpolate.interp1d(self.XZ[:, 0], self.XZ[:, 1], kind='linear', bounds_error=True)
        # # self.__is_propagated = True


    def closest_step_to_point (self, target: np.ndarray):  # Calculate closest approach
        d_target = np.sqrt(np.power(target[0]-self.XZ[:, 0], 2) + np.power(target[1]-self.XZ[:, 1], 2))  # TODO: self.XZ
        return np.argmin(d_target)

    def plot (self, fig, **kwargs):  # TODO: Plot gain: self.L, self.G[freq]
        if True:  # self.__is_propagated:
            plt.plot(self.XZ[:,0], self.XZ[:,1], figure=fig, **kwargs)

# TODO: move out of Ray2D
# # Calculate time at target  # calculate step closest to target (on the x-axis)  # TODO: DISTANCE ON BOTH AXES (linalg.norm)
# self.step_target = np.argmin(np.abs(self.XZ[:, 0] - self.target[0])) if self.target is not None else None  # Step closest to target on the x-axis
# self.T_target = self.T[self.step_target] if self.target is not None else None
        
# Calculate distance to target  # TODO: better check interpolation range (if range_x < target.x) and find target step
# if self.target is not None:
# try: self.dist_to_target = abs(self.target[1] - self.calc_z(self.target[0]))
# except: self.dist_to_target = np.nan