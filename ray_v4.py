import numpy as np
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
        self.angle: Union[float, int] = angle

        self.data = dict()  # TODO: store G, T, etc in a big dictionary
        # self.XZ = np.expand_dims(self.source.copy(), axis=0)  # TODO: rename to self.xz?
        # self.T = np.array([0.,])

    def propagate (self, **kwargs):  # TODO: propagate with conditions (no backpropagation, …)
        # Parse kwargs
        n_steps_max = kwargs.get('n_steps_max', N_STEPS_MAX)
        n_rebounds_max = kwargs.get('n_rebounds_max', N_REBOUNDS_MAX)
        allow_backprop = kwargs.get('allow_backprop', False)
        dx_max = kwargs.get('dx_max', DX_MAX_DEFAULT)
        dz_max = kwargs.get('dz_max', DZ_MAX_DEFAULT)
        # TODO: n_rebounds

        # Calculate decision slope between min resolutions on both axes
        dx_z_max = np.abs(self.dz_max / self.dx_max)  # Decision slope between both min resolutions

        # Calculate initial parameters
        P = self.source
        angle = -1 * self.angle + (np.pi / 2)  # convert angle notation  # TODO: don't? cf equation
        dx_z = 1 / np.tan(angle)  # Initial slope
        dxdx_z = 0  # Initial curvature
        x_dir = 1.  # forwards propagation  # TODO: needed?
        k = np.array([x_dir, dx_z])




        # Initialise solver parameters
        c = self.env.penv.calc_c(0)  # TODO: 0? or z0?
        mult = -1 * np.power(c / np.sin(angle), 2)  # differential equation multiplier

        # Propagation loop
        for i in range(n_steps_max):

            # Enforce minimum resolution (k becomes this step's [dx, dz])
            if abs(dx_z) > dx_z_max: k *= self.dz_max / abs(dx_z)
            else: k *= self.dx_max

            # Unpack coordinates
            x, z = P
            x_new, z_new = P_new = P + k
            
            # Check horizontal simulation bounds
            # TODO

            # Check for rebounds
            # TODO

            # Check all bounds
            # TODO

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
            # TODO

        
        # Data aggregation
        self.n_steps = self.XZ.shape[0]



        #! Data aggregation
        # Populate ray's effective range
        self.range_min = np.array([np.min(self.XZ[:, 0]), np.min(self.XZ[:, 1])])
        self.range_max = np.array([np.max(self.XZ[:, 0]), np.max(self.XZ[:, 1])])

        # Calculate integration segments
        # TODO: self.data['dl', 'l', 'c', 'dt', 't']
        self.dL = np.linalg.norm(np.diff(self.XZ, axis=0), axis=1)  # dL at each arrival point (excluding initial point => shape of n_steps-1)
        self.L = np.cumsum(np.insert(self.dL, 0, 0.))
        self.C = self.env.penv.calc_c(self.XZ[:, 1])  # Velocity at each point
        self.dT = self.dL / self.C[:-1]  # dT at each arrival point (excluding initial point)
        self.T = np.cumsum(np.insert(self.dT, 0, 0.))  # Time from ray launch at each point

        # Generate interpolated path function
        self.calc_z = interpolate.interp1d(self.XZ[:, 0], self.XZ[:, 1], kind='linear', bounds_error=True)
        self.__is_propagated = True

    def closest_step_to_point (self, target: np.ndarray):  # Calculate closest approach
        d_target = np.sqrt(np.power(target[0]-self.XZ[:, 0], 2) + np.power(target[1]-self.XZ[:, 1], 2))
        return np.argmin(d_target)

# TODO: move out of Ray2D
# # Calculate time at target  # calculate step closest to target (on the x-axis)  # TODO: DISTANCE ON BOTH AXES (linalg.norm)
# self.step_target = np.argmin(np.abs(self.XZ[:, 0] - self.target[0])) if self.target is not None else None  # Step closest to target on the x-axis
# self.T_target = self.T[self.step_target] if self.target is not None else None
        
# Calculate distance to target  # TODO: better check interpolation range (if range_x < target.x) and find target step
# if self.target is not None:
# try: self.dist_to_target = abs(self.target[1] - self.calc_z(self.target[0]))
# except: self.dist_to_target = np.nan