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
        self.source: np.ndarray = source  # TODO: cast to ndarray
        self.angle: Union[float, int] = angle

        self.data = dict()  # TODO: store G, T, etc in a big dictionary

    def propagate (self, **kwargs):  # TODO: propagate with conditions (no backpropagation, …)
        # TODO: parse kwargs

        # Calculate initial parameters
        P = self.source
        angle = -1 * self.angle + (np.pi / 2)  # convert angle notation
        dx_z = 1 / np.tan(angle)  # Initial slope
        dxdx_z = 0  # Initial curvature
        x_dir = 1.  # forwards propagation  # TODO: needed?
        k = np.array([x_dir, dx_z])

        # TODO: loop for calculation

