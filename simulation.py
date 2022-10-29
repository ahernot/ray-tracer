# 2D Ray Simulation (v1.0)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np

from environment import Environment2D
from ray import Ray2D



class Simulation2D:

    def __init__ (self, env: Environment2D):
        self.env = env
        self.rays = list()  # Raw list of rays (stored as references to Ray2D objects)
        self.angles = dict()  # Dictionary of Ray2D objects grouped by angle
        # can group rays by frequency or some other key if needed
        self.n_rays = 0
        self.n_angles = 0

        self.range_min = np.zeros(2)
        self.range_max = np.zeros(2)
        self.size = np.zeros(2)

    def __repr__ (self):
        return f'2D simulation containing {self.n_rays} rays'  # TODO: to improve later

    def add_rays (self, source, freq, *angles, **kwargs):
        # Supports two rays for the same angle (for example two different frequencies)

        for angle in angles:
            # Generate ray
            ray = Ray2D (self.env, source, freq, angle)
            ray.propagate(**kwargs)

            # Add ray to simulation
            self.rays .append(ray)
            self.angles[angle] = self.angles[angle] + [ray] if angle in self.angles else [ray]

            # Update simulation range
            if ray.range_min[0] < self.range_min[0]: self.range_min[0] = ray.range_min[0]
            if ray.range_min[1] < self.range_min[1]: self.range_min[1] = ray.range_min[1]
            if ray.range_max[0] > self.range_max[0]: self.range_max[0] = ray.range_max[0]
            if ray.range_max[1] > self.range_max[1]: self.range_max[1] = ray.range_max[1]

        # Update simulation size
        self.size = self.range_max - self.range_min

        self.n_rays = len(self.rays)
        self.n_angles = len(self.angles)



    def plot (self):
        # Plot environment (floor & ceiling for now)
        # Plot rays (or selection of rays)
        pass

    def save (self):
        # save as a file
        pass

    def extract (self):
        # Return new Simulation with only the selected rays
        pass
