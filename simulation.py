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
        # self.n_freqs = 0  # TODO

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


    def heatmap (self, **kwargs):
        """
        Generate heatmap of ray power
        kwargs:
        :param res: Resolution (in meters)
        :param cutoff:
        """
        res = kwargs.get('resolution', np.array([50, 25]))
        cutoff = kwargs.get('cutoff', 150)

        # Initialise heatmap
        xdim, ydim = heatmap_shape = np.ceil(self.size/res).astype(int)
        heatmap_full = np.zeros(heatmap_shape)

        # Trace rays
        for ray in self.rays:
            rx = ray.XZ.copy().astype(int)

            # Homonegeneise bounds
            rx = np.insert(rx, 0, self.range_min, axis=0)
            rx = np.insert(rx, 0, self.range_max, axis=0)

            # Downsample
            x = rx[:, 0] // res[0]
            y = rx[:, 1] // -1 * res[1]

            # Add ray's 2D histogram
            heatmap, xedges, yedges = np.histogram2d(x, y, bins=heatmap_shape)
            heatmap_full += heatmap

        # Process heatmap
        heatmap_full[heatmap_full > cutoff] = cutoff
        heatmap_norm = heatmap_full / np.max(heatmap_full)
        heatmap_plot = np.log(heatmap_norm + 1.)

        # TODO: plot ground (convex assumption)
        # z_floor = self.env.floor(np.arange(0, xdim, 1) * res[0])
        # z_width = np.arange(0, ydim, 1) * -1 * res[1]
        # floor_mask = (np.tile(z_width, (xdim, 1)).T < z_floor).T

        return heatmap_plot

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
