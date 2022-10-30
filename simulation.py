# 2D Ray Simulation (v1.0)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np

from environment import Environment2D
from ray import Ray2D

# TODO: simulation to rename source, to include inside of a simulation

class Simulation2D:

    def __init__ (self, env: Environment2D, source, **kwargs):
        """
        :param env: Simulation environment
        :param source: Source point

        kwargs:
        :param energy: Combined source power => divided first into the different frequencies (using self.spectrum), and then among rays
        :param spectrum: Power spectrum
        """
        
        self.env = env
        self.source = source

        self.rays = list()  # Raw list of rays (stored as references to Ray2D objects)
        self.angles = dict()  # Dictionary of Ray2D objects grouped by angle
        self.freqs = dict()  # Dictionary of Ray2D objects grouped by frequency

        self.n_rays = 0
        self.n_angles = 0
        self.n_freqs = 0

        self.energy_total = kwargs.get('energy', 1.)
        self.spectrum = kwargs.get('spectrum', lambda x: 1.)  # spectral power distribution (only the frequencies generated make up the ray power)
        self.__spectrum_vals = dict()  # Samples of self.spectrum non-normalised distribution function
        self.__spectrum_distrib = dict()  # xf

        self.range_min = np.zeros(2)
        self.range_max = np.zeros(2)
        self.size = np.zeros(2)

    def __repr__ (self):
        return f'2D simulation containing {self.n_rays} rays'  # TODO: to improve later

    def add_rays (self, freq, *angles, **kwargs):
        # Supports two rays for the same angle (for example two different frequencies)

        for angle in angles:

            # Generate ray
            ray = Ray2D (self.env, self.source, freq, angle)
            ray.propagate(**kwargs)

            # Add ray to simulation
            self.rays .append(ray)
            self.angles[angle] = self.angles[angle] + [ray] if angle in self.angles else [ray]
            self.freqs[freq] = self.freqs[freq] + [ray] if freq in self.freqs else [ray]

            # Update simulation range
            if ray.range_min[0] < self.range_min[0]: self.range_min[0] = ray.range_min[0]
            if ray.range_min[1] < self.range_min[1]: self.range_min[1] = ray.range_min[1]
            if ray.range_max[0] > self.range_max[0]: self.range_max[0] = ray.range_max[0]
            if ray.range_max[1] > self.range_max[1]: self.range_max[1] = ray.range_max[1]

        # Update simulation size
        self.size = self.range_max - self.range_min

        # Redistribute spectral power (if new frequency added)
        if len(self.freqs) > self.n_freqs:
            self.__distribute_spectral_power()

        # Update counters
        self.n_rays = len(self.rays)
        self.n_angles = len(self.angles)
        self.n_freqs = len(self.freqs)

        # Regenerate normalised ray energy unit
        self.__energy_norm = self.energy_total / np.sum([len(self.freqs[freq]) * self.__spectrum_distrib[freq] for freq in self.freqs])
        self.ray_power = {freq: self.__spectrum_distrib[freq] * self.__energy_norm for freq in self.freqs}


    def __distribute_spectral_power (self):
        """
        [desc]
        Run when adding a new frequency
        """

        # Add samples of self.spectrum non-normalised distribution function
        for freq in self.freqs:
            if freq not in self.__spectrum_vals:
                self.__spectrum_vals [freq] = self.spectrum(freq)

        # Regenerate total of raw spectral values
        self.__spectrum_total = np.sum(list(self.__spectrum_vals.values()))

        # Regenerate dictionary of spectral energy distribution per frequency (xf)
        self.__spectrum_distrib = {freq: self.spectrum(freq) / self.__spectrum_total for freq in self.freqs}




    def heatmap (self, **kwargs):
        """
        # TODO: redo this function
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
