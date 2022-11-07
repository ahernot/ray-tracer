# 2D Ray Simulation (v1.1)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# TODO: Preset with self.env.range_{min/max} instead of extending for each ray? But range is useful for plotting too… so maybe add plot_range_{min/max}?

import numpy as np
# import matplotlib.pyplot as plt

# from preferences import *
from functions import coords_to_mask_2d

from environment import Environment2D
from ray import Ray2D



class Simulation2D:

    def __init__ (self, env: Environment2D, source, **kwargs):
        """
        :param env: Simulation environment
        :param source: Source point

        kwargs
        :param energy: Combined source power => divided first into the different frequencies (using self.spectrum), and then among rays
        :param spectrum: Power spectrum
        """
        
        self.env = env
        self.source = source

        self.rays = list()  # Raw list of rays (stored as references to Ray2D objects)
        self.angles = dict()  # Dictionary of Ray2D objects indexed by angle
        self.freqs = dict()  # Dictionary of Ray2D objects indexed by frequency
        self.stop_reasons = dict()  # Dictionary of Ray2D objects indexed by stop reason

        self.n_rays = 0
        self.n_angles = 0
        self.n_freqs = 0

        self.energy_total = kwargs.get('energy', 1.)
        self.spectrum = kwargs.get('spectrum', lambda x: 1.)  # spectral power distribution (only the frequencies generated make up the ray power)
        self.__spectrum_vals = dict()  # Samples of self.spectrum non-normalised distribution function
        self.__spectrum_distrib = dict()  # xf

        self.range_min = np.zeros(2)  # maximum extent is self.env.range_min
        self.range_max = np.zeros(2)  # maximum extent is self.env.range_min
        self.size = np.zeros(2)

    def __repr__ (self):  # TODO: to improve later
        stop_reasons_formatted = '\n'.join([f'\t{stop_reason}: {len(self.stop_reasons[stop_reason])}' for stop_reason in self.stop_reasons])
        return f'2D simulation containing {self.n_rays} rays\nStop reasons:\n{stop_reasons_formatted}'

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
            self.stop_reasons[ray.stop_reason] = self.stop_reasons[ray.stop_reason] + [ray] if ray.stop_reason in self.stop_reasons else [ray]

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
        self.ray_energy = {freq: self.__spectrum_distrib[freq] * self.__energy_norm for freq in self.freqs}


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
        Generate heatmap of ray power
        kwargs
        :param res: Resolution (in meters) - # TODO? must be less precise than the simulation's max dz and dx
        :param n_reductions:
        :param cutoff: Saturation percentage for normalised heatmap (pre-log scaling)
        """
        res = kwargs.get('resolution', np.array([50, 25]))
        reduction_power = kwargs.get('reduction_power', .5)
        cutoff = kwargs.get('cutoff', .02)

        # Initialise heatmap
        heatmap_shape = np.ceil(self.size/res).astype(int) + 1  # +1 because the downsampling of the coordinates is right bound inclusive
        heatmap_full = np.zeros((heatmap_shape))

        for ray in self.rays:

            # Generate array of downsampled coordinates of shape (N, 2) = ray.XZ.shape and convert from physical coordinates to matrix coordinates (inverting y-axis)
            rx = ray.XZ.copy().astype(int)
            rx[:, 0] //= res[0]
            rx[:, 1] //= -1 * res[1]

            # Plot ray heatmap according to ray energy
            vals = np.power(ray.Tmult, -0.1)
            heatmap_ray = coords_to_mask_2d(heatmap_shape, rx, vals) * self.ray_energy[ray.freq]
            heatmap_full += heatmap_ray

        # Normalise heatmap
        heatmap_norm = heatmap_full / np.max(heatmap_full)
        heatmap_norm[heatmap_norm > cutoff] = cutoff
        heatmap_plot = np.power(heatmap_norm, reduction_power) if reduction_power != 1 else heatmap_norm

        return heatmap_plot.T


    def plot (self, fig):

        for ray in self.rays:
            ray.plot(fig)
        self.env.plot(fig, c='red')


    def save (self):
        # save as a file
        pass


    def extract (self):
        # Return new Simulation with only the selected rays
        pass
