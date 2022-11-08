# 2D Ray Simulation (v1.1)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# TODO: Preset with self.env.range_{min/max} instead of extending for each ray? But range is useful for plotting too… so maybe add plot_range_{min/max}?

import numpy as np
# import matplotlib.pyplot as plt

from functions import coords_to_mask_2d

from environment import Environment2D
from ray import Ray2D, RayPack2D



class Simulation2D:

    def __init__ (self, env: Environment2D, source, **kwargs):
        """
        Simulation (total_energy=1.)
        :param env: Simulation environment
        :param source: Source point
        :param kwargs/spectrum: Power spectrum
        :param kwargs/pack_default: Default pack name
        """
        
        self.env = env
        self.source = source
        self.pack_default = kwargs.get('pack_default', 'main')
        self.raypacks = {self.pack_default: RayPack2D()}
        

        self.spectrum = kwargs.get('spectrum', lambda x: 1.)  # spectral power distribution (only the frequencies generated make up the ray power)

        # TODO: Generate a range_min and range_max for plotting which is updated with every ray added => for RayPack2D rather than Simulation2D
        self.range_min = self.env.range_min
        self.range_max = self.env.range_max
        self.size = self.env.size

    def __repr__ (self):  # TODO: to improve later
        stop_reasons_formatted = '\n'.join([f'\t{stop_reason}: {len(self.stop_reasons[stop_reason])}' for stop_reason in self.stop_reasons])
        return f'2D simulation containing {self.n_rays} rays\nStop reasons:\n{stop_reasons_formatted}'

    def __distribute_spectral_power (self, pack):

        # Get target raypack
        raypack = self.raypacks[pack]

        # Add samples of self.spectrum non-normalised distribution function
        for freq in raypack.freqs:
            if freq not in raypack.spectrum_vals:
                raypack.spectrum_vals [freq] = self.spectrum(freq)

        # Regenerate total of raw spectral values
        raypack.spectrum_total = np.sum(list(raypack.spectrum_vals.values()))

        # Regenerate dictionary of spectral energy distribution per frequency (xf)
        raypack.spectrum_distrib = {freq: self.spectrum(freq) / raypack.spectrum_total for freq in raypack.freqs}

    def cast (self, freq, *angles, **kwargs):
        """
        Propagate rays and add them to the simulation
        :param freq: Ray frequency (Hz)
        :param angles: Angles
        :param kwargs/pack: Raypack name
        """

        # Get target raypack
        pack = kwargs.get('pack', self.pack_default)
        if pack not in self.raypacks: self.raypacks[pack] = RayPack2D()
        raypack = self.raypacks[pack]

        # Power recalculation flag
        redistribute_power = False
        if freq not in raypack.freqs:
            redistribute_power = True

        for angle in angles:

            # Generate ray
            ray = Ray2D (self.env, self.source, freq, angle)
            ray.propagate(**kwargs)
            # Add ray to simulation
            raypack.add(ray)

        # Redistribute spectral power (if new frequency added)
        if redistribute_power:
            self.__distribute_spectral_power(pack)
        # Regenerate ray unit energy
        raypack.regen_energy()

    def heatmap (self, **kwargs):
        """
        Generate heatmap of ray power
        kwargs
        :param res: Resolution (in meters) - # TODO? must be less precise than the simulation's max dz and dx
        :param n_reductions:
        :param cutoff: Saturation percentage for normalised heatmap (pre-log scaling)
        :param kwargs/pack: Raypack name
        """

        # Get heatmap parameters
        res = kwargs.get('resolution', np.array([50, 25]))
        reduction_power = kwargs.get('reduction_power', .5)
        cutoff = kwargs.get('cutoff', .02)
        
        # Get target raypack
        pack = kwargs.get('pack', self.pack_default)
        raypack = self.raypacks[pack]

        # Initialise heatmap
        heatmap_shape = np.ceil(self.size/res).astype(int) + 1  # +1 because the downsampling of the coordinates is right bound inclusive
        heatmap_full = np.zeros((heatmap_shape))

        for ray in raypack.rays:

            # Generate array of downsampled coordinates of shape (N, 2) = ray.XZ.shape and convert from physical coordinates to matrix coordinates (inverting y-axis)
            rx = ray.XZ.copy().astype(int)
            rx[:, 0] //= res[0]
            rx[:, 1] //= -1 * res[1]

            # Plot ray heatmap according to ray energy
            vals = np.power(ray.Tmult, -0.1)
            heatmap_ray = coords_to_mask_2d(heatmap_shape, rx, vals) * raypack.ray_energy[ray.freq]
            heatmap_full += heatmap_ray

        # Normalise heatmap
        heatmap_norm = heatmap_full / np.max(heatmap_full)
        heatmap_norm[heatmap_norm > cutoff] = cutoff
        heatmap_plot = np.power(heatmap_norm, reduction_power) if reduction_power != 1 else heatmap_norm

        return heatmap_plot.T

    def plot (self, fig, **kwargs):
        """
        :param fig: Figure
        :param kwargs/pack: Raypack name
        """

        # Get target raypack
        pack = kwargs.get('pack', self.pack_default)
        raypack = self.raypacks[pack]

        for ray in raypack.rays:
            ray.plot(fig)
        self.env.plot(fig, c='red')

    def save (self):
        # save as a file
        pass

    def extract (self):
        # Return new Simulation with only the selected rays
        pass



class EigenraySim2D (Simulation2D):

    def __init__ (self, env: Environment2D, source, target, **kwargs):
        """
        Eigenray simulation
        :param env: Simulation environment
        :param source: Source point
        :param target: Target point
        :param kwargs/spectrum: Power spectrum
        :param kwargs/scan_freq: Scanning frequency (doesn't affect ray path)
        """
        super(EigenraySim2D, self).__init__(env, source, pack_default='scan', **kwargs)
        self.target = target
        self.scan_freq = kwargs.get('scan_freq', 1)
        
        ###
        self.init_kwargs = kwargs
        self.__scan()  # <= lacking angular range parameters, resolution, nb of rays to keep, distance threshold, etc.

        # 1. Scan (within reasonable angular range)
        # 2. Refine (iterations)
        # Generate one raypack per stage

        self.n_refines = 0


    def __scan (self):

        ###
        angle_min = -1 * np.pi / 2 + 0.01
        angle_max = np.pi / 2 - 0.01
        angles = np.linspace(angle_min, angle_max, 1000)

        self.cast (self.scan_freq, *angles, pack='scan', **self.init_kwargs)


    def refine (self, **kwargs):
        """
        Refine
        :param ?:
        :param kwargs/iterations: Number of back-to-back refine iterations
        """

        iterations = kwargs.get('iterations', 1)

        for i in range(iterations):

            # Initialise raypack
            self.n_refines += 1
            self.raypacks[f'refine_{self.n_refines}'] = RayPack2D()

            ###

            # self.cast (self.scan_freq, *angles, raypack=raypack, **self.init_kwargs)

    # Generate filter => requires > 1 refine iteration and will by default choose the final refine raypack