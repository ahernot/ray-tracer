# 2D Ray Simulation (v1.2)
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
        :param kwargs/pack_default: Default pack name (use None to manually create and assign raypacks)
        """
        
        # Get simulation objects
        self.env = env
        self.source = source
        self.raypacks = dict()

        # Generate default pack
        self.pack_default = kwargs.get('pack_default', 'main') 
        if self.pack_default: self.raypacks[self.pack_default] = RayPack2D()
        
        # Get spectral power distribution
        self.spectrum = kwargs.get('spectrum', lambda x: 1.)  # spectral power distribution (only the frequencies generated make up the ray power)

        # Initialise ranges
        self.range_min = self.env.range_min
        self.range_max = self.env.range_max
        self.size = self.env.size
        self.range_min_plot = np.zeros(2)
        self.range_max_plot = np.zeros(2)
        self.size_plot = np.zeros(2)

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
            raypack.add(ray, **kwargs)

            # Update plot range
            if ray.range_min[0] < self.range_min_plot[0]: self.range_min_plot[0] = ray.range_min[0]
            if ray.range_min[1] < self.range_min_plot[1]: self.range_min_plot[1] = ray.range_min[1]
            if ray.range_max[0] > self.range_max_plot[0]: self.range_max_plot[0] = ray.range_max[0]
            if ray.range_max[1] > self.range_max_plot[1]: self.range_max_plot[1] = ray.range_max[1]

        # Update plot size
        self.size_plot = self.range_max_plot - self.range_min_plot

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
        heatmap_shape = np.ceil(self.size_plot/res).astype(int) + 1  # +1 because the downsampling of the coordinates is right bound inclusive
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

        # Pack names
        self.__pack_temp_scan = '__iter0-scan'
        self.pack_scan = 'scan'
        self.__gen_pack_temp_refine = lambda n_refines, ray_id: f'__iter{n_refines}-refine-ray{ray_id}'  # iter3-refine-ray7
        self.__gen_pack_refine = lambda n_refines: f'refine-{n_refines}'

        super(EigenraySim2D, self).__init__(env, source, pack_default=None, **kwargs)
        self.target = target
        self.scan_freq = kwargs.get('scan_freq', 1)
        self.init_kwargs = kwargs
        
        # Initial scan
        self.n_rays = kwargs.get('n_rays', 100)
        n_rays_scan = self.n_rays * 10  # TODO: default multiplier
        self.__scan(n_rays_scan)  # <= lacking angular range parameters, resolution, nb of rays to keep, distance threshold, etc.

        # 1. Scan (within reasonable angular range)
        # 2. Refine (iterations)
        # Generate one raypack per stage

        # Refinement stages counter
        self.n_refines = 0


    def __scan (self, n_rays_scan):
        # TODO: scan(n_rays, angle_min, angle_max) with separate processing function taking average of env and nb of rebounds (+1)

        ###
        angle_min = -1 * np.pi / 2 + 0.01
        angle_max = np.pi / 2 - 0.01
        angles = np.linspace(angle_min, angle_max, n_rays_scan)

        # Generate temporary raypack
        self.raypacks[self.__pack_temp_scan] = RayPack2D()
        ray__pack_temp_scan = self.raypacks[self.__pack_temp_scan]

        # Cast n_rays_scan scanning rays
        self.cast (self.scan_freq, *angles, pack=self.__pack_temp_scan, target=self.target, **self.init_kwargs)
        # Keep rays
        dist_list_sorted = np.sort(list(ray__pack_temp_scan.dist.keys()))

        # Generate scan output raypack
        self.raypacks[self.pack_scan] = RayPack2D()
        raypack_scan = self.raypacks[self.pack_scan]

        # keep n rays for each nb of rebounds, up to max nb of rebounds



    def refine (self, **kwargs):
        """
        Refine
        :param ?:
        :param kwargs/iterations: Number of back-to-back refine iterations
        """

        iterations = kwargs.get('iterations', 1)

        # Load previous raypack
        pack_prev = self.pack_scan if self.n_refines == 0 else self.__gen_pack_refine(self.n_refines-1)
        raypack_prev = self.raypacks[pack_prev]

        # Generate new raypack
        pack_refine = self.__gen_pack_refine(self.n_refines)
        self.raypacks[pack_refine] = RayPack2D()
        raypack_refine = self.raypacks[pack_refine]

        for i in range(iterations):
            for ray_id, ray in enumerate():
                # TODO: IMPORTANT: ALWAYS KEEP RAYS IN SAME ORDER FOR EACH REFINE (CHILD OF RAY 1 IS NEW_RAY 1 ETC)
                # Increment refines counter (first refine has id=1)
                self.n_refines += 1

                # Generate temporary raypack
                pack_temp = self.__gen_pack_temp_refine(self.n_refines, ray_id)
                self.raypacks[pack_temp] = RayPack2D()
                raypack_temp = self.raypacks[pack_temp]
                # IF NO PROGRESS FOR A SPECIFIC RAY: DO SOMETHING ELSE (check dist_new - dist_prev)

                # Cast rays
                # self.cast (self.scan_freq, *angles, raypack=pack_temp)#, **self.init_kwargs)

                # DO SORTING STUFF (function?)

    # Generate filter => requires > 1 refine iteration and will by default choose the final refine raypack


    def get (self):
        # TODO: rename, etc
        pack = self.__gen_pack_refine(self.n_refines)
        raypack = self.raypacks[pack]
