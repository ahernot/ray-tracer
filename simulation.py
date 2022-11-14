# 2D Ray Simulation (v1.3)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# TODO: Preset with self.env.range_{min/max} instead of extending for each ray? But range is useful for plotting too… so maybe add plot_range_{min/max}?
# TODO: Calculate number of ray steps using rough max nb of rebounds to never stop midway
# TODO: # keep n rays for each nb of rebounds, up to max nb of rebounds? wise? idts
# TODO: always keep the rays in the same order after each refine: useful?

import numpy as np
# import matplotlib.pyplot as plt
from operator import itemgetter
import itertools

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
        :param kwargs/verbose: Verbose (default=False)
        :param kwargs/verbose_depth: Verbose depth
        :param kwargs/verbose_depth_max: Verbose max depth
        """

        # Verbose
        self.verbose = kwargs.get('verbose', False)
        self.verbose_depth_max = kwargs.get('verbose_depth_max', 0)
        self.verbose_depth =  kwargs.get('verbose_depth', 0)  # TODO: if verbose_depth > verbose_depth_max for a function (while applying verbose_depth++ for a deeper function each step) then set verbose=False regardless
        if self.verbose: self.verbose_next = (self.verbose_depth < self.verbose_depth_max) or (self.verbose_depth_max == -1)
        self.__vi = '\t' * self.verbose_depth + self.__class__.__name__ + ': '
        
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

        :param Ray2D.__init__/kwargs/target: Target point
        :param Ray2D.propagate/kwargs/backprop: Allow backpropagation (default=True)
        :param Ray2D.propagate/kwargs/dx_max:
        :param Ray2D.propagate/kwargs/dz_max:
        :param Ray2D.propagate/kwargs/n_steps_max:
        :param Ray2D.propagate/kwargs/n_rebounds_max: Maximum number of rebounds (default=infinite)
        :param Ray2D.propagate/kwargs/verbose: Verbose (default=False)
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
            # Generate ray and add to simulation
            kwargs.update(verbose=self.verbose_next, verbose_depth=self.verbose_depth+1, verbose_depth_max=self.verbose_depth_max)
            ray = Ray2D (self.env, self.source, freq, angle, **kwargs)
            ray.propagate(**kwargs)
            raypack.add(ray)

        # Redistribute spectral power (if new frequency added) and regenerate ray unit energy
        if redistribute_power: self.__distribute_spectral_power(pack)
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
        heatmap_shape = np.ceil(raypack.size_plot/res).astype(int) + 1  # +1 because the downsampling of the coordinates is right bound inclusive
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
        kwargs.pop('pack', None)

        for ray in raypack.rays:
            ray.plot(fig, **kwargs)
        self.env.plot(fig, c='red')

    def save (self):
        # save as a file
        pass

    def extract (self):
        # Return new Simulation with only the selected rays
        pass



class EigenraySim2D (Simulation2D):

    def __init__ (self, env: Environment2D, source: np.ndarray, target: np.ndarray, n_rebounds_max: int, **kwargs):
        """
        Eigenray simulation
        :param env: Simulation environment
        :param source: Source point
        :param target: Target point
        :param n_rebounds_max: Maximum number of rebounds

        :param kwargs/spectrum: Power spectrum
        :param kwargs/scan_freq: Scanning frequency (doesn't affect ray path)
        :param kwargs/n_rays: Number of rays (can be exceeded if rays have the same distance)
        :param n_rays_scan: Number of scanning rays (higher is better)
        :param kwargs/scan_angle_min: Override automatic min scan angle
        :param kwargs/scan_angle_max: Override automatic max scan angle
        :param kwargs/precision_cutoff: Focusing precision cutoff (in meters)

        :param Ray2D.propagate/kwargs/backprop: Allow backpropagation (default=True)
        :param Ray2D.propagate/kwargs/dx_max:
        :param Ray2D.propagate/kwargs/dz_max:
        :param Ray2D.propagate/kwargs/n_steps_max:
        :param Ray2D.propagate/kwargs/verbose: Verbose (default=False)
        """
        self.init_kwargs = kwargs

        # Pack names
        self.pack_temp_scan = '__iter0-scan'
        self.pack_scan = 'scan'
        self.gen_pack_temp_refine = lambda n_refines, ray_id: f'__iter{n_refines}-refine-ray{ray_id}'  # iter3-refine-ray7
        self.gen_pack_refine = lambda n_refines: f'refine-{n_refines}'

        # Initialise Simulation2D
        super(EigenraySim2D, self).__init__(env, source, pack_default=None, **kwargs)

        # Initialise EigenraySim2D
        self.target: np.ndarray = target
        self.n_rebounds_max = n_rebounds_max
        self.scan_freq = kwargs.get('scan_freq', 1)
        self.n_rays: int = kwargs.get('n_rays', 25)
        self.dist = dict()
        self.dist_avg = None
        self.n_refines = 0
        self.__angular_precision = None
        self.precision_cutoff = kwargs.get('precision_cutoff', 0.1)  # in meters

        # Initial scan        
        self.n_rays_scan = kwargs.get('n_rays_scan', self.n_rays * 10)  # TODO: default multiplier
        self.__scan(kwargs.get('scan_angle_min', None), kwargs.get('scan_angle_max', None))

    def __scan (self, angle_min = None, angle_max = None):
        """
        Initial scan (broad sweep)
        :param angle_min:
        :param angle_max:
        :param Ray2D.propagate/kwargs/backprop: Allow backpropagation (default=True)
        :param Ray2D.propagate/kwargs/dx_max:
        :param Ray2D.propagate/kwargs/dz_max:
        :param Ray2D.propagate/kwargs/n_steps_max:
        """

        # Calculate scanning angular range (cf README.md)
        Dz = abs(self.env.ceil_avg - self.env.floor_avg)  # Distance between avg_floor and avg_ceiling  # TODO: AVERAGES OVER THE COURSE OF THE RAY'S TRAJECTORY (BETWEEN SOURCE AND TARGET)
        Dx = abs(self.target[0] - self.source[0])  # Horizontal distance between source and target
        Dz_sc = abs(self.source[1])  # Vertical distance between source and avg_ceiling
        Dz_tc = abs(self.target[1])  # Vertical distance between target and avg_ceiling
        N = self.n_rebounds_max + 1  # Allow for one extra rebound to account for averaging imprecisions
        if N % 2 == 0:
            angle_min = angle_min if angle_min else -1 * np.arcsin( Dx / (N * Dz - Dz_sc + Dz_tc) )
            angle_max = angle_max if angle_max else np.arcsin( Dx / (N * Dz + Dz_sc - Dz_tc) )
        elif N % 2 == 1:
            angle_min = angle_min if angle_min else -1 * np.arcsin( Dx / ((N+1) * Dz - Dz_sc - Dz_tc) )
            angle_max = angle_max if angle_max else np.arcsin( Dx / ((N-1) * Dz + Dz_sc + Dz_tc) )

        # Calculate angular precision
        if self.verbose: print(f'{self._Simulation2D__vi}Scanning using {self.n_rays_scan} rays between angles {angle_min} and {angle_max}')
        self.__angular_precision = (angle_max - angle_min) / self.n_rays_scan# / 2  # Cone half-angle

        # Cast scanning rays
        self.raypacks[self.pack_temp_scan] = RayPack2D()
        angles = np.linspace(angle_min, angle_max, self.n_rays_scan)
        self.cast (self.scan_freq, *angles, pack=self.pack_temp_scan, target=self.target, n_rebounds_max=self.n_rebounds_max, **self.init_kwargs)
        
        # Sort and select rays
        raypack_temp_scan = self.raypacks[self.pack_temp_scan]
        dist_sorted = raypack_temp_scan.dist_sorted
        if np.isnan(dist_sorted[-1]): dist_sorted = dist_sorted[:-1]
        rays = list(itertools.chain.from_iterable( itemgetter(*dist_sorted[:self.n_rays])(raypack_temp_scan.dist) ))

        # Generate scan output raypack
        self.raypacks[self.pack_scan] = RayPack2D()
        raypack_scan = self.raypacks[self.pack_scan]
        raypack_scan.add(*rays)
        # Update scan raypack
        self._Simulation2D__distribute_spectral_power(self.pack_scan)
        raypack_scan.regen_energy()
        # Set default raypack
        self.pack_default = self.pack_scan

        self.dist_avg = np.mean(raypack_scan.dist_sorted)
        if self.verbose: print(f'{self._Simulation2D__vi}Scan mean distance: {self.dist_avg}')

    def refine (self, **kwargs):
        """
        Refine rays
        :param kwargs/cone_half_angle:
        :param kwargs/n_rays:
        :param kwargs/iterations: Number of back-to-back refine iterations
        :param kwargs/force_refine: Force refine past precision cutoff
        :param kwargs/verbose: Verbose (default=False)
        :param Ray2D.propagate/kwargs/backprop: Allow backpropagation (default=True)
        :param Ray2D.propagate/kwargs/dx_max:
        :param Ray2D.propagate/kwargs/dz_max:
        :param Ray2D.propagate/kwargs/n_steps_max:
        """
        
        cone_half_angle = kwargs.get('cone_half_angle', self.__angular_precision)
        n_rays = kwargs.get('n_rays', 5)
        iterations = kwargs.get('iterations', 1)
        force_refine = kwargs.get('force_refine', False)

        for i in range(iterations):

            # Precision check
            if self.dist_avg <= self.precision_cutoff and not force_refine:
                if self.verbose: print(f'{self._Simulation2D__vi}Target precision reached')
                break

            # Load previous raypack
            pack_prev = self.__get_last_pack()
            raypack_prev = self.raypacks[pack_prev]
            # Increment refines counter (first refine has id=1)
            self.n_refines += 1  # TODO: better log of refine steps (incl. parameters and number of cast and selected) for __repr__
            # Generate new raypack
            pack_refine = self.gen_pack_refine(self.n_refines)
            self.raypacks[pack_refine] = RayPack2D()
            raypack_refine = self.raypacks[pack_refine]

            for ray_id, ray in enumerate(raypack_prev.rays):
                # Generate temporary refine raypack
                pack_temp = self.gen_pack_temp_refine(self.n_refines, ray_id)
                self.raypacks[pack_temp] = RayPack2D()
                raypack_temp = self.raypacks[pack_temp]

                # Cast rays
                angle_min = ray.angle - cone_half_angle
                angle_max = ray.angle + cone_half_angle
                angles = np.linspace(angle_min, angle_max, n_rays)
                self.cast (self.scan_freq, *angles, pack=pack_temp, target=self.target, **self.init_kwargs)

                # Add initial ray (add after cast so that frequency is already generated)
                raypack_temp.add(ray)
                raypack_temp.regen_energy()
                
                # Select best ray
                dmin = raypack_temp.dist_sorted[0]
                ray_selected = raypack_temp.dist[dmin][0]  # Select best ray # TODO: only 1? or multiple if same distance?
                raypack_refine.add(ray_selected)

            # Update refine raypack
            self._Simulation2D__distribute_spectral_power(pack_refine)
            raypack_refine.regen_energy()

            # Update angular precision (precision cone half-angle)
            self.__angular_precision /= n_rays
            self.dist_avg = np.mean(raypack_refine.dist_sorted)

            if self.verbose: print(f'{self._Simulation2D__vi}Refine #{self.n_refines} mean distance: {self.dist_avg}')  # if self.verbose

        # Set default raypack
        self.pack_default = pack_refine

    def __get_last_pack (self):
        return self.pack_scan if self.n_refines == 0 else self.gen_pack_refine(self.n_refines)


    # TODO: Generate filter => requires > 1 refine iteration and will by default choose the final refine raypack
    # Can downsample
    # def downsample (self):
    # creates new downsampled pack

