# 2D Ray Simulation (v1.3)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# TODO: Calculate number of ray steps using rough max nb of rebounds to never stop midway
# TODO: Keep n rays for each nb of rebounds, up to max nb of rebounds? wise? idts
# TODO: Always keep the rays in the same order after each refine: useful?
# TODO: Simulation2D.__repr__()
# TODO: Make RayPack2D.add() work when adding a new ray without generating it using Simulation2D: all the power recalculations must be done from within the raypack
# TODO: add scanning and refining timers
# TODO: add freq_min, freq_max parameters to always load in [-24000, 24000] Hz
# TODO: Simulation2D.n_rays which counts all rays in all raypacks


import numpy as np
from scipy.fft import fft, fftfreq, ifft
import cmath

# import matplotlib.pyplot as plt
from operator import itemgetter
import itertools

# from functions import coords_to_mask_2d

from environment import Environment2D
from ray import Ray2D, RayPack2D



class Simulation2D:

    def __init__ (self, env: Environment2D, source: np.ndarray, **kwargs):
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
        
        # Initialise ranges
        self.range_min = self.env.range_min
        self.range_max = self.env.range_max
        self.size = self.env.size

    def __repr__ (self):  # TODO: __repr__
        # stop_reasons_formatted = '\n'.join([f'\t{stop_reason}: {len(self.stop_reasons[stop_reason])}' for stop_reason in self.stop_reasons])
        # return f'2D simulation containing {self.n_rays} rays\nStop reasons:\n{stop_reasons_formatted}'
        return self.__class__.__name__ + ' object'

    def cast (self, *angles, **kwargs):
        """
        Propagate rays and add them to the simulation
        :param angles: Angles
        :param kwargs/pack: Raypack name  # TODO: Rename to label (can either append to pack or create a new one)

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

        for angle in angles:
            # Generate ray and add to simulation
            kwargs.update(verbose=self.verbose_next, verbose_depth=self.verbose_depth+1, verbose_depth_max=self.verbose_depth_max)
            ray = Ray2D (self.env, self.source, angle, **kwargs)
            ray.propagate(**kwargs)
            raypack.add(ray)

    # def heatmap (self, **kwargs):
    #     # DEPRECATED FUNCTION
    #     """
    #     Generate heatmap of ray power
    #     kwargs
    #     :param resolution (array, x and y): Resolution (in meters) - # TODO? must be less precise than the simulation's max dz and dx
    #     :param n_reductions:
    #     :param cutoff: Saturation percentage for normalised heatmap (pre-log scaling)
    #     :param kwargs/pack: Raypack name
    #     """

    #     # Get heatmap parameters
    #     res = kwargs.get('resolution', np.array([50, 25]))
    #     reduction_power = kwargs.get('reduction_power', .5)
    #     cutoff = kwargs.get('cutoff', .02)
        
    #     # Get target raypack
    #     pack = kwargs.get('pack', self.pack_default)
    #     raypack = self.raypacks[pack]

    #     # Initialise heatmap
    #     heatmap_shape = np.ceil(raypack.size_plot/res).astype(int) + 1  # +1 because the downsampling of the coordinates is right bound inclusive
    #     heatmap_full = np.zeros((heatmap_shape))

    #     for ray in raypack.rays:

    #         # Generate array of downsampled coordinates of shape (N, 2) = ray.XZ.shape and convert from physical coordinates to matrix coordinates (inverting y-axis)
    #         rx = ray.XZ.copy().astype(int)
    #         rx[:, 0] //= res[0]
    #         rx[:, 1] //= -1 * res[1]

    #         # Plot ray heatmap according to ray energy
    #         ray.populate(1)
    #         vals = np.power(ray.Tmult[1], -0.1)  # TODO: need to populate first
    #         heatmap_ray = coords_to_mask_2d(heatmap_shape, rx, vals)
    #         heatmap_full += heatmap_ray

    #     # Normalise heatmap
    #     heatmap_norm = heatmap_full / np.max(heatmap_full)
    #     heatmap_norm[heatmap_norm > cutoff] = cutoff
    #     heatmap_plot = np.power(heatmap_norm, reduction_power) if reduction_power != 1 else heatmap_norm

    #     return heatmap_plot.T

    def plot (self, fig, **kwargs):
        """
        :param fig: Figure
        :param kwargs/pack: Raypack to plot (defaults to default raypack)  # TODO: rename to label
        """

        # Get target raypack
        pack = kwargs.get('pack', self.pack_default)
        raypack = self.raypacks[pack]
        kwargs.pop('pack', None)

        for ray in raypack.rays:
            ray.plot(fig, **kwargs)
        self.env.plot(fig, c='red')




class EigenraySim2D (Simulation2D):

    def __init__ (self, env: Environment2D, source: np.ndarray, target: np.ndarray, n_rays: int, n_rebounds_max: int, **kwargs):
        """
        Eigenray simulation
        :param env: Simulation environment
        :param source: Source point
        :param target: Target point
        :param n_rays: Number of rays to refine (can be exceeded if multiple rays have the same distance)
        :param n_rebounds_max: Maximum number of rebounds

        :param kwargs/spectrum: Power spectrum
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

        # Initialise Simulation2D & EigenraySim2D
        super(EigenraySim2D, self).__init__(env, source, pack_default=None, **kwargs)
        self.target = target
        self.n_rebounds_max = n_rebounds_max
        self.n_rays = n_rays
        self.precision_cutoff = kwargs.get('precision_cutoff', 0.1)  # in meters
        # Variables updated along the way
        self.dist = dict()
        self.dist_avg = None
        self.n_refines = 0
        self.__angular_precision = None
        
        # Initial scan        
        self.n_rays_scan = kwargs.get('n_rays_scan', self.n_rays * 10)  # TODO: default multiplier  # TODO: must be >= 2
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
        self.__angular_precision = (angle_max - angle_min) / (self.n_rays_scan - 1)  # / 2  # Cone half-angle  # TODO: ???

        # Cast scanning rays
        self.raypacks[self.pack_temp_scan] = RayPack2D()  # TODO: Redundant since raypack is initialized when running Simulation2D.cast
        angles = np.linspace(angle_min, angle_max, self.n_rays_scan)
        self.cast (*angles, pack=self.pack_temp_scan, target=self.target, n_rebounds_max=self.n_rebounds_max, **self.init_kwargs)
        
        # Sort and select rays
        raypack_temp_scan = self.raypacks[self.pack_temp_scan]
        dist_sorted = raypack_temp_scan.dist_sorted
        if np.isnan(dist_sorted[-1]): dist_sorted = dist_sorted[:-1]  # TODO: why?
        rays = list(itertools.chain.from_iterable( itemgetter(*dist_sorted[:self.n_rays])(raypack_temp_scan.dist) ))  # TODO: what?

        # Generate scan output raypack
        self.raypacks[self.pack_scan] = RayPack2D()
        raypack_scan = self.raypacks[self.pack_scan]
        raypack_scan.add(*rays)
        self.pack_default = self.pack_scan  # Set latest pack as default raypack

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
            pack_prev = self.pack_default
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

                # Cast rays
                angles = np.linspace(ray.angle - cone_half_angle, ray.angle + cone_half_angle, n_rays)
                self.cast (*angles, pack=pack_temp, target=self.target, **self.init_kwargs)

                # Add initial ray (add after cast so that frequency is already generated)
                raypack_temp = self.raypacks[pack_temp]
                raypack_temp.add(ray)
                
                # Select best ray and add to output raypack
                dmin = raypack_temp.dist_sorted[0]
                ray_selected = raypack_temp.dist[dmin][0]  # Select best ray # TODO: only 1? or multiple if same distance?
                raypack_refine.add(ray_selected)

            # Update angular precision (precision cone half-angle)
            self.__angular_precision /= n_rays
            self.dist_avg = np.mean(raypack_refine.dist_sorted)

            if self.verbose: print(f'{self._Simulation2D__vi}Refine #{self.n_refines} mean distance: {self.dist_avg}')  # if self.verbose
            self.pack_default = pack_refine  # Set latest pack as default raypack

    def get_filter_fourier (self, *freqs, **kwargs):
        # New version

        # Load raypack
        pack = kwargs.get('pack', self.pack_default)  # TODO: rename to label
        raypack = self.raypacks [pack]

        # Generate aperture percentages
        angles_sorted = np.sort(list(raypack.angles.keys()))
        diff = np.diff(angles_sorted) / 2
        aperture_max = np.max(angles_sorted) - np.min(angles_sorted)
        aperture_left = np.insert(diff, 0, diff[0])
        aperture_right = np.insert(diff, diff.shape[0], diff[-1])
        aperture_percents = (aperture_left + aperture_right) / aperture_max

        filter_data, offsets = list(), list()
        for i, angle in enumerate(angles_sorted):
            # Populate ray
            ray = raypack.angles[angle]
            ray.populate(*freqs)

            # Get ray data
            offset, filter_func = ray.T_target, ray.calc_Tmult_target  # TODO: move out of Ray2D
            filter_data.append({'offset': offset, 'filter_func': filter_func, 'aperture_mult': aperture_percents[i]})
            offsets.append(offset)
        
        # Process ray time offsets
        offset_min = np.min(offsets)
        offset_max = np.max(offsets)
        time_add = offset_max - offset_min


    def get_filter (self, *freqs, **kwargs):
        """
        :param kwargs/pack: Pack to generate the filter for
        # Signal must begin when fastest ray gets there
        """
        # Load raypack
        pack = kwargs.get('pack', self.pack_default)  # TODO: rename to label
        raypack = self.raypacks [pack]

        # Generate aperture percentages
        angles_sorted = np.sort(list(raypack.angles.keys()))
        diff = np.diff(angles_sorted) / 2
        aperture_max = np.max(angles_sorted) - np.min(angles_sorted)
        aperture_left = np.insert(diff, 0, diff[0])
        aperture_right = np.insert(diff, diff.shape[0], diff[-1])
        aperture_percents = (aperture_left + aperture_right) / aperture_max

        filter_data, offsets = list(), list()
        for i, angle in enumerate(angles_sorted):
            # Populate ray
            ray = raypack.angles[angle]
            ray.populate(*freqs)

            # Get ray data
            offset, filter_func = ray.T_target, ray.calc_Tmult_target  # TODO: move out of Ray2D
            filter_data.append({'offset': offset, 'filter_func': filter_func, 'aperture_mult': aperture_percents[i]})   # TODO: Add offset to phase
            offsets.append(offset)
        
        # Process ray time offsets
        offset_min = np.min(offsets)
        offset_max = np.max(offsets)
        time_add = offset_max - offset_min

        # Create filtering function  # TODO: multi-channel signals  # TODO: standardized filter format?
        def filter (sample_rate, signal):
            n_samples = signal.shape[0]
            data_add = np.ceil(time_add * sample_rate).astype(int)  # Calculate number of datapoints to add
            signal_filtered = np.zeros(n_samples + data_add, dtype=np.float64)

            # Apply Fourier transform to signal
            yf = fft(signal)
            xf = fftfreq(n_samples, 1/sample_rate)

            for ray_filter in filter_data:
                # Calculate offset in datapoints
                time_offset = ray_filter['offset'] - offset_min
                data_offset = np.ceil(time_offset * sample_rate).astype(int)

                # Apply filter on Fourier transform of signal
                filter = ray_filter['filter_func'] (xf)
                yf_filtered = np.multiply(filter, yf)

                # Add temporal filtered signal to result
                sig_filt_ray = ifft(yf_filtered).real
                mult = ray_filter['aperture_mult']
                signal_filtered [data_offset:data_offset+n_samples] += sig_filt_ray * mult  # TODO: not real? => TODO: return the FILTER as a complex numpy array... or as a numpy filter even, if these exist?
            
            return sample_rate, signal_filtered

        return filter        


    def get_filter_NEW (self, *freqs, **kwargs):
        # Load raypack
        label = kwargs.get('pack', self.pack_default)  # TODO: rename to label in kwargs
        raypack = self.raypacks [label]
        rf = np.array(freqs)

        # Generate aperture percentages  # TODO: replace with a fixed precision/resolution aperture (linked to the scanning density)

        # Pregenerate filter data
        angles_sorted = np.sort(list(raypack.angles.keys()))  # TODO: Sorting needed?

        T_target_dict = dict()
        Tmult_target_dict = dict()
        for angle in angles_sorted:
            # Populate ray
            ray = raypack.angles[angle]
            ray.populate_NEW(*freqs)

            # Generate ray filter data
            T_target, Tmult_target = ray.gen_filter(self.target)

            # Save ray filter data
            T_target_dict[angle] = T_target
            Tmult_target_dict[angle] = Tmult_target

        # Process ray time offsets
        T_target_min = np.min(list(T_target_dict.values()))
        T_target_max = np.max(list(T_target_dict.values()))
        T_stretch_total = T_target_max - T_target_min

        # Create filtering function  # TODO: multi-channel signals  # TODO: standardized filter format?
        def filter (sample_rate, signal):
            """_summary_

            Args:
                sample_rate (_type_): Signal sample rate
                signal (_type_): Discrete temporal signal sequence

            Returns:
                _type_: _description_
            """

            n_samples = signal.shape[0]
            data_add = np.ceil(T_stretch_total * sample_rate).astype(int)  # Calculate number of datapoints to add
            signal_filtered = np.zeros(n_samples + data_add, dtype=np.float64)

            # Apply Fourier transform to signal
            yf = fft(signal)
            xf = fftfreq(n_samples, 1/sample_rate)

            for angle in angles_sorted:
                # Unpack filter data
                T_target = T_target_dict[angle]
                Tmult_target = Tmult_target_dict[angle]

                # Calculate offset in datapoints
                time_offset = T_target - T_target_min
                data_offset = np.ceil(time_offset * sample_rate).astype(int)

                # Apply filter on Fourier transform of signal
                Tmult_target_interp = np.interp(xf, np.concatenate((-1 * rf[::-1], rf)), np.concatenate((Tmult_target[::-1], Tmult_target)))
                yf_filtered = np.multiply(Tmult_target_interp, yf)

                # Add temporal filtered signal to result
                sig_filt_ray = ifft(yf_filtered).real
                mult = 0.01  # ray_filter['aperture_mult']  # TODO: aperture multiplier
                signal_filtered [data_offset:data_offset+n_samples] += sig_filt_ray * mult
            
            return sample_rate, signal_filtered

        return filter
