# surface_reflection_model v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np


def calc_refcoef_surface (wavelength, angle, wave_height_rms):
    """
    Coherent rough surface reflection coefficient
    :param wavelength: Ray wavelength
    :param angle: grazing angle
    :param wave_height_rms: RMS wave height
    """
    return np.exp (-2 * np.power(2 * np.pi / wavelength * wave_height_rms * np.sin(angle), 2))
