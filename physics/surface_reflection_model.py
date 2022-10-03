# surface_reflection_model v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np


def calc_refcoef_surface (wavelength, wave_height, grazing_angle):
    """
    Coherent rough surface reflection coefficient
    """
    return np.exp (-2 * np.power(2 * np.pi / wavelength * wave_height * np.sin(grazing_angle), 2))
