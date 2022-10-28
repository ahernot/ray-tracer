# model_reflection v2.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.
# TODO: return complex impedance

import numpy as np

from physics.profile_impedance import Zp1_DEFAULT, Zp2_DEFAULT


def calc_refcoef_surface (wavelength, angle, wave_height_rms):
    """
    Coherent rough surface reflection coefficient
    :param wavelength: Ray wavelength
    :param angle: Grazing angle (in radians) - complementary to the incident angle
    :param wave_height_rms: RMS wave height
    """
    return np.exp (-2 * np.power(2 * np.pi / wavelength * wave_height_rms * np.sin(angle), 2))


def calc_refcoef_sediment (angle, Zp0, Zp1 = Zp1_DEFAULT, Zp2 = Zp2_DEFAULT, dz_sediment = 2):
    """
    
    :param angle: Grazing angle (in radians) - complementary to the incident angle
    :param Zp0: Specific acoustic impedance for compression in the water column (in kg/(s*m**2))
    :param Zp1: Specific acoustic impedance for compression in the sediment layer (in kg/(s*m**2))
    :param Zp2: Specific acoustic impedance for compression in the solid bottom (in kg/(s*m**2))
    :param dz_sediment: Sediment layer thickness (in m)
    """

    # Ignore shear waves
    Zs2 = 0
    theta_s2 = 0

    # Calculate vertical wave number for sediment layer
    gamma_p1 = np.sin(angle)  # sin of grazing angle (?)
    
    # Calculate reflection coefficients
    Zp2_cos2 = Zp2 if theta_s2 == 0             else Zp2 * np.power(np.cos(2 * theta_s2), 2)
    Zs2_sin2 = 0   if theta_s2 == 0 or Zs2 == 0 else Zs2 * np.power(np.sin(2 * theta_s2), 2)
    r01 = (Zp1 - Zp0) / (Zp1 + Zp0)
    r12 = (Zp2_cos2 + Zs2_sin2 - Zp1) / (Zp2_cos2 + Zs2_sin2 + Zp1)

    # Combine reflection coefficients
    r12_exp = r12 * np.exp(-2j * gamma_p1 * dz_sediment)
    refcoef = (r01 + r12_exp) / (1 + r01 * r12_exp)

    return np.abs(refcoef)
