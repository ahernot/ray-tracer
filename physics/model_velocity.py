# velocity_models v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np


def sound_velocity_simplified (S, T, z):
    """
    Sound velocity in seawater (simplified formula)
    http://lecalve.univ-tln.fr/oceano/fiches/fiche3F.htm

    :param S: Practical salinity (in ‰)
    :param T: Temperature (in ℃)
    :param z: Altitude (in m)
    :return: Sound velocity (in m.s^-1)
    """

    c = 1410 + 4.21 * T - 0.037 * np.power(T, 2) + 1.1 * S - 0.018 * z
    
    return c


def sound_velocity_leroy (S, T, z):
    """
    Sound velocity in seawater (Claude Leroy)
    http://lecalve.univ-tln.fr/oceano/fiches/fiche3F.htm

    :param S: Practical salinity (in ‰), range:     30‰ ≤ S ≤ 42‰
    :param T: Temperature (in ℃),        range:     -2℃ ≤ T ≤ 24.5℃
    :param z: Altitude (in m),           range: -1000m ≤ z ≤ 0m
    :return: Sound velocity (in m.s^-1)
    """

    c = 1492.9 + \
        3 * (T-10) - 0.006 * np.power(T-10, 2) - 0.04 * np.power(T-18, 2) + \
        1.2 * (S-35) - 0.01 * (T-18) * (S-35) - \
        z / 61
    
    return c
    

def sound_velocity_medwin (S, T, z):
    """
    Sound velocity in seawater (H. Medwin)
    http://lecalve.univ-tln.fr/oceano/fiches/fiche3F.htm

    :param S: Practical salinity (in ‰), range:     0‰ ≤ S ≤ 45‰
    :param T: Temperature (in ℃),        range:     0℃ ≤ T ≤ 35℃
    :param z: Altitude (in m),           range: -1000m ≤ z ≤ 0m
    :return: Sound velocity (in m.s^-1)
    """

    c = 1449.2 + \
        4.6 * T - 0.055 * np.power(T, 2) + 0.00029 * np.power(T, 3) + \
        (1.34 - 0.01 * T) * (S - 35) - \
        0.016 * z 
    
    return c


def sound_velocity_mackenzie (S, T, z):
    """
    Sound velocity in seawater (Mackenzie)
    http://lecalve.univ-tln.fr/oceano/fiches/fiche3F.htm

    :param S: Practical salinity (in ‰), range:     0‰ ≤ S ≤ 40‰
    :param T: Temperature (in ℃),        range:     0℃ ≤ T ≤ 30℃
    :param z: Altitude (in m),           range: -1000m ≤ z ≤ 0m
    :return: Sound velocity (in m.s^-1)
    """

    c = 1448.96 + \
        4.591 * T - 0.05304 * np.power(T, 2) + 0.0002374 * np.power(T, 3) + \
        1.34 * (S-35) - 0.01025 * T * (S-35) - \
        0.016 * z + 0.0000001675 * np.power(z, 2)
    
    return c
