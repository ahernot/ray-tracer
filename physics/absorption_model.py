# absorption_model v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np


def calc_absorption_dB (f: float, z: float, T: float, S: float, pH: float):
    """
    Ainslie and McColm 1998
    :param f: Signal frequency (in Hz)
    :param z: Altitude (in m)
    :param temperature: Temperature (in ℃)
    :param salinity: Salinity (in ppt)
    :param pH: Acidity (pH)
    """
    # http://resource.npl.co.uk/acoustics/techguides/seaabsorption/
    # http://resource.npl.co.uk/acoustics/techguides/seaabsorption/physics.html

    f /= 1000

    # Boric acid (relaxation absorption)
    f1 = 0.78 * np.sqrt(S / 35) * np.exp(T / 26)  # Boric acid relaxation frequency (kHz)
    A1 = 0.106
    alpha_1 = A1 * (f1 * np.power(f, 2)) / (np.power(f1, 2) + np.power(f, 2))

    # Magnesium sulfate (relaxation absorption)
    f2 = 42 * np.exp(T / 17)  # Magnesium sulfate relaxation frequency (kHz)
    A2 = 0.52 * (1 + T / 43) * (S / 35)
    alpha_2 =  A2 * (f2 * np.power(f, 2)) / (np.power(f2, 2) + np.power(f, 2))

    # Pure water (viscous absorption)
    A3 = 0.00049
    alpha_3 = A3 * np.power(f, 2)

    # Total seawater absorption (in dB/km)
    alpha = alpha_1 * np.exp((pH - 8) / 0.56) +\
            alpha_2 * np.exp(z / 6000)  +\
            alpha_3 * np.exp(z / 17000 - T / 27)


    return alpha

# output in dB/m or dB/km???
