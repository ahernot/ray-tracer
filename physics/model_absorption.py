# Absorption Model v2.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np

# Sound absorption in seawater
def calc_absorption_dB (f: float, z: float, T: float, S: float, pH: float):  # TODO: rename calc_dz_dG
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

    f /= 1000   # to get f in kHz for the formula

    # Boric acid (relaxation absorption)
    f1 = 0.78 * np.sqrt(S / 35) * np.exp(T / 26)  # Boric acid relaxation frequency (kHz)
    A1 = 0.106
    alpha_1 = A1 * (f1 * np.power(f, 2)) / (np.power(f1, 2) + np.power(f, 2))

    # Magnesium sulfate (relaxation absorption)
    f2 = 42 * np.exp(T / 17)  # Magnesium sulfate relaxation frequency (kHz)
    A2 = 0.52 * (1 + T / 43) * (S / 35)
    alpha_2 = A2 * (f2 * np.power(f, 2)) / (np.power(f2, 2) + np.power(f, 2))

    # Pure water (viscous absorption)
    A3 = 0.00049
    alpha_3 = A3 * np.power(f, 2)

    # Total seawater absorption (in dB/km)
    alpha = alpha_1 * np.exp((pH - 8) / 0.56) +\
            alpha_2 * np.exp(z / 6000) +\
            alpha_3 * np.exp(z / 17000 - T / 27)


    return alpha  # dB/km




def calc_dz_dG_coefs (z, T, S, pH):
        f1 = 0.78 * np.sqrt(S / 35) * np.exp(T / 26)  # Boric acid relaxation frequency (kHz)
        f2 = 42 * np.exp(T / 17)  # Magnesium sulfate relaxation frequency (kHz)
        A1 = 0.106
        A2 = 0.52 * (1 + T / 43) * (S / 35)
        A3 = 0.00049
        E1 = np.exp((pH - 8) / 0.56)
        E2 = np.exp(z / 6000)
        E3 = np.exp(z / 17000 - T / 27)
        return f1, f2, A1, A2, A3, E1, E2, E3


def calc_dz_dG (f: float, z: float, calc_dz_dG_coefs_interp):
    f1, f2, A1, A2, A3, E1, E2, E3 = calc_dz_dG_coefs_interp(z)
    f /= 1000   # to get f in kHz for the formula

    alpha_1 = A1 * (f1 * np.power(f, 2)) / (np.power(f1, 2) + np.power(f, 2))  # Boric acid (relaxation absorption)
    alpha_2 = A2 * (f2 * np.power(f, 2)) / (np.power(f2, 2) + np.power(f, 2))  # Magnesium sulfate (relaxation absorption)
    alpha_3 = A3 * np.power(f, 2)  # Pure water (viscous absorption)

    # Total seawater absorption (in dB/km)
    alpha = alpha_1 * E1 +\
            alpha_2 * E2 +\
            alpha_3 * E3

    return alpha
