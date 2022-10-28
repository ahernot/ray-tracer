# model_rho v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np


def calc_rho (S, z):
    """
    IES 80 – High Pressure International Equation of State of Seawater : https://unesdoc.unesco.org/ark:/48223/pf0000047363
    :param S: Practical salinity (‰), range:   0‰ ≤ S ≤ 42‰
    :param z: Depth (m), range: -1000 ≤ z ≤ 0
    :returns: Seawater density (in kg.m^-3)
    """

    T = calc_T(z)

    #pressure linearily increases with depth
    p = calc_p(z)  # pressure in bar
 
    # Density of the Standard Mean Ocean Water (SMOW) [Bigg, 1967]
    rho_W = 999.842594 + 6.793952e-2 * T - 9.09529e-3 * np.power(T, 2) + 1.001685e-4 * np.power(T, 3) - 1.120083e-6 * np.power(T, 4) + 6.536336e-9 * np.power(T, 5)
    
    # Coefficients
    A = 8.24493e-1 - 4.0899e-3 * T + 7.6438e-5 * np.power(T, 2) - 8.2467e-7 * np.power(T, 3) + 5.3875e-9 * np.power(T, 4)
    B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * np.power(T, 2)
    C = 4.8314e-4

    # One Atmosphere International Equation of State of Seawater (1980) (standard error = 3.6e-3 kg.m-3; validity 0‰≤S≤42‰, -2℃≤T≤40℃)
    rho_1atm = rho_W + A * S + B * np.power(S, 1.5) + C * np.power(S, 2)

    # Pure water terms
    K_W = 19652.21 + 148.4206 * T - 2.327105 * np.power(T, 2) + 1.360477e-2 * np.power(T, 3) - 5.155288e-5 * np.power(T, 4)
    A_W = 3.239908 + 1.43713e-3 * T + 1.16092e-4 * np.power(T, 2) - 5.79905e-7 * np.power(T, 3)
    B_W = 8.50935e-5 - 6.12293e-6 * T + 5.2787e-8 * np.power(T, 2)

    K_1atm = K_W + \
        (54.6746 - 0.603459 * T + 1.09987e-2 * np.power(T, 2) - 6.1670e-5 * np.power(T, 3)) * S + \
        (7.944e-2 + 1.6483e-2 * T - 5.3009e-4 * np.power(T, 2)) * np.power(S, 1.5)
    A = A_W + \
        (2.2838e-3 - 1.0981e-5 * T - 1.6078e-6 * np.power(T, 2)) * S + \
        1.91075e-3 * np.power(S, 1.5)
    B = B_W + \
        (-9.9348e-7 + 2.0816e-8 * T + 9.1697e-10 * np.power(T, 2)) * S

    # Secant Bulk Modulus
    K = K_1atm + A * p + B * np.power(p, 2)

    # High Pressure International Equation of State of Seawater
    rho = rho_1atm / (1 - p / K)
    
    return rho
