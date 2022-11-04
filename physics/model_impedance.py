# model_impedance v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

from physics.model_rho import calc_rho
from physics.profile_salinity import calc_S
from physics.profile_velocity import calc_c
from physics.profile_salinity import calc_S
from physics.profile_temperature import calc_T
from physics.profile_pressure import calc_P


# def calc_Z (z):
#     c = calc_c(z)
#     rho = calc_rho(z, calc_S(z), calc_T(z), calc_P(z))  # TODO: put all the calc_() functions in a class so they can all access each other's z etc
#     return c * rho

def calc_Z (rho, c) :
    return rho * c
