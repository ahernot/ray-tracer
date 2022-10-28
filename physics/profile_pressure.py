# profile_pressure v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np


def calc_P (z):
    """
    Calculate applied pressure
    """
    return np.abs(z) * 0.1
