# profile_salinity v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

# NO MODEL

import numpy as np


S = 35

def calc_S (z):
    return np.ones_like(z) * S
