# profile_impedance v1.0
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np


# Values from https://www.researchgate.net/publication/321661766_Acoustic_impedance_properties_of_seafloor_sediments_off_the_coast_of_Southeastern_Hainan
Z_REFERENCE_VALUES = {
    'clayey silt' : 2390000,
    'silt' : 2550000,
    'sand silt clay' : 2660000,
    'sandy silt' : 2680000,
    'silty clay' : 3160000,
    'rock' : 17000000
}

Zp1_keys = ['clayey silt', 'silt', 'sand silt clay', 'sandy silt', 'silty clay']
Zp1_DEFAULT = np.mean([Z_REFERENCE_VALUES[key] for key in Zp1_keys])
Zp2_DEFAULT = Z_REFERENCE_VALUES['rock']
