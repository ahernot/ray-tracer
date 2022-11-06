# Bathymetry (v1.0)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np


EARTH_RADIUS_AVG = 6378000  # in meters (https://en.wikipedia.org/wiki/Earth_radius)


class Map:

    def __init__ (self, map, xllcorner, yllcorner, cellsize):
        """
        :param map:
        :param xllcorner: Lower-left corner longitude  # TODO: rename to ll_long
        :param yllcorner: Lower-left corner latitude  # TODO: rename to ll_lat
        :param cellsize: Cell size (in degrees)
        """
        self.xllcorner = xllcorner
        self.yllcorner = yllcorner
        self.cellsize = cellsize
        self.res = EARTH_RADIUS_AVG * self.cellsize / 180 * np.pi  # Map resolution (in meters)
        self.map = map

    @classmethod
    def from_asc (cls, path: str):

        # Load header
        header_len = 6
        header = dict()
        with open(path, 'r') as map:
            for line in map.readlines()[:header_len]:
                key, val = line.strip().split()
                if key in ('NCOLS', 'NROWS'): val = int(val)
                else: val = float(val)
                header[key] = val
        # Load map
        map = np.loadtxt(path, skiprows=header_len)
        map[map==header['NODATA_VALUE']] = np.nan

        return cls(map, header['XLLCORNER'], header['YLLCORNER'], header['CELLSIZE'])
