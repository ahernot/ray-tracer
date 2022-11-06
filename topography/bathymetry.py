# Bathymetry (v1.0)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np
from scipy import interpolate


EARTH_RADIUS_AVG = 6378000  # in meters (https://en.wikipedia.org/wiki/Earth_radius)


class Map:

    def __init__ (self, map, xllcorner, yllcorner, cellsize):
        """
        :param map:
        :param xllcorner: Lower-left corner longitude  # TODO: rename to ll_long
        :param yllcorner: Lower-left corner latitude  # TODO: rename to ll_lat
        :param cellsize: Cell size (in degrees)
        """
        # self.xllcorner = xllcorner
        # self.yllcorner = yllcorner

        self.llcorner_deg = np.array([xllcorner, yllcorner])
        self.urcorner_deg = self.llcorner_deg + np.array([map.shape[1], map.shape[0]], dtype=float) * cellsize

        self.llcorner_rad = self.llcorner_deg / 180 * np.pi
        self.urcorner_rad = self.urcorner_deg / 180 * np.pi

        self.__x = np.linspace(self.llcorner_rad[0], self.urcorner_rad[0], map.shape[1])
        self.__y = np.linspace(self.llcorner_rad[1], self.urcorner_rad[1], map.shape[0])
        self.__z = map

        print(self.__x.shape)
        print(self.__y.shape)
        print(self.__z.shape)

        self.cellsize = cellsize
        self.res = EARTH_RADIUS_AVG * self.cellsize / 180 * np.pi  # Map resolution (in meters)

        self.z = interpolate.interp2d(self.__x, self.__y, self.__z, kind='linear', copy=False)

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

    def cut (self, start_lat, start_long, stop_lat, stop_long, res):

        # SLERP


        pass
