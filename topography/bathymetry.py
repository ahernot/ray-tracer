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

        # Resolution
        self.cellsize = cellsize
        self.res = EARTH_RADIUS_AVG * np.deg2rad(self.cellsize)  # Map resolution (in meters)

        # Lower-left corner
        self.llcorner_deg = np.array([xllcorner, yllcorner])
        self.llcorner_rad = np.deg2rad(self.llcorner_deg)
        # Upper-right corner
        self.urcorner_deg = self.llcorner_deg + np.array([map.shape[1], map.shape[0]], dtype=float) * cellsize
        self.urcorner_rad = np.deg2rad(self.urcorner_deg)
        
        # Map interpolation
        self.__x = np.linspace(self.llcorner_rad[0], self.urcorner_rad[0], map.shape[1])
        self.__y = np.linspace(self.llcorner_rad[1], self.urcorner_rad[1], map.shape[0])
        self.__z = map
        self.__z_interp = interpolate.interp2d(self.__x, self.__y, self.__z, kind='linear', copy=False)
        self.z_from_rad = lambda phi, theta: np.diagonal(self.__z_interp(phi, theta))
        self.z_from_deg = lambda long, lat: self.z_rad(np.deg2rad(long), np.deg2rad(lat))

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
        """
        Generate a cut from the interpolated map using a SLERP spherical linear interpolation
        """

        NPOINTS = 100  # from res (in meters)

        # TODO: check if start and end point in map range

        P0_rad = self.llcorner_rad
        P1_rad = self.urcorner_rad
        theta_rad = np.arccos(P0_rad.dot(P1_rad) / (EARTH_RADIUS_AVG**2))

        t = np.linspace(0, 1, NPOINTS)
        u = np.sin(theta_rad * (1 - t))
        v = np.sin(theta_rad * t)

        P_rad = (P0_rad * np.tile(u,(2,1)).T + P1_rad * np.tile(v,(2,1)).T) / np.sin(theta_rad)

        z = self.z_from_rad(P_rad[:,0], P_rad[:,1])

        # Generate x
        x = res * t

        return x, z


# # Plot map
# PLOT = np.empty((map.map.shape[0], map.map.shape[1], 3))

# array = np.copy(map.map)
# array[np.isnan(array)] = 0
# array = np.abs(array) / np.amax(np.abs(array))

# PLOT[:, :, 0] = np.ones(array.shape) - array
# PLOT[:, :, 1] = np.ones(array.shape) - array
# PLOT[:, :, 2] = np.ones(array.shape) - array

# PLOT[np.isnan(map.map)] = [0.6, 0.3, 0]
# plt.imshow(PLOT)
