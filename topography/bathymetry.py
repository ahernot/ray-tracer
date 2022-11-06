# Bathymetry (v1.0)
# Copyright Anatole Hernot (Mines Paris), 2022. All rights reserved.

import numpy as np
from scipy import interpolate

import matplotlib.pyplot as plt


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
        self.map_array = map
        self.__x = np.linspace(self.llcorner_rad[0], self.urcorner_rad[0], map.shape[1])
        self.__y = np.linspace(self.llcorner_rad[1], self.urcorner_rad[1], map.shape[0])
        self.__z = map.copy()
        self.__z_interp = interpolate.interp2d(self.__x, self.__y, self.__z, kind='linear', copy=False)
        self.z_from_rad = lambda phi, theta: np.diagonal(self.__z_interp(phi, theta))
        self.z_from_deg = lambda long, lat: self.z_rad(np.deg2rad(long), np.deg2rad(lat))

    def __repr__ (self):
        return f'map with llcorner: {self.llcorner_deg} and urcorner: {self.urcorner_deg}.'

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
        map[map==header['NODATA_VALUE']] = 0.  # np.nan  # TODO: keep np.nan?

        return cls(map, header['XLLCORNER'], header['YLLCORNER'], header['CELLSIZE'])

    def cut (self, start, stop, res):
        """
        Generate a cut from the interpolated map using a SLERP spherical linear interpolation
        # TODO: check if start and end point in map range
        :param start: start point long and lat (in degrees)
        :param stop: stop point long and lat (in degrees)
        """
        p0 = self.llcorner_rad  # theta, phi (spherical coords)
        p1 = self.urcorner_rad  # theta, phi (spherical coords)
        # p0 = np.deg2rad(start)
        # p1 = np.deg2rad(stop)

        # Find angle between the two points (using cartesian coordinates)
        x0 = np.sin(p0[1]) * np.cos(p0[0])
        y0 = np.sin(p0[1]) * np.sin(p0[0])
        z0 = np.cos(p0[1])
        x1 = np.sin(p1[1]) * np.cos(p1[0])
        y1 = np.sin(p1[1]) * np.sin(p1[0])
        z1 = np.cos(p1[1])
        v0 = np.array([x0, y0, z0])
        v1 = np.array([x1, y1, z1])
        theta = np.arccos(np.dot(v0, v1))

        NPOINTS = 100  # TODO: from res (in meters)
        t = np.linspace(0, 1, NPOINTS)
        u = np.sin(theta * (1 - t))
        v = np.sin(theta * t)

        p = (p0 * np.tile(u,(2,1)).T + p1 * np.tile(v,(2,1)).T) / np.sin(theta)
        z = self.__z_interp(p[:,0], p[:,1])
        x = res * t
        z = np.array([float(self.__z_interp(p[i,0], p[i,1])) for i in range(NPOINTS)])

        return x, z

        # TODO: plot cut line onto map to verify the depth data


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
