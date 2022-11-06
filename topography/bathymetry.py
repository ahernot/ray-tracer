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
        self.map_array = map.T[:, ::-1]
        self.__map_nan = np.isnan(self.map_array)  # nan mask
        self.cellsize = cellsize
        self.res = EARTH_RADIUS_AVG * np.deg2rad(self.cellsize)  # Map resolution (in meters)
        self.shape_y, self.shape_x = map.shape
        self.shape = self.map_array.shape

        # Lower-left corner
        self.llcorner_deg = np.array([xllcorner, yllcorner])
        self.llcorner_rad = np.deg2rad(self.llcorner_deg)
        # Upper-right corner
        self.urcorner_deg = self.llcorner_deg + np.array([self.shape_x, self.shape_y], dtype=float) * cellsize
        self.urcorner_rad = np.deg2rad(self.urcorner_deg)
        
        # Map interpolation
        self.__x = np.linspace(self.llcorner_rad[0], self.urcorner_rad[0], self.shape_x)  # Longitude angular values (horizontal axis)
        self.__y = np.linspace(self.llcorner_rad[1], self.urcorner_rad[1], self.shape_y)  # Latitude angular values (vertical axis)
        self.__z = self.map_array.copy()
        self.__z[self.__map_nan] = 0.
        self.__z_interp = interpolate.RectBivariateSpline(self.__x, self.__y, self.__z)  # 2D interpolation
        self.z_from_rad = lambda phi, theta: self.__z_interp(phi, theta)
        self.z_from_deg = lambda long, lat: self.__z_interp(np.deg2rad(long), np.deg2rad(lat))

    def __repr__ (self):
        return f'map with llcorner: {self.llcorner_deg} and urcorner: {self.urcorner_deg}.'  # \nshape: ({self.shape_x}, {self.shape_y}).'

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

    def cut (self, start, stop, npoints=100):
        """
        Generate a cut from the interpolated map using a SLERP spherical linear interpolation
        # TODO: check if start and end point in map range
        :param start: Start point long and lat (in degrees)
        :param stop: Stop point long and lat (in degrees)
        # :param res: Resolution (in meters)  # TODO: optional kwargs param res_min which specifies max dx
        :param npoints: Number of points
        """

        # Convert to spherical coordinates
        p0 = np.deg2rad(start)
        p1 = np.deg2rad(stop)

        # Find angle between start and stop points (using cartesian coordinates)
        x0 = np.sin(p0[1]) * np.cos(p0[0])
        y0 = np.sin(p0[1]) * np.sin(p0[0])
        z0 = np.cos(p0[1])
        x1 = np.sin(p1[1]) * np.cos(p1[0])
        y1 = np.sin(p1[1]) * np.sin(p1[0])
        z1 = np.cos(p1[1])
        v0 = np.array([x0, y0, z0])
        v1 = np.array([x1, y1, z1])
        theta = np.arccos(np.dot(v0, v1))

        # Calculate distances
        d = theta * EARTH_RADIUS_AVG
        res = d / (npoints - 1)
        
        # Map points (SLERP)
        t = np.linspace(0, 1, npoints)
        u = np.sin(theta * (1 - t))
        v = np.sin(theta * t)
        p = (p0 * np.tile(u,(2,1)).T + p1 * np.tile(v,(2,1)).T) / np.sin(theta)  # points

        # Generate interpolated points
        x = res * t
        z = self.__z_interp(p[:,0], p[:,1], grid=False)

        return x, z

    def display (self):
        depth = np.copy(self.map_array)[:,::-1]
        depth[self.__map_nan[:,::-1]] = 0.
        depth = 1. - (np.abs(depth) / np.amax(np.abs(depth)))

        img = np.zeros((self.shape_x, self.shape_y, 3))
        img[:, :, 0] = depth
        img[:, :, 1] = depth
        img[:, :, 2] = depth
        img[self.__map_nan[:,::-1]] = [0.6, 0.3, 0.]

        plt.imshow(np.transpose(img, (1, 0, 2)))
