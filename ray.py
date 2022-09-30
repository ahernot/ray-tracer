

class Point2D:

    def __init__ (self, x, z):
        self.__x = x
        self.__z = z

    def __repr__ (self):
        return f'2D point at coordinates ({self.__x}, {self.__z})'
    
    def coordinates (self):
        return self.__x, self.__z




class Environment2D:

    def __init__ (self, source: Point2D, floor, ceiling):
        self.__source = source

    def __generate_normals (self):
        pass

# class EnvironmentStatic (Environment):
#     def __init__ (self, source: Point2D, floor, ceiling):
#         pass





class Ray2D:

    def __init__ (self, env: Environment2D):
        self.env = env

    def propagate (self):
        pass


class Simulation2D:

    def __init__ (self, env: Environment2D):
        self.__env = env
        self.__rays = list()

    def add_ray (self):
        pass

