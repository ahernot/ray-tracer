



class Point2D:

    def __init__ (self, x, z):
        self.__x = x
        self.__z = z

    def __repr__ (self):
        return f'2D point at coordinates ({self.__x}, {self.__z})'
    
    def coordinates (self):
        return self.__x, self.__z




class Environment:

    def __init__ (self, source: Point2D, floor, ceiling):
        self.__source = source

    def __generate_normals (self):
        pass

# class EnvironmentStatic (Environment):
#     def __init__ (self, source: Point2D, floor, ceiling):
#         pass


class Ray:

    def __init__ (self, env: Environment):
        self.__env = env

    def propagate (self):
        pass


class Simulation:

    def __init__ (self, env: Environment, *rays: Ray):
        self.__rays = rays  # ordered by angle??
        self.__env = env

    pass