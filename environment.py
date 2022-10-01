from geometry import Point2D


class Environment2D:

    def __init__ (self, source: Point2D, floor, ceiling):
        self.__source = source

    def __generate_normals (self):
        pass

# class EnvironmentStatic (Environment):
#     def __init__ (self, source: Point2D, floor, ceiling):
#         pass
