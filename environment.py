from geometry import Point2D, Vector2D


class Environment2D:

    def __init__ (self, floor, ceiling):
        # self.__source = source
        # function to model ray speed as a function of x and z
        # current as additional displacement (vector field)

        #self.range = {'x': (0, 100000), 'z': (-10000, 10000)}
        self.range_min = Vector2D (0, -10000)
        self.range_max = Vector2D (100000, 10000)


        # reflection power coefficient as a function of distance

    def __generate_normals (self):
        pass

# class EnvironmentStatic (Environment):
#     def __init__ (self, source: Point2D, floor, ceiling):
#         pass
