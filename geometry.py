
class Point2D:

    def __init__ (self, x, z):
        self.__x = x
        self.__z = z

    def __repr__ (self):
        return f'2D point at coordinates ({self.__x}, {self.__z})'
    
    def coordinates (self):
        return self.__x, self.__z
