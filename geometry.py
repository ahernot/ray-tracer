
class Point2D:

    def __init__ (self, x, z):
        self.x = x
        self.z = z

    def __repr__ (self):
        return f'2D point at coordinates ({self.x}, {self.z})'
    
    def coordinates (self):
        return self.x, self.z


class Vector2D:

    def __init__ (self, x, z):
        self.x = x
        self.z = z
    
    def __repr__ (self):
        return f'2D vector of coordinates ({self.x}, {self.z})'
    
    def coordinates (self):
        return self.x, self.z
