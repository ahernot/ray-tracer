
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



def ccw (A: Point2D, B: Point2D, C: Point2D):
    return (C.z-A.z) * (B.x-A.x) > (B.z-A.z) * (C.x-A.x)

# Return true if line segments AB and CD intersect
def intersect (A: Point2D, B: Point2D, C: Point2D, D: Point2D):
    return ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)
