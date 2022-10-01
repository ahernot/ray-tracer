from environment import Environment2D

class Ray2D:

    def __init__ (self, env: Environment2D, angle):
        self.__env = env
        self.__angle = angle

    def __repr__ (self):
        return f'Ray object'  # TODO: improve repr

    def propagate (self):  # propagation params
        
        # Minimum resolutions
        dx_max = 0.1
        dz_max = 0.1
        # adaptive dz based on angle




        
        pass
