from ray import Ray2D

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








class Simulation2D:

    def __init__ (self, env: Environment2D):
        self.__env = env
        self.__rays = list()  # Raw list of rays (stored as references to Ray2D objects)
        self.__angles = dict()  # Dictionary of Ray2D objects grouped by angle
        # can group rays by frequency or some other key if needed
        self.n_rays = 0
        self.n_angles = 0

    def __repr__ (self):
        return f'2D simulation containing {self.n_rays} rays'  # TODO: to improve later

    def add_rays (self, *angles):
        # Supports two rays for the same angle (for example two different frequencies)

        for angle in angles:
            # Generate ray
            ray = Ray2D (env=self.__env, angle=angle)  # more params
            ray.propagate()  # feed propagation params for the simulation

            # Add ray to simulation
            self.__rays .append(ray)
            self.__angles[angle] = self.__angles[angle] + [ray] if angle in self.__angles else [ray]
        
        self.n_rays = len(self.__rays)
        self.n_angles = len(self.__angles)

    def plot (self):
        pass

    def save (self):
        # save as a file
        pass
