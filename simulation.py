from geometry import Point2D
from environment import Environment2D
from ray import Ray2D



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
        # Plot environment (floor & ceiling for now)
        # Plot rays (or selection of rays)
        pass

    def save (self):
        # save as a file
        pass

    def extract (self):
        # Return new Simulation with only the selected rays
        pass
