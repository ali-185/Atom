import numpy as np
import math
from point import Point

class Surface:
    def __init__(self, data, simplify=False):
        # Data is a 2d np.array of Points or None
        self.data = data
        self.height = self.data.shape[0]
        self.width  = self.data.shape[1]
        if simplify:
            self.simplify()
    
    def __getitem__(self, arg):
        return self.data[arg]
    
    def x(self):
        # Return 2d np.array of the x values
        return np.vectorize(lambda p: np.nan if p == None else p.x)(self.data)
    
    def y(self):
        # Return 2d np.array of the y values
        return np.vectorize(lambda p: np.nan if p == None else p.y)(self.data)
    
    def z(self):
        # Return 2d np.array of the z values
        return np.vectorize(lambda p: np.nan if p == None else p.z)(self.data)
    
    def join(self, surface):
        # Returns the joining Surface
        if (self.height, self.width) != (surface.height, surface.width):
            raise ValueError("Surfaces must have the same height and width.")
        edges1 = self.edges()
        edges2 = surface.edges()
        join = np.full((2, self.height*(self.width+1)), None, dtype=np.object)
        for i in xrange(0, self.height):
            join[0,i*(self.width+1):(i+1)*(self.width+1)-1] = edges1[i]
            join[1,i*(self.width+1):(i+1)*(self.width+1)-1] = edges2[i]
        return Surface(join, True)
    
    def edges(self):
        # Returns a Surface of only the edge Points
        nan_map = np.equal(self.data, None)
        edge_map = np.full((self.height, self.width), False, dtype=bool)
        for roll, axis in [(1, 1), (-1, 1), (1, 0), (-1, 0)]:
            edge_map += np.logical_and(~nan_map, np.roll(nan_map, roll, axis=axis))
        edge = np.full((self.height, self.width), None, dtype=np.object)
        edge[edge_map] = self.data.copy()[edge_map]
        return Surface(edge)
    
    def simplify(self):
        # Removes superfluous data (nans)
        nan_map = np.equal(self.data, None)
        for axis in [0, 1]:
            nan_arr = np.all(nan_map, axis=axis)
            nan_roll = np.logical_and(np.roll(nan_arr, 1), np.roll(nan_arr, -1))
            nan_edge = np.logical_and(nan_roll, nan_arr)
            if axis == 0:
                self.data = self.data[:, ~nan_edge]
            elif axis == 1:
                self.data = self.data[~nan_edge]
    
def create_surfaces(fn, max_r, res):
    # Generates a list of surfaces from a function. These surfaces will enclose the
    # True values. The function must take r, t, p as arguments.
    # Ray cast radially to find create surfaces
    surfaces = []
    for i, t in enumerate(np.linspace(0, math.pi, res)):
        for j, p in enumerate(np.linspace(0, 2*math.pi, 2*res)):
            prev = None
            lvl = 0
            for r in np.linspace(0, max_r, res):
                pres = fn(r, t, p)
                if prev != None and pres != prev:
                    for _0 in 1,2:
                        try:
                            surface = surfaces[lvl]
                            surface[i][j] = Point(r, t, p, spherical=True)
                            lvl += 1
                            break
                        except IndexError:
                            surface = Surface(np.full((res, 2*res), None, dtype=np.object))
                            surfaces.append(surface)
                prev = pres
    # Connect the surfaces
    for surface1, surface2 in zip(surfaces[-1::-2], surfaces[-2::-2]):
        surfaces.append(surface1.join(surface2))
    return surfaces