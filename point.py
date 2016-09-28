import math
import numpy as np

def args2floats(*args):
    return [float(x) for x in args]

def cart2sph(x, y, z):
    # Returned range: r: (0, inf), t: (0, pi), p: (0, 2*pi)
    x, y, z = args2floats(x, y, z)
    r = math.sqrt(x**2 + y**2 + z**2)
    t = 0 if r==0 else math.acos(z/r)
    p = math.atan2(y, x)%(2*math.pi)
    return r, t, p

def sph2cart(r, t, p):
    x = r * math.sin(t) * math.cos(p)
    y = r * math.sin(t) * math.sin(p)
    z = r * math.cos(t)
    return x, y, z

class Point:
    def __init__(self, c1, c2, c3, spherical=False):
        if spherical:
            self.r, self.t, self.p = c1, c2, c3
            self.x, self.y, self.z = sph2cart(c1, c2, c3)
        else:
            self.x, self.y, self.z = c1, c2, c3
            self.r, self.t, self.p = cart2sph(c1, c2, c3)
    
    def __repr__(self):
        return '({}, {}, {})'.format(self.x, self.y, self.z)