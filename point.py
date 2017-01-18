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
    
    def __getitem__(self, i):
        if i == 0: return self.x
        if i == 1: return self.y
        if i == 2: return self.z
        else: raise IndexError

    def __setitem__(self, i, a):
        if i == 0: self.x = a
        if i == 1: self.y = a
        if i == 2: self.z = a
        else: raise IndexError
    
    def _mathop(self, fun):
        mapping = map(fun, self)
        if NotImplemented in mapping:
            return NotImplemented
        return Point(*mapping)
    
    def __abs__(self):
        return self._mathop(lambda x: abs(x))
    
    def __add__(self, a):
        return self._mathop(lambda x: x + a)
    
    def __radd__(self, a):
        return self + a
        
    def __sub__(self, a):
        return self._mathop(lambda x: x - a)
    
    def __rsub__(self, a):
        return self - a
    
    def __mul__(self, a):
        return self._mathop(lambda x: x * a)
    
    def __rmul__(self, a):
        return self * a
    
    def __div__(self, a):
        return self._mathop(lambda x: x / a)
    
    def __floordiv__(self, a):
        return self._mathop(lambda x: x // a)
    
    def __mod__(self, a):
        return self._mathop(lambda x: x % a)
    
    def __pow__(self, a):
        return self._mathop(lambda x: x ** a)