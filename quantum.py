# Quantum equations library

# Z is the atomic number
# n, l, and m are quantum numbers
# r, t, and p are the spherical coordinates
#     with ranges r:[0,inf), t:[0,pi], p:[0,2*pi] 
import scipy.constants
import math
from sympy import symbols,  Integer, simplify, lambdify
from sympy import summation, factorial as fact, sqrt
from sympy import pi, cos, sin, exp, diff, integrate, I, conjugate
    
a0 = scipy.constants.physical_constants['Bohr radius'][0] # Bohr radius

def args2Integers(*args):
    return [Integer(x) for x in args]

def get_radial_wave_eq(r, n, l):
    # TODO arguments: Z
    n, l = args2Integers(n, l)
    def L(k, n, x):
        # Associated Laguerre polynomial
        m = symbols('m')
        return summation((-1)**m * fact(n+k) / (fact(n-m) * fact(k+m) * fact(m)) * x**m, (m, 0, n))
    a = symbols('a', real=True, positive=True) # Bohr radius
    N = sqrt((2/(n*a))**3 * fact(n-l-1) / (2*n*fact(n+l))) # Normalization
    q = 2 * r / (n * a)
    R = N * q**l * L(2*l+1, n-l-1, q) * exp(-q/2)
    return simplify(R), {a: a0}

def get_inclination_wave_eq(t, l, m):
    l, m = args2Integers(l, m)
    def P(m, l, x):
        # Associated Legendre polynomials
        if m < 0:
            return fact(l+m) / fact(l-m) * P(-1*m, l, x)
        else:
            return(-1)**m/(2**l * fact(l)) * (1-x**2)**(m/2) * diff((x**2-1)**l, x, l+m)
    N = (-1)**m * sqrt((2*l+1)/2*fact(l-m)/fact(l+m)) # Normalization
    return simplify(N * P(m, l, cos(t))), {}

def get_azimuthal_wave_eq(p, m):
    m, = args2Integers(m)
    return 1/sqrt(2*pi) * exp(I*m*p), {}

def get_radial_density_eq(r, n, l):
    r_eq, subs = get_radial_wave_eq(r, n, l)
    r_eq = r_eq.subs(subs)
    return simplify(r_eq * conjugate(r_eq) * r**2)
    
def get_inclination_density_eq(t, l, m):
    t_eq, subs = get_inclination_wave_eq(t, l, m)
    t_eq = t_eq.subs(subs)
    return simplify(t_eq * conjugate(t_eq) * sin(t))
    
def get_azimuthal_density_eq(p, m):
    p_eq, subs = get_azimuthal_wave_eq(p, m)
    p_eq = p_eq.subs(subs)
    return simplify(p_eq * conjugate(p_eq))
    
def get_density_fn(n, l, m):
    r, t, p = symbols('r t p', real=True)
    r_fn = lambdify(r, get_radial_density_eq(r, n, l))
    t_fn = lambdify(t, get_inclination_density_eq(t, l, m))
    p_fn = lambdify(p, get_azimuthal_density_eq(p, m))
    return lambda r, t, p: r_fn(r) * t_fn(t) * p_fn(p)

def get_probability_fn(n, l, m):
    r, t, p = symbols('r t p', real=True)
    r_eq = integrate(get_radial_density_eq(r, n, l), r)
    r_fn = lambdify(r, simplify(r_eq))
    t_eq = integrate(get_inclination_density_eq(t, l, m), t)
    t_fn = lambdify(t, simplify(t_eq))
    p_eq = integrate(get_azimuthal_density_eq(p, m), p)
    p_fn = lambdify(p, simplify(p_eq))
    return lambda r1, r2, t1, t2, p1, p2: (r_fn(r2)-r_fn(r1))*(t_fn(t2)-t_fn(t1))*(p_fn(p2)-p_fn(p1))

def get_radius(probability_fn, percentile=0.99):
    # Returns the radius (multiple of Bohr radius, a0) which encloses
    # the requested percentile
    for i in xrange(0, 10000):
        P = probability_fn(0, i*a0, 0, math.pi, 0, 2*math.pi)
        if P > percentile:
            break
    return i*a0