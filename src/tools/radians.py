#encoding: utf-8
"""
radians.py -- A radian numerical type and various angle-handling functions

Written by Joe Monaco
Center for Theoretical Neuroscience
Copyright (c) 2007-2008 Columbia Unversity. All Rights Reserved.  

This software is provided AS IS under the terms of the Open Source MIT License. 
See http://www.opensource.org/licenses/mit-license.php.
"""

# Library imports
from math import atan
from numpy import (pi, dot, cos, sin, fmod, absolute, logical_and, empty, 
    arange, histogram)
from enthought.traits.api import HasTraits, Trait, Property, Float, TraitError

# Constants
TWO_PI = 2*pi


# Validation function for radian values
def radian_domain(object, name, value):
    try:
        value = fmod(float(value), TWO_PI)
    except ValueError:
        raise TraitError
    if value < 0:
        value += TWO_PI
    return value


class radian(HasTraits):
    
    """
    Radian angle numeric type
    """
    
    _r = Trait(0.0, radian_domain)
    degree = Property(Float)
    
    def __init__(self, r0=0.0, **traits):
        HasTraits.__init__(self, **traits)
        self._r = r0
    
    def _get_degree(self):
        return (180/pi) * self._r
    
    def subtract(self, oper):
        return radian(circle_diff(self, oper))
    
    def __add__(self, y):   return radian(self._r + float(y))
    def __sub__(self, y):   return radian(self._r - float(y))
    def __mul__(self, y):   return radian(self._r * float(y))
    def __div__(self, y):   return radian(self._r / float(y))
    def __radd__(self, x):  return radian(float(x) + self._r)
    def __rsub__(self, x):  return radian(float(x) - self._r)
    def __rmul__(self, x):  return radian(float(x) * self._r)
    def __rdiv__(self, x):  return radian(float(x) / self._r)
    def __iadd__(self, o):  self._r += float(o); return self
    def __isub__(self, o):  self._r -= float(o); return self
    def __imul__(self, o):  self._r *= float(o); return self
    def __idiv__(self, o):  self._r /= float(o); return self
    def __eq__(self, y):    return self._r == float(radian(y))
    def __ne__(self, y):    return self._r != float(radian(y))
    def __ge__(self, y):    return self._r >= float(radian(y))
    def __gt__(self, y):    return self._r >  float(radian(y))
    def __le__(self, y):    return self._r <= float(radian(y))
    def __lt__(self, y):    return self._r <  float(radian(y))
    def __abs__(self):      return self
    def __repr__(self):     return str(self)
    def __str__(self):      return str(self._r) + 'r'
    def __float__(self):    return self._r


# Convenience functions for computations on radian angles

def rot2D_vec(v, psi):
    """v, psi => v' rotated in plane by psi radians
    """
    assert v.shape[0] == 2, 'first dimension of points array must have length 2'
    assert v.ndim <= 2, 'array must be a single point or array of points'
    return dot([[cos(psi), -sin(psi)], [sin(psi),  cos(psi)]], v)
    
def rot2D_pt(x, y, psi):
    """x, y, psi => x', y' rotated in plane by psi radians
    """
    return dot([[cos(psi), -sin(psi)], [sin(psi),  cos(psi)]], [x, y])

def circle_diff(a, b):
    """Difference on the semi-circle
    """
    delta = fmod(a, TWO_PI) - fmod(b, TWO_PI)
    mag = abs(delta)
    if mag > pi:
        if delta > 0:
            return delta - TWO_PI
        else:
            return TWO_PI - mag
    else:
        return delta

def circle_diff_vec(u, v):
    """Element-wise circle_diff for arrays of radian values
    """
    delta = fmod(u, TWO_PI) - fmod(v, TWO_PI)
    mag = absolute(delta)
    res = empty(delta.shape, 'd')
    # cond 1: mag > pi, delta > 0
    ix = logical_and(mag>pi, delta>0)
    res[ix] = delta[ix] - TWO_PI
    # cond 2: mag > pi, delta <= 0
    ix = logical_and(mag>pi, delta<=0)
    res[ix] = TWO_PI - mag[ix]
    # cond 3: mag <= pi
    ix = (mag<=pi)
    res[ix] = delta[ix]
    return res

def shortcut(a, b):
    """Smallest angle (<pi) between two radian values
    """
    a = float(radian(a))
    b = float(radian(b))
    delta = abs(a - b)
    if delta > pi:
        delta = TWO_PI - delta
    return delta

def rad_avg(a, b, weight=0.5):
    """Semi-circle average, optionally weighted
    """
    a = float(radian(a))
    b = float(radian(b))
    delta = a - b
    if abs(delta) < pi:
        sc = radian(weight*a + (1-weight)*b)
    elif a < b:
        sc = radian(weight*a + (1-weight)*(b-TWO_PI))
    else:
        sc = radian(weight*(a-TWO_PI) + (1-weight)*b)
    return sc

def xy_to_rad_vec(dx, dy):
    """Vectorized radian conversion, dx/dy are ndarrays
    """
    from numpy import empty
    if type(dx) is float or dx.size == 1:
        return xy_to_rad(dx, dy)
    rad = empty((len(dx),), 'd')
    for i in xrange(len(dx)):
        rad[i] = xy_to_rad(dx[i], dy[i])
    return rad

def xy_to_rad(dx, dy):
    """Radian angle for relative x, y pair
    """
    x, y = float(dx), float(dy)
    if x == 0:
        if y == 0:
            return 0.0
        if y > 0:
            return pi/2
        if y < 0:
            return 3*pi/2
    if x < 0:
        return atan(y/x) + pi
    if y < 0:
        return atan(y/x) + TWO_PI
    return atan(y/x)

def xy_to_deg(dx, dy):
    """Degree angle for relative x, y pair
    """
    return (180/pi)*xy_to_rad(dx, dy)

def get_angle_array(bins, degrees=False, zero_center=False):
    """Get array vector of radian angles for a binned angle range
    """
    amax = degrees and 360.0 or TWO_PI
    if zero_center:
        return arange(-amax/2, amax/2, amax/bins)
    return arange(0, amax, amax/bins)

def get_angle_histogram(x0, y0, bins):
    """Bin x and y offset arrays into a directional histogram on [0, 2PI)
    """
    return histogram(xy_to_rad_vec(x0, y0), bins=bins, 
        range=(0, TWO_PI))[0].astype('d')
