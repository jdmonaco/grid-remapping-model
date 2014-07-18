#encoding: utf-8
"""
tools.interp -- Several tools for doing 2D interpolation over data samples

These objects return interpolated values by being called with the query point
after object creation with the set of data samples.

Exported namespace: BilinearInterp2D, GSmoothInterp2D

Created by Joe Monaco, 09-17-2008
(c) Copyright 2008 Columbia University. All Rights Reserved. 
"""

import numpy as N
from enthought.traits.api import HasTraits, Array, Int, Float, Tuple


class BilinearInterp2D(HasTraits):
    
    """
    Grid-based 2D piece-wise bilinear interpolation
    
    Required keyword arguments:
    x -- 1D array of x-values (should be sorted ascending)
    y -- 1D array of y-values (should be sorted ascending)
    z -- flattened 1D array of data values f(x,y) for x,y pairs
    
    Optional keyword arguments:
    fill_value -- out of bounds queries will return this value (0.0)
    """
    
    x = Array
    y = Array
    z = Array
    _x = Array
    _y = Array
    _sx = Int
    _sy = Int
    fill_value = Float
    
    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        xgrid, ygrid = N.meshgrid(self.x, self.y)
        self._x, self._y = xgrid.flatten(), ygrid.flatten()
        self._sx, self._sy = self.x.size, self.y.size
        self._process_z_data()
    
    def __call__(self, xt, yt):
        
        # Handle out of bounds points
        if xt < self.x[0] or xt > self.x[-1]:
            return self.fill_value
        if yt < self.y[0] or yt > self.y[-1]:
            return self.fill_value
        
        # Indices into x and y for lower-left cell vertex
        x_ix = (xt >= self.x).nonzero()[0][-1]
        y_ix = (yt >= self.y).nonzero()[0][-1]
        
        # Force upper and right bounds to be part of existing cells
        if x_ix == self._sx - 1:
            x_ix -= 1
        if y_ix == self._sy - 1:
            y_ix -= 1
        
        # Index into z for lower-left cell vertex
        z_00_ix = x_ix + self._sx * y_ix
        
        # Function values at each corner of cell bounds
        z_00 = self.z[z_00_ix]
        z_10 = self.z[z_00_ix + 1]
        z_01 = self.z[z_00_ix + self._sx]
        z_11 = self.z[z_00_ix + self._sx + 1]
        
        # X- and Y-values of cell bounds 
        x0, x1 = self.x[x_ix], self.x[x_ix + 1]
        y0, y1 = self.y[y_ix], self.y[y_ix + 1]
        
        # Linear weights for query position within cell
        C_x = 1 - (xt - x0) / (x1 - x0)
        C_y = 1 - (yt - y0) / (y1 - y0)
        
        # Return bilinear weighting of grid vertices as interpolated value
        return    C_y  * (C_x * z_00 + (1-C_x) * z_10) + \
               (1-C_y) * (C_x * z_01 + (1-C_x) * z_11)
    
    def _process_z_data(self):
        """
        Take z data passed into constructor and properly flatten it for the
        interpolation algorithm.
        """
        self.z = N.squeeze(self.z)
        if self.z.ndim == 1:
            if self.z.size != self._sx * self._sy:
                raise AttributeError, "incorrect size for flattened z data"
        elif self.z.ndim == 2:
            if self.z.shape != (self._sx, self._sy):
                raise AttributeError, "size mismatch for z matrix data"
            else:
                self.z = self.z.flatten()
        elif self.z.ndim == 3:
            if self.z.shape[1:] != (self._sy, self._sx):
                raise AttributeError, "size mismatch for rank-3 z data"
            else:
                print "Flattening rank-3 array for vectorized interpolation..."
                _sz = self.z.shape[0]
                _z = N.empty((self._sx*self._sy, _sz), 'd')
                ix = 0
                for i in xrange(self._sy):
                    for j in xrange(self._sx):
                        _z[ix] = self.z[:,i,j]
                        ix += 1
                self.z = _z
        else:
            raise AttributeError, "z data must be rank-3 or less"
        return
                      

class GSmoothInterp2D(HasTraits):
    
    """
    A simple but robust 2D interpolation from random grid samples based on
    a Gaussian kernel averaging nearest neighbors.
    
    If the 2D range of interpolation is not 1:1, you may specify the *aspect* 
    keyword argument. This will correct for unequal x and y scaling. Make sure
    that *k* is in units relative to the y-range.
    
    If Voronoi tiling is apparent, k is too small. If details are washed over, 
    k is too large. Generally, k ~ 0.1 * (y.max()-y.min()).
    """ 
    
    x = Array
    y = Array
    z = Array
    k = Float(1)
    aspect = Float(1)
    neighbors = Int(10)
    
    def __call__(self, xt, yt):
        
        # Distances square
        d = N.sqrt(((self.x - xt)/self.aspect)**2 + (self.y - yt)**2)
        
        # Indices for nearest neighbors
        ix = N.argsort(d)[:self.neighbors]
        
        # Gaussian coefficients
        c = N.exp(-(d[ix]+self.k)**2/(2*self.k**2))

        # Weighted neighborhood average
        return (c*self.z[ix]).sum() / c.sum()
