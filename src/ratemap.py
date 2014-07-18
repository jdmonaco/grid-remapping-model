# encoding: utf-8
"""
ratemap.py -- Objects based on PlaceMap for constructing spatial ratemaps

Branched by Joe Monaco on 09-04-2008.
Copyright (c) 2008 Columbia University. All rights reserved.
"""

# Library imports
import numpy as N
from scipy import signal
from enthought.traits.api import Array

# Package imports
from .placemap import PlaceMap


class AbstractImpulseRatemap(PlaceMap):
    
    """
    Abstract base class with functionality for converting network output from 
    simulations based on impulse rasters (pixel-by-pixel dwell) into usable 
    ratemap objects.
    """
    
    PV = Array
    points = Array
    
    def _PV_default(self):
        return self.data.r
    
    def _points_default(self):
        return N.c_[self.data.x, self.data.y]
    

class CheckeredRatemap(AbstractImpulseRatemap):
    
    """
    CheckeredRatemap handles automatic creation of empirical ratemaps based on 
    network output as driven by a BipartiteRaster trajectory.
    """
    
    def _initialize(self):
        """
        Perform z-vector neighbor-averaging to compute full ratemap stacks
        """
        
        # Set map data for directly probed raster pixels
        self.out('Initializing...')
        self.Map[:] = 0.0
        for ix, p in enumerate(self.points):
            i, j = self.index(*p)
            self.Map[:, i, j] = self.PV[ix]
        
        # Raster scan stage to interpolate non-probed pixels
        for x in self._xrange:
            for y in self._yrange:
                i, j = self.index(x, y)
                
                # Skip directly probed pixels
                if N.fmod(j + N.fmod(i, 2) + 1, 2):
                    continue
                
                # Manually average on-axis neighbors
                n = 0
                S = self.Map[:, i, j]
                if i > 0:
                    S += self.Map[:, i-1, j]
                    n += 1
                if i < self.H - 1:
                    S += self.Map[:, i+1, j]
                    n += 1
                if j > 0:
                    S += self.Map[:, i, j-1]
                    n += 1
                if j < self.W - 1:
                    S += self.Map[:, i, j+1]
                    n += 1
                S /= n
        
        # Median-filter to despeckle and smooth
        self.Map = signal.medfilt(self.Map, kernel_size=[1,3,3])
        
        # All done!
        self.out('Done!')
