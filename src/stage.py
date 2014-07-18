# encoding: utf-8
"""
stage.py -- Handling 2D spatial maps in register with environmental stages

Created by Joe Monaco on 2007-10-23. Updated 2008-04-01; 2009-08-21.
Copyright (c) 2007-2009 Columbia University. All rights reserved.
"""

# Library imports
import numpy as N

# Package imports
from .tools.array_container import TraitedArrayContainer
from .tools.bash import CPrint
    
# Traits imports
from enthought.traits.api import Instance, Int, Array, Tuple, false


class StagingMap(TraitedArrayContainer):

    """
    Matrix-based maps in register with staging area environment
    
    Traited keyword arguments:
    num_maps -- number of computed (computable) maps to handle
    
    Instance attributes:
    Map -- the actual rank-3 ndarray storing the computed maps
    
    Public Methods:
    initialize -- construct the map array based on input data
    index -- get an array (i,j) index from a stage position (x,y)
    map_value -- get an array of each map's value at a stage position (x, y)
    inbounds -- whether an (x,y) position is within stage bounds
    masked_map -- get the ith map with wall pixels set to zero or other value
    
    Subclass override methods:
    _initialize -- construct map array from the input data
    """ 
    
    Map = Array
    W = Int(100)
    H = Int(100)
    num_maps = Int(1)
    x0 = Array
    out = Instance(CPrint)
    _xrange = Array
    _yrange = Array
    _initialized = false

    # Subclass override methods
    
    def _initialize(self):
        """Subclass overload; compute the map data"""
        raise NotImplementedError
        
    # Public methods
    
    def initialize(self):
        """Compute each map by calling the subclass _initialize method"""
        if self._initialized:
            self.out('Maps have already been initialized!', error=True)
            return
        
        # Run the map data initialization method
        self._initialize()
        
        # Set the flag
        self._initialized = True
    
    def index(self, x, y):
        """x, y in map -> index into matrix
        
        Returns (i, j) spatial indices into *Map* array attribute. NOTE: First
        dimension picks out spatial maps.
        
        Raises IndexError if the calculated index is outside the bounds of the 
        spatial dimensions of the staging map matrix.
        """
        i = int(N.ceil(self.H - 1 - y))
        j = int(x)
        if i < 0 or j < 0 or i > self.H-1 or j > self.W-1:
            raise IndexError
        return i, j
        
    def map_value(self, x, y):
        """Get an array of each map's value for a give stage position"""
        i, j = self.index(x, y)
        return self.Map[:, i, j]
    
    def inbounds(self, x, y): 
        """Whether the stage position (x, y) is a wall pixel or not"""
        try:
            test = self.index(x, y)
        except IndexError:
            return False
        else:
            return True
    
    # Traits defaults
    
    def _Map_default(self):
        return N.empty((self.num_maps, self.H, self.W), 'd')

    def _x0_default(self):
        """Start trajectories at mid-point position of square stage"""
        return N.array([0, 0], 'd')
    
    def __xrange_default(self):
        return N.arange(0, self.W)
    
    def __yrange_default(self):
        return N.arange(self.H, 0, -1) - 1
    
    def _out_default(self):
        return CPrint(prefix=self.__class__.__name__)
