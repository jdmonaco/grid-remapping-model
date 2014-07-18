# encoding: utf-8
"""
timeseries.py -- Data collection and time-series tracking

Example:
To manually track a set of evolving variables (e.g.):
   x, y = x0, y0
   ts = TimeSeries(['x', 'y'], dt=0.05, T=200)
   ...
   while ts:
       <modify x>
       <modify y>
       ts.set('x', x)
       ts.set('y', y)
       ts.advance()

Example: 
Automated storage of Traited attributes in a model class:
    x = Trait(..., track=True)
    y = Trait(..., track=True)
    
    def compute_something(self, ...):
        <modify self.x>
        <modify self.y>

    def run(self):
        while self.ts:
            self.compute_something(...)
            self.ts.advance()

The new values of x and y will get stored at the current timepoint in the 
TimeSeries object automagically.

Copyright (c) 2007 Columbia University. All rights reserved.

This software is provided AS IS under the terms of the Open Source MIT License. 
See http://www.opensource.org/licenses/mit-license.php.
"""

# Library imports
from copy import copy as _copy
import numpy as _numpy

# Package imports
from ..tools.array_container import ArrayContainer

# Traits imports
from enthought.traits.api import HasTraits, Range, List, Dict, Array, Int, \
    Float, false


class TSPostMortem(ArrayContainer):
    """Skeleton object for storing time-series data in named attributes
    """
    
    def __init__(self, **data):
        self.__dict__.update(data)
        self.tracking = filter(lambda k: k != 't', data.keys())


class TimeSeries(HasTraits):
    
    """
    Handle collection of time-series data streams, both automatically with
    Traited source objects and manually with any scalar or vector data you
    can come up with.
    
    Three options for required argument source:
       1) a list of data labels
         - set(...) method must be called manually for each timestep
         - only scalar values may be tracked this way
       2) a dictionary of data label keys
         - arrays can be tracked manually
         - source[label] must yield a tuple with the shape of the array
       3) a Traited object instance
         - any Traited attributes with track=True metadata will be tracked
    """ 
    
    traited = false
    data = Dict
    dt = Range(low=0.0, exclude_low=True)
    T = Range(low=10.0)
    i = Int
    t = Float
    time = Array
    _done = false
    _update = List
    
    def __init__(self, source, **traits):
        HasTraits.__init__(self, **traits)
        nTs = self.time.shape[0]
        if isinstance(source, HasTraits):
            self.traited = True
            self.tracking = source.traits(track=True).keys()
            for attr in self.tracking:
                shape = (nTs,)
                if type(getattr(source, attr)) is _numpy.ndarray:
                    shape += getattr(source, attr).shape
                self.data[attr] = _numpy.zeros(shape, 'd')
        else:
            if type(source) is type([]):
                self.tracking = source
                for var in source:
                    self.data[var] = _numpy.zeros((nTs,), 'd')
                self._update = _copy(source)
            elif type(source) is type({}):
                self.tracking = source.keys()
                for var in self.tracking:
                    shape = (nTs,) + tuple(source[var])
                    self.data[var] = _numpy.zeros(shape, 'd')
                self._update = _copy(self.tracking)
            else: 
                raise TypeError, self.__class__.__doc__
        self.source = source

    def __nonzero__(self): 
        return not self._done
    
    def __call__(self, key): 
        return self.get(key)
    
    def __len__(self): 
        return self.T
    
    def __repr__(self): 
        return str(self)
    
    def __str__(self):
        s = "TimeSeries: t=%.2fs to %.2fs @ %1.1fms\nTracking: %s"
        return s%(self.t, self.T, self.dt*1000, self.tracking)
    
    def add_tracker(self, name):
        if i != 0:
            return
        if name not in self.data:
            self.data[name] = _numpy.empty((self.time.shape[0],), 'd')
            self.tracking += [name]
            self._update = _copy(self.tracking)
    
    def copy(self): 
        """Get a fresh TimeSeries object"""
        return TimeSeries(self.source, dt=self.dt, T=self.T)

    def set(self, name, new):
        """Set current value for named data stream"""
        try:
            self.data[name][self.i] = new
        except KeyError: 
            pass
        else: 
            self._update.remove(name)
    
    def get(self, key):
        """Get series vector if done; current value otherwise
        """
        try:
            if self._done:
                retdata = self.data[key]
            elif self.i:
                retdata = self.data[key][self.i - 1]
            else:
                retdata = self.data[key][self.i]
        except KeyError: 
            pass
        else:
            return retdata

    def advance(self):
        """Store tracked data for this timestep and procede to the next
        """
        if self._done: 
            return
        self._advance()
    
    def _advance(self):
        if self.traited:
            for key in self.tracking:
                self.data[key][self.i] = getattr(self.source, key)
        else:
            if self.i:
                for key in self._update:
                    self.data[key][self.i] = self.data[key][self.i-1]
            self._update = _copy(self.tracking)
        if self.i + 1 == self.time.shape[0]:
            self.finish()
        else:
            self.i += 1
            self.t = self.time[self.i]

    def truncate(self): 
        """Synonym for finish()
        """
        self.finish()
    
    def finish(self):
        """Stop collecting series data at current timestep
        """
        self._done = True
        if self.i + 1 != self.time.shape[0]:
            if self.i == 0:
                self.i += 1
            for var in self.data:
                self.data[var] = self.data[var][:self.i]
            self.time = self.time[:self.i]
            self.t = self.T = self.time[-1]
            self.i -= 1

    def data_slice(self, keys, start=0.0, end=None):
        """
        Return time/data-slice
          - keys must be a tracked data label or tuple of labels
          - Use 't' as a label for time data
          - start/end are times in seconds
          - negative start time requests slice from that far back 
            into the past to the current time
        Returns tuple of array slices for each label passed in
        """
        if end == None or end < 0 or end > self.t: 
            end = self.t
        if start < 0.0:
            end = self.t
            if self.t > -start:
                start = self.t + start
            else:
                start = 0.0
        s, e = self.index(start), self.index(end)
        if type(keys) is type(''):
            keys = (keys,)
        slices = []
        for k in keys:
            if k == 't':
                slices.append(self.time[s:e])
            else:
                try:
                    slices.append(self.data[k][s:e])
                except KeyError:
                    pass
        return tuple(slices)
    
    def index(self, t=None):
        """Array index returned for time (t) in seconds
        """
        if t is None:
            return self.i
        return abs(self.time - t).argmin()

    def post_mortem(self):
        """Get a TSPortMortem object for this instance
        """
        if not self._done: 
            return None
        return TSPostMortem(t=self.time, **self.data)
    
    def _time_default(self):
        return _numpy.arange(0.0, self.T + self.dt, self.dt, 'd')
