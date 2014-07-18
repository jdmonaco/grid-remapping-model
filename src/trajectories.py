# encoding: utf-8
"""
grid.trajectories -- Spatiotemporal trajectories within staging environments.

Copyright (c) 2007, 2008 Columbia University. All rights reserved.
"""

# Library imports
import numpy as N, scipy as S
from scipy.interpolate import interp1d as _i1d

# Package imports
from .stage import StagingMap
from .core.timeseries import TimeSeries

# Traits imports
from enthought.traits.api import (HasTraits, Trait, Constant, Range, Array,
    Float, Int)


class BaseTrajectory(HasTraits):
    
    """
    Superclass for spatiotemporal trajectories. Subclasses should override
    the array constructor __full_traj_default() so that it returns the 
    full trajectory data in a (2, _time.size) matrix.
    
    Public methods:
    advance -- update x and y attributes to next position
    
    Keyword arguments:
    dt -- timestep between contiguous spatial samples
    T -- total duration for the trajectory
    """
    
    dt = Trait(TimeSeries.__class_traits__['dt'])
    T = Trait(TimeSeries.__class_traits__['T'])

    Map = Trait(StagingMap)
    x = Float
    y = Float
    _full_traj = Array
    _time = Array
    _i = Int
    
    def advance(self):
        """Move the current position one step along the trajectory"""
        self._i += 1
        try:
            self.x = self._full_traj[self._i,0]
            self.y = self._full_traj[self._i,1]
        except IndexError:
            pass
                
    def reset(self):
        """Return this trajectory to its initial position"""
        self._i = 0
        self.x = self._x_default()
        self.y = self._y_default()
    
    def _Map_default(self):
        return StagingMap(map_type='TrajectoryMap', quiet=True)
        
    def _x_default(self):
        return self._full_traj[0,0]
    
    def _y_default(self):
        return self._full_traj[0,1]
    
    def __time_default(self):
        return N.arange(0.0, self.T + 5*self.dt, self.dt, 'd')
    
    def __full_traj_default(self):
        """Construct 2 x nTimeSteps array containing trajectory"""
        return N.zeros((2, self.T / self.dt), 'd')


class RandomWalk(BaseTrajectory):

    """
    A smoothed random-walk trajectory.
    
    Keyword arguments:
    v_bar -- mean velocity for the trajectory
    step_freq -- average frequency for random turns
    """
    
    v_bar = Float(15.0, unit='cm/s')    
    step_freq = Float(2.0, unit='Hz')
    
    def __full_traj_default(self):
        step = int(1  / (self.step_freq * self.dt))
        tstep = N.arange(0, self.T + self.dt*step, self.dt*step)
        v_sigma = self.v_bar * self.dt * step
        X = N.empty((tstep.shape[0], 2), 'd')
        X[0] = self.Map.x0
        
        def random_step(x0, x1):
            _angle_ = 2*S.pi * N.random.random_sample()
            _x = x0 + (v_sigma * N.cos(_angle_), v_sigma * N.sin(_angle_))
            while not self.Map.inbounds(_x[0], _x[1]):
                _angle_ = 2*S.pi * N.random.random_sample()
                _x = x0 + (v_sigma * N.cos(_angle_), v_sigma * N.sin(_angle_))
            x1[:] = _x
        
        for t in xrange(1, tstep.shape[0]):
            random_step(X[t-1], X[t])
        
        return N.c_[
            _i1d(tstep, X[:,0], kind='cubic', bounds_error=False, fill_value=X[-1,0])(self._time),
            _i1d(tstep, X[:,1], kind='cubic', bounds_error=False, fill_value=X[-1,1])(self._time)]


class AbstractImpulseRaster(BaseTrajectory): 

    """
    Abstract base class provides functionality for creating x,y trajectories 
    through a StagingMap environment that sequentially 'clamp' on a set of 
    stage pixels for a predetermined 'dwell-time'. 
    
    Subclasses must implement the _get_sample_index method to specify the 
    subset of pixels to clamp.
    """

    dwell = Float(0.2) 
    sample_freq = Range(low=0.01, high=1.0, value=0.1)
    
    _full_raster = Array
    _nsamples = Int
    _transit = Float(10**-7)
    _nsteps = Int
    _req_time = Float
    _init_factor = Float(10)
    
    def get_points(self):
        """
        Convenience method to return a 2-row matrix containing the raster
        points scanned by this trajectory
        """
        return self._full_raster[:, self._get_sample_index()]
    
    def _get_sample_index(self):
        """Subclass provided; return column-index array into stage raster"""
        raise NotImplementedError
    
    def __full_raster_default(self):
        """
        A (2, H*W)-shaped matrix containing full stage raster locations
        """
        Xfull = N.empty((2, self._nsteps), 'd')
        
        # X-values
        Xfull[0] = N.repeat(self.Map._xrange, self.Map._yrange.shape[0])

        # Y-values
        _tmp = N.empty(
            (self.Map._xrange.shape[0], self.Map._yrange.shape[0]), 'd')
        _tmp[:] = self.Map._yrange[N.newaxis]
        Xfull[1] = _tmp.flatten()
        
        return Xfull
    
    def __full_traj_default(self):
        """
        A (2, ntimesteps)-shaped matrix containing the full temporal trajectory
        """
        # Down sample the stage raster according to sample index
        X = self._full_raster[:, N.repeat(self._get_sample_index(), 2)]
        
        # Dwell vector and full-series time vector
        _init = self._init_factor*self.dwell
        _dwell_t = N.linspace(0, self._req_time - _init, int(self._nsteps/2))
        t = N.repeat(_dwell_t, 2)
        t[1::2] += self.dwell - self._transit
        
        # Insert initial dwell time for transients
        t[1:] += _init
        
        return N.c_[
            _i1d(t, X[0], kind='linear', bounds_error=False, 
                fill_value=X[0,-1])(self._time),
            _i1d(t, X[1], kind='linear', bounds_error=False, 
                fill_value=X[1,-1])(self._time)]
    
    def _T_default(self):
        return self._req_time
    
    def __req_time_default(self):
        return (self._init_factor + self._nsamples) * self.dwell
    
    def __nsteps_default(self):
        return self.Map._xrange.shape[0] * self.Map._yrange.shape[0]
    
    def __nsamples_default(self):
        return int(self._nsteps * self.sample_freq)


class BipartiteRaster(AbstractImpulseRaster): 

    """
    Raster-scan stage with clamped input impulses on every other pixel.
    
    Specifically, a bisampled, even partitioning into sampled and non-sampled
    pixels in a checkered pattern.
    
    Keyword arguments:
    dwell -- residence time for each pixel in the scan
    """
        
    sample_freq = Constant(0.5)

    def _get_sample_index(self):
        """
        Pick out bipartite checkered pattern of pixel samples
        """
        odd_cols = self.Map._xrange[1::2]
        s_ix = N.arange(0, self._nsteps, 2)

        for i in xrange(s_ix.size):
            if self._full_raster[0, s_ix[i]] in odd_cols:
                s_ix[i] += 1
                if s_ix[i] >= self._nsteps:
                    s_ix[i] = self._nsteps - 2
        
        return s_ix
