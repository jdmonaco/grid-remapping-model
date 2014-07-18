#encoding: utf-8
"""
grid.place_network -- Place network model driven by grid-cell inputs, with 
    competition mediated by global inhibtion

Written by Joe Monaco
Copyright (c) 2007, 2008 Columbia University. All rights reserved.
"""

# Library imports
import numpy as N

# Package imports
from .dmec import GridCollection
from .core.api import AbstractModel, TSPostMortem
from .trajectories import BaseTrajectory, RandomWalk, BipartiteRaster
from .tools.integrate import integrator

# Traits imports
from enthought.traits.api import Delegate, Instance, Property, Trait, Array, \
    Range, Enum, Float, CInt, Int, true, false


class PlaceNetwork(AbstractModel):
    
    """
    Simple CA3 model with global inhibition; no excitatory recurrence
    
    Model parameters and other attributes are Traits, so values should be set
    as keywork argumnets to the class constructor.
    
    Trial setup computes a new weight matrix (self.W) only if refresh_weights
    is True. The MEC input is modified according to refresh_orientation, 
    refresh_phase, phi and psi. The phi and psi Property traits allow direct
    specification of cortical input alignment arrays (you should also then
    set refresh_orientation and/or refresh_phases to False).
    
    Model settings and parameters:
    EC -- GridCollection instance to be used as input (required)
    traj_type -- type of trajectory to traverse (default RandomWalk)
    phi_lambda -- field nonlinearity threshold (default 1.5)
    phi_sigma -- field nonlinearity smoothness (default 0.1)
    J0 -- gain of global inhibition (default 2.5)
    tau_r -- time constant for input integration (default 0.05)
    C_W -- fraction of afferent connectivity (default 0.5)
    mu_W -- mean of afferent weight distribution (default 0.5)
    dwell_factor -- duration of pixel dwell in tau (default 5.0)
    
    Simulation control settings (per trial):
    refresh_traj -- new trajectory (default False)
    refresh_weights -- new network weights (default True)
    refresh_orientation -- randomize cortical orientation array (default False)
    refresh_phase -- new random array of cortical input phases (default False)
    
    Cortical input alignment:
    phi -- directly set the cortical phase array; shape: (2, num_maps)
    psi -- directly set the cortical orientation array; shape: (num_maps,)
    """
    
    label = 'Grid Model'
    app_name = label
    evolve = Instance(integrator)
    
    # Trajectory traits
    traj_type = Trait(['randwalk', 'checker'], user=True)
    traj = Instance(BaseTrajectory)
    refresh_traj = false(user=True)
    dwell_factor = Range(low=0.0, value=5.0, exclude_low=True, user=True)
    x = Property(Float, track=True)
    y = Property(Float, track=True)
    
    # Network and input traits
    EC = Trait(None, Instance(GridCollection), user=True)
    get_afferent_input = Delegate('EC', prefix='map_value')
    N_EC = Delegate('EC', prefix='num_maps')
    N_CA = CInt(500, user=True)
    refresh_orientation = false(user=True)
    refresh_phase = false(user=True)
    
    # Nonlinearity definition
    phi_lambda = Float(0.2, user=True)
    phi_sigma = Range(low=0.0, value=0.015, exclude_low=True, user=True)
    
    # Weights and rates traits :-)
    W = Array
    refresh_weights = true(user=True)
    C_W = Range(low=0.0, high=1.0, value=0.5, user=True)
    mu_W = Float(0.5, user=True)
    r = Array(track=True)
    r_EC = Array
    i_aff = Array
    tau_r = Float(0.05, user=True)
    dt = 0.005
    
    # Synaptic gains
    J0 = Float(250, user=True)
    beta = Float
    
    # AbstractModel override methods
    
    def trial_setup(self):
        self.r = self._r_default()
        self.beta = self._beta_default()
        self.evolve = integrator(self.drdt, self.r, dt=self.dt)
        if self.refresh_traj:
            self.out('Creating new stage trajectory')
            self.traj = self.new_trajectory()
        else:
            self.out('Reseting current stage trajectory')
            self.traj.reset()
        if self.refresh_weights:
            self.out('Computing a new weight matrix')
            self.W = self.new_weights()
        if self.refresh_phase:
            self.EC.randomize_phase()
            self.out('Computed new spatial phase vector')
        if self.refresh_orientation:
            self.EC.randomize_orientation()
            self.out('Computed new orientation vector')

    def run_timestep(self):
        """Simulation time-step computation
        """
        self.r_EC = self.get_afferent_input(self.x, self.y) # MEC input
        self.i_aff = N.dot(self.W, self.r_EC)               # afferent current
        self.evolve()                                       # evolve rates
        self.traj.advance()                                 # move trajectory
    
    def drdt(self, r, t0):
        """Rate equation: dr/dt = (Phi[h - lambda] - r)/tau
        """
        return (self.phi_h( self.beta * self.i_aff -
                            self.J0 * r.mean() - 
                            self.phi_lambda) - r) /  self.tau_r
    
    # Field nonlinearity
    
    def phi_h(self, h):
        phi = N.tanh(h/self.phi_sigma)
        phi[phi<0] = 0
        return phi
        
    # Create new trajectory, weight matrix
    
    def new_trajectory(self):
        """Get a new trajectory instance
        """
        traj = None
        if self.traj_type == 'randwalk':
            traj = RandomWalk(dt=self.dt, T=self.T)
        elif self.traj_type == 'checker':
            traj = BipartiteRaster( dt=self.dt, 
                                    dwell=self.dwell_factor*self.tau_r)
            self.T = traj.T
        return traj
    
    def new_weights(self):
        """Get a new weight matrix
        """
        from scipy.stats import uniform
        from numpy.random import permutation
        W = N.empty((self.N_CA, self.N_EC), 'd')
        Wdist = uniform.rvs(size=self.N_EC, loc=0, scale=2*self.mu_W)
        Wdist[int(self.C_W*self.N_EC):] = 0
        for Wi in W:
            Wi[:] = permutation(Wdist)
        return W

    # Property getters for current position
    
    def _get_x(self): 
        return self.traj.x
    
    def _get_y(self): 
        return self.traj.y
    
    # Trajectory and rate vectors are available by default
    
    def _traj_default(self):
        return self.new_trajectory()
    
    def _W_default(self):
        return self.new_weights()
    
    def _r_default(self):
        return N.zeros(self.N_CA, 'd')
    
    def _i_aff_default(self):
        return N.zeros(self.N_CA, 'd')
    
    def _r_EC_default(self):
        return N.zeros(self.N_EC, 'd')
    
    # Normalizing factor for feedforward input current
        
    def _beta_default(self):
        return 1 / float(self.C_W * self.N_EC * self.mu_W)
    

class PlaceNetworkRaster(PlaceNetwork):
    
    """
    PlaceNetwork variant optimized for raster-based trajectories
    
    It provides a modified timestep kernel and turns off trajectory and rate 
    data tracking. This method uses significantly less memory and should serve
    as the superclass for PlaceNetwork subclasses primarily intended for use with
    raster trajectories.
    
    Setting dwell_factor determines pixel dwell-time (xTau).
    """
    
    # Redefine x, y, r traits to disable data tracking
    x = Float(track=False)
    y = Float(track=False)
    r = Array(track=False)

    # Default to BipartiteRaster trajectory
    traj_type = 'checker'
    
    # Matrix structure to hold results of raster scan
    scan = Array
    scan_ix = Int
    trial_result = 'scan'
    zero_start = false
    
    # Counter for dwell-time clamp
    dwell = Int
    dwell_count = Int
    
    def trial_setup(self):
        PlaceNetwork.trial_setup(self)
        self.scan = self.new_scan()
        self.dwell = int(self.traj._init_factor * self.dwell_count)
        self.scan_ix = 0
        self.set_input()
        
    def set_input(self):
        """Set cortical input according to current scan_ix
        """
        x, y = self.scan[self.scan_ix, :2]
        self.r_EC = self.get_afferent_input(x, y)
        self.i_aff = N.dot(self.W, self.r_EC)
        
    def run_timestep(self):
        """Non-tracking simulation timestep kernel
        """
        if self.dwell:
            self.dwell -= 1
        else:
            # Store network state
            self.scan[self.scan_ix, 2:] = self.r
            
            # Advance the scan and set inputs
            self.scan_ix += 1
            self.set_input()
            
            # Reset dwell-time counter
            self.dwell = self.dwell_count
            
            # Reset rate vector
            if self.zero_start:
                self.r[:] = 0.0
        self.evolve()
    
    def new_trajectory(self):
        if self.traj_type is 'randwalk':
            self.out('Forcing bipartite raster trajectory')
            self.traj_type = 'checker'
        return PlaceNetwork.new_trajectory(self)

    def new_scan(self):
        points = self.traj.get_points().T
        return N.c_[points, N.zeros((points.shape[0], self.N_CA), 'd')]
    
    def post_mortem(self, trial=1):
        """Modified post_mortem to properly return single-trial scan data
        """
        assert self.done, 'model simulation must be finished running'
        assert trial <= self.num_trials, 'invalid trial number specified'
        scan = self.results[trial-1]
        return TSPostMortem(x=scan[:,0], y=scan[:,1], r=scan[:,2:], ntrials=1)

    def _run_trial(self):
        """Overload AbstractModel._run_trial to use a scan-based while condition
        """
        while self.scan_ix < self.scan.shape[0] - 1:
            self.run_timestep()
            self._handle_monitor_msg()
            self.ts.advance()

    def _dwell_count_default(self):
        return int(N.round(self.dwell_factor * (self.tau_r / self.dt)))


class PlaceNetworkStd(PlaceNetworkRaster):
    
    """
    Resets default dynamic parameter values to the maximum fit achieved during 
    the PlaceNetworkSearch genetic search.
    """
    
    J0 = 45.0
    N_CA = 500
    C_W = 0.33
    phi_lambda = 0.04
    phi_sigma = 0.02
