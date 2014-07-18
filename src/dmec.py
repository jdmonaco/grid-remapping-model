# encoding: utf-8
"""
dmec.py -- Real-time (non-caching) procedural model of MEC grid cell responses 
    capable of translational and rotational realignment. 

Created by Joe Monaco on 2009-07-29. 
Completed/frozen on 2009-12-15.
Copyright (c) 2009 Johns Hopkins University. All rights reserved.
"""

# Library imports
import numpy as np
from scipy import pi, rand, sqrt, sin, cos

# Package imports
from .tools.filters import halfwave
from .tools.array_container import TraitedArrayContainer
from .tools.radians import xy_to_rad, xy_to_rad_vec, shortcut

# Traits imports 
from enthought.traits.api import Float, Int, Tuple, Array, List, false

# Constant values
GRID_SPACING_RANGE = (30.0, 90.0)
ENVIRONMENT_SIZE = (100.0, 100.0)


class GridCollection(TraitedArrayContainer):
    
    """
    Procedural model of a collection of grid cell spatial response maps
    """
    
    num_maps = Int(1000)
    spacing_bounds = Tuple(GRID_SPACING_RANGE)
    mid = Tuple((ENVIRONMENT_SIZE[0] / 2.0, ENVIRONMENT_SIZE[1] / 2.0))
    peak_rate = Float(1)
    spacing = Array
    k = Array
    ellipticity = false
    ell_mag = Array
    ell_angle = Array
    zoom = false
    zoom_scale = Array
    _ellipticity = false(desc='cache')
    _ell_mag = Array(desc='cache')
    _ell_angle = Array(desc='cache')
    _zoom = false(desc='cache')
    _zoom_scale = Array(desc='cache')
    _phi = Array
    _psi = Array
    _phi0 = Array(desc='cache')
    _psi0 = Array(desc='cache')
    _phi_radius = Array
    _thetas = List([0.0, 2*pi/3, 4*pi/3])
    _norm = Float
                                        
    def __init__(self, **traits):
        TraitedArrayContainer.__init__(self, **traits)
        self.store()
    
    def map_value(self, x, y):
        """Get population rate vector of this grid collection at position (x,y)
        """
        x, y = self.map_transforms(x, y)
        return self._norm * self.__g(
            reduce(np.add, 
                [cos(
                    (sin(t-self._psi)*(x-self._phi[:,0]-self.mid[0]) + 
                     cos(t-self._psi)*(y-self._phi[:,1]-self.mid[0]))/self.k
                    ) 
                 for t in self._thetas]))
    
    def __g(self, x):
        """Monotonic gain function for grid responses
        """
        return halfwave(np.exp(0.25*x) - 0.75)
    
    # Ellipticity and zoom (scaling) transforms 
    
    def map_transforms(self, x, y):
        if self.ellipticity:
            # Get polar coordinates from midpoint
            dx = x - self.mid[0]
            dy = y - self.mid[1]
            r = sqrt(dx**2 + dy**2)
            theta = xy_to_rad_vec(dx, dy)
            
            # Rotational coordinate transform, back to Cartesian 
            theta_prime = theta - self.ell_angle
            dx_prime = r*cos(theta_prime)
            dy_prime = r*sin(theta_prime)
            
            # Do the elliptical transform, back to polar
            dx_ell = dx_prime / (1+self.ell_mag)
            dy_ell = dy_prime * (1+self.ell_mag)
            r_ell = sqrt(dx_ell**2 + dy_ell**2)
            theta_ell = xy_to_rad_vec(dx_ell, dy_ell) + self.ell_angle
            
            # Revert to absolute Cartesian coordinate frame
            x = self.mid[0] + r_ell*cos(theta_ell)
            y = self.mid[1] + r_ell*sin(theta_ell)
        
        if self.zoom:
            # Get polar coordinates from midpoint
            dx = x - self.mid[0]
            dy = y - self.mid[1]
            
            # Compute scaled radius and center-angles
            r_zoom = sqrt(dx**2 + dy**2) / self.zoom_scale
            theta = xy_to_rad_vec(dx, dy)
            
            # Project back to absolute Cartesian coordinates
            x = self.mid[0] + r_zoom*cos(theta)
            y = self.mid[1] + r_zoom*sin(theta)

        return x, y
    
    # Traits default values
    
    def _spacing_default(self):
        return self.spacing_bounds[0] + \
            (self.spacing_bounds[1] - self.spacing_bounds[0]) * \
            rand(self.num_maps)
    
    def _k_default(self):
        return (sqrt(3)/(4*pi)) * self.spacing
    
    def _ell_mag_default(self):
        return np.zeros(self.num_maps, 'd')
    
    def _ell_angle_default(self):
        return np.zeros(self.num_maps, 'd')
    
    def _zoom_scale_default(self):
        return np.ones(self.num_maps, 'd')    
    
    def __psi_default(self):
        return self.new_orientations()
    
    def __phi_default(self):
        return self.new_spatial_phases()
    
    def __norm_default(self):
        return self.peak_rate / self.__g(3)
    
    def __phi_radius_default(self):
        return (self.spacing/2) / cos(pi/6)

    # Rotate/shift remapping methods
    
    def shift(self, shift, mask=None):
        """Shift the grids
        
        The phase shift value can be a 2-element array to be applied to all
        grid phases (subject to the binary/index *mask* array) or a *phi*-shaped 
        array specifying per-grid phase shifts.
        
        The phases are wrapped on the half-spacing circle.
        """
        # Add the delta shift value to grid phases
        shift = np.squeeze(np.array(shift))
        try:
            if mask is not None:
                self._phi[mask] += shift
            else:
                self._phi += shift
        except ValueError:
            raise ValueError, 'mask and shift arrays must match'
            
        # Wrap the phase values on the half-spacing circle
        hex_angles = np.arange(0, 2*pi, pi/3)
        for i in xrange(self.num_maps):
            vertices = hex_angles + self._psi[i]
            while sqrt((self._phi[i]**2).sum()) > self._phi_radius[i]:
                orig = xy_to_rad(self._phi[i,0], self._phi[i,1]) - pi
                proj = vertices[np.argmin([shortcut(v, orig) for v in vertices])]
                self._phi[i,0] += self.spacing[i] * np.cos(proj)
                self._phi[i,1] += self.spacing[i] * np.sin(proj)
                                    
    def rotate(self, angle, mask=None):
        """Rotate the grids (arena centered)
        
        Grids to be rotated can be optionally specified by bool/index array
        *mask*, otherwise population is rotated. Specified *angle* can be a
        scalar value to be applied to the population or a population- or
        mask-sized array depending on whether *mask* is specified.
        """
        rot2D = lambda psi: [[cos(psi), sin(psi)], [-sin(psi),  cos(psi)]]
        if mask is not None and type(mask) is np.ndarray:
            if mask.dtype.kind == 'b':
                mask = mask.nonzero()[0]
            if type(angle) is np.ndarray and angle.size == mask.size:
                for i,ix in enumerate(mask):
                    self._phi[ix] = np.dot(self._phi[ix], rot2D(angle[i]))
            elif type(angle) in (int, float, np.float64):
                angle = float(angle)
                self._phi[mask] = np.dot(self._phi[mask], rot2D(angle))
            else:
                raise TypeError, 'angle must be mask-sized array or float'
            self._psi[mask] = np.fmod(self._psi[mask]+angle, 2*pi)
        elif mask is None:
            if type(angle) is np.ndarray and angle.size == self.num_maps:
                for i in xrange(self.num_maps):
                    self._phi[i] = np.dot(self._phi[i], rot2D(angle[i]))
            elif type(angle) in (int, float, np.float64):
                angle = float(angle)
                self._phi = np.dot(self._phi, rot2D(angle))
            else:
                raise TypeError, 'angle must be num_maps array or float'
            self._psi = np.fmod(self._psi+angle, 2*pi)
        else:
            raise TypeError, 'mask must be bool/index array'

    # Store/reset alignment
    
    def store(self):
        """Save the current grid configuration to be restored later
        """
        self._phi0 = self._phi.copy()
        self._psi0 = self._psi.copy()
        self._ellipticity = self.ellipticity
        self._ell_mag = self.ell_mag.copy()
        self._ell_angle = self.ell_angle.copy()
        self._zoom = self.zoom
        self._zoom_scale = self.zoom_scale.copy()
        
    def reset(self):
        """Reset the grid configuration to the stored configuration
        """
        self._phi[:] = self._phi0
        self._psi[:] = self._psi0
        self.ellipticity = self._ellipticity
        self.ell_mag[:] = self._ell_mag
        self.ell_angle[:] = self._ell_angle
        self.zoom = self._zoom
        self.zoom_scale[:] = self._zoom_scale

    # Convenience methoda

    def randomize_phase(self):
        """Randomize grid spatial phases noncoherently
        """
        self._phi = self.new_spatial_phases()
    
    def randomize_orientation(self):
        """Set grid orientations coherently to a random value
        """
        self._psi = self.new_orientations()
    
    def new_orientations(self):
        """Get a new coherent array of grid orientations
        """
        return (pi/3) * rand() + np.zeros(self.num_maps)
    
    def new_spatial_phases(self):
        """Get x,y array of random spatial phases on the half-spacing circle
        """
        p0 = 2*rand(self.num_maps, 2) - 1
        for m in xrange(self.num_maps):
            while (p0[m]**2).sum() > 1:
                p0[m] = 2*rand(2) - 1
        return p0 * self._phi_radius[:,np.newaxis]

    def get_modules(self, nmodules, freq_sort=False):
        """Get a list of index arrays for a modular partition of the grids
        
        Arguments:
        nmodules -- the number of equal-sized modular partitions
        freq_sort -- whether to partition based on spatial frequency
        """
        if freq_sort:
            grid_ix = np.argsort(self.spacing)
        else:
            grid_ix = np.arange(self.num_maps)
        return np.array_split(grid_ix, nmodules)
    
    def get_z_stack(self, size=ENVIRONMENT_SIZE):
        """Get a z-stack matrix of the population responses
        
        Convenience method to get a matrix array with the spatial responses
        of each grid-unit in this GridCollection object. Pixels get value from 
        the middle of the area represented by the pixel, and the origin is the
        lower left corner of the individual spatial maps (index (size[1]-1,0)).
        
        Keyword arguments:
        size -- (H,W)-tuple specifying the area in cm-pixels
        """
        M = np.squeeze(np.empty((self.num_maps, size[0], size[1]), 'd'))
        for i in xrange(int(size[0])):
            for j in xrange(int(size[1])):
                M[...,i,j] = self.map_value(j+0.5, size[1]-i-0.5)
        return M

    # Realignment helper functions

    @classmethod
    def get_delta_phi(cls, scale=None):
        """Generate a random spatial phase displacement
    
        Keyword arguments:
        scale -- set grid scale that determines range of possible phase shifts
        """
        if scale is None:
            scale = max(GRID_SPACING_RANGE)
        outer_bound = 0.5 * scale
        lower_bound = 0.2 * outer_bound
        
        # Generate and return random displacement
        r = (outer_bound - lower_bound) * rand() + lower_bound
        theta = 2 * pi * rand()
        return r * np.array([cos(theta), sin(theta)])
    
    @classmethod
    def get_delta_psi(cls):
        """Generate a random orientation realignment (-30 to +30 degrees)
        """
        return (pi/6) * (2 * rand() - 1)

    @classmethod
    def get_ellipticity(cls, ecc_range=(0.0, 0.2)):
        """Generate a random magnitude for the ellipticity transform
        """
        return (ecc_range[1] - ecc_range[0]) * rand() + ecc_range[0]
    
    @classmethod
    def get_elliptic_angle(cls):
        """Generate a random angle for the semimajor axis of ellipticity
        """
        return pi * (rand() - 0.5)

    @classmethod
    def get_zoom_scale(cls, zoom_range=(1.0, 1.2)):
        """Generate a random rescaling factor
        """
        return (zoom_range[1] - zoom_range[0]) * rand() + zoom_range[0]
