# encoding: utf-8
"""
placemap.py -- Constructing spatial ratemaps and basic placemap analysis

Created by Joe Monaco on 02-07-2008.
Copyright (c) 2007, 2008 Columbia University. All rights reserved.
"""

# Library imports
import numpy as N, scipy as S, os

# Package imports
from .stage import StagingMap
from .core.model import AbstractModel
from .tools.bash import CPrint
from .tools.interp import GSmoothInterp2D
from .tools.images import tile2D, array_to_image
from .tools.misc import Null
from .tools.path import unique_path

# Traits imports
from enthought.traits.api import Trait, List, Array, Property, Instance, \
    Float, Int, false


# Place field determinative criteria

NOISE_FLOOR = 0.20              # Noise floor relative to map peak
FIELD_CUTOFF = 0.20             # Unit-based cutoff for defining fields relative 
                                # to peak estimated activity rate 
                                # --> MullerKubie89: 20% of peak activity
MIN_FIELD_SIZE = 50             # Putative fields with areas smaller than this 
                                # area are ignored (cm^2).
                                # --> MullerKubie89: 200 cm^2 contiguous area


class AbstractPlaceMap(StagingMap):

    """
    Superclass with general functionality for ratemap classes
    
    Required constructor argument: the AbstractModel subclass object that has the
    trajectory ('x', 'y') and rate ('r') timeseries OR a ModelPostMortem
    object containing the data ('x', 'y', 'r').
    
    To get maps for network model simulations containing multiple trials, use 
    the classmethod trial_split() to get a python list object of ratemap 
    instances for each trial.
    """

    # Network output data
    data = Instance(object)
    
    # Per-map normalized ratemap array
    Norm = Array

    def __init__(self, data, **traits):
        StagingMap.__init__(self, **traits)
        if isinstance(data, AbstractModel):
            data = data.post_mortem()
        self._validate_data(data)
        if data.ntrials > 1:
            self.out('Use trial_split classmethod for multiple trials!',
                error=True)
            return
        self.data = data
    
    @classmethod
    def _validate_data(cls, data):
        """
        Return a valid post mortem object of the constructor data argument
        """
        if not (hasattr(data, 'r') and hasattr(data, 'x') and
                hasattr(data, 'y')):
            raise ValueError, 'Data must have x, y, and r attributes'

    # Public methods
    
    def initialize(self):
        """Delete raw data after map initialization"""
        super(AbstractPlaceMap, self).initialize()
        self.data = None
    
    def initialize_norm(self):
        """Compute the normalized maps after creating the raw ratemaps"""
        if not self._initialized:
            self.initialize()
        
        self.Norm = N.empty((self.num_maps, self.H, self.W), 'd')
        for m in xrange(self.num_maps): 
            if self.Map[m].max():
                self.Norm[m] = self.Map[m] / self.Map[m].max()
            else:
                self.Norm[m] = 0
        self.out('Normalized ratemaps computed!')
    
    def save_map_image(self, desc=None, imgdir='spatial_maps', which='Map'):
        """
        Save PNG images of tiled maps stored in this ratemap object
        
        Parameters:
        desc -- short description of maps (default to model description)
        subdir -- save images in this subdirectory (default 'spatial_maps')
        which -- string name of attribute map array to tile and save (any
            (num_maps, H, W) shaped array) (default 'Map')
        """
        if not self._initialized:
            self.out('Ratemaps are not initialized!', error=True)
            return
            
        # Validate attribute map object specified by which keyword
        themap = getattr(self, which, None)
        if themap is None or type(themap) is not N.ndarray:
            self.out('Nonexistent or non-array attribute \'%s\''%which, 
                error=True)
            return
        if themap.shape[0] != self.num_maps or themap.ndim != 3:
            self.out('Map object %s does not have shape %dxHxW'%(which, 
                self.num_maps), error=True)
            return
        
        # Set map description
        if desc is None: 
            desc = self.desc
        
        # Unique image file path
        img_stem = '_'.join(desc.split() + [which.lower()])
        img_path = unique_path(os.path.join(imgdir, img_stem), 
            ext='png', fmt='%s_%02d')
        
        # Save the tiled maps to disk 
        array_to_image(tile2D(themap, gridvalue=1.03*themap.max()), img_path)
        self.out('Saved %s as a tiled image:\n-->%s'%(which.lower(), img_path), 
            prefix='Maps saved')
    
    @classmethod
    def trial_split(cls, data):
        """Get a list of ratemap objects for a multi-trial simulation"""
        if isinstance(data, AbstractModel):
            data = data.post_mortem()
        cls._validate_data(data)
        
        maps_list = []
        for t in xrange(1, data.ntrials+1):
            maps_list.append(cls(data.get_trial_data(t), desc='trial %d'%t))
        
        return maps_list
    
    # Traits defaults
    
    def _num_maps_default(self):
        return self.data.r.shape[1]


class PlaceMap(AbstractPlaceMap):
    
    """
    Baseclass with functionality for dealing with ratemaps of place fields
    
    Public methods:
    find_maxima -- save the location and magnitude of the ratemap maxima 
    find_peaks -- locate all peaks within each ratemap
    store_fields -- create coverage maps for all individual place fields
    get_unit_data -- return a records array of per-unit data*
    get_field_data -- return a recrods array of per-field data*
    
    * Records data can be accessed as:
    rec_data[i] -- entire record for ith entry
    rec_data['col_name'] -- entire data column
    rec_data[i]['col_name'] -- value of column col_name for ith entry
    
    Data attributes:
    maxima -- [find_maxima()] N_CA x 3 array of position (x, y) and rate of 
        each map's maxima; rows are 0's for silent units
    peaks -- [find_peaks()] object-array of arrays of peak postions (x, y) and 
        rates found in each ratemap; element is None for silent units
    fields -- [store_fields()] object-array of N_fields[m] x H x W arrays of 
        boolean coverage maps for each field in the mth map; element is None 
        for silent units
    
    Stage representation data:
    (These attributes are set by calling compute_coverage())
    coverage_maps -- maps of total field coverage for each place cell
    stage_coverage_map -- a stage map showing coverage across cells
    stage_coverage -- fractional coverage across cells
    stage_repr_map -- a stage map showing relative representation
        of each pixel in the stage (0.0 means no cells have overlapping
        fields, 1.0 means all cells have overlapping fields)
    stage_repr -- average number of fields representing any stage-pixel
    num_active -- the number of active place cells in this map
    sparsity -- the proportion of silent cells in this map
    """
    
    # Primary data sources
    maxima = Array
    peaks = Array
    fields = Array
    
    # Per-PC data
    coverage_maps = Array
    single_maps = Array
    
    # Across-PC data
    stage_coverage_map = Array
    stage_coverage = Float
    stage_repr_map = Array
    stage_repr = Float
    peak_rate = Float
    num_active = Int
    sparsity = Float
    
    # Private attributes
    _surround = Array
    _peaks_found = false
    _fields_stored = false
    _coverage_computed = false
    
    # Pubic methods for computing data sources
    
    def reset(self, reinitialize=False):
        """Reset state of place map computations
        """
        if reinitialize:
            self._initialized = False
        self._peaks_found = False
        self._fields_stored = False
        self._coverage_computed = False
        self.maxima = self._maxima_default()
        self.peaks = self._peaks_default()
        self.fields = self._fields_default()
        self.coverage_maps = self._coverage_maps_default()
        self.single_maps = self._single_maps_default()
        self.stage_coverage_map = self._stage_coverage_map_default()
        self.stage_repr_map = self._stage_repr_map_default()
    
    def find_maxima(self):
        """Find maxima of each ratemap, storing locations and magnitudes
        
        The following attribute are set:
        maxima -- array data is set such that maxima[m] = [x_m, y_m, z_m]
        num_active -- number of active place units
        sparsity -- proportion of inactive units in the network
        peak_rate -- maximum firing rate across all maps
        """
        if not self._initialized:
            self.initialize()
        self.out('Determining activity maxima...')
        
        x, y = N.mgrid[0:self.W, 0:self.H]
        xf = x.flatten() + 0.5
        yf = y.flatten() + 0.5
        
        # Scan the maps and store the maxima
        for m, M in enumerate(self.Map):
            zf = N.flipud(M).T.flatten()
            maxix = zf.argmax()
            if zf[maxix]:
                self.maxima[m] = xf[maxix], yf[maxix], zf[maxix]
        self.peak_rate = self.maxima[:,2].max()
    
    def find_peaks(self):
        """Find all field peaks in each ratemap
        
        Location and magnitudes are stored as 3-column arrays in *peaks* list:
            [x, y, z]
        """
        
        # Locate maxima so that cutoffs are computed
        self.find_maxima()
        self.out('Scanning maps for peaks...')
        
        # For each map, scan each pixel (x, y)
        for m, M in enumerate(self.Map):

            # Progress bar and skip empty maps
            self.out.printf('.', color='purple')
            if M.sum() == 0.0:
                continue
            
            # Define minimum rate for determining a new field peak as maximum
            # of per-cell peak or population noise floor
            cutoff = max(FIELD_CUTOFF*self.maxima[m,2], 
                NOISE_FLOOR*self.peak_rate)
                
            # Raster scan to find peaks to add to list
            peak_list = []
            for x in self._xrange:
                for y in self._yrange:
                    
                    # Get z value at the corresponding index
                    i, j = self.index(x, y)
                    z = self.Map[m, i, j]
                    if z < cutoff:
                        continue
                        
                    # Add x, y iff z is greater than all surrounding pixels
                    is_peak = True
                    for dx in self._surround:
                        try:
                            di, dj = self.index(x+dx[0], y+dx[1])
                        except IndexError:
                            continue
                        else:
                            if self.Map[m, di, dj] > z:
                                is_peak = False
                                break
                    if is_peak:
                        peak_list.append([x, y, z])
                        
            self.peaks[m] = N.array(peak_list)
            
        self._peaks_found = True
        self.out.printf('\n')
        self.out('Done!')
    
    def store_fields(self):
        """
        Store boolean arrays representing individual firing fields in each map
        
        The number of fields will be stored in *num_fields*.
        """        
        if not self._peaks_found:
            self.find_peaks()

        self.out('Scanning maps to store individual fields...')
        
        # Scan each rough-cut map using stored peaks
        for m, M in enumerate(self.Map):
            
            # Get the peaks and rough-cut for this map
            peaks = self.peaks[m]
            if peaks is None:
                self.out.printf('.', color='lightred')
                continue
                
            # Master cut of cell activity based on field cutoff
            map_cut = M > FIELD_CUTOFF*self.maxima[m, 2]
            
            # Handle nonspecific responses
            if map_cut.sum() / float(self.H*self.W) > .4:
                field_list = [map_cut]
            else:
                field_list = []
            
                # Scan peaks for unique fields
                for p, peak in enumerate(peaks):
                    field = N.zeros((self.H, self.W), '?')
                    self._mark_field(peak[0], peak[1], map_cut, field)
                
                    # Enforce field size minimum and kill dupes
                    if field.sum() > MIN_FIELD_SIZE:
                        duplicate = False
                        for f in field_list:
                            if (f*field).sum():
                                duplicate = True
                                break
                        if not duplicate:
                            field_list.append(field)
            
            # If valid fields were found, store them
            if len(field_list):
                self.fields[m] = N.array(field_list)
                if len(field_list) == 1:
                    self.out.printf('.', color='yellow')
                elif len(field_list) == 2:
                    self.out.printf('.', color='lightblue')
                else:
                    self.out.printf('.', color='purple')
            else:
                self.out.printf('.', color='lightred')

        
        self._fields_stored = True
        self.out.printf('\n')
        self.out('Done!')
            
    # Field marking recursion
    
    def _mark_field(self, x, y, master, field):
        """
        Recursively mark off complete field coverage
        
        Required parameters:
        x, y -- location within the field to mark
        master -- array (HxW) of all valid field pixels
        field -- array (HxW) containing the field picked by (x, y) in master
        """
        # Stop if out of bounds
        try:
            i, j = self.index(x, y)
        except IndexError:
            return
        
        # Stop if already marked
        if field[i, j]:
            return
        
        # Set field pixel, otherwise stop
        if master[i, j]:
            field[i, j] = True
        else:
            return
        
        # Probe surrounding pixels
        for dx in self._surround:
            try:
                self._mark_field(x+dx[0], y+dx[1], master, field)
            except RuntimeError: # max recursion depth
                break
    
    # Field coverage and stage representation methods
    
    def compute_coverage(self):
        """
        Compute coverage maps and single-field ratemaps for each place cell
        
        This method computes the following attributes:
        coverage_maps, single_maps, stage_coverage_map, stage_coverage, 
        stage_repr_map, stage_repr, num_active and sparsity
        """
        if not self._fields_stored:
            self.store_fields()
        
        self.out('Collapsing fields for coverage maps...')
        
        # Sum each set of individual field maps and store data
        npixels = self.H * self.W
        self.num_active = 0
        for m, fields in enumerate(self.fields):
            if fields is None:
                continue
            self.coverage_maps[m] = (fields.sum(axis=0) != 0)
            if len(fields) == 1:
                self.single_maps[m] = fields[0] * self.Map[m]
            else:
                for f in fields:
                    tmp = f * self.Map[m]
                    if tmp.max() == self.maxima[m, 2]:
                        self.single_maps[m] = tmp
                        break
            self.num_active += 1
        self.sparsity = 1 - (self.num_active / float(self.num_maps))
        
        # Compute full coverage of stage
        self.stage_coverage_map[:] = (self.coverage_maps.sum(axis=0) != 0)
        self.stage_coverage = float(self.stage_coverage_map.sum()) / npixels
        
        # Compute relative representation of stage
        self.stage_repr_map[:] = \
            self.coverage_maps.sum(axis=0) / float(self.num_maps)
        self.stage_repr = self.stage_repr_map.mean() * self.num_maps
        
        self._coverage_computed = True
        self.out('Done!')
    
    # Create record arrays of per-unit and per-field data 
    
    def get_unit_data(self):
        """
        Return a recarray object with data records for each place unit
        
        Record columns (i.e., data fields):
        unit -- unique integer id for each place unit
        max_x, max_y -- position where rate maximum occurs each map
        max_r -- value of the peak rate for each map
        num_fields -- number of fields
        avg_area -- average area of all fields in this map
        avg_diameter -- average diameter of all fields in this map
        coverage -- fractional coverage by area of all fields
        """
        if not self._coverage_computed:
            self.compute_coverage()
        
        # Create some useful data fields
        num_fields = N.array([(f is not None) and f.shape[0] or 0
            for f in self.fields], 'h')
        active = num_fields.astype(bool)
        num_fields = num_fields[active]
        unit = N.arange(self.num_maps)[active]
        max_x, max_y = self.maxima[active, 0], self.maxima[active, 1]
        max_rate = self.maxima[active, 2]
        avg_area = N.array([N.array([fld.sum() for fld in f]).mean() 
            for f in self.fields[active]])
        avg_diameter = N.array([N.array([2*N.sqrt(fld.sum()/N.pi) 
            for fld in f]).mean()
            for f in self.fields[active]])
        coverage = self.coverage_maps[active].sum(axis=-1).sum(axis=-1) / \
            float(self.H*self.W)
        
        # Create records array
        unit_data = N.rec.fromarrays(
            [unit, max_x, max_y, max_rate, num_fields, avg_area, avg_diameter, 
                coverage], 
            names='unit, max_x, max_y, max_r, num_fields, avg_area, '
                'avg_diameter, coverage',
            formats='l, d, d, d, h, d, d, d')
        
        return unit_data
    
    def get_field_data(self):
        """
        Return a recarray object with data records for each place field
        
        Record columns (i.e., data fields):
        id -- unique integer id for each place field
        unit -- place unit whose response this field is a part
        area -- place field area in sq-cm
        diameter -- approximate field diameter in cm
        radius -- approximate field radius in cm
        peak -- peak activity rate within each field
        average -- average activity rate across each field
        x, y -- position of the rate-weighted centroid 
        """
        if not self._coverage_computed:
            self.compute_coverage()
        
        # Find the total number of fields
        num_fields = 0
        for m in xrange(self.num_maps):
            if self.fields[m] is not None:
                num_fields += self.fields[m].shape[0]
        
        # Initialize data fields
        field_id = N.arange(num_fields)
        unit_id = N.empty(num_fields, 'h')
        area = N.empty(num_fields, 'd')
        diameter = N.empty(num_fields, 'd')
        radius = N.empty(num_fields, 'd')
        maximum = N.empty(num_fields, 'd')
        average = N.empty(num_fields, 'd')
        center_x = N.empty(num_fields, 'd')
        center_y = N.empty(num_fields, 'd')
        
        # Quantify place field characteristics
        f_id = 0
        for m in xrange(self.num_maps):            
            if self.fields[m] is None:
                continue
                
            for field in self.fields[m]:
                
                # Single field-masked ratemap and sum
                rates = field * self.Map[m]
                rates_sum = float(rates.sum())
                
                # Place unit identification
                unit_id[f_id] = m
                
                # Coverage geometry
                area[f_id] = field.sum()
                diameter[f_id] = 2*N.sqrt(area[f_id]/N.pi)
                radius[f_id] = diameter[f_id] / 2
                
                # Rate-dependent quantities
                maximum[f_id] = rates.max()
                average[f_id] = rates_sum / area[f_id]
                center_x[f_id] = (self._xrange[N.newaxis,:] * rates).sum() \
                    / rates_sum
                center_y[f_id] = (self._yrange[:,N.newaxis] * rates).sum() \
                    / rates_sum
                
                f_id += 1
            
        # Create records array
        field_data = N.rec.fromarrays(
            [field_id, unit_id, area, diameter, radius, maximum, average,
                    center_x, center_y],
            names='id, unit, area, diameter, radius, peak, average, x, y',
            formats='l, l, l, d, d, d, d, d, d')
        
        return field_data
    
    # Traits properties and defaults
    
    def _maxima_default(self):
        return N.zeros((self.num_maps, 3), 'd')
        
    def _peaks_default(self):
        return [None] * self.num_maps
    
    def _fields_default(self):
        return [None] * self.num_maps
    
    def __surround_default(self):
        """Relative indices for surrounding pixels
        """
        full = N.array([z.flatten() for z in N.mgrid[-1:2, -1:2]]).T
        return N.r_[full[:4], full[-4:]]
    
    def _coverage_maps_default(self):
        return N.zeros((self.num_maps, self.H, self.W), '?')
    
    def _single_maps_default(self):
        return N.zeros((self.num_maps, self.H, self.W), 'd')
    
    def _stage_coverage_map_default(self):
        return N.empty((self.H, self.W), '?')
    
    def _stage_repr_map_default(self):
        return N.empty((self.H, self.W), 'd')
