#encoding: utf-8
"""
grid.analysis.scan -- AbstractAnalysis subclass for exploring statistics of spatial 
    map properties by scanning a single parameter with even sampling.

Written by Joe Monaco, 05/14/2008.
Copyright (c) 2008 Columbia University. All rights reserved.
"""

# Library imports
import numpy as N, scipy as S, os

# Package imports
from .. import PlaceNetworkStd, CheckeredRatemap, GridCollection
from ..core.analysis import AbstractAnalysis
from ..tools.string import snake2title

# Traits imports
from enthought.traits.api import Enum, Button
from enthought.traits.ui.api import View, Group, Item, Include

# Chaco imports for custom analysis view
from enthought.chaco.api import (ArrayPlotData, ArrayDataSource, Plot, 
    GridContainer)
from enthought.enable.component_editor import ComponentEditor

# Tuples of available placemap data
STAGE_DATA = ('sparsity', 'stage_coverage', 'stage_repr', 'peak_rate')
UNIT_DATA = ('max_rate', 'num_fields', 'coverage')
FIELD_DATA = ('area', 'diameter', 'peak', 'average')
DATA = STAGE_DATA + UNIT_DATA + FIELD_DATA


class BaseScan(AbstractAnalysis):

    label = 'Stat Scan'
    save_current_plot = Button
    
    traits_view = \
        View(
            Item('figure', show_label=False, editor=ComponentEditor()),
            title='Network Population Scan',
            kind='live',
            resizable=True,
            width=0.75,
            height=0.75,
            buttons=['Cancel', 'OK'])
    
    def collect_data(self, ntrials=10, npoints=4, param='J0', bounds=(1, 4), 
        **kwargs):
        raise NotImplementedError 

    def create_plots(self):
        """Create a simple 2D image plot of the parameter sweep"""
        
        # Figure is horizontal container for main plot + colorbar
        self.figure = \
            container = GridContainer(fill_padding=True, spacing=(5,5),
                padding=[20, 20, 40, 10], bgcolor='linen', shape=(3,4))
        
        # Create datasource for means and confidence intervals
        data_dict = {}
        for d in DATA:
            # The means
            data_dict[d] = self.results[d]
            data_dict[d][N.isnan(data_dict[d])] = 0.0
            
            # The 95% confidence intervals
            conf_interval = 1.96*self.results[d + '_err']
            conf_interval[N.isnan(conf_interval)] = 0.0
            data_dict[d + '_err'] = \
                N.r_[self.results[d] + conf_interval,
                    (self.results[d] - conf_interval)[::-1]]
        
        data = ArrayPlotData(
            index=self.results['samples'], 
            err_ix=N.r_[self.results['samples'], self.results['samples'][::-1]], 
            **data_dict)
        
        # Create individual plots and add to grid container
        pad_factor = 0.08
        for d in DATA:
            p = Plot(data, padding=[25, 25, 25, 40])
            styles = {'line_width':1.5, 'color':'darkcyan'}
            
            # Plot the error, line, and scatter plots
            p.plot(('err_ix', d + '_err'), name=d+'_p', type='polygon',
                edge_color='transparent', face_color='silver')
            p.plot(('index', d), name=d+'_l', type='line', **styles)
            p.plot(('index', d), name=d+'_s', type='scatter', marker='circle', 
                marker_size=int(2*styles['line_width']), color=styles['color'],
                line_width=0)
            
            # Y-axis padding
            low = (self.results[d] - self.results[d + '_err']).min()
            high = (self.results[d] + self.results[d + '_err']).max()
            padding = pad_factor * (high - low)
            low -= padding
            high += padding
            if low == high:
                low -= 1
                high += 1
            p.value_range.set_bounds(low, high)
            
            # X-axis padding
            padding = pad_factor * self.results['samples'].ptp()
            p.index_range.set_bounds(
                self.results['samples'][0] - padding, 
                self.results['samples'][-1] + padding)
            
            # Labels, grids and ticks
            p.title = snake2title(d)
            p.x_axis.title = snake2title(self.results['param'])
            p.y_grid.visible = p.x_grid.visible = False
            p.x_axis.tick_in = p.y_axis.tick_in = 0
            
            # Add the plot to the grid container
            container.add(p)


class MultiNetworkScan(BaseScan):
    
    """
    Analyze a 1D scan of regularly-spaced distributions of network simulations
    across parameter space. 
    
    See core.analysis.AbstractAnalysis documentation and collect_data method
    signature and docstring for usage.
    """
    
    label = 'Network Scan'
    
    def collect_data(self, ntrials=10, npoints=5, param='J0', bounds=(0.5, 8), 
        **kwargs): 
        """
        Store statistics about placemap data from multiple trials along a 1D 
        parameter scan
        
        Keyword arguments:
        ntrials -- number of network trials to run per sample point
        npoints -- number of sample points, inclusive of the bounds
        param -- string name of PlaceNetwork parameter to scan
        bounds -- bounds for the parameter scan
        """
        
        # Store bounds and scan parameter
        self.results['bounds'] = N.array(bounds)
        self.results['param']  = param
        
        # Load cortex
        self.out('Creating grid collection object...')
        EC = GridCollection()
        os.chdir(self.datadir)

        # Set default model parameters
        pdict = dict(   EC=EC, 
                        growl=False, 
                        desc='scan', 
                        projdir=self.datadir, 
                        refresh_weights=True, 
                        refresh_orientation=False,
                        refresh_phases=False,
                        refresh_traj=False,
                        traj_type='checker',
                        num_trials=ntrials,
                        monitoring=True)
        pdict.update(kwargs)
        
        # Update with keyword arguments
        if param not in PlaceNetworkStd().traits(user=True).keys():
            raise ValueError, 'param (%s) is not a user parameter'%param
        
        # Create the list of sample points to scan
        self.out('Creating %s scan vector from %.2f to %.2f'%((param,)+bounds))
        if bounds[0] > bounds[1]:
            bounds = bounds[::-1]
        pts = N.linspace(bounds[0], bounds[1], num=npoints)
        self.results['samples'] = pts
        
        # Initialize stage map sample data arrays
        sparsity = N.empty(npoints, 'd')
        sparsity_err = N.empty(npoints, 'd')
        stage_coverage = N.empty(npoints, 'd')
        stage_coverage_err = N.empty(npoints, 'd')
        stage_repr = N.empty(npoints, 'd')
        stage_repr_err = N.empty(npoints, 'd')
        peak_rate = N.empty(npoints, 'd')
        peak_rate_err = N.empty(npoints, 'd')
        
        # Initialize per-unit sample data arrays
        max_rate = N.zeros(npoints, 'd')
        max_rate_err = N.zeros(npoints, 'd')
        num_fields = N.zeros(npoints, 'd')
        num_fields_err = N.zeros(npoints, 'd')
        coverage = N.zeros(npoints, 'd')
        coverage_err = N.zeros(npoints, 'd')

        # Initialize per-field sample data arrays
        area = N.zeros(npoints, 'd')
        area_err = N.zeros(npoints, 'd')
        diameter = N.zeros(npoints, 'd')
        diameter_err = N.zeros(npoints, 'd')
        peak = N.zeros(npoints, 'd')
        peak_err = N.zeros(npoints, 'd')
        average = N.zeros(npoints, 'd')
        average_err = N.zeros(npoints, 'd')
        
        # Error calculation
        def error(values):
            return N.std(values) / N.sqrt(len(values))
        
        # Per-sample data collection method
        def run_sample_point(i, model):
            self.out('Running (%d): %s = %.4f'%(i, param, getattr(model, param)))
            
            # Run the model simulation and save the results
            model.advance_all()

            # Create ratemap objects
            ir_list = [None] * ntrials
            fdata_list = [None] * ntrials
            udata_list = [None] * ntrials
            for trial in xrange(ntrials):
                ir = CheckeredRatemap(model.post_mortem(trial=trial+1))
                ir.compute_coverage()
                ir_list[trial] = ir
                fdata_list[trial] = ir.get_field_data()
                udata_list[trial] = ir.get_unit_data()
            
            # Collate the stage map data
            sparsity[i] = N.mean([ir.sparsity for ir in ir_list])
            sparsity_err[i] = error([ir.sparsity for ir in ir_list])
            stage_coverage[i] = N.mean([ir.stage_coverage for ir in ir_list])
            stage_coverage_err[i] = error([ir.stage_coverage for ir in ir_list])
            stage_repr[i] = N.mean([ir.stage_repr for ir in ir_list])
            stage_repr_err[i] = error([ir.stage_repr for ir in ir_list])
            peak_rate[i] = N.mean([ir.peak_rate for ir in ir_list])
            peak_rate_err[i] = error([ir.peak_rate for ir in ir_list])
            
            # Collate the per-unit data
            _max_rate = N.array([], 'd')
            _num_fields = N.array([], 'd')
            _coverage = N.array([], 'd')

            for udata in udata_list:
                _max_rate = N.r_[_max_rate, udata['max_r']]
                _num_fields = N.r_[_num_fields, udata['num_fields']]
                _coverage = N.r_[_coverage, udata['coverage']]

            max_rate[i] = _max_rate.mean()
            max_rate_err[i] = error(_max_rate)
            num_fields[i] = _num_fields.mean()
            num_fields_err[i] = error(_num_fields)
            coverage[i] = _coverage.mean()
            coverage_err[i] = error(_coverage)
            
            # Collate the per-field data
            _area = N.array([], 'd')
            _diameter = N.array([], 'd')
            _peak = N.array([], 'd')
            _average = N.array([], 'd')
            
            for fdata in fdata_list:
                _area = N.r_[_area, fdata['area']]
                _diameter = N.r_[_diameter, fdata['diameter']]
                _peak = N.r_[_peak, fdata['peak']]
                _average = N.r_[_average, fdata['average']]
            
            area[i] = _area.mean()
            area_err[i] = error(_area)
            diameter[i] = _diameter.mean()
            diameter_err[i] = error(_diameter)
            peak[i] = _peak.mean()
            peak_err[i] = error(_peak)
            average[i] = _average.mean()
            average_err[i] = error(_average)
                
        # Execute data collection process for each sample point
        self.out('Beginning data collection process')
        for i, p in enumerate(pts):
            pdict[param] = p
            self.execute(run_sample_point, i, PlaceNetworkStd(**pdict))
        
        # Store the mean data results
        self.results['sparsity'] = sparsity
        self.results['stage_coverage'] = stage_coverage
        self.results['stage_repr'] = stage_repr
        self.results['peak_rate'] = peak_rate
        self.results['max_rate'] = max_rate
        self.results['num_fields'] = num_fields
        self.results['coverage'] = coverage
        self.results['area'] = area
        self.results['diameter'] = diameter
        self.results['peak'] = peak
        self.results['average'] = average
        
        # ... and the error data
        self.results['sparsity_err'] = sparsity_err
        self.results['stage_coverage_err'] = stage_coverage_err
        self.results['stage_repr_err'] = stage_repr_err
        self.results['peak_rate_err'] = peak_rate_err
        self.results['max_rate_err'] = max_rate_err
        self.results['num_fields_err'] = num_fields_err
        self.results['coverage_err'] = coverage_err
        self.results['area_err'] = area_err
        self.results['diameter_err'] = diameter_err
        self.results['peak_err'] = peak_err
        self.results['average_err'] = average_err
        
        # Good-bye!
        self.out('All done!')
