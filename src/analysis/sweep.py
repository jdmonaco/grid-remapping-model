#encoding: utf-8
"""
grid.analysis.sweep -- AbstractAnalysis subclass for exploring spatial map properties 
    by sweeping a 2D region of parameter space with random sampling.

Written by Joe Monaco, 04/30/2008.
Copyright (c) 2008 Columbia University. All rights reserved.
"""

# Library imports
from IPython.kernel import client as IPclient
import numpy as N, scipy as S, os

# Package imports
from .. import PlaceNetworkStd, CheckeredRatemap, GridCollection
from ..core.analysis import AbstractAnalysis
from ..tools.interp import BilinearInterp2D
from ..tools.cmap_ui import ColormapControl
from ..tools.string import snake2title

# Enthought imports
from enthought.traits.api import Enum, Bool, Button
from enthought.traits.ui.api import View, Group, Item, Include
from enthought.chaco.api import ArrayPlotData, HPlotContainer, Plot
from enthought.enable.component_editor import ComponentEditor


def run_sample_point(**kwargs):
    gc.collect()
    
    # Handle file save
    do_save = False
    if 'save_file' in kwargs:
        do_save = True
        save_file = kwargs['save_file']
        del kwargs['save_file']
    
    # Check for pre-existing data to load
    if do_save and os.path.exists(save_file):
        self.out('Loading found data:\n%s'%save_file)
        pmap = CheckeredRatemap.fromfile(save_file)
    else:
        # Run the simulation and save the results
        model = PlaceNetworkStd(W=W, EC=EC, **kwargs)
        model.advance()
        pmap = CheckeredRatemap(model)
        pmap.compute_coverage()
        if do_save:
            pmap.tofile(save_file)

    # Get field and unit data record arrays
    fdata = pmap.get_field_data()
    udata = pmap.get_unit_data()
    
    # Collate the place map sample data
    sample = {}
    sample['sparsity'] = pmap.sparsity
    sample['stage_coverage'] = pmap.stage_coverage
    sample['stage_repr'] = pmap.stage_repr
    sample['peak_rate'] = pmap.peak_rate
    
    # Collate the per-unit data
    if udata.shape[0] != 0:
        sample['max_rate'] = udata['max_r'].mean()
        sample['num_fields'] = udata['num_fields'].mean()
        sample['coverage'] = udata['coverage'].mean()
    else:
        sample['max_rate'] = sample['num_fields'] = \
            sample['coverage'] = 0.0
    
    # Collate the per-field data
    if fdata.shape[0] != 0:
        sample['area'] = fdata['area'].mean()
        sample['diameter'] = fdata['diameter'].mean()
        sample['peak'] = fdata['peak'].mean()
        sample['average'] = fdata['average'].mean()
    else:
        sample['area'] = sample['diameter'] = sample['peak'] = \
            sample['average'] = 0.0
        
    return sample
        

class SingleNetworkSweep(AbstractAnalysis, ColormapControl):
    
    """
    Analyze a 2D random sample of single-trial network simulations across
    parameter space. 
    
    See core.analysis.AbstractAnalysis documentation and collect_data method
    signature and docstring for usage.
    """
    
    label = 'Single Sweep'
    save_current_plot = Button
    show_sample_points = Bool(True)
    
    #
    # These traits must be kept up-to-date with the data made available in
    # the field and unit record arrays of PlaceMap:
    #
    # display_data -- the actual data to display in the figure plot
    # map_data -- the subset of per-map data
    # unit_data -- the subset of unit-averaged data
    # field_data -- the subset of field-averaged data
    #
    display_data = Enum('sparsity', 'stage_coverage', 'stage_repr', 
        'peak_rate', 'max_rate', 'num_fields', 'coverage', 'area', 
        'diameter', 'peak', 'average')
    map_data = Enum('sparsity', 'stage_coverage', 'stage_repr', 
        'peak_rate', 'none')
    unit_data = Enum('max_rate', 'num_fields', 'coverage', 'none')
    field_data = Enum('area', 'diameter', 'peak', 'average', 'none')
    
    traits_view = \
        View(
            Group(
                Item('figure', label='Data Map', height=450, 
                    editor=ComponentEditor()),
                Group(
                    Group(
                        Item('map_data', style='custom'),
                        Item('unit_data', style='custom'),
                        Item('field_data', style='custom'),
                        label='Data to Display',
                        show_border=True),
                    Group(
                        Include('colormap_group'),
                        Group(
                            Item('show_sample_points'),
                            label='Samples',
                            show_border=True),
                        Item('save_current_plot', show_label=False),
                        show_border=False),
                    show_border=False,
                    orientation='horizontal'),
                layout='split',
                orientation='vertical',
                show_border=False),
            title='Single Network Sweep',
            kind='live',
            resizable=True,
            width=0.6,
            height=0.8,
            buttons=['Cancel', 'OK'])
    
    def collect_data(self, x_density=10, x_bounds=(0.5,8), x_param='J0', 
        y_density=10, y_bounds=(0,2.5), y_param='phi_lambda', save_maps=True, 
        **kwargs): 
        """Store placemap data from a grid-sampled 2D region of parameter space
        
        The same network and inputs are used for the simulation at each point.
        
        Keyword arguments:
        nsamples -- the number of random samples to collect
        x_param -- string name of PlaceNetwork parameter to sweep along the x-axis
        y_param -- ibid for y-axis
        x_bounds -- bounds on sampling the parameter specified by x_param
        y_bounds -- ibid for y_param
        """
        # Store bounds and sweep parameters
        self.results['x_bounds'] = N.array(x_bounds)
        self.results['y_bounds'] = N.array(y_bounds)
        self.results['x_param']  = x_param
        self.results['y_param']  = y_param
        
        # Get ipcontroller clients
        mec = self.get_multiengine_client()
        tc = self.get_task_client()
                
        # Setup namespace on ipengine instances
        self.out('Setting up ipengines for task-farming...')
        mec.clear_queue()
        mec.reset()
        mec.execute('import gc, os')
        mec.execute('from grid_remap.place_network import PlaceNetworkStd')
        mec.execute('from grid_remap.dmec import GridCollection')
        mec.execute('from grid_remap.ratemap import CheckeredRatemap')
        
        # Set default model parameters
        pdict = dict(   growl=False, 
                        refresh_weights=False, 
                        refresh_orientation=False,
                        refresh_phase=False
                        )
        pdict.update(kwargs)
        
        # Update with keyword arguments
        all_params = PlaceNetworkStd().traits(user=True).keys()
        if x_param not in all_params:
            raise KeyError, 'x_param (%s) is not a PlaceNetwork parameter'%x_param
        if y_param not in all_params:
            raise KeyError, 'y_param (%s) is not a PlaceNetwork parameter'%y_param
        
        # Send some network weights and a grid cell object
        self.out('Pushing network weights and grid configuration...')
        EC = GridCollection()
        mec.push(dict(W=PlaceNetworkStd(EC=EC, **pdict).W, 
            spacing=EC.spacing, phi=EC._phi, psi=EC._psi))
        self.out('...and reconstructing grid collection...')
        mec.execute(
            'EC = GridCollection(spacing=spacing, _phi=phi, _psi=psi)')

        # Build the sample grid according to specifications
        pts_x = N.linspace(x_bounds[0], x_bounds[1], x_density)
        pts_y = N.linspace(y_bounds[0], y_bounds[1], y_density)
        x_grid, y_grid = N.meshgrid(pts_x, pts_y)
        pts = N.c_[x_grid.flatten(), y_grid.flatten()]
        self.results['samples'] = pts
        
        def interpolate_data(z, density=256):
            """Interpolate value z across sample points with *density* points
            """
            M = N.empty((density, density), 'd')
            x_range = N.linspace(x_bounds[0], x_bounds[1], num=density)
            y_range = N.linspace(y_bounds[1], y_bounds[0], num=density)
            
            f = BilinearInterp2D(x=pts_x, y=pts_y, z=z)
            
            for j, x in enumerate(x_range):
                for i, y in enumerate(y_range):
                    M[i,j] = f(x, y)
                    
            return M
        
        # Execute data collection process for each sample point
        self.out('Initiating task farming...')
        save_dir = os.path.join(self.datadir, 'data')
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        tasks = []
        for i, p in enumerate(pts):
            self.out('Point %d: %s = %.4f, %s = %.4f'%(i, x_param, p[0], 
                y_param, p[1]))
            pdict[x_param], pdict[y_param] = p
            if save_maps:
                pdict['save_file'] = \
                    os.path.join(save_dir, 'point_%03d.tar.gz'%i)
            tasks.append(
                tc.run(
                    IPclient.MapTask(run_sample_point, kwargs=pdict)))
        tc.barrier(tasks)
        
        # Collate sample data returned from map task
        samples = [tc.get_task_result(t_id) for t_id in tasks]
        
        # Populate result arrays for interpolation
        sparsity = N.array([pt['sparsity'] for pt in samples])
        stage_coverage = N.array([pt['stage_coverage'] for pt in samples])
        stage_repr = N.array([pt['stage_repr'] for pt in samples])
        peak_rate = N.array([pt['peak_rate'] for pt in samples])
        max_rate = N.array([pt['max_rate'] for pt in samples])
        num_fields = N.array([pt['num_fields'] for pt in samples])
        coverage = N.array([pt['coverage'] for pt in samples])
        area = N.array([pt['area'] for pt in samples])
        diameter = N.array([pt['diameter'] for pt in samples])
        peak = N.array([pt['peak'] for pt in samples])
        average = N.array([pt['average'] for pt in samples])

        # Create interpolated maps for the collated data
        def dot(): 
            self.out.printf('.', color='purple')
        self.out('Creating interpolated parameter maps for collected data'); dot()
        self.results['sparsity'] = interpolate_data(sparsity); dot()
        self.results['stage_coverage'] = interpolate_data(stage_coverage); dot()
        self.results['stage_repr'] = interpolate_data(stage_repr); dot()
        self.results['peak_rate'] = interpolate_data(peak_rate); dot()
        self.results['max_rate'] = interpolate_data(max_rate); dot()
        self.results['num_fields'] = interpolate_data(num_fields); dot()
        self.results['coverage'] = interpolate_data(coverage); dot()
        self.results['area'] = interpolate_data(area); dot()
        self.results['diameter'] = interpolate_data(diameter); dot()
        self.results['peak'] = interpolate_data(peak); dot()
        self.results['average'] = interpolate_data(average); dot()
        self.out.printf('\n')
        
        # Good-bye!
        self.out('All done!')
    
    def create_plots(self):
        """Create a simple 2D image plot of the parameter sweep
        """
        
        # Figure is horizontal container for main plot + colorbar
        self.figure = \
            container = HPlotContainer(fill_padding=True, padding=25, 
                bgcolor='linen')
        
        # Convert old data sets to the new generalized style
        if 'J0_bounds' in self.results:
            self.results['x_bounds'] = self.results['J0_bounds']
            self.results['x_param'] = 'J0'
        if 'lambda_bounds' in self.results:
            self.results['y_bounds'] = self.results['lambda_bounds']
            self.results['y_param'] = 'phi_lambda'
        
        # Data and bounds for main plot
        raw_data = self.results[self.display_data]
        data = ArrayPlotData(image=self.get_rgba_data(raw_data), raw=raw_data, 
            x=self.results['samples'][:,0], y=self.results['samples'][:,1])
        x_range = tuple(self.results['x_bounds'])
        y_range = tuple(self.results['y_bounds'])
        bounds = dict(xbounds=x_range, ybounds=y_range)

        # Create main plot
        p = Plot(data)
        p.img_plot('image', name='sweep', origin='top left', **bounds)
        p.contour_plot('raw', name='contour', type='line', origin='top left', **bounds)
        p.plot(('x', 'y'), name='samples', type='scatter', marker='circle', 
            color=(0.5, 0.6, 0.7, 0.4), marker_size=2)
        
        # Tweak main plot
        p.title = snake2title(self.display_data)
        p.x_axis.orientation = 'bottom'
        p.x_axis.title = snake2title(self.results['x_param'])
        p.y_axis.title = snake2title(self.results['y_param'])
        p.plots['samples'][0].visible = self.show_sample_points
    
        # Add main plot and colorbar to figure
        container.add(p)
        container.add(
            self.get_colorbar_plot(bounds=(raw_data.min(), raw_data.max())))
        
        # Set radio buttons
        self.unit_data = self.field_data = 'none'
        
    # Traits notifications for the interactive GUI
    
    def _cmap_notify_changed(self):
        """Respond to changes to the colormap specification by updating
        """
        self._update_figure_plot()
    
    def _save_current_plot_fired(self):
        self.save_plots(fmt='png')
    
    def _show_sample_points_changed(self, new):
        self.figure.components[0].plots['samples'][0].visible = new
        self.figure.request_redraw()

    def _update_figure_plot(self):
        if self.figure is None:
            return
        
        # Update data for the main plot
        raw_data = self.results[self.display_data]
        main_plot = self.figure.components[0]
        main_plot.data.set_data('image', self.get_rgba_data(raw_data))
        main_plot.data.set_data('raw', raw_data)
        main_plot.title = snake2title(self.display_data)
        
        # Remove old colorbar and add new one
        del self.figure.components[1]
        self.figure.add(
            self.get_colorbar_plot(bounds=(raw_data.min(), raw_data.max())))
        
        self.figure.request_redraw()
    
    def _display_data_changed(self, old, new):
        if new in self.results:
            self._update_figure_plot()
        else:
            self.display_data = old
            self.out('This analysis does not contain \'%s\' data'%new, 
                error=True)
    
    def _map_data_changed(self, new):
        if new != 'none':
            self.unit_data = self.field_data = 'none'
            self.display_data = new

    def _unit_data_changed(self, new):
        if new != 'none':
            self.map_data = self.field_data = 'none'
            self.display_data = new

    def _field_data_changed(self, new):
        if new != 'none':
            self.unit_data = self.map_data = 'none'
            self.display_data = new
        
