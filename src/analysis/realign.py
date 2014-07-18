#encoding: utf-8
"""
realign.py -- A two-dimensional sweep of realignment magnitudes and 
    within-module variances to explore the parameters under which positional
    and rate remapping occur. 

Written by Joe Monaco, 08/04/2008. Updated 01/21/09.
Copyright (c) 2008 Columbia University. All rights reserved.
"""

# Library imports
from IPython.kernel import client as IPclient
import numpy as N, scipy as S, os, gc

# Package imports
from ..place_network import PlaceNetworkStd
from ..ratemap import CheckeredRatemap
from ..dmec import GridCollection
from ..tools.interp import BilinearInterp2D
from ..tools.string import snake2title
from .sweep import SingleNetworkSweep
from .compare import compare_AB

# Enthought imports
from enthought.traits.api import Enum
from enthought.chaco.api import ArrayPlotData, HPlotContainer, Plot


def run_sample_point(save_file, d_x, d_y):
    gc.collect()

    # Create the modules index arrays
    mods = nmodules
    if x_type == 'modules':
        mods = int(d_x)
    elif y_type == 'modules':
        mods = int(d_y)
    modules = EC.get_modules(mods, freq_sort=freq_modules)

    # Reset grids and activate transforms if necessary
    EC.reset()
    if 'ellipticity' in (x_type, y_type):
        EC.ellipticity = True
    if 'zoom' in (x_type, y_type):
        EC.zoom = True        

    # Modulate grid responses according to realignment parameters
    for m, m_ix in enumerate(modules):
        # Handle x-axis realignment
        if x_type == 'shift':
            EC.shift(d_x * delta_phi[m], mask=m_ix)
        elif x_type == 'rotate':
            EC.rotate(d_x * delta_psi[m], mask=m_ix)
        elif x_type == 'ellipticity':
            EC.ell_mag[m_ix] = d_x * ell_mags[m]
            EC.ell_angle[m_ix] = d_x * ell_angles[m]
        elif x_type == 'zoom':
            EC.zoom_scale[m_ix] = 1 + d_x * (zoom_scales[m] - 1)

        # Handle y-axis realignment
        if y_type == 'shift':
            EC.shift(d_y * delta_phi[m], mask=m_ix)
        elif y_type == 'rotate':
            EC.rotate(d_y * delta_psi[m], mask=m_ix)
        elif y_type == 'ellipticity':
            EC.ell_mag[m_ix] = d_y * ell_mags[m]
            EC.ell_angle[m_ix] = d_y * ell_angles[m]
        elif y_type == 'zoom':
            EC.zoom_scale[m_ix] = 1 + d_y * (zoom_scales[m] - 1)

    # Simulate and save the realigned spatial map
    model = PlaceNetworkStd(EC=EC, W=W, **pdict)
    model.advance()
    B = CheckeredRatemap(model)
    B.compute_coverage()
    B.tofile(save_file)
    return 


class RealignmentSweep(SingleNetworkSweep):
    
    """
    Analyze a 2D random sample of single-trial network simulations across
    realignment magnitudes or variances in A-B environment comparisons. 
    
    See core.analysis.AbstractAnalysis documentation and collect_data method
    signature and docstring for usage.
    """
    
    label = 'Realign Sweep'
    
    display_data = Enum('remapping', 'rate_remapping', 'turnover', 'sparsity', 
        'stage_coverage', 'stage_repr', 'peak_rate', 'max_rate', 'num_fields', 
        'coverage', 'area', 'diameter', 'peak', 'average')
    map_data = Enum('remapping', 'rate_remapping', 'turnover', 'sparsity', 
        'stage_coverage', 'stage_repr', 'peak_rate', 'none')


    def collect_data(self, x_type='shift', y_type='rotate', x_density=10, y_density=10, 
        nmodules=1, freq_modules=False, x_max=None, y_max=None, **kwargs): 
        """
        Store placemap data from a randomly sampled 2D region of parameter space
        for realignment magnitudes or variances (spatial phase vs. orientation).
        
        The same network is used for the simulation at each point, and each sample
        is compared to a reference (A) spatial map. 
        
        Keyword arguments:
        x_type -- realignment type along x axis; must be one of 'shift', 'rotate', 
            'ellipticity', 'zoom', or 'modules' (default 'shift')
        y_type -- realignment type along y axis (default 'rotate)
        x_density -- number of x_type samples along the defined x_bounds (10)
        y_density -- number of y_type samples along the defined y_bounds (10)
        nmodules -- number of independent alignment modules; used as max number
            of modules if x_type or y_type is set to 'modules'
        freq_modules -- whether modules are spatial frequency partitions
        x_max -- set upper bound for extent of x_type realignment along x axis;
             (shift should be a 2-tuple value)
        y_max -- set upper bound for extent of y_type realignment along y axis
        """
        # Parse the realignment types
        realignment_types = ('shift', 'rotate', 'ellipticity', 'zoom', 'modules')
        if x_type not in realignment_types:
            raise ValueError, 'invalid realignment type specification (x_type)'
        if y_type not in realignment_types:
            raise ValueError, 'invalid realignment type specification (y_type)'
        
        # Split cortical population into modules
        self.results['nmodules'] = nmodules = int(nmodules)
        self.results['freq_modules'] = freq_modules
        self.results['x_type'] = x_type
        self.results['y_type'] = y_type

        # Make data directory
        map_dir = os.path.join(self.datadir, 'data')
        if not os.path.exists(map_dir):
            os.makedirs(map_dir)
        
        # Set default model parameters
        pdict = dict(   refresh_weights=False,
                        refresh_phase=False,
                        refresh_orientation=False
                        )
        pdict.update(kwargs)
        
        # Simulate reference spatial map for environment A
        self.out('Simulating reference spatial map...')
        EC = GridCollection()
        model = PlaceNetworkStd(EC=EC, **pdict)
        model.advance()
        A = CheckeredRatemap(model)
        A.compute_coverage()
        A.tofile(os.path.join(map_dir, 'map_A'))
        
        # Setup namespace on ipengine instances
        self.out('Setting up ipengines for task-farming...')
        mec = self.get_multiengine_client()
        tc = self.get_task_client()
        mec.clear_queue()
        mec.reset()
        mec.execute('import gc')
        mec.execute('from grid_remap.place_network import PlaceNetworkStd')
        mec.execute('from grid_remap.dmec import GridCollection')
        mec.execute('from grid_remap.ratemap import CheckeredRatemap')
        
        # Send some network weights, grid configuration and sweep info
        self.out('Pushing network weights and grid configuration...')
        W = model.W
        mec.push(dict(  W=model.W, 
                        pdict=pdict,
                        spacing=EC.spacing, 
                        phi=EC._phi, 
                        psi=EC._psi,
                        nmodules=nmodules,
                        freq_modules=freq_modules, 
                        x_type=x_type,
                        y_type=y_type
                        ))
        mec.execute('EC = GridCollection(spacing=spacing, _phi=phi, _psi=psi)')
        
        # Set up modular realignment parameters, pushing data out to engines
        self.results['bounds'] = bounds = N.array([[0, 1]]*2, 'd')
        density = [x_density, y_density]
        r_max = (x_max, y_max)
        r_type = (x_type, y_type)
        for i in 0, 1:
            if r_type[i] == 'shift':
                if nmodules == 1 and r_max[i] is not None:
                    delta_phi = N.array([r_max[i]], 'd')
                elif nmodules > 1 and r_max[i] is not None:
                    delta_phi = N.array(r_max[i], 'd')
                else:
                    grid_scale = None
                    if freq_modules and r_type[1-i] == 'modules':
                        grid_scale = 60.0 # cf. lab notebook @ p.147
                    delta_phi = \
                        N.array([GridCollection.get_delta_phi(scale=grid_scale)
                            for m in xrange(nmodules)])
                mec.push(dict(delta_phi=delta_phi))
                self.results[r_type[i] + '_params'] = delta_phi
                self.out('Pushed shift parameters:\n%s'%str(delta_phi))
                
            elif r_type[i] == 'rotate':
                if nmodules == 1 and r_max[i] is not None:
                    delta_psi = N.array([r_max[i]], 'd')
                elif nmodules > 1 and r_max[i] is not None:
                    delta_psi = N.array(r_max[i], 'd')
                else:
                    delta_psi = N.array([GridCollection.get_delta_psi() 
                        for m in xrange(nmodules)])
                mec.push(dict(delta_psi=delta_psi))                
                self.results[r_type[i] + '_params'] = delta_psi
                self.out('Pushed rotate parameters:\n%s'%str(delta_psi))
                
            elif r_type[i] == 'ellipticity':
                if nmodules == 1 and r_max[i] is not None:
                    ell_mags = N.array([r_max[i]], 'd')
                    ell_angles = N.array([0.0])
                else:
                    ell_mags = N.array([GridCollection.get_ellipticity() 
                        for m in xrange(nmodules)])
                    ell_angles = N.array([GridCollection.get_elliptic_angle() 
                        for m in xrange(nmodules)])
                mec.push(dict(ell_mags=ell_mags, ell_angles=ell_angles))
                self.results[r_type[i] + '_params'] = \
                    N.c_[ell_mags, ell_angles]
                self.out('Pushed ellipticity parameters:\n' + 
                    'Flattening: %s\nAngles: %s'%(str(ell_mags), 
                        str(ell_angles)))

            elif r_type[i] == 'zoom':
                if nmodules == 1 and r_max[i] is not None:
                    zoom_scales = N.array([r_max[i]], 'd')
                else:
                    zoom_scales = N.array([GridCollection.get_zoom_scale() 
                        for m in xrange(nmodules)])
                mec.push(dict(zoom_scales=zoom_scales))
                self.results[r_type[i] + '_params'] = zoom_scales
                self.out('Pushed zoom parameters:\n%s'%str(zoom_scales))
                
            elif r_type[i] == 'modules':
                density[i] = nmodules
                bounds[i] = 1, nmodules
                self.out('Setting up modularity sweep for %d modules'%nmodules)
        
        # Build the sample grid according to specifications
        pts_x = N.linspace(bounds[0,0], bounds[0,1], density[0])
        pts_y = N.linspace(bounds[1,0], bounds[1,1], density[1])
        x_grid, y_grid = N.meshgrid(pts_x, pts_y)
        pts = N.c_[x_grid.flatten(), y_grid.flatten()]
        self.results['samples'] = pts

        # Initialize stage map sample data arrays
        nsamples = density[0] * density[1]
        self.results['remapping_samples'] = remapping = N.empty(nsamples, 'd')
        self.results['rate_remapping_samples'] = rate_remapping = N.empty(nsamples, 'd')
        self.results['turnover_samples'] = turnover = N.empty(nsamples, 'd')
        self.results['sparsity_samples'] = sparsity = N.empty(nsamples, 'd')
        self.results['stage_coverage_samples'] = stage_coverage = N.empty(nsamples, 'd')
        self.results['stage_repr_samples'] = stage_repr = N.empty(nsamples, 'd')
        self.results['peak_rate_samples'] = peak_rate = N.empty(nsamples, 'd')
        self.results['max_rate_samples'] = max_rate = N.zeros(nsamples, 'd')
        self.results['num_fields_samples'] = num_fields = N.zeros(nsamples, 'd')
        self.results['coverage_samples'] = coverage = N.zeros(nsamples, 'd')
        self.results['area_samples'] = area = N.zeros(nsamples, 'd')
        self.results['diameter_samples'] = diameter = N.zeros(nsamples, 'd')
        self.results['peak_samples'] = peak = N.zeros(nsamples, 'd')
        self.results['average_samples'] = average = N.zeros(nsamples, 'd')
        
        # Method for creating interpolated maps of collated data
        def interpolate_data(z, pixels=256):
            """Interpolate value z across sample points with *density* points
            """
            M = N.empty((pixels,)*2, 'd')
            f = BilinearInterp2D(x=pts_x, y=pts_y, z=z)
            x_range = N.linspace(bounds[0,0], bounds[0,1], num=pixels)
            y_range = N.linspace(bounds[1,1], bounds[1,0], num=pixels)
            for j, x in enumerate(x_range):
                for i, y in enumerate(y_range):
                    M[i,j] = f(x, y)
            return M
        
        # Execute data collection process for each sample point
        tasks = []
        for i, p in enumerate(pts):
            self.out('Submitting: d_%s = %.2f, d_%s = %.2f'%
                (x_type, p[0], y_type, p[1]))
            save_file = os.path.join(map_dir, 'map_%03d.tar.gz'%i)
            tasks.append(
                tc.run(
                    IPclient.MapTask(run_sample_point, 
                        args=(save_file, float(p[0]), float(p[1])))))
        tc.barrier(tasks)
        
        # Collate data return from task farming
        for i in xrange(nsamples):
            self.out('Loading data from map %d for analysis...'%i)
            B = CheckeredRatemap.fromfile(os.path.join(map_dir, 'map_%03d.tar.gz'%i))
            
            # Get field and unit data record arrays
            fdata = B.get_field_data()
            udata = B.get_unit_data()
        
            # Collate the stage map data
            sparsity[i] = B.sparsity
            stage_coverage[i] = B.stage_coverage
            stage_repr[i] = B.stage_repr
            peak_rate[i] = B.peak_rate
        
            # Collate the per-unit data
            if udata.shape[0] != 0:
                max_rate[i] = udata['max_r'].mean()
                num_fields[i] = udata['num_fields'].mean()
                coverage[i] = udata['coverage'].mean()
        
            # Collate the per-field data
            if fdata.shape[0] != 0:
                area[i] = fdata['area'].mean()
                diameter[i] = fdata['diameter'].mean()
                peak[i] = fdata['peak'].mean()
                average[i] = fdata['average'].mean()
        
            # Compute remapping strength from map A
            cmp_AB = compare_AB(A, B)
            remapping[i] = cmp_AB['remapping']
            rate_remapping[i] = cmp_AB['rate_remapping']
            turnover[i] = cmp_AB['turnover']
                    
        # Create interpolated maps for the collated data
        def dot(): 
            self.out.printf('.', color='purple')
        self.out('Creating interpolated parameter maps for collected data'); dot()
        self.results['remapping'] = interpolate_data(remapping); dot()
        self.results['rate_remapping'] = interpolate_data(rate_remapping); dot()
        self.results['turnover'] = interpolate_data(turnover); dot()
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
        """Create a simple 2D image plot of the parameter sweep"""
        
        # Figure is horizontal container for main plot + colorbar
        self.figure = \
            container = HPlotContainer(fill_padding=True, padding=25, 
                bgcolor='linen')
        
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
            color=(0.5, 0.6, 0.7, 0.4), marker_size=4)
        
        # Tweak main plot
        p.title = snake2title(self.display_data)
        p.x_axis.orientation = 'bottom'
        p.x_axis.title = 'Spatial Phase (cm)'
        p.y_axis.title = 'Orientation (rads)'
        p.plots['samples'][0].visible = self.show_sample_points
    
        # Add main plot and colorbar to figure
        container.add(p)
        container.add(
            self.get_colorbar_plot(bounds=(raw_data.min(), raw_data.max())))
        
        # Set radio buttons
        self.unit_data = self.field_data = 'none'
        

# Convenience function to reorganize results data

def get_module_columns(res, module_dim='y', which='remapping'):
    """Get matrix of columns of line data from results samples to plot
    
    Arguments:
    res -- results dict from a completed RealignmentSweep analysis object
    module_dim -- set to 'x' or 'y' to specify modularity axis
    which -- which data to retrieve ('remapping', 'turnover', etc.)
    
    Returns modules array, sweep (realignment) array, and column data matrix.
    """
    pts, data = res['samples'], res[which+'_samples']
    
    # Get the module and sweep information 
    mod_dim = int(module_dim == 'y')
    modules = N.unique(pts[:,mod_dim]).astype('i')
    sweep = pts[pts[:,mod_dim] == modules[0], 1-mod_dim]
    
    # Fill the column matrix
    lines = N.empty((len(modules), len(sweep)), 'd')
    for m,module in enumerate(modules):
        pts_ix = (pts[:,mod_dim] == module).nonzero()[0]
        lines[:,m] = data[pts_ix]
    
    return modules, sweep, lines
