#encoding: utf-8
"""
grid.analysis.point -- AbstractAnalysis subclass for significant random sampling 
    of a single point in parameter space.

Written by Joe Monaco, 06/08/2008.
Updated for task-farming, 01/22/2010.
Copyright (c) 2008 Columbia University. All rights reserved.
Copyright (c) 2010 Johns Hopkins University. All rights reserved.
"""

# Library imports
import numpy as N, scipy as S, os
from IPython.kernel import client as IPclient

# Package imports
from .. import PlaceNetworkStd, CheckeredRatemap, GridCollection
from ..core.analysis import AbstractAnalysis


# Task function to simulate an individual spatial map
def simulate_map(**kwargs):
    """Simulate a spatial map based on provided model parameters
    """
    gc.collect()
    
    # Handle file save
    do_save = False
    if 'save_file' in kwargs:
        do_save = True
        save_file = kwargs['save_file']
        del kwargs['save_file']
    
    # Create new input set for this simulation
    kwargs['EC'] = GridCollection()
            
    # Compute and return the spatial map
    model = PlaceNetworkStd(**kwargs)
    model.advance()
    pmap = CheckeredRatemap(model)
    pmap.compute_coverage()
    
    # Do file save if specified
    if do_save:
        pmap.tofile(save_file)
    
    # Reduce spatial map data for output
    data = dict(
        sparsity = pmap.sparsity,
        stage_coverage = pmap.stage_coverage,
        stage_repr = pmap.stage_repr,
        peak_rate = pmap.peak_rate,
        udata = pmap.get_unit_data(),
        fdata = pmap.get_field_data()
        )    
    return data
    

class PointSample(AbstractAnalysis):
    
    """
    Perform a random-environment, random-network sampling of a single point in
    parameter space. The results are simply a single set of means and 
    confidence intervals for the available spatial map statistics.
    
    See core.analysis.AbstractAnalysis documentation and collect_data method
    signature and docstring for usage.
    """
    
    label = "point sample"
    
    def collect_data(self, nsamples=25, save_maps=False, **kwargs):
        """Collect statistics from an environmentally and network-wise random 
        sample of the parameter point specified by keyword arguments.
        
        Set save_maps to True to save the spatial maps for the sample.
        """
        # Get ipcontroller clients
        mec = self.get_multiengine_client()
        tc = self.get_task_client()
                
        # Setup namespace on ipengine instances
        self.out('Setting up ipengines for task-farming...')
        mec.clear_queue()
        mec.reset()
        mec.execute('import gc')
        mec.execute('from grid_remap.place_network import PlaceNetworkStd')
        mec.execute('from grid_remap.dmec import GridCollection')
        mec.execute('from grid_remap.ratemap import CheckeredRatemap')
        
        # Setup default model parameters
        pdict = dict(   growl=False, 
                        desc='point', 
                        projdir=self.datadir, 
                        monitoring=True, 
                        dwell_factor=6.5, 
                        num_trials=1)
        pdict.update(kwargs)
        
        # Create save directory if specified
        if save_maps:
            save_dir = os.path.join(self.datadir, 'data')
            if not os.path.isdir(save_dir):
                os.makedirs(save_dir)
        
        # Create simulation tasks, run them and wait for completion
        self.out('Creating and running simulation tasks...')
        tasks = []
        for sample in xrange(nsamples):
            if save_maps:
                pdict['save_file'] = \
                    os.path.join(save_dir, 'map_%03d.tar.gz'%sample)
            tasks.append(
                tc.run(
                    IPclient.MapTask(simulate_map, kwargs=pdict)))
        tc.barrier(tasks)
            
        # Gather the spatial map results
        self.out('Collecting results from sample simulations...')
        sample_list = [tc.get_task_result(t_id) for t_id in tasks]
        udata_list = [data['udata'] for data in sample_list]
        fdata_list = [data['fdata'] for data in sample_list]

        # Collate the spatial map data
        self.out('Collating results data and statistics...')
        sparsity_arr = N.array([data['sparsity'] for data in sample_list])
        stage_coverage_arr = N.array([data['stage_coverage'] for data in sample_list])
        stage_repr_arr = N.array([data['stage_repr'] for data in sample_list])
        peak_rate_arr = N.array([data['peak_rate'] for data in sample_list])
        
        # Confidence interval computation
        def interval(values):
            return 1.96 * N.std(values) / N.sqrt(len(values))

        # Compute basic statistics on collated spatial map data
        self.results['sparsity'] = N.mean(sparsity_arr)
        self.results['sparsity_std'] = N.std(sparsity_arr)
        self.results['sparsity_err'] = interval(sparsity_arr)
        self.results['stage_coverage'] = N.mean(stage_coverage_arr)
        self.results['stage_coverage_std'] = N.std(stage_coverage_arr)
        self.results['stage_coverage_err'] = interval(stage_coverage_arr)
        self.results['stage_repr'] = N.mean(stage_repr_arr)
        self.results['stage_repr_std'] = N.std(stage_repr_arr)
        self.results['stage_repr_err'] = interval(stage_repr_arr)
        self.results['peak_rate'] = N.mean(peak_rate_arr)
        self.results['peak_rate_std'] = N.std(peak_rate_arr)
        self.results['peak_rate_err'] = interval(peak_rate_arr)
        
        # Collate the per-unit data
        _max_rate = N.array([], 'd')
        _num_fields = N.array([], 'd')
        _coverage = N.array([], 'd')
        for udata in udata_list:
            _max_rate = N.r_[_max_rate, udata['max_r']]
            _num_fields = N.r_[_num_fields, udata['num_fields']]
            _coverage = N.r_[_coverage, udata['coverage']]

        # Compute per-unit statistics
        self.results['max_rate'] = _max_rate.mean()
        self.results['max_rate_std'] = _max_rate.std()
        self.results['max_rate_err'] = interval(_max_rate)
        self.results['num_fields'] = _num_fields.mean()
        self.results['num_fields_std'] = _num_fields.std()
        self.results['num_fields_err'] = interval(_num_fields)
        self.results['coverage'] = _coverage.mean()
        self.results['coverage_std'] = _coverage.std()
        self.results['coverage_err'] = interval(_coverage)
        
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
        
        # Compute per-field statistics
        self.results['area'] = _area.mean()
        self.results['area_std'] = _area.std()
        self.results['area_err'] = interval(_area)
        self.results['diameter'] = _diameter.mean()
        self.results['diameter_std'] = _diameter.std()
        self.results['diameter_err'] = interval(_diameter)
        self.results['peak'] = _peak.mean()
        self.results['peak_std'] = _peak.std()
        self.results['peak_err'] = interval(_peak)
        self.results['average'] = _average.mean()
        self.results['average_std'] = _average.std()
        self.results['average_err'] = interval(_average)
            
        # Good-bye!
        self.out('All done!')


def write_sample_file(sample, fn=None, sep='\t'):
    """
    Write out a simple tab-separated file of the sample statistics
    """
    if fn is None:
        fn = 'sample.dat'
    fn = os.path.realpath(fn)
    if os.path.exists(fn):
        os.sys.stderr.write('File path already exists:\n%s\n'%fn)
        return
    fd = file(fn, 'w')
    hdr = ['Value', 'Mean', '95-CI', 'SD']
    fd.write(sep.join(hdr) + '\n\n')
    if type(sample) is type({}):
        data = sample
    elif type(sample) is PointSample:
        data = sample.results
    else:
        os.sys.stderr.write('Bad sample type: %s'%type(sample))
    keys = data.keys()
    m_values = ['sparsity', 'stage_coverage', 'stage_repr', 'peak_rate']
    u_values = ['max_rate', 'num_fields', 'coverage']
    f_values = ['area', 'diameter', 'peak', 'average']
    for v_list in (m_values, u_values, f_values):
        for v in v_list:
            rec = [v, '%.6e'%data[v], '%.6e'%data[v+'_err'], '%.6e'%data[v+'_std']]
            fd.write(sep.join(rec) + '\n')
        fd.write('\n')
    fd.close()
