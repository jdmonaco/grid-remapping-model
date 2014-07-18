# encoding: utf-8
"""
search.py -- Genetic algorithm parameter space for the PlaceNetwork model

Created by Joe Monaco on 2009-07-30.
Refactored by Joe Monaco on 2010-01-27.
Copyright (c) 2009,2010 Johns Hopkins University. All rights reserved.
"""

from numpy import median
from ..core.search import GeneticSearch


class PlaceNetworkSearch(GeneticSearch):
    
    """
    Perform a basic genetic algorithm search of the parameter space for the
    competitive grid-to-place model using spatial map properties as components
    of the fitness function.
        
    See core.search.GeneticSearch documentation and collect_data method
    signature and docstring for usage.
    """
    
    label = "grid search"    
    min_fitness = 2.4887892376681613

    def set_targets(self):
        """Set objective targets for basic place coding characteristics
        """
        self.targets = dict(
            sparsity = (0.65, 1.0),
            coverage = (1.0, 0.8),
            iscvr = (16.0, 0.8),
            peak_rate = (0.99, 0.9),
            num_fields = (1.0, 1.0),
            area = (200.0, 0.6)
            )

    def setup_engines(self, mec):
        """Import grid modules for simulating spatial maps
        """
        mec.execute('import gc')
        mec.execute('from grid_remap.place_network import PlaceNetworkStd')
        mec.execute('from grid_remap.dmec import GridCollection')
        mec.execute('from grid_remap.ratemap import CheckeredRatemap')
    
    def _sample_func_default(self):
        def sample_point(**kwargs):
            """PlaceNetwork spatial map simulation sampling function
            """
            gc.collect()

            # Set up model parameters
            pdict = dict(   
                growl=False, 
                desc='spawn', 
                dwell_factor=5.0, 
                num_trials=1
                )
            pdict.update(kwargs)

            # Create new input set for this simulation
            pdict['EC'] = GridCollection()
        
            # Simulate and compute the spatial map
            model = PlaceNetworkStd(**pdict)
            model.advance()
            pmap = CheckeredRatemap(model)
            pmap.compute_coverage()
            pmap.data = None 
            
            # Pare down the output data
            sample = dict(
                udata = pmap.get_unit_data(),
                fdata = pmap.get_field_data(),
                sparsity = pmap.sparsity,
                coverage = pmap.stage_coverage,
                peak_rate = pmap.peak_rate,
                stage_repr_map = pmap.stage_repr_map
                )
            return sample
        return sample_point

    def calculate_outputs(self, sample):
        """Provides a scalar fitness value given a spatial map
        """
        # Get data from spatial map
        udata = sample['udata']
        fdata = sample['fdata']
    
        # Extract fitness components
        outputs = {}
        outputs['sparsity'] = sample['sparsity']
        outputs['coverage'] = sample['coverage']
        outputs['peak_rate'] = sample['peak_rate']
        if udata.shape[0]:
            outputs['num_fields'] = udata.num_fields.mean()
        else:
            outputs['num_fields'] = 0.0
        if fdata.peak.shape[0]:
            outputs['area'] = median(fdata.area)
        else:
            outputs['area'] = 0.0
    
        # ISCVR measure for even representation:
        # Inverse square coefficient of variation of representation
        srep = sample['stage_repr_map'].flatten()
        if srep.sum():
            outputs['iscvr'] = (srep.std()/srep.mean())**-2
        else:
            outputs['iscvr'] = 0.0
        return outputs
        