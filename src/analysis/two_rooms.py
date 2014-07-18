#encoding: utf-8
"""
grid.analysis.two_rooms -- AbstractAnalysis subclass for examining how spatial maps
    remap between two orthogonal environments (i.e., different "rooms").

Written by Joe Monaco, 06/16/2008.
Copyright (c) 2008 Columbia University. All rights reserved.
"""

# Library imports
from IPython.kernel.client import MapTask
import numpy as N, os

# Package imports
from ..place_network import PlaceNetworkStd
from ..ratemap import CheckeredRatemap
from ..dmec import GridCollection
from ..core.analysis import AbstractAnalysis
from ..tools.array_container import archive_load
from .compare import compare_AB


def smooth_remap_task(savefile, r_step):
    gc.collect()

    # Load realignment parameters
    delta_phi = rdict['delta_phi']
    delta_psi = rdict['delta_psi']
    ell_mag = rdict['ell_mag']
    ell_angle = rdict['ell_angle']
    zoom = rdict['zoom']
    
    # Realign the inputs according to parameters
    EC.reset()
    for m, m_ix in enumerate(modules):
        if shift:
            EC.shift(r_step * delta_phi[m], mask=m_ix)
        if rotate:
            EC.rotate(r_step * delta_psi[m], mask=m_ix)
        if ellipticity:
            EC.ell_mag[m_ix] = r_step * ell_mag[m]
            EC.ell_angle[m_ix] = r_step * ell_angle[m]
        if rescaling:
            EC.zoom_scale[m_ix] = 1 + r_step * (zoom[m] - 1)
    
    # Simulate the spatial map and save to the specified path
    model = PlaceNetworkStd(EC=EC, W=W, **pdict)
    model.advance()
    pmap = CheckeredRatemap(model)
    pmap.compute_coverage()
    pmap.tofile(savefile)
    return                

class SmoothRemap(AbstractAnalysis):
    
    """
    Assess remapping for smooth transition between two environments
    
    See core.analysis.AbstractAnalysis documentation and collect_data method
    signature and docstring for usage.
    """
    
    label = 'smooth realign'
    
    def collect_data(self, nsteps=10, nmodules=1, freq_modules=False, 
        shift=False, rotate=False, delta_phi=None, delta_psi=None, 
        ellipticity=False, rescaling=False, **kwargs):
        """
        A random environment A will be computed and then a random realignment 
        will be determined to achieve the end-point environment B. A specified 
        number of intermediary environments (nsteps) will be computed by scaling 
        the magnitude of the phase and orientation realignment vectors.
        
        Optional keyword arguments:
        load_dir -- if specified, will look for saved data to re-analyze
        nsteps -- end-point inclusive number of environments in the sweep
        nmodules -- number of independent alignment modules; if set to >1, then
            random alignments are sampled for each module ignoring keyword
            specification (default 1)
        shift -- whether to perform phase realignment 
        rotate -- whether to perform orientation realignment
        delta_phi -- tuple (x,y) specifying phase realignment
        delta_psi -- radian value specificing orientation realignment
        ellipticity -- whether to include elliptical modulation of grids
        rescaling -- whether to include rescaling (enlargement) of grids    
        
        Extra keywords will be passed to the model instantiation.
    
        The pair-wise distance correlation measure of remapping strength as well
        as total spatial correlations to A and B will be stored as the output.        
        """
        
        map_dir = os.path.join(self.datadir, 'data')
        
        if os.path.exists(map_dir):
            
            # Set data directory
            self.out('Using saved data from:\n%s'%map_dir)
            
            # Load pickled results for realignment information
            old_data = N.load(os.path.join(self.datadir, 'analysis.pickle'))
            self.results['realignment'] = old_data['realignment']
            self.results['delta_phi'] = old_data['delta_phi']
            self.results['delta_psi'] = old_data['delta_psi']
            self.results['ellipticity'] = old_data['ellipticity']
            self.results['elliptical_angles'] = old_data['elliptical_angles']
            self.results['rescaling'] = old_data['rescaling']
            self.results['nmodules'] = old_data['nmodules']
            self.results['freq_modules'] = old_data['freq_modules']

        else:
            
            # Save module information to results
            self.results['nmodules'] = nmodules = int(nmodules)
            self.results['freq_modules'] = freq_modules
            
            # Make the map directory
            os.mkdir(map_dir)
        
            # Get ipcontroller clients
            mec = self.get_multiengine_client()
            tc = self.get_task_client()
            
            # Setup namespace on the ipengines
            self.out('Setting up ipengines for task-farming...')
            mec.clear_queue()
            mec.reset()
            mec.execute('import gc')
            mec.execute('from grid_remap.place_network import PlaceNetworkStd')
            mec.execute('from grid_remap.dmec import GridCollection')
            mec.execute('from grid_remap.ratemap import CheckeredRatemap')
            
            # Set default model parameters and push to ipengines
            pdict = dict(   growl=False, 
                            refresh_weights=False, 
                            refresh_orientation=False,
                            refresh_phase=False
                            )
            pdict.update(kwargs)
            mec.push(dict(pdict=pdict))
            
            # Send some network weights and a grid cell object
            self.out('Pushing network weights and grid configuration...')
            EC = GridCollection()
            modules = EC.get_modules(nmodules, freq_sort=freq_modules)
            mec.push(dict(  W=PlaceNetworkStd(EC=EC, **pdict).W, 
                            spacing=EC.spacing, 
                            phi=EC._phi, 
                            psi=EC._psi, 
                            modules=modules,
                            shift=shift,
                            rotate=rotate,
                            ellipticity=ellipticity,
                            rescaling=rescaling
                            ))
            mec.execute('EC = GridCollection(' +
                            'ellipticity=ellipticity, ' +
                            'zoom=rescaling, ' + 
                            'spacing=spacing, ' + 
                            '_phi=phi, _psi=psi)')
            
            # Set phase realignment
            d_phi_list = []
            if shift:
                if nmodules == 1 and delta_phi is not None:
                    d_phi_list.append(N.array(delta_phi))
                else:
                    grid_scale = None
                    for i in xrange(nmodules):
                        if freq_modules:
                            grid_scale = EC.spacing[modules[i]].max()
                        d_phi_list.append(
                            GridCollection.get_delta_phi(scale=grid_scale))
            d_phi = N.array(d_phi_list, 'd')
            self.results['delta_phi'] = d_phi

            # Set orientation realignment
            d_psi_list = []
            if rotate:
                if nmodules == 1 and delta_psi is not None:
                    d_psi_list.append(float(delta_psi))
                else:
                    for i in xrange(nmodules):
                        d_psi_list.append(GridCollection.get_delta_psi())
            d_psi = N.array(d_psi_list, 'd')
            self.results['delta_psi'] = d_psi
            
            # Set ellipticity of response
            ell_mag_list = []
            ell_angle_list = []
            if ellipticity:
                for m in xrange(nmodules):
                    ell_mag_list.append(GridCollection.get_ellipticity())
                    ell_angle_list.append(GridCollection.get_elliptic_angle())
            ell_mags = N.array(ell_mag_list)
            ell_angles = N.array(ell_angle_list)
            self.results['ellipticity'] = ell_mags
            self.results['elliptical_angles'] = ell_angles
            
            # Set rescaling factors for grid response
            zoom_list = []
            if rescaling:
                for m in xrange(nmodules):
                    zoom_list.append(GridCollection.get_zoom_scale())
            zooms = N.array(zoom_list)
            self.results['rescaling'] = zooms
            
            # Push a dictionary of pre-computed realignment parameters
            mec.push(dict(rdict=dict(   delta_phi=d_phi,
                                        delta_psi=d_psi,
                                        ell_mag=ell_mags,
                                        ell_angle=ell_angles,
                                        zoom=zooms
                                        )))
            
            # Compute spatial maps and coverage data
            self.results['realignment'] = r_step = N.linspace(0, 1, nsteps)
            tasks = []
            for i, c in enumerate(r_step):
                self.out('Setting up %.1f%% realignment task...'%(100*c))
                savefile = os.path.join(map_dir, 'map_%03d.tar.gz'%i)
                tasks.append(
                    tc.run(MapTask(smooth_remap_task, args=(savefile, c))))
            self.out('Waiting for tasks to complete...')
            tc.barrier(tasks)
            self.out('...done!')
        
        # Load the saved spatial map data
        map_list = []
        nsteps = 0
        while os.path.isfile(os.path.join(map_dir, 'map_%03d.tar.gz'%nsteps)):
            self.out('Loading map %d data...'%nsteps)
            map_list.append(
                archive_load(os.path.join(map_dir, 'map_%03d.tar.gz'%nsteps)))
            nsteps += 1
        map_A, map_B = map_list[0], map_list[-1]
        
        # Calculate statistics on changes between the maps
        self.out('Computing remapping and correlation results...')
        self.results['remapping'] = remapping = N.empty(nsteps, 'd')
        self.results['rate_remapping'] = rate_remapping = N.empty(nsteps, 'd')
        self.results['turnover'] = turnover = N.empty(nsteps, 'd')
        self.results['A_corr'] = A_corr = N.empty(nsteps, 'd')
        self.results['B_corr'] = B_corr = N.empty(nsteps, 'd')
        self.results['remapping_matrix'] = P = N.empty((nsteps, nsteps), 'd')
        self.results['rate_remapping_matrix'] = R = N.empty((nsteps, nsteps), 'd')
        self.results['turnover_matrix'] = T = N.empty((nsteps, nsteps), 'd')
        self.results['ratecorr_matrix'] = C = N.empty((nsteps, nsteps), 'd')
        for i in xrange(nsteps):
            self.out.printf('.', color='green')
            cmpAi = compare_AB(map_A, map_list[i])
            cmpBi = compare_AB(map_B, map_list[i])
            remapping[i] = cmpAi['remapping']
            rate_remapping[i] = cmpAi['rate_remapping']
            turnover[i] = cmpAi['turnover']
            A_corr[i] = cmpAi['ratecorr']
            B_corr[i] = cmpBi['ratecorr']
            for j in xrange(i, nsteps):
                cmpij = compare_AB(map_list[i], map_list[j])
                P[i,j] = P[j,i] = cmpij['remapping']
                R[i,j] = R[j,i] = cmpij['rate_remapping']
                T[i,j] = T[j,i] = cmpij['turnover']
                C[i,j] = C[j,i] = cmpij['ratecorr']
        self.out.printf('\n')
        
        # Good-bye
        self.out('All done!')
    
    def create_plots(self):
        from matplotlib.pyplot import figure, axis, draw
        from matplotlib import cm
        
        # Alias to results data
        r = self.results
        
        # Create the figure
        self.figure = f = figure(figsize=(12,9))
        
        # Add the remapping lineplot
        ax = f.add_subplot(211)
        ax.plot(100*r['realignment'], r['remapping'], 'k-', lw=2, label='position')
        ax.plot(100*r['realignment'], r['rate_remapping'], 'k--', lw=2, label='rate')
        ax.plot(100*r['realignment'], r['A_corr'], 'b-', label='corr(A)')
        ax.plot(100*r['realignment'], r['B_corr'], 'r-', label='corr(B)')
        ax.set_xlabel('% Realigned')
        ax.set_ylabel('Remapping [Correlation]')
        ax.set_title('Smooth Realignment')
        ax.legend(loc=9)
        ax.set_ylim((0,1))
        
        # Pcolor plot function
        def add_matrix_plot(subp, M, name):
            """Add the specified matrix pcolor"""
            ax = f.add_subplot(subp)
            ax.pcolor(M, cmap=cm.hot)
            axis('image')
            if N.fmod(N.fmod(subp, 230), 3) == 1:
                ax.set_ylabel('A --> B')
            ax.set_xlabel('A --> B [%s]'%name)
        
        # Plot the matrix results data
        add_matrix_plot(234, r['remapping_matrix'], 'position')
        add_matrix_plot(235, r['rate_remapping_matrix'], 'rate')
        add_matrix_plot(236, r['ratecorr_matrix'], 'corr')
        
        # Redraw the figure
        draw()


def sample_remap_task(rdict, pdict):
    gc.collect()

    do_random = rdict['do_random']
    shift = rdict['shift']
    rotate = rdict['rotate']
    ellipticity = rdict['ellipticity']
    rescaling = rdict['rescaling']
    
    # Run A-B remapping experiment according to parameters
    EC = GridCollection()
    modules = EC.get_modules(nmodules, freq_sort=freq_modules)
    model = PlaceNetworkStd(EC=EC, **pdict)
    model.advance()
    A = CheckeredRatemap(model)
    A.compute_coverage()
    
    # Activate map transforms
    EC.ellipticity = ellipticity
    EC.zoom = rescaling
    
    if do_random:
        model.EC = GridCollection()
    else:
        grid_scale = None
        for m_ix in modules:
            if shift:
                if freq_modules:
                    grid_scale = EC.spacing[m_ix].max()
                EC.shift(GridCollection.get_delta_phi(scale=grid_scale), mask=m_ix)
            if rotate:
                EC.rotate(GridCollection.get_delta_psi(), mask=m_ix)
            if ellipticity:
                EC.ell_mag[m_ix] = GridCollection.get_ellipticity()
                EC.ell_angle[m_ix] = GridCollection.get_elliptic_angle()
            if rescaling:
                EC.zoom_scale[m_ix] = GridCollection.get_zoom_scale()
    
    model.reset()
    model.advance()
    B = CheckeredRatemap(model)
    B.compute_coverage()

    # Collate results dictionary for return
    AB = compare_AB(A, B)
    res = dict( remapping=AB['remapping'],
                rate_remapping=AB['rate_remapping'],
                turnover=AB['turnover']                
                )
        
    return res
        
class SampleRemap(AbstractAnalysis):
    
    """
    Assess remapping for a random sampling of A-B realignment experiments.
    
    See core.analysis.AbstractAnalysis documentation and collect_data method
    signature and docstring for usage.
    """
    
    label = 'Remap Sample'
    
    def collect_data(self, nsamples=10, nmodules=1, freq_modules=False, 
        shift=False, rotate=False, ellipticity=False, rescaling=False, do_random=False, 
        **kwargs):
        """
        Perform a series of A-B remapping experiments under observed parameters
        of cortical realignment. Save statistics about rate and positional
        remapping strength.
        
        All realignment is off by default, you must set the realignment types to
        True that you want to enable.
        
        Keyword arguments:
        nsamples -- number of random samples to compute for the analysis
        nmodules -- number of cortical alignment modules (default 1)
        rotate -- include orientation-based realignment (default False)
        shift -- include spatial phase-based realignment (default False)
        ellipticity -- whether to include elliptical modulation of grids
        rescaling -- whether to include rescaling (enlargement) of grids    
        do_random -- map B is another random environment (default False)
        
        Stores *remapping*, *rate_remapping* and *turnover* means, confidence
        intervals and SDs.
        """
        
        # Set default model parameters
        pdict = dict(   growl=False, 
                        monitoring=True, 
                        refresh_weights=False, 
                        refresh_orientation=False,
                        refresh_phase=False
                        )
        pdict.update(kwargs)
        self.results['nmodules'] = nmodules = int(nmodules)
        self.results['freq_modules'] = freq_modules
        
        # Get ipcontroller clients
        mec = self.get_multiengine_client()
        tc = self.get_task_client()
        
        # Setup namespace on the ipengines
        self.out('Setting up ipengines for task-farming...')
        mec.clear_queue()
        mec.reset()
        mec.execute('import numpy as N, gc')
        mec.execute('from grid_remap.place_network import PlaceNetworkStd')
        mec.execute('from grid_remap.dmec import GridCollection')
        mec.execute('from grid_remap.ratemap import CheckeredRatemap')
        mec.execute('from grid_remap.analysis.compare import compare_AB')
        mec.push(dict(nmodules=nmodules, freq_modules=freq_modules))
        
        # Confidence interval computation
        def interval(values):
            return 1.96 * N.std(values) / N.sqrt(len(values))
        
        # Saved data structures
        self.results['remapping_samples'] = remapping = N.empty(nsamples, 'd')
        self.results['rate_remapping_samples'] = rate_remapping = N.empty(nsamples, 'd')
        self.results['turnover_samples'] = turnover = N.empty(nsamples, 'd')
        
        # Set up dictionary of remapping parameters
        rdict = dict(   do_random=do_random,
                        shift=shift,
                        rotate=rotate,
                        ellipticity=ellipticity,
                        rescaling=rescaling
                        )
        
        # Main loop generating random remapping experiments
        tasks = []
        for i in xrange(nsamples):
            self.out('Setting up remapping task %d...'%i)
            tasks.append(
                tc.run(
                    MapTask(sample_remap_task, args=(rdict, pdict))))
        tc.barrier(tasks)
        
        # Collate results
        for i,t_id in enumerate(tasks):
            res = tc.get_task_result(t_id)
            remapping[i] = res['remapping']
            rate_remapping[i] = res['rate_remapping']
            turnover[i] = res['turnover']
        
        # Compute confidence intervals and S.D.'s
        self.results['remapping'] = remapping.mean()
        self.results['remapping_int'] = interval(remapping)
        self.results['remapping_std'] = remapping.std()
        self.results['rate_remapping'] = rate_remapping.mean()
        self.results['rate_remapping_int'] = interval(rate_remapping)
        self.results['rate_remapping_std'] = rate_remapping.std()
        self.results['turnover'] = turnover.mean()
        self.results['turnover_int'] = interval(turnover)
        self.results['turnover_std'] = turnover.std()
        
        # Good-bye
        self.out("All done!")
