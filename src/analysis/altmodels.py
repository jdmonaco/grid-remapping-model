#encoding: utf-8
"""
grid.analysis.altmodels -- Analysis simulating model variants for comparison

Exports: ModelComparison

Written by Joe Monaco, 02/05/2011.
Copyright (c) 2011 Johns Hopkins University. All rights reserved.
"""

# Library imports
from scipy.stats import sem
import os, numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
    
# Package imports
from ..place_network import PlaceNetworkStd
from ..core.analysis import AbstractAnalysis
from ..tools.images import array_to_image
from ..ratemap import CheckeredRatemap
from ..dmec import GridCollection
from .compare import compare_AB
from .map_funcs import get_tuned_weights


class ModelComparison(AbstractAnalysis):
    
    """
    Load a standard simulation from pre-existing data (or simulate a new map) 
    and then simulate several model variants to compare place fields size
    and location differences.
    
    See core.analysis.AbstractAnalysis documentation and collect_data method
    signature and docstring for usage.
    """
    
    label = "alt models"
    
    def collect_data(self, load_dir=None, alpha=0.3, gamma=1.0, rec_tuned=False):
        """Run a standard simulation and then variants using the same network
        
        Keyword arguments:
        load_dir -- if loading pre-existing network, set directory here
        alpha -- learning parameter for tuned weights (get_tuned_weights)
        gamma -- gain of recurrent excitation (based on overlap)
        rec_tuned -- whether recurrent variant is based on tuned output (True)
            or the standard output (False)
        
        Set save_maps to True to save the spatial maps for the sample.
        """
        self.results['model_types'] = ('std', 'fwd', 'tuned', 'rec')        
        if load_dir is not None:
            if not os.path.isdir(load_dir):
                raise ValueError, 'invalid load directory'
            self.results['load_dir'] = os.path.abspath(load_dir)
            self.out('Loading network from\n%s...'%self.results['load_dir'])
            os.chdir(load_dir)
            l = np.load
            EC = GridCollection(
                _phi=l('phi.npy'), _psi=l('psi.npy'), spacing=l('spacing.npy'))
            model = PlaceNetworkStd(EC=EC, W=l('W.npy'), refresh_weights=False)
            os.chdir(self.datadir)
        else:
            self.out('Creating new grid inputs and place network...')
            EC = GridCollection()
            model = PlaceNetworkStd(EC=EC)
        W = model.W
        
        def get_norms(M):
            return np.sqrt((M**2).sum(axis=0))
        
        def store_data(prefix, pmap):
            udata = pmap.get_unit_data()
            fdata = pmap.get_field_data()
            self.results['%s_sparsity'%prefix] = pmap.sparsity
            self.results['%s_num_fields'%prefix] = udata['num_fields']
            self.results['%s_area'%prefix] = fdata['area']
            self.results['%s_diameter'%prefix] = fdata['diameter']
            self.results['%s_x'%prefix] = fdata['x']
            self.results['%s_y'%prefix] = fdata['y']
            if not os.path.exists('%s_map.tar.gz'%prefix):
                pmap.tofile('%s_map'%prefix)
            return
            
        # Get input strength map
        self.out('Computing grid input strengths...')
        EC_R = EC.get_z_stack()
        EC_norms = get_norms(EC_R)
        np.save('EC_norms.npy', EC_norms)
        array_to_image(EC_norms, 'EC_norms.png', cmap=mpl.cm.gray_r)
        array_to_image(EC_norms, 'EC_norms_jet.png', cmap=mpl.cm.jet)

        # Run the standard simulation
        if not os.path.exists('std_map.tar.gz'):
            self.out('Running standard simulation...')
            model.advance()
            pmap = CheckeredRatemap(model)
        else:
            self.out('Loading standard simulation data...')
            pmap = CheckeredRatemap.fromfile('std_map.tar.gz')
        store_data('std', pmap)
        std_num_active = pmap.num_active
        self.out('Standard active units = %d'%std_num_active)
        R = pmap.Map
        array_to_image(get_norms(R), 'std_norms.png', cmap=mpl.cm.gray_r)
        array_to_image(get_norms(R), 'std_norms_jet.png', cmap=mpl.cm.jet)
        
        def sparsity_match_threshold(Map):
            self.out('Searching for sparsity-matching threshold...')
            N, H, W = Map.shape
            I = np.empty((N,), 'd')
            for i in xrange(N):
                I[i] = Map[i].max()
            
            # Test activity peaks as thresholds to find sparsity-matching threshold
            I.sort()
            R_ = np.empty(Map.shape, 'd') # probe workspace
            thresh = 0
            for i in xrange(N):
                R_[:] = Map # reset
                Rmax = R_.max()
                num_active = 0

                for j in xrange(N):
                    if (R_[j].max()>0.2*Rmax):
                        if (R_[j]>0.2*R_[j].max()).sum() > 50:
                            num_active += 1
                self.out.printf('%d '%num_active)
                if num_active < std_num_active:
                    self.out.printf('\n')
                    self.out('... sparsity match at %.4f ...'%thresh)
                    break

                thresh = I[i] # get next peak 
                R_ -= thresh # and apply test threshold
                R_[R_<0] = 0
            del R_
            if num_active >= std_num_active:
                self.out.printf('\n')
            if thresh:
                Map -= thresh
                Map[Map<0] = 0
            return

        # Run feedforward inhibition simulation
        if not os.path.exists('fwd_map.tar.gz'):
            self.out('Computing feedforward model variant...')
            R[:] = 0 # using R matrix as a spatial map workspace
            for i in xrange(model.N_CA):
                R[i] = model.beta * (W[i].reshape(model.N_EC, 1, 1) * EC_R).sum(axis=0)
            
            # Feedforward inhibition as sparsity-matching threshold
            sparsity_match_threshold(R)
            
            pmap.reset()
            pmap.compute_coverage()
            self.out('Feedforward active units = %d'%pmap.num_active)
        else:
            self.out('Loading feedforward model data...')
            pmap = CheckeredRatemap.fromfile('fwd_map.tar.gz')
            R = pmap.Map
        array_to_image(get_norms(R), 'fwd_norms.png', cmap=mpl.cm.gray_r)
        store_data('fwd', pmap)
        
        # Run associatively tuned simulation
        if not os.path.exists('tuned_map.tar.gz'):
            self.out('Running input tuned simulation (alpha = %.2f)...'%alpha)
            model.W = get_tuned_weights(
                CheckeredRatemap.fromfile('std_map.tar.gz'), W, EC, alpha,
                    grow_synapses=True)
            model.reset()
            model.advance()
            pmap = CheckeredRatemap(model)
            pmap.compute_coverage()
            self.out('Tuned active units = %d'%pmap.num_active)
        else:
            self.out('Loading input tuned model data...')
            pmap = CheckeredRatemap.fromfile('tuned_map.tar.gz')
        R = pmap.Map
        array_to_image(get_norms(R), 'tuned_norms.png', cmap=mpl.cm.gray_r)
        store_data('tuned', pmap)
        
        # Run recurrent excitation simulation
        if not os.path.exists('rec_map.tar.gz'):
            # Construct the E-E weight matrix
            self.out('Constructing E-E weight matrix...')
            if rec_tuned:
                self.out('--> Using input-tuned output as base')
            else:
                self.out('--> Using standard output as base')
                pmap = CheckeredRatemap.fromfile('std_map.tar.gz')
            R = pmap.Map
            N, H, W = R.shape
            J = np.zeros((N, N), 'd')
            for i in xrange(N):
                for j in xrange(i+1, N):
                    J[i,j] = J[j,i] = gamma * \
                        (pmap.single_maps[i] * pmap.single_maps[j]).sum()
                    if J[i,j] > 0:
                        J[i,j] = J[j,i] = J[i,j] / \
                            min(pmap.single_maps[i].sum(), 
                                pmap.single_maps[j].sum())
            
            # Add in first-order recurrent excitation across the map
            self.out('Adding first-order recurrent excitation to map...')
            for i in xrange(H):
                for j in xrange(W):
                    R[:,i,j] += np.dot(R[:,i,j], J) # feedforward
                    R[:,i,j] += np.dot(R[:,i,j], J) # feedback
            
            # Feedforward threshold to maintain activity level
            sparsity_match_threshold(R)
            
            pmap.reset()
            pmap.compute_coverage()
            self.out('Recurrent active units = %d'%pmap.num_active)
        else:
            self.out('Loading recurrent model data...')
            pmap = CheckeredRatemap.fromfile('rec_map.tar.gz')
        R = pmap.Map
        array_to_image(get_norms(R), 'rec_norms.png', cmap=mpl.cm.gray_r)
        store_data('rec', pmap)
        
        # Good-bye!
        self.out('All done!')
        
    def create_plots(self, legend=False):
        # Move into data directoary and start logging
        os.chdir(self.datadir)
        self.out.outfd = file('figure.log', 'w')
        
        # Set up main figure for plotting
        self.figure = {}
        figsize = 8, 10
        plt.rcParams['figure.figsize'] = figsize
        self.figure['altmodels'] = f = plt.figure(figsize=figsize)
        f.suptitle(self.label.title())
        
        # Load data
        data = self.results
        models = data['model_types']
        getval = lambda pre, k: data[pre + '_' + k]
        
        # Log some data
        def print_mean_sem(value, arr):
            if type(arr) is float:
                self.out('%s = %.4f'%(value, arr))
            else:
                self.out('%s = %.4f +/- %.4f'%(value, arr.mean(), sem(arr)))
            
        for prefix in models:
            for val in ('sparsity', 'num_fields', 'area', 'diameter'):
                key = prefix + '_' + val
                print_mean_sem(key, data[key])
        
        # Draw place fields as circles
        def draw_circle_field_plots(ax, prefix):
            x = getval(prefix, 'x')
            y = getval(prefix, 'y')
            d = getval(prefix, 'diameter')
            nfields = len(x)
            ax.plot(x, y, 'k+', ms=6, aa=False)
            for i in xrange(nfields):
                ell = mpl.patches.Ellipse((x[i], y[i]), d[i], d[i], 
                    fill=False, lw=1, ec='k')
                ell.clip_box = ax.bbox
                ax.add_artist(ell)
            ax.axis("image")
            ax.set_xlim(0, 100)
            ax.set_ylim(0, 100)
            ax.set_title(prefix)
            return ax
        
        # Render place field plots
        rows = 3
        cols = 2
        for i,prefix in enumerate(models):
            draw_circle_field_plots(plt.subplot(rows, cols, i+1), prefix)
            
        # Statistics plot
        ax = plt.subplot(rows, cols, 5)
        markers = "ods^"
        for i,prefix in enumerate(models):
            a = getval(prefix, 'area')
            nf = getval(prefix, 'num_fields')
            ax.errorbar(a.mean(), nf.mean(), xerr=sem(a), yerr=sem(nf),
                fmt=markers[i], ecolor='k', elinewidth=1, capsize=4, 
                ms=6, mfc='k', mec='k', mew=1)
        # ax.set_ylim(1, 2)
        # ax.set_xlim(xmax=245)
        ax.set_xlabel('area')
        ax.set_ylabel('num. fields')
        
        # Remapping data
        if os.path.exists('remapping.npy'):
            self.out('Loading remapping/turnover values...')
            remapping, turnover = np.load('remapping.npy') 
        else:
            self.out('Computing remapping/turnover measures...')
            pmaps = [CheckeredRatemap.fromfile('%s_map.tar.gz'%p) for p in models]
            remapping = []
            turnover = []
            for pm in pmaps[1:]:
                cmpAB = compare_AB(pmaps[0], pm)
                remapping.append(cmpAB['remapping'])
                turnover.append(cmpAB['turnover'])
            np.save('remapping.npy', np.array([remapping, turnover]))
        self.out('Remapping: %s'%str(remapping))
        self.out('Turnover: %s'%str(turnover))
        
        # Set up bar plot data
        ax = plt.subplot(rows, cols, 6)
        left = []
        height = []
        xticklabels = models[1:]
        bar_w = 1/float(len(xticklabels))
        c = 0
        for i in xrange(len(xticklabels)):
            left.extend([c-bar_w, c])
            height.extend([remapping[i], turnover[i]])
            c += 1
            
        # Render the bar chart and legend
        bar_cols = mpl.cm.gray(([0.25, 0.6])*c)
        bar_h = ax.bar(left, height, width=bar_w, 
            ec='k', color=bar_cols, linewidth=0, ecolor='k', aa=False)
        if legend:
            ax.legend(bar_h[:2], ['Remapping', 'Turnover'], loc=1)
        ax.hlines(1.0, xmin=-0.5, xmax=c-0.5, linestyle=':', color='k')
        ax.set_xlim(-0.5, c-0.5)
        ax.set_ylim(0.0, 1.1)
        ax.set_xticks(np.arange(c))
        ax.set_xticklabels(xticklabels)

        plt.draw()
        plt.rcParams['figure.figsize'] = plt.rcParamsDefault['figure.figsize']
        self.out.outfd.close()
