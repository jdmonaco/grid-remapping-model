#encoding: utf-8
"""
grid.analysis.compare -- Omnibus spatial map comparison function and associated
    figure plot function

Exported namespace: compare_AB, compare_AB_figure
    
Written by Joe Monaco, 12/19/2008
Copyright (c) 2008 Columbia Unversity. All Rights Reserved.  
"""

# Library imports
import numpy
from scipy.stats import pearsonr

# Package imports
from .map_funcs import remap_quiver_plot
from ..tools.setops import intersection, difference, union, symmetric_difference


def compare_AB(A, B, sparsity=0.614):
    """
    Perform several analyses on spatial maps A and B, returning a dict of the 
    results for visualization and analysis.
    
    Arguments:
    A,B -- PlaceMap subclass instances of spatial maps to be compared
    sparsity -- expected spatial map sparsity for computing turnover
    """    
    results = {}
    
    # Get data for active units
    udata_A = results['udata_A'] = A.get_unit_data()
    udata_B = results['udata_B'] = B.get_unit_data()
    
    # Indices of units active in both rooms
    AB_active = intersection(udata_A['unit'], udata_B['unit'])
    results['num_active'] = num_active = AB_active.shape[0]
    results['frac_active'] = num_active / float(A.num_maps)
    results['A_active'] = A.num_active
    results['B_active'] = B.num_active
    
    # Allocate paired distance arrays
    num_pairs = num_active*(num_active-1)/2
    D_A = results['D_A'] = numpy.empty(num_pairs, 'd')
    D_B = results['D_B'] = numpy.empty(num_pairs, 'd')
    R_A = results['R_A'] = numpy.empty(num_pairs, 'd')
    R_B = results['R_B'] = numpy.empty(num_pairs, 'd')
    
    # Compute pair-wise positional and rate distances in both rooms
    ix = 0
    for i in xrange(num_active):
        x_iA, y_iA, r_iA = A.maxima[AB_active[i]]
        x_iB, y_iB, r_iB = B.maxima[AB_active[i]]
        for j in xrange(i+1, num_active):
            x_jA, y_jA, r_jA = A.maxima[AB_active[j]]
            x_jB, y_jB, r_jB = B.maxima[AB_active[j]]
            D_A[ix] = numpy.sqrt((x_iA-x_jA)**2 + (y_iA-y_jA)**2)
            D_B[ix] = numpy.sqrt((x_iB-x_jB)**2 + (y_iB-y_jB)**2)
            R_A[ix] = (r_iA - r_jA) / (r_iA + r_jA)
            R_B[ix] = (r_iB - r_jB) / (r_iB + r_jB)
            ix += 1
    
    # Distribution of remapped distances for active units
    D_AB = results['D_AB'] = numpy.empty(num_active, 'd')
    D_AB[:] = numpy.sqrt(
        ((A.maxima[AB_active] - B.maxima[AB_active])**2).sum(axis=1))
    
    # Store peak locations for active units in both maps
    results['A_xy'] = A.maxima[AB_active, :2].T
    results['B_xy'] = B.maxima[AB_active, :2].T
    
    # Distribution of active unit rate remapping strength: (max-min)/max
    R_AB = results['R_AB'] = numpy.empty(num_active, 'd')
    peak_rates = numpy.c_[A.maxima[AB_active, 2], B.maxima[AB_active, 2]]
    r_max, r_min = peak_rates.max(axis=1), peak_rates.min(axis=1)
    R_AB[:] = (r_max - r_min) / r_max
        
    # Active environment counts
    counts = results['env_counts'] = numpy.zeros((2, 3), 'i')
    counts[0] = 0, 1, 2
    counts[1,0] = difference(numpy.arange(A.num_maps), 
        union(udata_A['unit'], udata_B['unit'])).shape[0]
    counts[1,1] = symmetric_difference(udata_A['unit'], 
        udata_B['unit']).shape[0]
    counts[1,2] = num_active
    
    # Compute independent "turnover" as 1-RMSD from expected random turnover
    E_rand = numpy.array([sparsity**2, 2*sparsity*(1-sparsity), (1-sparsity)**2])
    E0 = numpy.array([sparsity, 0, 1-sparsity])
    RMSD = lambda cbar: numpy.sqrt(((cbar - E_rand)**2).mean())
    results['turnover'] = 1 - RMSD(counts[1]/float(A.num_maps)) / RMSD(E0)
    
    # Population ratemap correlation coefficient
    results['ratecorr'] = pearsonr(A.Map.flatten(), B.Map.flatten())[0]
    
    # Active pair-wise distance correlation: positional remapping strength
    results['remapping'] = 1.0 - pearsonr(D_A, D_B)[0]
    
    # Active unit peak-rate correlation: rate remapping strength
    results['rate_remapping'] = 1.0 - pearsonr(R_A, R_B)[0]
    
    return results

def compare_AB_figure(r, f=None):
    """
    Visualize some A-B spatial map comparison data in a figure
    """
    from matplotlib.pyplot import figure, draw
    figsize = (13, 10)
    if f is None:
        f = figure(figsize=figsize)
    else:
        f.clf()
        f.set_size_inches(figsize)
    
    # Plot the inter-environmnet paired distance scatter plot
    remap_plot = f.add_subplot(321)
    remap_plot.plot(r['D_A'], r['D_B'], 'b.', ms=1)
    remap_plot.plot([0, numpy.sqrt(2)*100], [0, numpy.sqrt(2)*100], 'k:')
    remap_plot.axis([0, numpy.sqrt(2)*100, 0, numpy.sqrt(2)*100])
    remap_plot.set_xlabel('D(A) (cm)', size='smaller')
    remap_plot.set_ylabel('D(B) (cm)', size='smaller')
    remap_plot.text(3, numpy.sqrt(2)*90, '1 - r = %.2f'%r['remapping'])
    
    # Plot the inter-environmnet paired rate difference scatter plot
    rate_plot = f.add_subplot(323)
    rate_plot.plot(r['R_A'], r['R_B'], 'b.', ms=1)
    rate_plot.axis('tight')
    rate_plot.set_xlabel('R(A)', size='smaller')
    rate_plot.set_ylabel('R(B)', size='smaller')
    rate_plot.text(0.05, 0.9, '1 - r = %.2f'%r['rate_remapping'], 
        transform=rate_plot.transAxes)
    
    # Plot a histogram of inter-environment remapping distances
    dist_hist = f.add_subplot(3, 4, 9)
    dist_hist.hist(r['D_AB'], bins=15, histtype='step', 
        edgecolor='g', lw=2)
    dist_hist.set_xlabel('Remapped Distance', size='smaller')
    dist_hist.set_ylabel('Count', size='smaller')
    v = dist_hist.axis()
    dist_hist.set_ylim(ymax=v[3]+3)
    
    rate_hist = f.add_subplot(3, 4, 10)
    rate_hist.hist(r['R_AB'], bins=15, histtype='step', 
        edgecolor='r', lw=2)
    rate_hist.set_xlabel('Rate Remapping', size='smaller')
    
    # Bar chart of environment counts (# envs where ith cell active)
    # env_plot = f.add_subplot(427)
    # env_plot.plot(r['env_counts'][0], r['env_counts'][1], 'kd', ms=12, mew=3,
    #     mec='k', mfc='w')
    # env_plot.set_xticks(r['env_counts'][0])
    # env_plot.set_xticklabels(['None', 'Single', 'Both'])
    # env_plot.set_xlabel('Environmental Activity', size='smaller')
    # env_plot.set_ylabel('# Cells', size='smaller')
    # env_plot.axis([-0.5, 2.5, 0, 1.1*max(r['env_counts'][1])])
    # env_plot.grid(True)

    # Remapping quiver plot on the right column
    quiver_plot = f.add_subplot(122)
    remap_quiver_plot(r, ax=quiver_plot, rate_colors=True, border_style=False)
    quiver_plot.set_xlabel('X (cm)', size='smaller')
    quiver_plot.set_ylabel('Y (cm)', size='smaller')
    quiver_plot.set_title('Positional/Rate Remapping Vectors')
    draw()
    
    return f    
