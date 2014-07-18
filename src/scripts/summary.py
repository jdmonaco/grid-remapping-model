#!/usr/bin/env python
# encoding: utf-8
"""
new_summary.py -- Run the sample analyses for the remapping summary figure.

Created by Joe Monaco on 2010-06-21.
Copyright (c) 2010 Johns Hopkins University. All rights reserved.
"""

from grid_remap.analysis.two_rooms import SampleRemap
from scipy.stats import ks_2samp
from numpy import pi, empty, arange, array
import os, cPickle

DATA_DIR = '/Users/joe/projects/grid_model/remapping/sample_sets/'

def main():
    """
    Perform random sampling under various realignment conditions and collate
    the statistical results for display.
    """
    old_dir = os.getcwd()
    os.chdir(DATA_DIR)
    
    # Load data from various sampling conditions
    rnd = SampleRemap.load_data('random/analysis.pickle').results    
    s1 = SampleRemap.load_data('shift/shift_1/analysis.pickle').results
    s1f = SampleRemap.load_data('shift/shift_1_freq/analysis.pickle').results
    s2 = SampleRemap.load_data('shift/shift_2/analysis.pickle').results
    s2f = SampleRemap.load_data('shift/shift_2_freq/analysis.pickle').results
    s4 = SampleRemap.load_data('shift/shift_4/analysis.pickle').results
    s4f = SampleRemap.load_data('shift/shift_4_freq/analysis.pickle').results
    s8 = SampleRemap.load_data('shift/shift_8/analysis.pickle').results
    s8f = SampleRemap.load_data('shift/shift_8_freq/analysis.pickle').results
    s16 = SampleRemap.load_data('shift/shift_16/analysis.pickle').results
    s16f = SampleRemap.load_data('shift/shift_16_freq/analysis.pickle').results
    srnd = SampleRemap.load_data('shift/shift_rnd/analysis.pickle').results
    e1 = SampleRemap.load_data('ellipticity/ellipticity_1/analysis.pickle').results
    e1f = SampleRemap.load_data('ellipticity/ellipticity_1_freq/analysis.pickle').results
    e2 = SampleRemap.load_data('ellipticity/ellipticity_2/analysis.pickle').results
    e2f = SampleRemap.load_data('ellipticity/ellipticity_2_freq/analysis.pickle').results
    e4 = SampleRemap.load_data('ellipticity/ellipticity_4/analysis.pickle').results
    e4f = SampleRemap.load_data('ellipticity/ellipticity_4_freq/analysis.pickle').results
    e8 = SampleRemap.load_data('ellipticity/ellipticity_8/analysis.pickle').results
    e8f = SampleRemap.load_data('ellipticity/ellipticity_8_freq/analysis.pickle').results
    e16 = SampleRemap.load_data('ellipticity/ellipticity_16/analysis.pickle').results
    e16f = SampleRemap.load_data('ellipticity/ellipticity_16_freq/analysis.pickle').results
    ernd = SampleRemap.load_data('ellipticity/ellipticity_rnd/analysis.pickle').results
    z1 = SampleRemap.load_data('rescaling/rescaling_1/analysis.pickle').results
    z1f = SampleRemap.load_data('rescaling/rescaling_1_freq/analysis.pickle').results
    z2 = SampleRemap.load_data('rescaling/rescaling_2/analysis.pickle').results
    z2f = SampleRemap.load_data('rescaling/rescaling_2_freq/analysis.pickle').results
    z4 = SampleRemap.load_data('rescaling/rescaling_4/analysis.pickle').results
    z4f = SampleRemap.load_data('rescaling/rescaling_4_freq/analysis.pickle').results
    z8 = SampleRemap.load_data('rescaling/rescaling_8/analysis.pickle').results
    z8f = SampleRemap.load_data('rescaling/rescaling_8_freq/analysis.pickle').results
    z16 = SampleRemap.load_data('rescaling/rescaling_16/analysis.pickle').results
    z16f = SampleRemap.load_data('rescaling/rescaling_16_freq/analysis.pickle').results
    zrnd = SampleRemap.load_data('rescaling/rescaling_rnd/analysis.pickle').results
    
    # Data collation
    labels = ('rnd', 's1', 's2', 's4', 's8', 's16', 's1f', 's2f', 's4f', 's8f', 's16f', 'srnd',
        'e1', 'e2', 'e4', 'e8', 'e16', 'e1f', 'e2f', 'e4f', 'e8f', 'e16f', 'ernd',
        'z1', 'z2', 'z4', 'z8', 'z16', 'z1f', 'z2f', 'z4f', 'z8f', 'z16f', 'zrnd')
    results = (rnd, s1, s2, s4, s8, s16, s1f, s2f, s4f, s8f, s16f, srnd,
        e1, e2, e4, e8, e16, e1f, e2f, e4f, e8f, e16f, ernd,
        z1, z2, z4, z8, z16, z1f, z2f, z4f, z8f, z16f, zrnd)

    # Create column-wise summary data matrix
    N_conds = len(labels)
    data = empty((N_conds, 9), 'd')
    for i in xrange(N_conds):
        data[i,0] = results[i]['remapping']
        data[i,1] = results[i]['rate_remapping']
        data[i,2] = results[i]['turnover']
        data[i,3] = results[i]['remapping_int']
        data[i,4] = results[i]['rate_remapping_int']
        data[i,5] = results[i]['turnover_int']
        data[i,6] = results[i]['remapping_std']
        data[i,7] = results[i]['rate_remapping_std']
        data[i,8] = results[i]['turnover_std']
    
    # Create K-S p-value matrices across all conditions
    ks = empty((3, N_conds, N_conds), 'd')
    for c1 in xrange(N_conds):
        for c2 in xrange(c1, N_conds):
            ks[0,c1,c2] = ks[0,c2,c1] = \
                ks_2samp(results[c1]['remapping_samples'], results[c2]['remapping_samples'])[1]
            ks[1,c1,c2] = ks[1,c2,c1] = \
                ks_2samp(results[c1]['rate_remapping_samples'], results[c2]['rate_remapping_samples'])[1]
            ks[2,c1,c2] = ks[2,c2,c1] = \
                ks_2samp(results[c1]['turnover_samples'], results[c2]['turnover_samples'])[1]

    # Save the data
    data_fd = file('summary.pickle', 'w')
    data_dict = dict(   data=data, 
                        labels=array(labels), 
                        visible=array([True]*N_conds), 
                        color=array(['k']*N_conds),
                        pvals=ks
                        )
    cPickle.dump(data_dict, data_fd)
    data_fd.close()
    
    os.chdir(old_dir)
    return data_dict


def summary_plot(summary, data_x=0, data_y=1, annotate=True):
    """Create a scatter ellipse plot of means and confidence intervals for 
    summary data
    
    Keyword arguments:
    data_x/y -- which data to display on the x- and y- axes; 0=remapping,
        1=rate_remapping, and 2=turnover.
    """
    from pylab import figure, axes, rcParams
    from matplotlib.patches import Ellipse
    
    # Plot parameters
    fig_size = 6.5, 7.5
    anno_x = 150 # pts to right of left axis border
    anno_y = -50 # pts above the bottom axis border
    anno_dy = 25 # annotate row height
    padding = 1.5 # factor of ellipse radius
    plot_data_x = data_x
    plot_data_y = data_y
    ellipse_kw = dict(  alpha=0.75, 
                        aa=True, 
                        lw=2, 
                        fill=True, 
                        fc='0.6'
                        )
    
    # Create the figure and axes
    old_fig_size = rcParams['figure.figsize']
    rcParams['figure.figsize'] = fig_size
    f = figure()
    f.set_size_inches(fig_size)
    f.suptitle('Remapping Sample Summary Data')
    ellipse_kw['axes'] = ax = axes()
    ellipse_kw['clip_box'] = ax.bbox
    
    # Draw ellipses and annotations to the axes
    box = [1]*4
    c = 0
    for i,row in enumerate(summary['data']):
        if not summary['visible'][i]:
            continue
        xy = row[plot_data_x], row[plot_data_y]
        ellipse_kw['ec'] = summary['color'][i]
        ax.add_artist(
            Ellipse(xy, 2*row[plot_data_x+3], 2*row[plot_data_y+3], **ellipse_kw))
        if annotate:
            ax.annotate(
                summary['labels'][i], xy, 
                xytext=(anno_x, anno_y-c*anno_dy),
                textcoords='axes points',
                arrowprops=dict(width=0.5, frac=0.0, headwidth=0.0, shrink=0.0))
        if xy[0] - padding*row[plot_data_x+3] < box[0] or c == 0:
            box[0] = xy[0] - padding*row[plot_data_x+3]
        if xy[0] + padding*row[plot_data_x+3] > box[2] or c == 0:
            box[2] = xy[0] + padding*row[plot_data_x+3]
        if xy[1] - padding*row[plot_data_y+3] < box[1] or c == 0:
            box[1] = xy[1] - padding*row[plot_data_y+3]
        if xy[1] + padding*row[plot_data_y+3] > box[3] or c == 0:
            box[3] = xy[1] + padding*row[plot_data_y+3]
        c += 1
            
    # Set 1.0 lines and axis extent
    ax.hlines(1.0, xmin=box[0], xmax=box[2], linestyle=':', color='k')
    ax.vlines(1.0, ymin=box[1], ymax=box[3], linestyle=':', color='k')
    ax.set_xlim(box[0], box[2])
    ax.set_ylim(box[1], box[3])
    
    # Axis labels
    if plot_data_x == 0:
        ax.set_xlabel('Place Remapping')
    elif plot_data_x == 1:
        ax.set_xlabel('Rate Remapping')
    elif plot_data_x == 2:
        ax.set_xlabel('Turnover')
    if plot_data_y == 0:
        ax.set_ylabel('Place Remapping')
    elif plot_data_y == 1:
        ax.set_ylabel('Rate Remapping')
    elif plot_data_y == 2:
        ax.set_ylabel('Turnover')
    
    rcParams['figure.figsize'] = old_fig_size
    return f


def summary_bar_plot(summary, sig_ref=None, sig_thresh='>0.05', ymin=0.5, 
    show_data=[0,1,2], legend=False):
    """Bar plot of summary statistics with SEM errorbars
    
    Arguments:
    summary -- summary results object
    ymin -- y-axis minimum
    show_data -- which data to show: list of integers, where 0=remapping,
        1=rate_remapping, and 2=turnover
    sig_ref -- index number of condition that will be reference for determining
        statistical significance
    sig_thresh -- string indicating threshold relationship for placing a 
        significance mark (star/asterisk): either '<' or '>' followed by float
        threshold value
    legend -- whether to display a legend
    """
    from pylab import figure, axes, rcParams
    from matplotlib import cm
    
    # Plot parameters
    fig_size = 9, 6
    show_data = array(show_data)
    num_data = len(show_data)
    bar_w = 1/float(num_data+1)
    
    # Significance parameters
    sig_less = True
    if sig_thresh[0] == '>':
        sig_less = False
    sig_value = float(sig_thresh[1:])
    if sig_ref is not None:
        print '* Significance indicates p%s%.5f.'%(sig_thresh[0], sig_value)

    # Create the figure and axes
    old_fig_size = rcParams['figure.figsize']
    rcParams['figure.figsize'] = fig_size
    f = figure()
    f.set_size_inches(fig_size)
    f.suptitle('Remapping Sample Summary Data')
    ax = axes()
    
    # Create the bar data
    left = []
    height = []
    yerr = []
    xticklabels = []
    sig_x = []
    sig_y = []
    c = 0
    for i,row in enumerate(summary['data']):
        if not summary['visible'][i]:
            continue
        if num_data == 1:
            left.extend([c-0.5*bar_w])
        elif num_data == 2:
            left.extend([c-bar_w, c])
        elif num_data == 3:
            left.extend([c-1.5*bar_w, c-0.5*bar_w, c+0.5*bar_w])
        height.extend(list(row[show_data]))
        yerr.extend(list(row[3+show_data]/1.96))
        xticklabels.append(summary['labels'][i])
        
        if sig_ref is not None and sig_ref <> i:
            for n in xrange(num_data):
                sig_mark = False
                if sig_less:
                    if summary['pvals'][show_data[n], sig_ref, i] < sig_value:
                        sig_mark = True
                else:
                    if summary['pvals'][show_data[n], sig_ref, i] >= sig_value:
                        sig_mark = True
                if sig_mark:
                    sig_x.append(left[n-num_data]+0.5*bar_w)
                    sig_y.append(height[n-num_data]+yerr[n-num_data]+0.02)
        
        c += 1
    
    # Create the bar chart and legend
    bar_cols = cm.gray(([0.25, 0.6, 0.8][:num_data])*c)
    bar_h = ax.bar(left, height, width=bar_w, yerr=yerr, 
        ec='k', color=bar_cols, linewidth=0, ecolor='k', capsize=3, aa=False)
    if legend:
        labels = ['Remapping', 'Rate Remapping', 'Turnover']
        data_labels = [labels[i] for i in show_data]
        ax.legend(bar_h[:num_data], data_labels)
    if sig_ref is not None:
        ax.plot(sig_x, sig_y, 'k*', ms=8)
    ax.hlines(1.0, xmin=-0.5, xmax=c-0.5, linestyle=':', color='k')
    ax.set_xlim(-0.5, c-0.5)
    ax.set_ylim(ymin, 1.1)
    ax.set_xticks(arange(c))
    ax.set_xticklabels(xticklabels)
    
    rcParams['figure.figsize'] = old_fig_size
    return f


if __name__ == '__main__':
    main()
