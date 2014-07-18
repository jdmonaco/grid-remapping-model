#encoding: utf-8
"""
grid.analysis.map_funcs -- A repository for standalone functions useful for 
    performing various computations on spatial map arrays.

Exported namespace: remap_quiver_plot, spatial_corr, linearize_spatial_corr,
    peak_vs_neighbors, input_vs_output_norms, haff_vs_r_peaks, 
    field_comparison_matrix, linear_rate_corr_matrix
    
Written by Joe Monaco, 10/21/2008
Copyright (c) 2008 Columbia Unversity. All Rights Reserved.  
"""

# Library imports
import scipy.signal, numpy
from scipy.stats import pearsonr


def get_tuned_weights(pmap, W, EC, alpha=0.5, grow_synapses=False):
    """
    Perform afferent tuning on the weight matrix and return new weights
    
    Required parameters:
    pmap -- a PlaceMap object resulting from the spatial map simulation
    W0 -- the afferent weight matrix used in the simulation
    EC -- the GridCollection instance used as input in the simulation
    
    Keyword arguments:
    alpha -- 0.0 to 1.0 value of how much tuning to (default 0.5)
    """
    norm = numpy.sqrt((W[0]**2).sum(axis=0))
    W0 = W / norm
    W1 = numpy.empty((pmap.num_maps, EC.num_maps), 'd')
    for i in xrange(pmap.num_maps):
        W1[i] = numpy.tanh(3*(pmap.maxima[i,2]-0.5)) * \
                    EC.map_value(pmap.maxima[i,0], pmap.maxima[i,1])
        if not grow_synapses:
            W1[i] *= W0[i] > 0.0
        W1[i] /= numpy.sqrt((W1[i]**2).sum(axis=0)) # normalize
    W2 = (1.0-alpha)*W0 + alpha*W1 # mixed old and tuned matrices
    for i in xrange(pmap.num_maps):
        W2[i] *= norm / numpy.sqrt((W2[i]**2).sum(axis=0)) # hetersynaptic LTD
    return W2

def remap_quiver_plot(cmp_AB, ax=None, rate_colors=False, 
    border_style=True, arrow_width=None, **kwextra):
    """
    Draw a remapping quiver plot for spatial map comparison data
    
    Requires a compare_AB dictionary as first argument.
    
    Keyword arguments:
    ax -- if specified, quiver plot is drawn to the given axes, otherwise
        a new figure and axes are created
    rate_colors -- whether to color the arrows based on rate remapping
    border_style -- if *rate_colors* is True, whether to use a black-bordered 
        arrow or not (if so, the Reds colormap is used; otherwise, a RKB 
        diffmap is used)
    
    Additional keywords are passed to the quiver call.
    """
    from matplotlib.pyplot import figure, axes, draw
    if ax is None:
        f = figure()
        ax = axes()
    
    # Set vector components for drawing arrows
    X, Y = cmp_AB['A_xy']
    U, V = cmp_AB['B_xy'] - cmp_AB['A_xy']
    args = (X, Y, U, V)
    
    # Calculate rate remapping vector for colors: (max-min)/max
    if rate_colors:
        C = cmp_AB['R_AB']
        args += (C,)
    
    # Set keyword arguments to format the quiver field
    if arrow_width is None:
        set_width = 0.5                 # set width here
    else:
        set_width = arrow_width
    kwargs = {  'units':'x',            # scale based on data range
                'scale':1,              # data per arrow unit
                'width':set_width,      # arrow units
                'headwidth':4,          # width units
                'headlength':5,         # width units
                'headaxislength':4,     # width units
                'minshaft':1,           # headlength units, scaling threshold
                'minlength':2.5/set_width }  # width units, display threshold
    if rate_colors:
        color_lims = numpy.array([0.0, 1.0])
        if border_style:
            from matplotlib import cm
            kwargs.update({
                    'cmap':cm.Reds,         # colormap for arrows
                    'clim':color_lims,      # colors on a (0,1) scale
                    'edgecolor':'k',        # arrow outline color
                    'lw':0.5 })             # arrow outline line-width                
        else:
            from ..tools.colormaps import diffmap
            kwargs.update({
                    'headwidth':4.0,        # scale up head with no borders
                    'headlength':5.0,       # 
                    'headaxislength':3.8,   #
                    'cmap':diffmap(use_black=True),
                    'clim':color_lims,      # colors on a (0,1) scale
                    'lw':0.0 })             # arrow outline line-width                
    kwargs.update(kwextra)
    
    # Execute the quiver command and draw the plot
    ax.cla()
    q = ax.quiver(*args, **kwargs)
    ax.axis('image')
    ax.axis([0, 100, 0, 100])
    draw()
    return q

def scatter_linreg_plot(x, y, ax=None, label='data', d_fmt='b.', l_fmt='k-',
    d_kw={}, l_kw={}):
    """Draw a scatter plot with linear regression fit line
    
    Keyword arguments:
    ax -- if specified, scatter plot is drawn to the given axes, otherwise
        a new figure and axes are created
    label -- label for this scatter data if a legend is created
    d_fmt/l_fmt -- format specifier for data and line, respectively
    d_kw/l_kw -- additional keyword dictionaries for the plot calls
    
    Prints Pearson r and corresponding p-value to console.
    
    Returns the Pearson r coefficient.
    """
    assert len(x) == len(y), 'scatter data must be same length'
    from matplotlib.pyplot import figure, axes, draw
    from scipy.stats import linregress
    if ax is None:
        f = figure()
        ax = axes()
    
    # Draw the scatter data
    ax.plot(x, y, d_fmt, label=label, **d_kw)
    
    # Get the linear regression
    m, b, r, p, sem = linregress(x, y)
    print '(r = %.4f, p = %.4e)' % (r, p)
    x0 = numpy.array([x.min(), x.max()], 'd')
    y0 = m * x0 + b
    
    # Plot the regression line
    ax.plot(x0, y0, l_fmt, zorder=-1, label='_nolegend_', **l_kw)
    draw()
    return r

def spatial_corr(*args):
    """2D spatial autocorrelation of rank-3 population arrays
    
    Pass in a single z-stack [(num_maps, H, W)-shaped rank-3] array to compute
    and return its spatial autocorrelogram.
    
    Pass in two z-stack maps (e.g., A and B) to compute the cross-correlogram
    of A with respect to B. 
    """
    # Handle input arguments
    if len(args) == 0 or len(args) > 2:
        raise ValueError, 'requires one or two arguments'
    if len(args) == 1:
        A = B = args[0]
    else:
        A, B = args    
    assert A.shape == B.shape, 'shape mismatch between input arrays'
    assert A.ndim == 3, 'input arrays must be rank-3'
    
    # Map and correlogram dimensions
    num_maps, H, W = A.shape  
    corr_shape = 2*H-1, 2*W-1 
    
    # Fourier transforms
    A_ = scipy.signal.fft2(A, shape=corr_shape)
    B_ = scipy.signal.fft2(B[:, ::-1, ::-1], shape=corr_shape)
    AB_conv = (A_ * B_).sum(axis=0)
    
    return scipy.signal.real(scipy.signal.ifft2(AB_conv))/num_maps

def linearize_spatial_corr(Mcorr):
    """Perform a radial collapse of a 2D spatial autocorrelation to get a 
    linear autocorrelation
    
    NOTE: This should not be used for cross-correlations!
    
    Mcorr should be a 199x199 autocorrelogram of a 100x100 map.
    
    Returns 2-row autocorrelation (lag, corr) array.
    """
    assert type(Mcorr) is numpy.ndarray, 'bad type for matrix argument'
    assert Mcorr.shape == (199,199), 'invalid shape for autocorrelation matrix'
    
    # Scan the correlogram and compute the radius from the midpoint
    n = numpy.zeros(101, 'h')
    c = numpy.zeros(101, 'd')
    mid_x, mid_y = 99.5, 99.5
    for i in xrange(199):
        for j in xrange(199):
            d = numpy.sqrt((mid_y - i)**2 + (mid_x - j)**2)
            if d > 100:
                d = 100
            n[int(d)] += 1
            c[int(d)] += Mcorr[i, j]
    c /= n # get the sample means
    
    # Create the return array: reflect 0->Max correlations to -Max->0 
    Lcorr = numpy.zeros((2,201), 'd')
    Lcorr[0] = numpy.arange(-100, 101)
    Lcorr[1] = numpy.r_[c[::-1], c[1:]]
    
    return Lcorr

def peak_vs_neighbors(pmap, k=4, median_dist=True, use_primary=False):
    """Compute scatter data for looking at the relationship between field peak 
    rates and a measure of nearest neighbor distance.
    
    A PlaceMap (or subclass) instance must be passed in.

    Keyword arguments:
    k -- number of nearest neighbors to factor into the measure
    median_dist -- use the median neighbor distance (if False, the maximum
        distance of the k-neighbors is used)
    use_primary -- only use primary place fields (most active field per unit)
    
    Returns 2-row concatenation of peaks and neighbor-distance arrays.
    """
    # Get field centroids and peak rates from the spatial map    
    if use_primary:
        udata = pmap.get_unit_data()
        x, y, peaks = udata['max_x'], udata['max_y'], udata['max_r']
        nfields = len(udata)
    else:
        fdata = pmap.get_field_data()
        x, y, peaks = fdata['x'], fdata['y'], fdata['peak']
        nfields = len(fdata)
    
    # Main loop through place fields
    neighbor_dists = numpy.empty(nfields, 'd')
    for f in xrange(nfields):
        d = numpy.sqrt((x - x[f])**2 + (y - y[f])**2)
        nearest_k = numpy.argsort(d)[1:k+1]
        if median_dist:
            neighbor_dists[f] = numpy.median(d[nearest_k])
        else:
            neighbor_dists[f] = d[nearest_k[-1]]
    
    return numpy.c_[peaks, neighbor_dists].T

def peaks_vs_area(pmap):
    """Get scatter data for field peak rates vs field area in cm^2
    
    A PlaceMap (or subclass) instance must be passed in.
    
    Returns 2-row (peak, area) array.
    """
    fdata = pmap.get_field_data()
    return numpy.c_[fdata['peak'], fdata['area']].T

def secondary_field_data(pmap):
    """Get scatter data for normalized secondary peak vs. distance from primary
    
    A PlaceMap (or subclass) instance must be passed in.
    
    Returns 2-row (primary normed rate, primary distance) array.
    """
    # Get place field data from spatial map
    fdata = pmap.get_field_data()
    units = numpy.unique(fdata['unit'])
    
    # Find dual fields and store data
    norm_peak = []
    primary_dist = []
    for u in units:
        ix = (fdata['unit'] == u).nonzero()[0]
        if len(ix) <= 1:
            continue
        fields = fdata[ix]
        sort_ix = numpy.argsort(fields['peak'])
        P = fields[sort_ix[-1]]
        for S in fields[sort_ix[:-1]]:
            norm_peak.append(S['peak']/P['peak'])
            primary_dist.append(
                numpy.sqrt((P['x']-S['x'])**2 + (P['y']-S['y'])**2))
    
    # Return array data
    return numpy.c_[numpy.array(norm_peak), numpy.array(primary_dist)].T

def input_vs_output_norms(EC, R):
    """Get scatter data for input and output population vector norms
    
    Arguments:
    EC -- the GridCollection cortical object used in the simulation
    R -- the PlaceMap object containing the output spatial map
    
    Returns 2-row (|EC|, |R|) scatter array.
    """
    return numpy.c_[numpy.sqrt((EC.Map * EC.Map).sum(axis=0)).flatten(), 
                    numpy.sqrt((R.Map * R.Map).sum(axis=0)).flatten()].T

def haff_vs_r_peaks(ca3, pmap=None):
    """Get input-output per-field scatter data to show effects of competition
    
    Arguments:
    ca3 -- PlaceNetwork model instance to run comparison data 
    pmap -- precomputed PlaceMap for ca3 [optional: if omitted, the spatial
        map is computed]
        
    Returns 2-row concatenation of field h_aff^i vs. r_i scatter points.
    """
    # Deprecate norm keyword
    if norm:
        import warnings
        warnings.warn('The \'norm\' keyword argument is deprecated.')
    
    # Compute the spatial map if not passed in
    if pmap is None:
        from ..ratemap import CheckeredRatemap
        pmap = CheckeredRatemap(ca3)
        pmap.compute_coverage()
        
    # Get field data
    fdata = pmap.get_field_data()
    x, y, peaks, unit = fdata['x'], fdata['y'], fdata['peak'], fdata['unit']
    nfields = len(fdata)
    
    # Main loop through place fields
    h_aff = numpy.empty(nfields, 'd')
    beta = ca3.gamma * ca3.beta_EC
    for f in xrange(nfields):
        r_EC = ca3.get_afferent_input(x[f], y[f])
        h_aff[f] = beta * numpy.dot(ca3.W[unit[f]], r_EC)
    
    return numpy.c_[h_aff, peaks].T    

def field_comparison_matrix(pmap, which='overlap'):
    """Get a matrix of pair-wise comparisons of single-max fields
    
    A PlaceMap (or subclass) instance must be passed in. The units are sorted
    by the quadrant into which their respective peaks fall. 
    
    Keyword arguments:
    which -- what sort of comparison to perform: 'overlap' (default) for
        a pixel overlap count, 'sim' for cosine similarity
    
    Returns a NxN matrix where N is the number of active place-units.
    """
    if which not in ('overlap', 'sim'):
        raise ValueError, 'invalid comparison type specified by which keyword'
    
    # Get indices of active units
    udata = pmap.get_unit_data()
    x, y, units = udata['max_x'], udata['max_y'], udata['unit']
    
    # Spatial sort of units based on peak location
    AND = numpy.logical_and
    mid_x, mid_y = 50.0, 50.0
    ll = AND(x<mid_x, y<mid_y).nonzero()[0]
    lr = AND(x>=mid_x, y<mid_y).nonzero()[0]
    ul = AND(x<mid_x, y>=mid_y).nonzero()[0]
    ur = AND(x>=mid_x, y>=mid_y).nonzero()[0]
    units = units[numpy.r_[ll, lr, ul, ur]]
    
    # Set up the matrix
    nfields = pmap.num_active
    M = numpy.empty((nfields,)*2, 'd')
    
    # Main loop pair-wise for fields
    if which is 'overlap':
        print 'Pixel overlap matrix...'
        for i in xrange(nfields):
            i_map = pmap.single_maps[units[i]].astype(bool)
            for j in xrange(i, nfields):
                j_map = pmap.single_maps[units[j]].astype(bool)
                M[i, j] = M[j, i] = (i_map * j_map).sum()
    elif which is 'sim':
        print 'Field vector cosine matrix...'
        for i in xrange(nfields):
            i_map = pmap.single_maps[units[i]].flatten()
            i_norm = numpy.sqrt(numpy.dot(i_map, i_map))
            for j in xrange(i, nfields):
                j_map = pmap.single_maps[units[j]].flatten()
                j_norm = numpy.sqrt(numpy.dot(j_map, j_map))
                M[i, j] = M[j, i] = numpy.dot(i_map, j_map) / (i_norm * j_norm)
    
    return M

def linear_rate_corr_matrix(R, which='corrcoef'):
    """Get a correlation matrix of the population rate vector for a line 
    scanned through the environment (a diagonal by default)
    
    Arguments:
    R -- the 3-index population rate matrix of responses
    which -- specify 'corrcoef' for Pearson correlations or 'sim' for cosine
        vector similarities
    
    Returns a NxN matrix where N is the number of pixel in a diagonal scan.
    """
    if which not in ('corrcoef', 'sim'):
        raise ValueError, 'invalid comparison type specified by which keyword'
    
    # Set us up the matrix
    npixels = 100
    M = numpy.empty((npixels,)*2, 'd')
    
    # Scan the diagonal from (0,0) to (100,100)
    if which is 'corrcoef':
        print 'Pearson correlation matrix...'
        for i in xrange(npixels):
            r_i = R[:,npixels-i-1, i]
            for j in xrange(npixels):
                r_j = R[:,npixels-j-1, j]
                r_corr = pearsonr(r_i, r_j)[0]
                if numpy.isnan(r_corr) or r_corr < 0:
                    M[i, j] = M[j, i] = 0.0
                else:
                    M[i, j] = M[j, i] = r_corr
    elif which is 'sim':
        print 'Cosine similarity matrix...'
        for i in xrange(npixels):
            r_i = R[:,npixels-i-1, i]
            r_i_norm = numpy.sqrt(numpy.dot(r_i, r_i))
            for j in xrange(npixels):
                r_j = R[:,npixels-j-1, j]
                r_j_norm = numpy.sqrt(numpy.dot(r_j, r_j))
                r_sim = numpy.dot(r_i, r_j) / (r_i_norm * r_j_norm)
                if numpy.isnan(r_sim):
                    M[i, j] = M[j, i] = 0.0
                else:
                    M[i, j] = M[j, i] = r_sim
    else:
        raise ValueError, 'invalid correlation measure specified: \'%s\''%which

    return M
