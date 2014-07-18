#encoding: utf-8
"""
tools.stats -- Toolbox functions for computing statistics 

Exported namespace: bootstrap, smooth_pdf, integer_hist

Written by Joe Monaco
Center for Theoretical Neuroscience
Copyright (c) 2007-2008 Columbia Unversity. All Rights Reserved.  
"""

import os as _os
import numpy as _N
from sys import platform as _plat, stderr as _err

def bootstrap(X, N, H0, *args):
    """Get a one-tail p-value for an algorithmic sampling process
    
    H0(*args) must return a scalar null sample value.
    
    The sign of the returned p-value indicates whether X is less than (-) or
    greater than (+) the median of the sample distribution. 
    
    Arguments:
    X -- the value for which to return a p-value
    N -- sampling size of the empirical distribution (beware O(n))
    H0 -- function that implements sampling process for the null result
    args -- additional arguments will be passed to H0
    """
    assert callable(H0), 'H0 must be a callable that returns a scalar sample'
    tail = 0
    for i in xrange(N):
        tail += int(H0(*args) >= X)
    if tail > float(N)/2:
        tail = tail - N # negative p-value for X less than median
    if not tail:
        _err.write('warning: bootstrap needs N > %d; returning upper bound\n'%N)
        tail = 1
    return tail / float(N)
    
def smooth_pdf(a, sd=None):
    """Get a smoothed pdf of an array of data for visualization
    
    Keyword arguments:
    sd -- S.D. of the gaussian kernel used to perform the smoothing (default is
        1/20 of the data range)
    
    Return 2-row (x, pdf(x)) smoothed probability density estimate.
    """
    from scipy.signal import gaussian, convolve
    from numpy import array, arange, cumsum, trapz, histogram, diff, r_, c_
    if sd is None:
        sd = 0.05 * a.ptp()
    data = a.copy().flatten() # get 1D copy of array data
    nbins = len(data) > 1000 and len(data) or 1000 # num bins >~ O(len(data))
    f, l = histogram(data, bins=nbins, normed=True) # fine pdf
    sd_bins = sd * (float(nbins) / a.ptp()) # convert sd to bin units
    kern_size = int(10*sd_bins) # sufficient convolution kernel size
    g = gaussian(kern_size, sd_bins) # generate smoothing kernel
    c = cumsum(f, dtype='d') # raw fine-grained cdf
    cext = r_[array((0,)*(2*kern_size), 'd'), c, 
        array((c[-1],)*(2*kern_size), 'd')] # wrap data to kill boundary effect
    cs = convolve(cext, g, mode='same') # smooth the extended cdf
    ps = diff(cs) # differentiate smooth cdf to get smooth pdf
    dl = l[1]-l[0] # get bin delta
    l = r_[arange(l[0]-kern_size*dl, l[0], dl), l, 
        arange(l[-1]+dl, l[-1]+kern_size*dl, dl)] # pad index to match bounds
    ps = ps[kern_size:kern_size+len(l)] # crop pdf to same length as index
    ps /= trapz(ps, x=l) # normalize pdf integral to unity
    return c_[l, ps].T # return 2-row concatenation of x and pdf(x)

def integer_hist(a, relative=False):
    """Get histogram data using integer bins
    
    Parameters:
    a -- the data to be histogrammed (ndim > 1 is flattened)
    relative -- whether count should be relative frequency or raw counts
    
    Returns (center, count):
    center -- integer bin centers for the histogram
    count -- bin frequencies, whether relative frequency or raw count        
    """
    data = a.round().flatten()
    center = _N.arange(int(a.min()), int(a.max())+1)
    if relative:
        count = _N.empty(center.shape[0], 'd')
    else:
        count = _N.empty(center.shape[0], 'l')
    for bin, c in enumerate(center):
        count[bin] = (data == c).sum()
    if relative:
        count /= count.sum()
    return center, count

