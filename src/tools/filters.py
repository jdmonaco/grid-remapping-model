# encoding: utf-8
"""
filters.py -- Custom signal filtering functions

Exported namespace: halfwave, quick_boxcar, circular_blur

Created by Joe Monaco on 2007-11-15. Updated 2009-09-11.
Copyright (c) 2007-2009 Columbia University. All rights reserved.
Copyright (c) 2009-2011 Johns Hopkins University. All rights reserved.

This software is provided AS IS under the terms of the Open Source MIT License. 
See http://www.opensource.org/licenses/mit-license.php.
"""

from numpy import r_, empty, zeros, ceil, trapz, ndarray
from scipy.signal import gaussian, convolve


def halfwave(x, copy=False):
    """Half-wave rectifier for arrays or scalars
    
    NOTE: Specify copy=True if array data should be copied before performing
    halwave rectification.
    """
    if type(x) is ndarray and x.ndim:
        if copy:
            x = x.copy()
        x[x<0.0] = 0.0
    else:
        x = float(x)
        if x < 0.0:
            x = 0.0
    return x

def quick_boxcar(s, M=4, centered=False):
    """Returns a boxcar-filtered version of the input signal
    
    Keyword arguments:
    M -- number of averaged samples (default 4)
    centered -- recenter the filtered signal to reduce lag (default False)
    """
    # Sanity check on signal and filter window
    length = s.shape[0]
    if length <= 2*M:
        raise ValueError, 'signal too short for specified filter window'
    
    # Set up staggered arrays for vectorized average
    z = empty((M, length+M-1), 'd')
    for i in xrange(M):
        z[i] = r_[zeros(i)+s[0], s, zeros(M-i-1)+s[-1]]
    
    # Center the average if specified
    start_ix = 0
    end_ix = length
    if centered:
        start_ix += int(M/2)
        end_ix += int(M/2)
    
    return z.mean(axis=0)[start_ix:end_ix]
    
def circular_blur(s, blur_width):
    """Return a wrapped gaussian smoothed (blur_width in degrees) signal for 
    data binned on a full circle range [0, 2PI/360).
    """
    bins = s.shape[0]
    width = blur_width / (360.0/bins)
    size = ceil(8*width)
    if size > bins:
        size = bins
    wrapped = r_[s[-size:], s, s[:size]]
    G = gaussian(size, width)
    G /= trapz(G)
    S = convolve(wrapped, G, mode='same')
    return S[size:-size]
