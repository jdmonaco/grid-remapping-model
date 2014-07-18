#encoding: utf-8
"""
tools.colormaps -- Custom colormap definitions

Exported namespace: diffmap

Written by Joe Monaco
Center for Theoretical Neuroscience
Copyright (c) 2007-2008 Columbia Unversity. All Rights Reserved.  
"""

from matplotlib.colors import LinearSegmentedColormap as LSC


def diffmap(use_black=False):
    """Conventional differencing map with graded red and blue for values less 
    than and greater than, respectively, the mean of the data. Values approaching 
    the mean are increasingly whitened, and the mean value is white.
    
    Keyword arguments:
    use_black -- if True, the mean value is black instead of white by default
    """
    m = int(not use_black)
    segmentdata = { 'red':   [(0, 1, 1), (.5, m, m), (1, 0, 0)],
                    'green': [(0, 0, 0), (.5, m, m), (1, 0, 0)],
                    'blue':  [(0, 0, 0), (.5, m, m), (1, 1, 1)] }
    return LSC('RdWhBu', segmentdata)   
