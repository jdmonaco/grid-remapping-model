#encoding: utf-8
"""
tools.images -- Toolbox functions for creating image output

Exported namespace: image_blast, array_to_rgba, array_to_image, diffmap

Written by Joe Monaco
Center for Theoretical Neuroscience
Copyright (c) 2007-2008 Columbia Unversity. All Rights Reserved.  
"""

import os as _os
import numpy as _N
from sys import platform as _plat

if _plat == "win32":
    import Image
else:
    from PIL import Image


def image_blast(M, savedir, stem='image', fmt='%s_%03d', rev=False, **kwargs):
    """Save a rank-3 stacked intensity matrix *M* to a set of individual PNG 
    image files in the directory *savedir*.
    
    If *savedir* does not exist it will be created. Set **stem** to specify 
    the filename suffix.
    
    Keyword arguments:
    stem -- file name stem to be used for output images
    fmt -- a unique_path fmt specification (need an %s followed by a %d)
    rev -- indicate use of a reversed fmt specification (%d followed by a %s)
    
    Extra keyword arguments will get passed through to array_to_rgba. See its 
    doc string for details.
    """
    assert M.ndim == 3, 'requires rank-3 array of intensity values'
    d = _os.path.realpath(str(savedir))
    if not _os.path.exists(d):
        _os.makedirs(d)
    stem = _os.path.join(d, stem)
    N = M.shape[0]
    first, middle, last = "", "", ""
    for i,m in enumerate(M):
        image_fn = unique_path(stem, fmt=fmt, ext="png", reverse_fmt=rev)
        if i == 0:
            first = image_fn
        elif i == N-1:
            last = image_fn
        array_to_image(m, image_fn, **kwargs)
    if N == 2:
        middle += '\n'
    elif N > 2:
        middle += '\n\t...\n'
    print first, middle, last
    return    

def array_to_rgba(mat, cmap=None, norm=True, cmin=0, cmax=1):
    """Intensity matrix (float64) -> RGBA colormapped matrix (uint8)
    
    Keyword arguments:
    cmap -- a matplotlib.cm colormap object
    norm -- whether the color range is normalized to values in M
    
    If *norm* is set to False:
    cmin -- minimum clipping bound of the color range (default 0)
    cmax -- maximum clipping bound of the color range (default 1)
    """
    if cmap is None:
        from matplotlib import cm
        cmap = cm.hot
    M = mat.copy()
    data_min, data_max = M.min(), M.max()
    if norm:
        cmin, cmax = data_min, data_max
    else:
        if cmin > data_min:
            M[M < cmin] = cmin # clip lower bound
        if cmax < data_max:
            M[M > cmax] = cmax # clip uppder bound
    return cmap((M-cmin)/float(cmax-cmin), bytes=True)

def array_to_image(M, filename, **kwargs):
    """Save matrix, autoscaled, to image file (use PIL fmts)
    
    Keyword arguments are passed to array_to_rgba.
    """
    if M.ndim != 2:
        raise ValueError, 'requires rank-2 matrix'
    img = Image.fromarray(array_to_rgba(M, **kwargs), 'RGBA')
    img.save(filename)
    return

def tile2D(M, mask=None, gridvalue=0.5, shape=None):
    """
    Construct a tiled 2D matrix from a 3D matrix
    
    Keyword arguments:
    mask -- an (H,W)-shaped binary masking array for each cell
    gridvalue -- the intensity value for the grid
    shape -- a (rows, columns) tuple specifying the shape of the tiling to use
    
    If shape is specified, rows+columns should equal M.shape[0].
    """
    if len(M.shape) != 3: 
        return
    N, H, W = M.shape
    if mask is not None and (H,W) != mask.shape:
        mask = None
    if shape and (type(shape) is type(()) and len(shape) == 2):
        rows, cols = shape
    else:
        rows, cols = tiling_dims(N)
    Mtiled = _N.zeros((rows*H, cols*W), 'd')
    for i in xrange(N):
        r, c = int(i/cols), _N.fmod(i, cols)
        if mask is None:
            Mtiled[r*H:(r+1)*H, c*W:(c+1)*W] = M[i]
        else:
            Mtiled[r*H:(r+1)*H, c*W:(c+1)*W] = mask * M[i]
    Mtiled[H::H,:] = gridvalue
    Mtiled[:,W::W] = gridvalue
    return Mtiled

def tiling_dims(N):
    """Square-ish (rows, columns) for tiling N things
    """
    d = _N.ceil(_N.sqrt(N))
    return int(_N.ceil(N / d)), int(d)
