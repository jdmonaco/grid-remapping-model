#encoding: utf-8
"""
tools.path -- Toolbox functions for creating or handling paths

Exported namespace: unique_path

Written by Joe Monaco
Center for Theoretical Neuroscience
Copyright (c) 2007-2008 Columbia Unversity. All Rights Reserved.  
"""

import os.path as _path


def unique_path(stem, ext="", fmt="%s%02d", reverse_fmt=False):
    """Insert a unique identifier into a file or directory path
    
    Required argument:
    stem -- the full path up until a unique identifier is needed
    
    Optional keyword arguments:
    ext -- this extension will be added to the path (default none)
    fmt -- format specification with one string (%s) and one integer (%d) for
        the stem and unique identifier; an extension may be included, but then 
        **ext** should not be specified (default '%s%02d')
    reverse_fmt -- specify whether the stem string and the unique id integer 
        are reversed in the *fmt* specification (i.e., %d before %s)
        
    Returns a modified path based on **stem**, a unique identifier and **ext**, 
    as formatted by **fmt**. 
    """
    if ext:
        if ext[0] != '.':
            ext = '.' + ext
    if reverse_fmt:
        head, tail = _path.split(stem)
        filename = lambda i: _path.join(head, fmt%(i, tail)) + ext
    else:
        filename = lambda i: fmt%(stem, i) + ext
    i = 0
    while _path.exists(filename(i)):
        i += 1
    return filename(i)
