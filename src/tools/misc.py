#encoding: utf-8
"""
tools.misc -- Miscellaneous toolbox functions

Exported namespace: fit_exp_linear, halfwave, set_max_recursion, Null

Written by Joe Monaco
Center for Theoretical Neuroscience
Copyright (c) 2007-2008 Columbia Unversity. All Rights Reserved.  
"""

import numpy as _N
import sys


def set_max_recursion():
    """Set platform-dependent maximum recursion depth
    
    NOTE: These values were obtained using the find_recursionlimit.py
    script that comes with python.
    """
    if sys.platform == 'darwin':
        sys.setrecursionlimit(4400)
    elif sys.platform == 'linux2':
        sys.setrecursionlimit(6700)
    return        
    
class Null(object):
    
    """
    Null object design pattern
    
    Python Cookbook, Second Edition: Recipe 6.17
    """
    
    def __new__(cls, *p, **kw):
        "force only one instance"
        if '_inst' not in vars(cls):
            cls._inst = object.__new__(cls, *p, **kw)
        return cls._inst
    
    def __init__(self, *p, **kw): pass
    def __call__(self, *p, **kw): return self
    def __str__(self): return "Null()"
    def __repr__(self): return "Null()"
    def __nonzero__(self): return False
    def __getattr__(self, name): return self
    def __delattr__(self, name): return self
    def __setattr__(self, name): return self
    def __getitem__(self, i): return self
    def __setitem__(self, *p): pass
