#encoding: utf-8
"""
tools.string -- Toolbox functions for handling strings

Exported namespace: snake2title, float_or_string, DumbStringIO

Written by Joe Monaco
Center for Theoretical Neuroscience
Copyright (c) 2007-2008 Columbia Unversity. All Rights Reserved.  
"""

from sys import stdout as _out


def snake2title(var_name):
    """Title-ize a variable name in snake_case
    """
    return ' '.join(var_name.split('_')).title()

def float_or_string(arg):
    """Force float or string
    """
    try:
        res = float(arg)
    except ValueError:
        res = arg
    return res

class DumbStringIO(object):
    
    """
    Dumb version of StringIO
    """
    
    def __init__(self, s=''):
        self._str = str(s)
        
    def __str__(self): return str(self._str)

    def write(self, new_str):
        self._str += str(new_str)
    
    def println():
        _out.write(self._str)
        if not _str.endswith('\n'):
            _out.write('\n')
        _out.flush()
