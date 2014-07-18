# encoding: utf-8
"""
bash.py -- Functions relating to bash scripting and output

Created by Joe Monaco on 2007-11-14.
Copyright (c) 2007 Columbia University. All rights reserved.
"""

# Library imports
from enthought.traits.api import HasTraits, Instance, Trait, String
from sys import stdout as _out, stderr as _err, platform as _plat
import time as _time

# Bash color codes
COLOR_CODE = {
    'WHI':"\033[1;37m", 
    'LGY':"\033[0;37m", 
    'GRY':"\033[1;30m", 
    'BLK':"\033[0;30m", 
    'RED':"\033[0;31m", 
    'LRD':"\033[1;31m", 
    'GRN':"\033[0;32m", 
    'LGN':"\033[1;32m", 
    'BRN':"\033[0;33m", 
    'YLW':"\033[1;33m", 
    'BLU':"\033[0;34m", 
    'LBL':"\033[1;34m", 
    'PUR':"\033[0;35m", 
    'PNK':"\033[1;35m", 
    'CYN':"\033[0;36m", 
    'LCN':"\033[1;36m", 
    'UNC':"\033[0m"
}

# Get bash-colored text
def _color_format(c, text):
    try:
        code = COLOR_CODE[c]
    except KeyError:
        code = COLOR_CODE['BLU']
    return code + text + COLOR_CODE['UNC']

# Individual coloring functions
def white(s):       return _color_format('WHI', s)
def lightgray(s):   return _color_format('LGY', s)
def gray(s):        return _color_format('GRY', s)
def black(s):       return _color_format('BLK', s)
def red(s):         return _color_format('RED', s)
def lightred(s):    return _color_format('LRD', s)
def green(s):       return _color_format('GRN', s)
def lightgreen(s):  return _color_format('LGN', s)
def brown(s):       return _color_format('BRN', s)
def yellow(s):      return _color_format('YLW', s)
def blue(s):        return _color_format('BLU', s)
def lightblue(s):   return _color_format('LBL', s)
def purple(s):      return _color_format('PUR', s)
def pink(s):        return _color_format('PNK', s)
def cyan(s):        return _color_format('CYN', s)
def lightcyan(s):   return _color_format('LCN', s)

# Mapping between color names and functions
COL_FUNC = {
    'white': white,
    'lightgray': lightgray,
    'gray': gray,
    'black': black,
    'red': red,
    'lightred': lightred,
    'green': green,
    'lightgreen': lightgreen,
    'brown': brown,
    'yellow': yellow,
    'blue': blue,
    'lightblue': lightblue,
    'purple': purple,
    'pink': pink,
    'cyan': cyan,
    'lightcyan': lightcyan,
}
COLORS = COL_FUNC.keys()
COLORS.sort()

class CPrint(HasTraits):
    
    """
    A simple callable printing object
    
    Constructor arguments:
    prefix -- set the default prefix string for all messages
    outfd -- optional file descriptor for timestamped log output
    color -- optional color specification for the prefix (default 'cyan')
    
    Call parameters:
    s -- required string message to output
    prefix -- optionally override default prefix string
    error -- specify whether this is an error message
    """
    
    prefix = String('Message')
    outfd = Instance(file)
    color = Trait('cyan', COLORS)
    
    def __call__(self, s, prefix=None, error=False):

        """
        Pretty message output: colored and multi-line indented
        """
        
        first_cap = lambda s: s[0].upper() + s[1:]
        
        # Set the prefix
        if prefix is None:
            pre = first_cap(self.prefix)
        else:
            pre = first_cap(str(prefix))
        pre += ': '
        
        # Message colors and error output
        if _plat == 'win32':
            prec = msgc = str
        else:
            prec = COL_FUNC[self.color]
            msgc = error and lightred or lightgray
        console = error and _err or _out
        
        # Console print with indentation
        first = True
        pre_len = len(pre)
        for l in s.split('\n'):
            if first:
                console.write(prec(pre) + msgc(first_cap(l.strip())) + '\n')
                first = False
            elif l:
                console.write(' '*pre_len + msgc(l.strip()) + '\n')
        console.flush()
        
        # Timestamped file print if available
        if self.outfd and not self.outfd.closed:
            self.outfd.write('[%s]\n'%_time.strftime('%H:%M:%S on %m/%d/%y'))
            if error:
                self.outfd.write('!---> ERROR MESSAGE <---!\n\a')
            self.outfd.write(first_cap(s.strip()) + '\n\n')
            self.outfd.flush()
    
    def printf(self, s, color='lightgreen'):
        if _plat == 'win32' or color not in COLORS:
            _out.write(s)
        else:
            _out.write(COL_FUNC[color](s))
        _out.flush()
