#encoding: utf-8
"""
messager.py -- Messager class handles console messaging for simulations

Copyright (c) 2011 Johns Hopkins University. All rights reserved.

This software is provided AS IS under the terms of the Open Source MIT License. 
See http://www.opensource.org/licenses/mit-license.php.
"""

import time
from sys import stdout, stderr

DEFAULT_NOTIFICATION = 'info'
DEFAULT_TITLE = "Message"


class Messager(object):
    
    """
    Simple time-stamped message output via console
    """
        
    def __init__(self, title=DEFAULT_TITLE):
        self.default_title = title
        self.console_active = True
    
    def __call__(self, *args, **kwargs):
        self.notify(*args, **kwargs)
    
    def notify(self, msg, title=None, notification=None):
        """Send time-stamped pretty format message out to console
        """        
        # Handle defaults
        if not title:
            title = self.default_title
        if not notification:
            notification = DEFAULT_NOTIFICATION
            
        # Format string outputs
        nprefix = notification.upper() + ': '
        title = str(title).title()
        msg = str(msg)
        msg = msg[0].upper() + msg[1:]
        
        # Pretty echo out to the console
        if notification == 'error':
            console = stderr.write
        else:
            console = stdout.write
        header = title.center(58, '-') + time.strftime('[%H:%M:%S]') + '--'
        console('%s\n%s %s\n'%(header, nprefix, msg))
