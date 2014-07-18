# encoding: utf-8
"""
Core simulation package -- Functionality for scientific models and simulations

Author: Joe Monaco, jmonaco@jhu.edu

Copyright (c) 2009-2011 Johns Hopkins University. All rights reserved.

This software is provided AS IS under the terms of the Open Source MIT License. 
See http://www.opensource.org/licenses/mit-license.php.
"""

import os as _os
import sys as _sys


# Time formatting

TS_FMT = "%A %b %d %Y [%I:%M:%S%p]"


# Determine and create (if necessary) analysis output directorgrid_modely

ANA_SUBDIR = 'grid_model'
RC_FILE = '.grid_model'

if _sys.platform == 'win32':
    HOME_DIR = _os.path.join(_os.getenv('SystemDrive'), _os.getenv('HOMEPATH'))
else:
    HOME_DIR = _os.getenv('HOME')

ANA_DIR = HOME_DIR
found = False
rc_path = _os.path.join(HOME_DIR, RC_FILE)
if _os.path.exists(rc_path):
    try:
        fd = file(rc_path, 'r')
        ANA_DIR = _os.path.abspath(fd.readline().strip())
    except IOError:
        pass
    else:
        if _os.path.isdir(ANA_DIR):
            found = True
        elif (raw_input('Output directory not found: %s\nRecreate? (y/n) '%ANA_DIR).lower()=='y'):
                _os.makedirs(ANA_DIR)
                found = True
    finally:
        fd.close()
if not found:
    ANA_DIR = _os.path.join(HOME_DIR, ANA_SUBDIR)
    if not _os.path.isdir(ANA_DIR):
        if (raw_input('Create default output directory? (%s)\n(y/n) '%ANA_DIR).lower()=='y'):
            _os.makedirs(ANA_DIR)
        else:
            _tempdir = raw_input('Specify output directory [then hit Enter]:\n')
            ANA_DIR = _os.path.abspath(_os.path.realpath(_tempdir))
            if not _os.path.isdir(ANA_DIR):
                _os.makedirs(ANA_DIR)
    fd = file(rc_path, 'w')
    fd.write(ANA_DIR + '\n')
    fd.close()
    _sys.stdout.write('Analysis output directory:\n\tSettings file: %s\n\tDirectory: %s\n'%(rc_path, ANA_DIR))


# Message strings used in Model class

MODEL_MESSAGES = {
	'TRIALSTOP': ('Enter \'stop\' to end the program. Anything else will cancel and rerun the current trial.', 
		'Trial Stopped', 'info'),
	'UNHEXC': ('Trial %d at t=%.2f seconds\n\n%s: %s', 'Unhandled exception', 'error'),
	'TRIALRST': ('Resetting to trial %d due to missing trial data!', 'Missing data', 'monitor'),
	'SAVEFAIL': ('KeyError: %s\n\nRerunning trial...', 'Save failure', 'error'),
	'DPATH': ('Failed to create data directory: %s', '', 'error'),
	'DPFILE': ('Regular file exists and is not a directory: %s', 'Path error', 'error'),
	'DPFAIL': ('Failed to create data directory: %s', 'Path error', 'error'),
	'DPCREATE': ('Data directory created: %s', 'Data path', 'info')
}
