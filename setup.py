#!/usr/bin/env python

from distutils.core import setup
    
setup(      
            name="grid_remap",
            description="Hippocampal Remapping and Modular Grid Cells",
            version="1.0",
            author="Joe Monaco",
            author_email="jmonaco@jhu.edu",
            url="http://jdmonaco.com/grid-remap/",
            platforms=['Mac OSX', 'Linux', 'Windows'],
            license='The MIT License',
            package_dir={   'grid_remap': 'src'},
            packages=[  'grid_remap', 
                        'grid_remap.core', 
                        'grid_remap.analysis', 
                        'grid_remap.tools'  ]
)
