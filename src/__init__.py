"""
Grid Remap Package -- Data structures, models, and classes for simulating 
    grid cells and characterizing place cell responses

Author: Joe Monaco
Created: 02/01/2007
Maintainer: joe@neurotheory.columbia.edu

Copyright (c) 2007-2009 Columbia University. All rights reserved.
"""


# Environment maps

from .stage import StagingMap
from .placemap import AbstractPlaceMap, PlaceMap
from .ratemap import AbstractImpulseRatemap, CheckeredRatemap
from .placemap_viewer import PlaceMapViewer


# Environment trajectories

from .trajectories import (  BaseTrajectory, 
                             RandomWalk, 
                             AbstractImpulseRaster, 
                             BipartiteRaster    )


# Grid cell model

from .dmec import GridCollection


# Place network model

from .place_network import PlaceNetwork, PlaceNetworkRaster, PlaceNetworkStd
from .place_network_ui import PlaceNetworkUI


# Analysis classes

from .analysis import *
