#encoding: utf-8
"""
tools.cmap_ui -- Mix-in functionality for Chaco interactive colormaps

Classes implementing a Chaco plot in need of interactive colormapping 
should inherit from ColormapControl.

See ColormapControl docstring for subclassing usage.

Written by Joe Monaco, 05/13/2008.
Copyright (c) 2008 Columbia University. All rights reserved.
"""

# Library imports
import numpy as N

# Import the matplotlib colormaps and sort the names
from matplotlib import cm
cm._cmapnames.sort()

# Package imports
from .images import array_to_rgba
from .colormaps import diffmap

# Traits imports
from enthought.traits.api import HasTraits, Property, Trait, Bool, Int
from enthought.traits.ui.api import Group, Item

# Chaco imports
from enthought.chaco.api import Plot, ArrayPlotData


class ColormapControl(HasTraits):
    
    """
    Provides instant colormap handling for interactive Chaco plots
    
    Public traits:
    colormap -- string name of current colormap; can be set to the name of any
        valid matplotlib colormap or 'difference'
    reversed_colormap -- whether the map is reversed or not
    
    Public methods for subclass use:
    get_rgba_data -- convert an intensity matrix into an RGBA image array
    get_colorbar_plot -- get a Plot object containing a new colorbar
    get_colormap_object -- get the current MPL colormap object
    
    Update notification:
    _cmap_notify_changed -- if defined, this method will listen for changes to the 
        colormap; it should update and redraw the plot
    
    NOTE: You can include *colormap_group* in a Traits UI View:    
    from enthought.traits.ui.api import Include
    ...
    class MyFigure(ColormapControl):
        traits_view = View(...
            Include('colormap_group'),
            ...)
    """
    
    colormap_group = \
        Group(
            Item('colormap', label='Color Map'),
            Item('reversed_colormap', label='Reversed'),
            label='Colors',
            show_border=True)
    
    # Public traits
    colormap = Trait('hot', cm._cmapnames, 'difference')
    reversed_colormap = Bool(False)
    
    # Colormap change notifier
    cmap_notify = Bool
    
    # Colorbar width and orientation
    _cbar_width = Int(30)
    _cbar_orientation = Trait('v', 'h')
    
    def get_rgba_data(self, M):
        """Get an RGBA image array based on an array M of intensity values"""
        return array_to_rgba(M, cmap=self.get_colormap_object())
        
    def get_colorbar_plot(self, bounds=(0,1)):
        """
        Create a colorbar plot to be added to a plot-container
        
        Arguments:
        bounds -- (min, max) tuple sets the intensity range for the colormap
        
        Returns a Chaco2 Plot object containing the colorbar.        
        """
        cb_rgba = array_to_rgba(
            N.repeat(N.linspace(1, 0, num=1024)[:,N.newaxis], 20, axis=1), 
            cmap=self.get_colormap_object())
        if self._cbar_orientation is 'h':
             cb_rgba = cb_rgba.T[::-1]
        data = ArrayPlotData(colorbar=cb_rgba)
        
        # Create the plot object
        cb = Plot(data, width=self._cbar_width, resizable=self._cbar_orientation, 
            padding_left=0, padding_top=0)
        cb.img_plot('colorbar', name='colorbar', xbounds=bounds, ybounds=bounds,
            origin='top left')
        
        # Plot tweaks
        if self._cbar_orientation is 'v':
            cb.x_axis.visible = False
            cb.y_axis.orientation = 'right'
        else:
            cb.y_axis.visible = False
            cb.x_axis.orientation = 'bottom'
    
        return cb

    def get_colormap_object(self):
        """
        Return the colormap specified by *colormap* and *reversed_colormap*
        """
        cmap_name = self.colormap
        if cmap_name == 'difference':
            return diffmap()
        if self.reversed_colormap:
            cmap_name += '_r'
        return getattr(cm, cmap_name)

    def _anytrait_changed(self, name):
        if name in ('colormap', 'reversed_colormap'):
            self.cmap_notify = not self.cmap_notify
