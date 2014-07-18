# encoding: utf-8
"""
placemap_viewer.py -- An interactive GUI interface for individual spatial maps

Created by Joe Monaco on 04-30-2008.
Copyright (c) 2008 Columbia University. All rights reserved.
"""

# Library imports
import numpy as N, scipy as S
from matplotlib import cm

# Package imports
from .ratemap import PlaceMap
from .tools.images import array_to_rgba
from .tools.stats import integer_hist
from .tools.bash import CPrint

# Traits imports
from enthought.traits.api import HasTraits, Instance, Trait, TraitError, \
    Property, Enum, Int, Float, Range, Delegate
from enthought.traits.ui.api import View, Group, Item, Heading

# Chaco imports
from enthought.chaco.api import ArrayPlotData, Plot, BasePlotContainer, VPlotContainer, copper
from enthought.enable.component_editor import ComponentEditor


class PlaceMapViewer(HasTraits):
    
    """
    Chaco viewer for placemap data
    
    Constructor arguments:
    pmap -- PlaceMap (or subclass) object to view
    
    Public methods:
    view -- Bring up the Chaco View window for looking at data
    """
    
    # Console output
    out = Instance(CPrint)
    
    # Reference to PlaceMap object
    PMap = Trait(PlaceMap)
    
    # Stage map traits
    stage_map = Instance(Plot)
    stage_map_type = Enum('representation', 'coverage', 'field_centers')
    sparsity = Delegate('PMap')
    num_active = Delegate('PMap')
    stage_coverage = Delegate('PMap')
    stage_repr = Delegate('PMap')
    peak_rate = Delegate('PMap')
    
    # Unit map traits
    _unit = Int
    unit_map = Instance(Plot)
    unit_map_type = Enum('ratemap', 'single', 'fields')
    num_fields = Int
    coverage = Float
    avg_area = Float
    avg_diameter = Float
    max_rate = Float
    
    # Unit data traits
    unit_data_plots = Instance(BasePlotContainer)
    unit_bins = Range(low=5, high=50, value=20)
    
    # Field data traits 
    field_data_plots = Instance(BasePlotContainer)
    field_bins = Range(low=5, high=50, value=20)
    
    # Chaco view definition
    traits_view = \
        View(
            Group(
                Group(
                    Item('stage_map_type'),
                    Item('stage_map', editor=ComponentEditor(), show_label=False),
                    Group(
                        Item('sparsity', style='readonly'),
                        Item('num_active', style='readonly'),
                        Item('stage_coverage', label='Coverage', style='readonly'),
                        Item('stage_repr', label='Representation', style='readonly'),
                        Item('peak_rate', style='readonly'),
                        label='Stage Coding',
                        show_border=True),
                    label='Stage Maps',
                    orientation='v'), 
                Group(
                    Item('unit_map_type'),
                    Item('unit', style='custom'),
                    Item('unit_map', editor=ComponentEditor(), show_label=False),
                    Group(
                        Item('max_rate', style='readonly'),
                        Item('num_fields', style='readonly'),
                        Item('coverage', style='readonly'),
                        Item('avg_area', label='Mean Field Area', style='readonly'),
                        Item('avg_diameter', label='Mean Field Diameter', style='readonly'),
                        label='Place Unit',
                        show_border=True),
                    label='Unit Maps',
                    orientation='v'),
                Group(
                    Heading('Distributions of Single-Unit Properties'),
                    Item('unit_data_plots', editor=ComponentEditor(), show_label=False),
                    Item('unit_bins', label='Bins'),
                    label='Unit Data'),
                Group(
                    Heading('Distributions of Single-Field Properties'),
                    Item('field_data_plots', editor=ComponentEditor(), show_label=False),
                    Item('field_bins', label='Bins'),
                    label='Field Data'),
                layout='tabbed'),
            title='Placemap Viewer',
            resizable=True,
            height=800,
            width=700,
            kind='live',
            buttons=['Cancel', 'OK'])    

    def __init__(self, pmap, **traits):
        HasTraits.__init__(self, **traits)
        try:
            self.PMap = pmap
        except TraitError:
            self.out('PlaceMap subclass instance required', error=True)
            return
        self.fdata = self.PMap.get_field_data()
        self.udata = self.PMap.get_unit_data()
        self.add_trait('unit', Range(low=0, high=self.PMap.num_maps-1))
        self._update_unit_values()
        self.out('Bringing up place-map visualization...')
        self.view()
        self.out('Done!')
    
    def view(self):
        self.configure_traits()
    
    # Plot creation methods
    
    def _stage_map_default(self):
        
        # RGBA maps
        rep_map = array_to_rgba(self.PMap.stage_repr_map, cmap=cm.hot)
        cov_map = array_to_rgba(self.PMap.stage_coverage_map, cmap=cm.gray)
        
        # Data sources and plot object
        data = ArrayPlotData(fields_x=self.fdata['x'], fields_y=self.fdata['y'], 
            fields_z=self.fdata['peak'], rep=rep_map, cov=cov_map)
        p = Plot(data)
        
        # Plot the field centers
        p.plot(('fields_x', 'fields_y', 'fields_z'), name='centers', type='cmap_scatter', 
            marker='dot', marker_size=5, color_mapper=copper, line_width=1, fill_alpha=0.6)
        
        # Plot the representation and coverage maps
        p.img_plot('rep', name='rep', xbounds=(0, self.PMap.W), ybounds=(0, self.PMap.H),
            origin='top left')
        p.img_plot('cov', name='cov', xbounds=(0, self.PMap.W), ybounds=(0, self.PMap.H),
            origin='top left')
        
        # Start with only the representation map visible
        p.plots['cov'][0].visible = False
        p.plots['centers'][0].visible = False
        
        # Plot tweaks
        p.aspect_ratio = 1.0
        p.y_axis.title = 'Y (cm)'
        p.x_axis.title = 'X (cm)'
        p.x_axis.orientation = 'bottom'
        p.title = 'Stage Maps'
        
        return p
    
    def _unit_map_default(self):
        
        # Set the initial unit map
        data = ArrayPlotData(unit_map=self._get_unit_map_data())
        p = Plot(data)
        
        # Plot the map
        p.img_plot('unit_map', name='unit', xbounds=(0, self.PMap.W), ybounds=(0, self.PMap.H),
            origin='top left')
        
        # Plot tweaks
        p.aspect_ratio = 1.0
        p.y_axis.title = 'Y (cm)'
        p.x_axis.title = 'X (cm)'
        p.x_axis.orientation = 'bottom'
        p.title = 'Single Unit Maps'
        
        return p
    
    def _unit_data_plots_default(self):
        
        # Plot data and vertical container object
        data = ArrayPlotData(**self._get_unit_plots_data())
        container = VPlotContainer()
        
        # Add individual distribution plots to container
        for key in ('avg_diameter', 'avg_area', 'coverage', 'max_r', 'num_fields'):
            p = Plot(data)
            p.plot((key+'_bins', key), name=key, type='polygon', edge_width=2, 
                edge_color='mediumblue', face_color='lightsteelblue')
            p.x_axis.title = key
            p.y_axis.title = 'count'
            p.padding = [50, 30, 20, 40]
            if key == 'num_fields':
                p.x_axis.tick_interval = 1
            container.add(p)
        
        return container
        
    def _field_data_plots_default(self):
        
        # Plot data and vertical container object
        data = ArrayPlotData(**self._get_field_plots_data())
        container = VPlotContainer()
        
        # Add individual distributions plots to container
        for key in ('area', 'diameter', 'average', 'peak'):
            p = Plot(data)
            p.plot((key+'_bins', key), name=key, type='polygon', edge_width=2, 
                edge_color='red', face_color='salmon')
            p.x_axis.title = key
            p.y_axis.title = 'count'
            p.padding = [50, 30, 20, 40]
            container.add(p)
        
        return container
        
    # Plot update methods
    
    def _update_stage_map(self):
        """Handle switching between different stage maps"""
        
        # Update and equalize bounds for all subplots
        self.stage_map.plots['rep'][0].bounds = self.stage_map.bounds
        self.stage_map.plots['cov'][0].bounds = self.stage_map.bounds
        self.stage_map.plots['centers'][0].bounds = self.stage_map.bounds
        
        # Set visibility flags
        if self.stage_map_type is 'representation':
            self.stage_map.title = 'Relative Representation'
            vis_plots = (True, False, False)
        elif self.stage_map_type is 'coverage':
            self.stage_map.title = 'Total Stage Coverage'
            vis_plots = (False, True, False)
        elif self.stage_map_type is 'field_centers':
            self.stage_map.title = 'Place Field Centroids'
            vis_plots = (False, False, True)
        
        # Toggle plot visibility and redraw
        self.stage_map.plots['rep'][0].visible, \
            self.stage_map.plots['cov'][0].visible, \
            self.stage_map.plots['centers'][0].visible = vis_plots
        self.stage_map.request_redraw()
    
    def _update_unit_map(self):
        """Update current image source and title; then redraw the plot"""
        self.unit_map.data.set_data('unit_map', self._get_unit_map_data())
        self.unit_map.title = '%s of Unit %d'%(self.unit_map_type.capitalize(), self.unit)
        self.unit_map.request_redraw()
    
    def _update_unit_values(self):
        """Update the scalar readonly values"""
        if self._unit == -1:
            self.num_fields = 0
            self.coverage = self.avg_area = self.avg_diameter = 0.0
            self.max_rate = self.PMap.maxima[self.unit, 2]
        else:
            self.num_fields = int(self.udata[self._unit]['num_fields'])
            self.coverage = float(self.udata[self._unit]['coverage'])
            self.avg_area = float(self.udata[self._unit]['avg_area'])
            self.avg_diameter = float(self.udata[self._unit]['avg_diameter'])
            self.max_rate = float(self.udata[self._unit]['max_r'])
    
    def _get_unit_map_data(self):
        """Helper function to get RGBA array for current unit and map type"""
        if self.unit_map_type is 'ratemap':
            map_data = array_to_rgba(self.PMap.Map[self.unit], cmap=cm.jet, 
                norm=False, cmax=self.peak_rate)
        elif self.unit_map_type is 'single':
            map_data = array_to_rgba(self.PMap.single_maps[self.unit], cmap=cm.hot)
        elif self.unit_map_type is 'fields':
            map_data = array_to_rgba(self.PMap.coverage_maps[self.unit], cmap=cm.gray)
        return map_data
    
    def _get_unit_plots_data(self):
        """Helper function for getting unit data distributions"""        
        
        # Integer distribution for number of fields
        data = {}
        data['num_fields_bins'], data['num_fields'] = integer_hist(self.udata['num_fields'])
        
        # Continuous distributions of other unit statistics
        for key in ('avg_area', 'avg_diameter', 'coverage', 'max_r'):
            keyb = key + '_bins'
            data[key], data[keyb] = S.histogram(self.udata[key], bins=self.unit_bins)
            data[keyb] += (data[keyb][1] - data[keyb][0]) / 2
            data[keyb] = data[keyb][:-1]
        
        # Add 0-value end-points for polygon display
        for key in data:
            if key[-4:] == 'bins':
                data[key] = N.r_[data[key][0], data[key], data[key][-1]]
            else:
                data[key] = N.r_[0, data[key], 0]
        
        return data
    
    def _get_field_plots_data(self):
        """Helper function for getting field data distributions"""
        
        # Continuous distributions of place field properties
        data = {}
        for key in ('area', 'diameter', 'average', 'peak'):
            keyb = key + '_bins'
            data[key], data[keyb] = S.histogram(self.fdata[key], bins=self.field_bins)
            data[keyb] += (data[keyb][1] - data[keyb][0]) / 2
            data[keyb] = data[keyb][:-1]
            
        # Add 0-value end-points for polygon display
        for key in data:
            if key[-4:] == 'bins':
                data[key] = N.r_[data[key][0], data[key], data[key][-1]]
            else:
                data[key] = N.r_[0, data[key], 0]
        
        return data

    # Map traits notifications
    
    def _unit_bins_changed(self):
        """Update plot data for unit distributions"""
        data = self._get_unit_plots_data()
        plot_data = self.unit_data_plots.components[0].data
        for key in data:
            plot_data.set_data(key, data[key])
    
    def _field_bins_changed(self):
        data = self._get_field_plots_data()
        plot_data = self.field_data_plots.components[0].data
        for key in data:
            plot_data.set_data(key, data[key])
    
    def _stage_map_type_changed(self):
        self._update_stage_map()
    
    def _unit_map_type_changed(self):
        self._update_unit_map()
        
    def _unit_changed(self):
        """Update the unit map and scalar values"""
        find_unit = (self.udata['unit'] == self.unit).nonzero()[0]
        if find_unit.shape[0]:
            self._unit = find_unit[0]
        else:
            self._unit = -1
        self._update_unit_map()
        self._update_unit_values()
    
    # Output object default
    
    def _out_default(self):
        return CPrint(prefix=self.__class__.__name__, color='purple')
            