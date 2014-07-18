#!/usr/bin/env python
#encoding: utf-8
"""
grid.place_network_ui.py -- Chaco2 and Traits UI interface for PlaceNetwork model class

Subclass of grid.place_network.PlaceNetwork which enables a graphical interface.

Copyright (c) 2007 Columbia University. All rights reserved.
"""

# Library imports
import numpy, time, wx

# Package imports
from .place_network import PlaceNetwork
from .core.view_model import ViewModel, ViewModelHandler
from .tools.images import array_to_rgba

# Traits imports
from enthought.traits.api import Property, Array, Instance, Range, Float, Bool
from enthought.traits.ui.api import View, VSplit, Group, Item, Heading, Include

# Chaco2 imports
from enthought.chaco.api import Plot, PlotLabel, ArrayPlotData, DataRange1D, hot
from enthought.enable.component_editor import ComponentEditor

# Global constants
COL_WIDTH = 240
DISP_UNITS = 200

# View for model initialization
initialization_view = View(
                            Heading('Set values prior to running simulation'),
                            Group(
                                Group(
                                    Item(name='desc', label='Description', width=COL_WIDTH),
                                    Item(name='num_trials', label='Trials'),
                                    Item(name='traj_type', label='Trajectory'),
                                    Item(name='dwell_factor', enabled_when="traj_type!='randwalk'"),
                                    Item(name='T', label='Duration (s)', enabled_when="traj_type=='randwalk'"),
                                    Item(name='dt', label='dt (s)'),
                                    Item(name='monitor_dt', label='Display dt (s)'),
                                    Item(name='N_CA', label='Output Units', style='readonly'),
                                    Item(name='N_EC', label='Grid Inputs', style='readonly'),
                                    Item(name='C_W', label='Connectivity'),
                                    label='Initialization', 
                                    show_border=True),
                                Group(
                                    Item(name='J0', width=COL_WIDTH),
                                    Item(name='tau_r', label='Tau (s)'),
                                    Group(
                                        Item(name='phi_lambda', label='Threshold'),
                                        Item(name='phi_sigma', label='Smoothness'),
                                        label='Field Nonlinearity:'),
                                    label='Parameters', 
                                    show_border=True),
                                orientation='horizontal'),
                            Item(name='growl', label='Growl On/Off', style='simple'),
                            Item(name='projdir', label='Project Directory'),
                            buttons=['Revert', 'Cancel', 'OK'],
                            title='Model Setup', 
                            kind='livemodal',
                            resizable=False)

# View for main simulation visualization
simulation_view = View(
                    VSplit(
                        Group(
                            Item(name='field_plots', editor=ComponentEditor()), 
                            Item(name='units_plot', editor=ComponentEditor()), 
                            Item(name='traj_plot', editor=ComponentEditor()), 
                            show_border=False,
                            show_labels=False,
                            orientation='horizontal'),
                        Group(
                            Group(
                                Group(
                                    Item(name='J0'), 
                                    label='Input Gain', 
                                    show_border=True),
                                Group(  
                                    Item(name='phi_lambda', label='Lambda'), 
                                    Item(name='phi_sigma', label='Sigma'), 
                                    label='Field Nonlinearity', 
                                    show_border=True),
                                Group(  
                                    Item(name='trail_length'), 
                                    Item(name='_program_flow', show_label=False, enabled_when='done==False'),
                                    Item(name='_reset_simulation', show_label=False),
                                    show_border=False),
                                springy=True,
                                show_border=False),
                            Item(name='phi_plot', editor=ComponentEditor(), width=0.4), 
                            show_border=False,
                            show_labels=False,
                            orientation='horizontal'),
                        show_labels=False,
                        show_border=False),
                    title='PlaceNetwork Simulation',
                    resizable=True,
                    height=1.0,
                    width=1.0,
                    buttons=['Revert', 'Cancel', 'OK'],
                    kind='live')
                    

class PlaceNetworkUI(ViewModel, PlaceNetwork):
    
    """
    PlaceNetwork model with Traits UI graphical real-time interface
    """
    
    pause = True
    
    # View traits
    sim_view = simulation_view
    init_view = initialization_view
    
    # Plot instances
    field_plots = Instance(Plot)
    units_plot = Instance(Plot)
    traj_plot = Instance(Plot)
    phi_plot = Instance(Plot)
    
    # Redefine user parameters as Range traits for sliders
    J0 = Range(low=0.0, high=100.0, value=45)
    
    # Control nonlinearity variables with sliders
    phi_lambda = Range(low=0.0, high=0.5, value=0.04)
    phi_sigma = Range(low=0.001, high=1.0, value=0.02)
    
    # Field plots tracking data
    h_aff = Property(Float, track=True)
    h_rec = Property(Float, track=True)
    h_sum = Property(Float, track=True)
    
    # Phi plot data
    h_range = Property(Array)
    phi_sample = Property(Array)
    _phi_updated = Bool(False)
    
    t0 = Float
    
    # Add to the simulation timestep
    
    def run_timestep(self):
        PlaceNetwork.run_timestep(self)
        
        # NOP to allow GUI to process
        time.sleep(.001)

    # Creating Plot instances as trait default functions
    
    def _field_plots_default(self):
        zero = numpy.array([0], 'd')
        data = ArrayPlotData(t=zero, h_rec=zero, h_aff=zero, h_sum=zero)
        p = Plot(data)
        p.plot(('t', 'h_aff'), name='Afferent', type='line', line_width=1, color='royalblue')
        p.plot(('t', 'h_rec'), name='Recurrent', type='line', line_width=1, color='tomato')
        p.plot(('t', 'h_sum'), name='Total', type='line', line_width=1, color='sienna')
        
        p.legend.visible = True
        p.legend.border_visible = False
        p.legend.align = 'ur'
        p.legend.bgcolor = (0.8, 0.8, 1.0, 0.4)
        p.legend.border_padding = 6
        p.legend.labels = ['Afferent', 'Recurrent', 'Total']
        
        p.y_grid.visible = p.x_grid.visible = False
        p.title = 'Synaptic Fields'
        p.x_axis.title = 'Time (s)'
        p.y_axis.title = 'Field Strength'
        p.bgcolor = 'mintcream'
        return p
    
    def _units_plot_default(self):
        N = min([self.N_CA, DISP_UNITS])
        data = ArrayPlotData(i=numpy.arange(N), r=self.r[:N], i_aff=self.i_aff[:N])
        p = Plot(data)
        p.plot(('i', 'r', 'i_aff'), type='cmap_scatter', color_mapper=hot, 
            marker='circle', marker_size=3, line_width=0)
        p.title = 'Place Cell Output'
        p.x_axis.title = 'Output Units'
        p.y_axis.title = 'Rate / Iaff'
        p.value_range.set_bounds(0.0, 1.0)
        p.x_grid.visible = p.y_grid.visible = False
        p.bgcolor = 'slategray'
        return p
    
    def _traj_plot_default(self):
        """Trajectory plot based on TrajectoryView.t_plot in chaco_threading_demo"""
        zero = numpy.array([0], 'd')
        data = ArrayPlotData(x=zero, y=zero)

        h, w = self.traj.Map.H, self.traj.Map.W
        data.set_data('x0', zero + self.traj.Map.x0[0])
        data.set_data('y0', zero + self.traj.Map.x0[1])
        
        p = Plot(data)
        p.plot(('x', 'y'), name='trail', color='red')
        p.plot(('x0', 'y0'), name='head', type='scatter', marker='circle', color='red')
        p.y_axis.visible = p.x_axis.visible = False
        p.y_grid.visible = p.x_grid.visible = False
        p.border_visible = True
        p.border_width = 2
        p.title = 'Rat Trajectory'
        
        p.index_range.set_bounds(0, w)
        p.value_range.set_bounds(0, h)
        
        p.overlays.append(PlotLabel('X (%d cm)'%w, component=p, overlay_position='bottom'))
        p.overlays.append(PlotLabel('Y (%d cm)'%h, component=p, overlay_position='left', angle=90))
        return p
        
    def _phi_plot_default(self):
        data = ArrayPlotData(h=self.h_range, phi=self.phi_sample)
        p = Plot(data)
        p.plot(('h', 'phi'), type='line', name='phi', color='slateblue', line_width=2.7)
        p.x_axis.title = 'h'
        p.y_axis.title = 'Phi[h]'
        p.x_grid.line_color = p.y_grid.line_color = 'slategray'
        p.bgcolor = 'khaki'
        p.title = 'Nonlinearity'
        return p
        
    # Callback for updating plot data
    
    def _update_plots(self):
        # Field plots data trails
        t, h_aff, h_rec, h_sum = self._trails('t', 'h_aff', 'h_rec', 'h_sum')
        if self.t > self.dt:
            self.field_plots.data.set_data('t', t)
            self.field_plots.data.set_data('h_aff', h_aff)
            self.field_plots.data.set_data('h_rec', h_rec)
            self.field_plots.data.set_data('h_sum', h_sum)
            
            # Trajectory trails
            new_x, new_y = self._trails('x', 'y')
            self.traj_plot.data.set_data('x', new_x)
            self.traj_plot.data.set_data('y', new_y)
            self.traj_plot.data.set_data('x0', numpy.array([new_x[-1]]))
            self.traj_plot.data.set_data('y0', numpy.array([new_y[-1]]))
        
        # Units plot update
        N = min([self.N_CA, DISP_UNITS])
        self.units_plot.data.set_data('r', self.r[:N])
        self.units_plot.data.set_data('i_aff', self.i_aff[:N])
        self.units_plot.value_range.high_setting = max([1, 1.05*self.r[:N].max()])
            
        # Phi data update
        self._update_phi_plot()
    
    def _update_phi_plot(self):
        if self._phi_updated:
            self.phi_plot.data.set_data('h', self.h_range)
            self.phi_plot.data.set_data('phi', self.phi_sample)
            self._phi_updated = False       
    
    # Trajectory changes refresh the plot
    
    def _stage_changed(self):
        self._refresh_traj_plot()
    
    def _traj_type_changed(self):
        self._refresh_traj_plot()
    
    def _refresh_traj_plot(self):
        self.traj = self.new_trajectory()
        self.traj_plot = self._traj_plot_default()
    
    # Field tracking properties and trait notifications
    
    def _get_h_aff(self):
        return self.i_aff.mean()
    
    def _get_h_rec(self):
        return -self.J0 * self.r.sum()
    
    def _get_h_sum(self):
        return self.h_aff + self.h_rec
    
    # Nonlinearity plot automation and line data
    
    def _get_phi_sample(self):
        return self.phi_h(self.h_range - self.phi_lambda)
    
    def _get_h_range(self):
        return numpy.arange(0, max([2.5, 2.5*self.phi_lambda]), 0.02)
     
    def _phi_lambda_changed(self):
        self._phi_pause_update()
        
    def _phi_sigma_changed(self):
        self._phi_pause_update()

    def _phi_pause_update(self):
        self._phi_updated = True
        if self.pause:
            self._update_phi_plot()
    
    # Convenience functions for calling views
    
    def setup(self):
        self.configure_traits(view='init_view')
    
    def simulation(self):
        self.configure_traits(view='sim_view', handler=ViewModelHandler())


if __name__ == "__main__":
    import os
    from .dmec import GridCollection
    EC = GridCollection()
    ca3 = PlaceNetworkUI(EC=EC, C_W=0.33, growl=False, T=300, desc='demo run')
    ca3.setup()
    ca3.simulation()
    