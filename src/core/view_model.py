# encoding: utf-8
"""
model.view_model -- Functionality for Model-based GUIs

This is complementary to the AbstractModel class and not a subclass of it. 
Subclasses should inherit from both Model (or a Model subclass) and this class
to take advantage of the visualization capability:

class MyModelGUI(ViewModel, Model):
    pass

ViewModel -- Super class with useful functionality for real-time data views
ViewModelHandler -- Handler subclass that starts a wx.Timer for refresh

Copyright (c) 2008 Columbia University. All rights reserved.
"""

# Library imports
import wx, numpy, threading

# Traits imports
from enthought.traits.api import HasTraits, Property, Range, Trait, Instance, \
    Float, Bool, Button
from enthought.traits.ui.api import Handler


class ViewModel(HasTraits):

    """
    Visualization helper class for AbstractModel subclasses
    
    The ViewModel base class complements AbstractModel by providing capabilities 
    for visualizing data traces in real-time as a AbstractModel simulation 
    progresses.
    
    Subclasses must override:
    _update_plots() -- update data sources for plots
    
    Traits and methods for subclass use:
    trail_length -- length of time into the past for data trails
    refresh_rate -- rate (Hz) at which plots are refreshed
    _trail(key) -- get data trail for a time-series (or time for key='t')
    _program_flow -- Button trait for pause/unpause button
    _reset_model -- Button trait for a model reset button
    
    Subclasses should inherit from ViewModel and Model (in that order).
    """
    
    # Traits for Play/Pause and Reset buttons
    _program_flow = Button
    _reset_simulation = Button

    # Data trail length in seconds
    trail_length = Range(low=0.0, value=5.0, exclude_low=True)
    
    # Thread for running the simulation loop
    _thread_ = Instance(threading.Thread)
    _view_open_ = Bool(False)
    
    # wx.Timer for refreshing plots
    refresh_rate = Range(low=2, high=30, value=12)
    _timer_ = Instance(wx.Timer)
    _t = Float(-1)
    
    # Public methods
    
    def _update_plots(self):
        """Subclasses override this to update data sources"""
        raise NotImplementedError
    
    def _trails(self, *keys):
        """
        Get the data trails for time-series indicated by keys tuple:        
        args -- string names of a trait attributes being tracked; passing in 
                't' returns the current time slice.
        
        Returns the trail data as a tuple of rank-1 arrays.
        """
        return self.ts.data_slice(keys, -self.trail_length)
    
    # Private methods
    
    def _update_(self, event):
        """Update plot data as model advances its timestep"""
        if self.t != self._t:
            self._update_plots()
            self._t = self.t
    
    def _pause_changed(self, paused):
        """Overriding Model to call advance() in its own thread"""
        if paused:
            if self._timer_.IsRunning():
                self._timer_.Stop()
        elif not (self._thread_ is not None and self._thread_.isAlive()):
            self._thread_ = threading.Thread(target=self.advance)
            self._thread_.start()
            self._timer_.Start(1000/float(self.refresh_rate), wx.TIMER_CONTINUOUS)
    
    def _done_changed(self, finished):
        if finished:
            self.pause = True
    
    # Actions for pause and reset button presses

    def __program_flow_fired(self):
        self.pause = not self.pause
    
    def __reset_simulation_fired(self):
        self.done = True
        self.reset()
        

class ViewModelHandler(Handler):
    """Handler subclass that initiates a timer for a ViewModel subclass"""

    def init(self, info):
        """Start timer after window creation and store in model object"""
        obj = info.object
        obj._view_open_ = True
        frame = info.ui.owner.control
        timer_id = wx.NewId()
        timer = wx.Timer(frame, timer_id)
        frame.Bind(wx.EVT_TIMER, obj._update_, id=timer_id)
        obj._timer_ = timer
    
    def closed(self, info, is_ok):
        """Join simulation thread for OK if still running, otherwise cleanup"""
        obj = info.object
        obj._view_open_ = False
        if is_ok:
            if obj.pause and not obj.done:
                obj.growl = True
                obj.out('Simulation is still paused. Start by setting pause ' +
                    'attribute to False.', 'Reminder', 'info')
            elif obj._thread_ is not None and obj._thread_.isAlive():
                obj._thread_.join()
        else:
            obj.done = True
