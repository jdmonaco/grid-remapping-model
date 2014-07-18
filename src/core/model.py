# encoding: utf-8
"""
model.py -- AbstractModel is a framework for writing time-series models

AbstractModel provides a Traits-based foundation for models. It provides a flexible
messaging service with Messager, time-series data collection and tracking with 
TimeSeries, exception handling and saves results and complete information 
about each simulation. 

Model subclasses simply define all the Traits and support methods needed for 
the scientific problem and then implement per-trial setup (trial_setup()) and
the per-timestep computation (run_timestep()). AbstractModel handles the rest.

AbstractModel -- Abstract modeling framework
ModelPostMortem -- Skeleton object storing completed simulation results

Copyright (c) 2007, 2008 Columbia University. All rights reserved.

This software is provided AS IS under the terms of the Open Source MIT License. 
See http://www.opensource.org/licenses/mit-license.php.
"""

# Library imports 
import time, os

# Package imports
from . import TS_FMT, MODEL_MESSAGES as MSG
from .timeseries import TimeSeries, TSPostMortem
from ..tools.array_container import ArrayContainer
from ..tools.messager import Messager
from ..tools.path import unique_path

# Traits imports
from enthought.traits.api import (HasTraits, Property, Instance, Directory, 
    File, Trait, Range, List, String, Int, Float, true, false)
        

class ModelPostMortem(ArrayContainer):
    
    """
    Converts AbstractModel.results data into attributes of skeleton object
    
    The shape of the array attributes depends on whether one or more trials were
    run in the model simulation:
    - num_trials == 1: time is the first dimension
    - num_trials > 1: trials x time are the first two dimensions 
    - in both cases: the shape of the tracked data follows time. 
    
    Public methods:
    get_trial_data -- single-trial post-mortem for a particular trial (trials are 
        indexed starting from 1)
    
    The **tracking** and **ntrials** attributes are set as a convenience.
    """ 
    
    def __init__(self, data):
        """Transform AbstractModel.results into array attributes on this object
        
        For multi-trial data, each data attribute is a numpy object array where
        each element is a trial data array. This avoids copying data since non-
        contiguous ndarray creation is iffy at this point.
        """
        from numpy import array
        ntrials = len(data)
        self.ntrials = ntrials
        if ntrials == 0:
            return
        self.tracking = data[0].tracking
        self.t = data[0].t
        _d = {}
        for key in self.tracking:
            if ntrials == 1:
                _d[key] = getattr(data[0], key)
            else:
                _d[key] = array([None]*ntrials, dtype='O')
                for trial in xrange(ntrials):
                    _d[key][trial] = getattr(data[trial], key)
        self.__dict__.update(_d)
    
    def get_trial_data(self, trial):
        """Extract a single trial of data from this object
        
        Required argument:
        trial -- trial number, indexed from 1
        
        Returns a single-trial data object.
        """
        if self.ntrials == 1:
            return self
        if trial < 1 or trial > self.ntrials:
            raise ValueError, 'bad trial number (%d)'%trial

        trial = int(trial)
        _d = {}
        for key in self.tracking:
            _d[key] = getattr(self, key)[trial-1]
        _d['t'] = self.t
        
        return ModelPostMortem([TSPostMortem(**_d)])


class AbstractModel(HasTraits):

    """
    Generalized base model class for constructing a time-series model
    
    Model parameters and settings should be specified as keyword arguments to 
    the subclass constructor.
    
    Attributes are categorized by Trait metadata:
    track -- setting to True will cause attribute to be tracked as a time-series
    user -- specify whether this attribute is a user-settable parameter
    
    Model attributes:
    label -- short name for the model (default module name)
       * this will be used as the project subdirectory
    desc -- descriptive phrase to be added to filenames (default label)
    app_name -- full name for this model (default label)
    icon_path -- full path to the model icon (default mandle_joe.png)
    console -- bool value toggles console messages (default True)
    num_trials -- number of trials to run (default 1)
    
    Important attributes:
    out -- the Messager object that handles messaging
    ts -- the TimeSeries object that handles all data collection
    trial -- current trial number
    trial_result -- if no variables are being tracked by *ts*, then the
        attribute whose name is stored in trial_result will be saved as the
        trial result data for each trial
    t -- current time in seconds of the current trial
    done -- whether all trials have completed
    pause -- setting to True pauses simulation (advance() restarts)
    monitoring -- (bool) toggles progress monitoring messages

    Public methods:
    info -- display a nice complete description of this model
    post_mortem -- get a ModelPostMortem object for the completed data
    reset -- trash any progress and reset the simulation
    advance -- run the next trial (or resume from pausing)
    advance_all -- run through all trials
    save_data -- save simulation data and a model information file
    parameter_dict -- get a dictionary of all user-settable parameters

    Subclasses override these methods:
    trial_setup -- setup to be performed prior to each trial
    run_timestep -- the per-timestep computation for this model; this is 
        wrapped by a while loop within advance(), e.g.:
            while self.ts:
                self.run_timestep()
                self.ts.advance()
    """
    
    # Output control
    out = Instance(Messager)
    
    # Program flow
    pause = false(user=True, label='Pause', desc='Toggle active state of model')
    done = false(user=True, label='Stop', desc='End model simulation now')
    
    # Metadata
    label = String(__name__)
    desc = String
    app_name = label
    
    # Progress monitoring
    monitoring = true(user=True)
    monitor = Property(false)
    monitor_dt = Float(1)
    
    # General model machinery
    trial_result = String
    timestamp = Trait(time.strftime(TS_FMT))
    success = true
    num_trials = Range(low=1)
    results = List
    trial = Int(1)
    ts = Instance(TimeSeries)
    dt = Trait(TimeSeries.__class_traits__['dt'])
    T = Trait(TimeSeries.__class_traits__['T'])
    t = Property(Float)
    
    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self._renew_timeseries()
        
    def __str__(self): 
        if self:
            progress = '\nAt: Trial %d of %d @ %.2f/%.2fs\n'%(self.trial, 
                self.num_trials, self.t, self.T)
        else:
            progress = '\nDone: Completed %d trials out of %d\n'% \
                (len(self.results), self.num_trials)
        user_params = self.traits(user=True).keys()
        user_params.sort()
        hdr = self.__class__.__name__ + '(Model) object\n' + '-'*32
        p_str = \
            '\nParameters:\n' + '\n'.join(['\t%s : %s'%(k, repr(getattr(self, k))) \
                for k in user_params])
        if not self.traits(track=True):
            t_str = ''
        else:
            t_str = \
                '\nTracking:\n\t' + ', '.join(self.traits(track=True).keys())
        return hdr + p_str + t_str
            
    def __repr__(self): 
        return self.__str__().split('\n')[0]
    
    def __nonzero__(self): 
        return not self.done 
    
    # Subclass override methods
    
    def trial_setup(self):
        """Subclass override: perform trial setup
        """
        raise NotImplementedError
        
    def run_timestep(self):
        """Subclass override: simulation time-step
        """
        raise NotImplementedError
    
    # Public methods
    
    def info(self): 
        self.out(self, 'Model information')
    
    def post_mortem(self): 
        """Get a PostMortem object for this instance
        """
        return ModelPostMortem(self.results)

    def reset(self):
        """Wipe out all results and go back to initial state
        """
        self.done = False
        self.pause = False
        self.trial = 1
        self.results = []
        self._renew_timeseries()
        self.out('Model has been reset!', 'Reset', 'init')
    
    def advance(self):
        """Run the current trial, save data, and handle exceptions
        """
        if self.done or self.pause:
            return time.sleep(.2)
        if self.t and self.ts:
            self.out('Unpausing trial %d...'%self.trial, notification='monitor')
        else:
            self.out('Starting trial %d...'%self.trial, notification='init')
            self.trial_setup()
            self._renew_timeseries()
        try:
            self._run_trial()
        except KeyboardInterrupt:
            self.out(*MSG['TRIALSTOP'])
            self.done = raw_input('> ').lower() == 'stop'
        except Exception, e:
            import pdb
            self.out(MSG['UNHEXC'][0]%(self.trial, self.t, repr(e).split('(')[0], 
                e.message), MSG['UNHEXC'][1], MSG['UNHEXC'][2])
            pdb.post_mortem(os.sys.exc_info()[2])
            self.success = False
            self.done = True
        else:
            if not self.pause:
                self._store_trial()
                if len(self.results) < self.trial:
                    self.trial = len(self.results)
                    self.out(MSG['TRIALRST'][0]%(self.trial+1), MSG['TRIALRST'][1], 
                        MSG['TRIALRST'][2])
                self.done = self.trial == self.num_trials
                self.trial += not self.done
        finally:
            if self.done:
                self.out('Finished simulation!', notification='complete')
    
    def advance_all(self):
        """Run through all trials
        """
        while not self.done:
            self.advance()

    def save_data(self, dpath='.'):
        """Save pickle files of model info and results data
        
        Saved file details:
        X.info -- a text file describing various properties of the model, its 
            state, and its parameters.
        X.tar.gz -- an ArrayContainer format archive of the post-mortem data of
            this model's results
        """     
        import cPickle
        if not self._create_datapath(dpath):
            self.out(MSG['DPATH'][0]%dpath, MSG['DPATH'][1], MSG['DPATH'][2])
            return

        fn_title = 'Save model data'
        stem_fn = '_'.join(self.label.lower().split())
        if self.desc:
            stem_fn += '-' + '_'.join(self.desc.lower().split())
        info_fn = unique_path(os.path.join(dpath, stem_fn), ext='info', 
            fmt='%s-%03d')
        stem = info_fn[:-4]
        data_fn = stem
        try:
            info_fd = open(info_fn, 'w')
        except IOError:
            self.out('Could not open file(s) for writing:\n->' + stem + 
                '{info,data}', fn_title, 'error')
        else:            
            # Write out the info file
            self._write_info_file(info_fd, stem)
            info_fd.close()
            
            # Save the post-mortem data in a compressed archive
            self.post_mortem().tofile(data_fn)
            
            # All done!
            self.out('Saved files:\n->%s\n->%s'%(info_fn, data_fn+'tar.gz'), 
                fn_title, sticky=True)
        
    def parameter_dict(self):
        """Get a dictionary with all user-settable parameters for this model
        """
        params = self.traits(user=True).keys()
        pdict = {}
        for p in params:
            pdict[p] = getattr(self, p)
        return pdict

    # Traits notifications and properties

    def _get_monitor(self): 
        return bool(self.t % self.monitor_dt < 0.5001*self.dt)
    
    def _get_t(self):
        return self.ts.t
    
    def _pause_changed(self, paused):
        if not paused:
            self.advance()
    
    def _out_default(self):
        return Messager(title=self.label + ' Simulation')
    
    # Private support methods
    
    def _run_trial(self):
        """Run a trial by wrapping subclass timestep method
        """
        while self.ts:
            
            # Enable pausing functionality
            if self.pause: 
                break
            
            # Subclass-provided timestep computation
            self.run_timestep()
            
            # Progress messages
            self._handle_monitor_msg()        
            
            # Store current state and advance the timeseries
            self.ts.advance()

    def _handle_monitor_msg(self):
        """Handle progress monitoring messages
        """
        if self.monitoring and self.monitor:
            self.out('t = %.2fsec (%.1f%%)'%(self.t, 100*self.t/self.T), 
                notification='monitor')
    
    def _renew_timeseries(self):
        """Construct a fresh data tracker
        """
        self.ts = TimeSeries(self, dt=self.dt, T=self.T)
        if self.monitor_dt < self.T / 10.0:
            self.monitor_dt = int(self.T) / 10.0
    
    def _store_trial(self):
        try:
            self.ts.finish()
        except Exception, e:
            self.out(MSG['SAVEFAIL'][0]%e.message, MSG['SAVEFAIL'][1], 
                MSG['SAVEFAIL'][2])
        else:
            if self.ts.tracking == [] and self.trial_result != '':
                self.results.append(getattr(self, self.trial_result))
                self.out('Trial %d data (\'%s\') saved!'%(self.trial, 
                    self.trial_result), 'Trial Saved')
            else:
                self.results.append(self.ts.post_mortem())
                self.out('Trial %d data saved!'%self.trial, 'Trial Saved')
    
    def _create_datapath(self, dpath):
        if not os.path.isdir(dpath):
            if os.path.isfile(dpath):
                self.out(MSG['DPFILE'][0]%dpath, MSG['DPFILE'][1], 
                    MSG['DPFILE'][2])
                return False
            else:
                try:
                    os.makedirs(dpath)
                except OSError:
                    self.out(MSG['DPFAIL'][0]%dpath, MSG['DPFAIL'][1], 
                        MSG['DPFAIL'][2])
                    return False
                else:
                    self.out(MSG['DPCREATE'][0]%dpath, MSG['DPCREATE'][1], 
                        MSG['DPCREATE'][2])
        return True

    def _write_info_file(self, fd, stem):
        from numpy import asarray
        tstamp = time.strftime(TS_FMT)
        div = '='*70
        info = [div, 'Simulation Information File', div]
        div = '-'*70
        info += ['', 'Model subclass    : ' + self.__class__.__name__]
        info += ['Subclass module   : ' + self.__module__, '']
        info += ['-> Instantiated : ' + self.timestamp, 
            '-> Saved        : ' + tstamp, '']
        if not self.success:
            info += ['', '  *** NOTE: This simulation encountered errors! ***', '']
        info += [div, self.__class__.__name__+' Simulation:', div]
        info += ['', div, self.__class__.__name__+' Parameters:', div]
        info += [' * ' + str(k).ljust(15) + '= ' + str(getattr(self, k)) 
            for k in self.traits(user=True)]
        info += ['', div, 'Collected Time-Series:', div]
        info += [' * ' + str(self.ts).split('\n')[0]]
        info += [' * Progress: trial %d @ %.2f seconds'%(self.trial, self.t)]
        if self.success:
            info += [' * Completed successfully!']
        else:
            info += [' * Errors were detected, check results!']
        info += ['', ' * Tracked Variables:']
        info += ['\t- %s : %s'%(k, str(asarray(getattr(self, k)).shape)) 
            for k in self.traits(track=True)]
        info += ['', div, 'Results Data Archive:', div]
        info += ['Post-mortem data:\n\t' + stem + 'tar.gz']
        info += ['', 'Access data by:']
        info += ['>>> from core.api import ModelPostMortem']
        info += ['>>> pm = ModelPostMortem.fromfile(\'%s\')'%(stem + 'tar.gz')]
        info += ['>>> print pm.tracking']
        info += ['', div, self.__class__.__name__ + ' Docstring:', div]
        info += ['', '\n'.join([s.strip() for s in self.__doc__.split('\n')])]
        fd.write('\n'.join(info))
        fd.close()
