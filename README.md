# Modular Remapping Model

This repository contains a Python package of the source code for the model described in this paper:

#### [Modular Realignment of Entorhinal Grid Cell Activity as a Basis for Hippocampal Remapping](http://www.jneurosci.org/content/31/25/9414)

Joseph D. Monaco [1] and L. F. Abbott [2]

[1] Zanvyl Krieger Mind/Brain Institute, Department of Neuroscience, Johns Hopkins University, Baltimore, MD, USA; [2] Departments of Neuroscience, and Physiology and Cellular Biophysics, Columbia University College of Physicians and Surgeons, New York, New York 10032

Hippocampal place fields, the local regions of activity recorded from place cells in exploring rodents, can undergo
large changes in relative location during remapping. This process would appear to require some form of modulated global
input. Grid-cell responses recorded from layer II of medial entorhinal cortex in rats have been observed to realign
concurrently with hippocampal remapping, making them a candidate input source. However, this realignment occurs
coherently across colocalized ensembles of grid cells (Fyhn et al., 2007). The hypothesized entorhinal contribution to
remapping depends on whether this coherence extends to all grid cells, which is currently unknown. We study whether
dividing grid cells into small numbers of independently realigning modules can both account for this localized
coherence and allow for hippocampal remapping. To do this, we construct a model in which place-cell responses arise
from network competition mediated by global inhibition. We show that these simulated responses approximate the sparsity
and spatial specificity of hippocampal activity while fully representing a virtual environment without learning. Place
field locations and the set of active place cells in one environment can be independently rearranged by changes to the
underlying grid-cell inputs. We introduce new measures of remapping to assess the effectiveness of grid-cell modularity
and to compare shift realignments with other geometric transformations of grid-cell responses. Complete hippocampal
remapping is possible with a small number of shifting grid modules, indicating that entorhinal realignment may be able
to generate place-field randomization despite substantial coherence.

## Installation

You will need a scientific Python distribution such as [Enthought Canopy](https://www.enthought.com/products/canopy/)
or [Continuum Anaconda](https://store.continuum.io/cshop/anaconda/). In the cloned repository, you can run the standard
distutils installation with something like ``python setup.py install``. The model can then be run interactively in an
IPython session.

## Libraries

Here is a brief description of the main modules and classes:

### Top-level Modules

- ``dmec``
    - ``GridCollection``: grid cell population model
- ``place_network``
    - ``PlaceNetwork``: model simulation class
    - ``PlaceNetworkStd``: model simulation class, search-optimized parameters
- ``place_network_ui``
    - ``PlaceNetworkUI``: Chaco graphical frontend for model simulation
- ``placemap``
    - ``PlaceMap``: spatial map class that computes place fields
- ``placemap_viewer``
    - ``PlaceMapViewer``: Chaco graphical interface for PlaceMap objects
- ``ratemap``
    - ``CheckeredRatemap``: PlaceMap subclass for rasterized simulation output
- ``stage``
    - ``StagingMap``: simple handler for defining and indexing the environment
- ``trajectories``
    - ``RandomWalk``: naturalistic random walk trajectory definition
    - ``BipartiteRaster``: checkered rasterization defintion
    
### Subpackages

- ``core`` 
    - Base classes for models, analyses, parameter searches, and time-series data
- ``analysis``
    - ``altmodels``: extensions to inhibitory model
        - ``ModelComparison``: analysis class for running model extensions
    - ``compare``: 
        - ``compare_AB``: function that computes remapping measures
    - ``map_funcs``: functions operating on spatial maps
    - ``movie``: 
        - ``SweepMovie``: analysis class for creating remapping videos
    - ``point``: 
        - ``PointSample``: analysis class for gathering statistics
    - ``realign``: 
        - ``RealignmentSweep``: analysis class for remapping sweeps
    - ``scan``: 
        - ``MultiNetworkScan``: analysis class for sampling parameter sweeps
    - ``search``: 
        - ``PlaceNetworkSearch``: model parameter search definition
    - ``sweep``: 
        - ``SingleNetworkSweep``: two-dimensional parameter sweeps
    - ``two_rooms``: 
        - ``SmoothRemap``: analysis class for progressive remapping simulations
        - ``SampleRemap``: analysis class for random sampling of remapping        
- ``tools``
    - A collection of scientific and utility support functions
    
Some of the ``analysis`` classes farm simulations out to IPython ipengine instances running on your machine. You must
first start them in another terminal:

      $ ipcluster local -n C

Set ``C`` to the number of cores available on your machine.
    

### Example Usage

You can run the model itself, specifying various parameters, or you can run pre-cooked analyses that were used as the
basis of figures in the paper.

#### Running the model

Start IPython in ``-pylab`` mode:

    $ ipython -pylab

Then, import the libraries and create a model instance:

    In [0]: from grid_remap import *
    In [1]: model = PlaceNetworkStd()

To see all the user-settable parameters, you can print the model:

    In [2]: print model
    PlaceNetworkStd(Model) object
    --------------------------------
    Parameters:
            C_W : 0.33000000000000002
            EC : None
            J0 : 45.0
            N_CA : 500
            done : False
            dwell_factor : 5.0
            monitoring : True
            mu_W : 0.5
            pause : False
            phi_lambda : 0.040000000000000001
            phi_sigma : 0.02
            refresh_orientation : False
            refresh_phase : False
            refresh_traj : False
            refresh_weights : True
            tau_r : 0.050000000000000003
            traj_type : 'checker'

Important model parameter definitions:

    C_W            feedforward connectivity 
    EC             the GridCollection to use as input
    J0             gain of global inhibition
    N_CA           the number of output units; each receives input from 
                        C_W*N_EC grid cells
    dwell_factor   multiple of tau_r that defines raster pixel dwell time
    mu_W           average weight of feedforward synapses
    phi_lambda     nonlinearity threshold
    phi_sigma      nonlinearity smoothness (gain)
    refresh_*      orientation/phase reset per trial; new random weight
                        matrix per trial
    tau_r          time constant of place-unit integration

Parameters can be changed by passing them as keyword arguments to the
constructor. To simulate only 100 place units, you would call
``PlaceNetworkstd(N_CA=100)``.

Run the simulation:

    In [3]: model.advance()

Look at the tracked data:

    In [4]: pmap = CheckeredRatemap(model)


#### Running analyses

To run the figure analyses, you simply create an analysis object and run it by
calling it with analysis parameters. To run progressive realignment
experiments using the ``RealignmentSweep`` analysis class, you would run:

    In [20]: fig = RealignmentSweep(desc='test')
    In [21]: fig.collect_data?

The first command creates an analysis object with the description 'test'. The
second command (with the ``?``) tells IPython to print out meta-data about the
``collect_data`` method. This is the method that actually performs the analysis
when you call the object, so this tells you the available parameters along with
their descriptions. We could run the analysis with modularity on the y-axis:

    In [22]: fig(y_type='modules')

This performs the simulations, collects data for the figures, and stores data,
statistics, and an *analysis.log* file in the analysis directory. When that
completes, you can bring up the resulting figure and save it:

    In [23]: fig.view()
    In [24]: fig.save_plots()

Running the ``view`` method renders the figures, outputs RGB image files, and
saves a *figure.log* file in the analysis directory. Some of the figures have
parameter arguments to change the figure. You will have to use the
``create_plots`` method, as this is what the ``view`` method actually calls. To 
see the figure parameters and make changes:

    In [25]: fig.create_plots?
    In [26]: fig.create_plots(...)

The same process can be used for the other figure analysis classes. You can
create your own analyses by subclassing from ``core.analysis.BaseAnalysis`` and
implementing the ``collect_data`` and ``create_plots`` methods.
