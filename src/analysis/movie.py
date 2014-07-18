#encoding: utf-8
"""
grid.analysis.movie -- Analysis subclass that creates parameter sweep movies
    of how the spatial map from a single network/environment pair is affected
    by changes to that parameter.
    
Written by Joe Monaco, 10/14/2008
Copyright (c) 2008 Columbia Unversity. All Rights Reserved.  
"""

# Library imports
import numpy as N, os
from pylab import figure
from matplotlib import cm

# Package imports
from .. import PlaceNetworkStd, CheckeredRatemap, GridCollection
from ..core.analysis import AbstractAnalysis
from ..tools.images import array_to_image, tile2D
from .map_funcs import spatial_corr, get_tuned_weights


class SweepMovie(AbstractAnalysis):
    
    """
    Using a single network/environment pair, run a specified parameter sweep 
    and save a series of spatial map images as frames of a movie
    
    The output of the analysis call (method collect_data) is several movie
    files in MPEG format in the data directory.
    """
    
    label = 'movie sweep'
    
    def collect_data(self, nframes=24, param='J0', bounds=(0.5, 6), nmodules=1, 
        freq_modules=False, shift=False, rotate=False, delta_phi=None, delta_psi=None, 
        ellipticity=False, rescaling=False, start_tuned=False, alpha=0.5, 
        no_frames=False, no_encode=False, **kwargs):
        """
        Run the sweep, storing spatial map matrices as frames of movies
        
        Keyword arguments:
        nframes -- the number of points (frames) traversing the specified
            parameter range (default 24)
        param -- string name of a model parameter to sweep (special cases:
            'realign' to simulate a realignment movie of a default-parameter 
            network (bounds are set to (0,1)); 'tuning' to simulate afferent
            weight matrix tuning, bounds specifies alpha bounds)
        bounds -- (min, max) tuple specifying the parameter range
        
        If param == 'realign', then these parameters can be specified:
        nmodules -- number of realignment modules if param == 'realign'
        freq_modules -- whether modules are spatial frequency partitions
        shift -- whether to perform phase realignment 
        rotate -- whether to perform orientation realignment
        delta_phi/delta_psi -- phase/orientation displacements for the modules;
            if nmodules > 1, must be a list of tuples (phi) or angles (psi)
        ellipticity -- whether to enable elliptical modulation of grids
        rescaling -- whether to enable rescaling (enlargement) of grids    
        start_tuned -- start realignment from a tuned network
        alpha -- if start_tuned is True, then alpha specifies tuning
        """
        
        # Store bounds and scan parameter
        self.results['bounds'] = N.array(bounds)
        self.results['param'] = param
        
        # Load cortex
        self.out('Creating grid collection object...')
        EC = GridCollection()
        
        ############################################
        # Realignment setup
        ############################################
        if param == 'realign':
            # Split the cortical population into modules
            self.results['nmodules'] = nmodules = int(nmodules)
            modules = EC.get_modules(nmodules, freq_sort=freq_modules)
            
            # Reset the bounds for full range of realignment
            bounds = (0.0, 1.0)
            self.results['bounds'] = N.array(bounds)

            # Set phase realignment
            d_phi_list = []
            if shift:
                if delta_phi is not None:
                    if nmodules > 1:
                        d_phi_list = list(delta_phi)
                        if len(d_phi_list) != nmodules:
                            raise ValueError, 'invalid shift specification'
                    elif nmodules == 1:
                        d_phi_list = [delta_phi]
                    else:
                        raise ValueError, 'invalid number of modules'
                else:
                    grid_scale = None
                    for i in xrange(nmodules):
                        if freq_modules:
                            grid_scale = EC.spacing[modules[i]].max()
                        d_phi_list.append(
                            GridCollection.get_delta_phi(scale=grid_scale))
            d_phi = N.array(d_phi_list, 'd')
            self.results['delta_phi'] = d_phi

            # Set orientation realignment
            d_psi_list = []
            if rotate:
                if delta_psi is not None:
                    if nmodules > 1:
                        d_psi_list = list(delta_psi)
                        if len(delta_psi) != nmodules:
                            raise ValueError, 'invalid rotate specification'
                    elif nmodules == 1:
                        d_psi_list = [delta_psi]
                    else:
                        raise ValueError, 'invalid number of modules'
                else:
                    for i in xrange(nmodules):
                        d_psi_list.append(GridCollection.get_delta_psi())
            d_psi = N.array(d_psi_list, 'd')
            self.results['delta_psi'] = d_psi
            
            # Set ellipticity of response
            ell_mag_list = []
            ell_angle_list = []
            if ellipticity:
                EC.ellipticity = ellipticity
                for m in xrange(nmodules):
                    ell_mag_list.append(GridCollection.get_ellipticity())
                    ell_angle_list.append(GridCollection.get_elliptic_angle())
            ell_mags = N.array(ell_mag_list)
            ell_angles = N.array(ell_angle_list)
            
            # Set rescaling factors for grid response
            zoom_list = []
            if rescaling:
                EC.zoom = rescaling
                for m in xrange(nmodules):
                    zoom_list.append(GridCollection.get_zoom_scale())
            zoom_scales = N.array(zoom_list)
        ############################################
        # End realignment setup
        ############################################
                
        # Set default model parameters
        pdict = dict(   EC=EC, 
                        desc='movie', 
                        projdir=self.datadir, 
                        refresh_weights=False, 
                        refresh_orientation=False,
                        refresh_phase=False,
                        refresh_traj=False,
                        monitoring=False)
        pdict.update(kwargs)
        
        # Update with keyword arguments
        if param not in PlaceNetworkStd().traits(user=True).keys():
            if param not in  ('dwell_factor', 'realign', 'tuning'):
                raise ValueError, 'param (%s) is not a valid parameter'%param
        
        # Create the list of sample points to scan
        self.out('Creating %s scan vector from %.2f to %.2f'%((param,)+bounds))
        pts = N.linspace(bounds[0], bounds[1], num=nframes)
        self.results['samples'] = pts
        
        # Initialize arrays for storing frame data
        self.results['stage_repr_map'] = N.empty((nframes, 100, 100), 'd')
        self.results['norms'] = N.empty((nframes, 100, 100), 'd')
        self.results['autocorr'] = N.empty((nframes, 199, 199), 'd')
        self.results['ratemaps'] = N.empty((nframes, 600, 800), 'd')
        self.results['activation'] = N.empty((nframes, 240, 320), 'd')
        
        # Method for running simulations and creating the frame images
        def run_sample_point(i, model):

            # Check for pre-existing data to load
            map_file = os.path.join(model.projdir, 'data', 'frame_%03d.tar.gz'%i)
            if os.path.exists(map_file):
                self.out('Loading found data:\n%s'%map_file)
                ratemap = CheckeredRatemap.fromfile(map_file)
            else:
                # Run the model simulation and save the results
                model.advance()

                # Create ratemap object
                ratemap = CheckeredRatemap(model.post_mortem(trial=1))
                ratemap.compute_coverage()
                ratemap.data = None
                ratemap.tofile(map_file)
            
            # Store ratemap data for movie frames
            if not no_frames:
                self.results['stage_repr_map'][i] = ratemap.stage_repr_map
                self.results['norms'][i] = N.sqrt((ratemap.Map**2).sum(axis=0))
                self.results['autocorr'][i] = spatial_corr(ratemap.Map)
                self.results['ratemaps'][i] = tile2D(ratemap.Map[:48], shape=(6,8), 
                    gridvalue=1)
            
                # Create the activation pattern image
                active_units = N.array([u for u in xrange(300) 
                    if ratemap.fields[u] is not None])
                pattern = N.zeros(300, 'd')
                pattern[::2] = 0.2 # dark gray checkerboard in the background
                pattern[active_units] = 1.0
                pattern = pattern.reshape(20,15).T # Stride over rows
                self.results['activation'][i] = \
                    pattern.repeat(16, axis=0).repeat(16, axis=1)                
            return
            
        # Run through the frames to calculate
        the_model = PlaceNetworkStd(**pdict)
        dpath = os.path.join(self.datadir, 'data')
        fpath = os.path.join(self.datadir, 'frames')
        the_model._create_datapath(dpath)
        the_model._create_datapath(fpath)
        pdict['W'] = the_model.W
        if param == 'tuning' or (param == 'realign' and start_tuned):
            self.out('Computing initial spatial map for tuning...')
            W0 = the_model.W.copy()
            the_model.advance()
            the_map = CheckeredRatemap(the_model)
            the_map.compute_coverage()
            if param == 'realign':
                pdict['W'] = get_tuned_weights(the_map, W0, EC, alpha=alpha)
        for i, p in enumerate(pts):
            if param == 'realign':
                for m, m_ix in enumerate(modules):
                    if shift:
                        EC.shift(d_phi[m] / nframes, mask=m_ix)
                    if rotate:
                        EC.rotate(d_psi[m] / nframes, mask=m_ix)
                    if ellipticity:
                        EC.ell_mag[m_ix] = p * ell_mags[m]
                        EC.ell_angle[m_ix] = p * ell_angles[m]
                    if rescaling:
                        EC.zoom_scale[m_ix] = 1 + p * (zoom_scales[m] - 1)
            elif param == 'tuning':
                pdict['W'] = get_tuned_weights(the_map, W0, EC, alpha=p)
            else:
                pdict[param] = p
            probe_model = PlaceNetworkStd(**pdict)
            self.out('Frame %d: %s = %.2f'%(i, param, p))
            self.execute(run_sample_point, i, probe_model)
        
        # Create a subdirectory for movie frames
        frame_dir = os.path.join(self.datadir, 'frames')
        if not os.path.isdir(frame_dir):
            os.mkdir(frame_dir)
        if no_frames:
            self.out('No frames or movies are being made. Done!')
            return
        
        # Save the frame data out as annotated still frames
        #
        # The ImageMagick tools 'convert' and 'montage' are used to create the
        # individual movie frames based on tools.array_to_image output.
        #
        smap_list_file = file(os.path.join(frame_dir, 'stage_repr.txt'), 'w')
        norms_list_file = file(os.path.join(frame_dir, 'norms.txt'), 'w')
        autocorr_file = file(os.path.join(frame_dir, 'autocorr.txt'), 'w')
        ratemaps_file = file(os.path.join(frame_dir, 'ratemaps.txt'), 'w')
        activation_file = file(os.path.join(frame_dir, 'activation.txt'), 'w')
        label_bgcolor = 'skyblue'
        
        # Some frame-constant data for making movies
        norm_min, norm_max = self.results['norms'].min(), self.results['norms'].max()

        for i in xrange(nframes):

            # Save the stage representation map, border/annotate/resize it
            s_file = os.path.abspath(
                os.path.join(frame_dir, 'stage_repr_%03d.png'%i))
            array_to_image(self.results['stage_repr_map'][i], s_file,
                norm=True, cmap=cm.hot)
            cmd = "convert -bordercolor white -border 27x0 %s %s"%(s_file, s_file)
            os.system(cmd)
            cmd = "montage -geometry +0+0 -background %s "%label_bgcolor
            cmd += "-label '%s = %.2f' "%(param, pts[i])
            cmd += "%s %s"%(s_file, s_file)
            os.system(cmd)
            cmd = "convert -geometry 320x240! %s %s"%(s_file, s_file)
            os.system(cmd)
            cmd = "convert %s %s"%(s_file, s_file[:-3]+'jpg')
            os.system(cmd)
            smap_list_file.write(s_file[:-3]+'jpg' + '\n')
            
            # Save the reduced map, border/annotate/resize it
            n_file = os.path.abspath(
                os.path.join(frame_dir, 'norms_%03d.png'%i))
            array_to_image(self.results['norms'][i], n_file,
                norm=False, cmin=norm_min, cmax=norm_max, cmap=cm.gist_heat)
            cmd = "convert -bordercolor black -border 27x0 %s %s"%(n_file, n_file)
            os.system(cmd)
            cmd = "montage -geometry +0+0 -background %s "%label_bgcolor
            cmd += "-label '%s = %.2f' "%(param, pts[i])
            cmd += "%s %s"%(n_file, n_file)
            os.system(cmd)
            cmd = "convert -geometry 320x240! %s %s"%(n_file, n_file)
            os.system(cmd)
            cmd = "convert %s %s"%(n_file, n_file[:-3]+'jpg')
            os.system(cmd)
            norms_list_file.write(n_file[:-3]+'jpg' + '\n')
            
            # Make a dual-flavor autocorrelation frame
            a_file = os.path.abspath(
                os.path.join(frame_dir, 'autocorr_%03d.png'%i))
            a_png0 = os.path.abspath(
                os.path.join(frame_dir, 'autocorr_%03d_0.png'%i))
            a_png1 = os.path.abspath(
                os.path.join(frame_dir, 'autocorr_%03d_1.png'%i))
            array_to_image(self.results['autocorr'][i], a_png0, 
                norm=True, cmap=cm.jet)
            array_to_image(self.results['autocorr'][i], a_png1, 
                norm=True, cmap=cm.flag)
            cmd = "convert %s %s +append -bordercolor black -border 1x1 %s"%(
                a_png0, a_png1, a_file)
            os.system(cmd)
            cmd = "montage -geometry +0+0 -background %s "%label_bgcolor
            cmd += "-label '%s = %.2f' %s %s"%(param, pts[i], a_file, a_file)
            os.system(cmd)
            cmd = "convert -bordercolor black -border 0x41 -background black "
            cmd += "-extent 400x300 %s %s"%(a_file, a_file)
            os.system(cmd)
            cmd = "convert %s %s"%(a_file, a_file[:-3]+'jpg')
            os.system(cmd)
            cmd = "rm -rf %s %s"%(a_png0, a_png1)
            os.system(cmd)
            autocorr_file.write(a_file[:-3]+'jpg' + '\n')
            
            # Create a tiled grid of raw ratemaps, annotate/resize it
            m_file = os.path.abspath(
                os.path.join(frame_dir, 'ratemaps_%03d.png'%i))
            array_to_image(self.results['ratemaps'][i], m_file,
                norm=False, cmin=0, cmax=1, cmap=cm.jet)
            cmd = "convert -fill white -undercolor '#00000080' -gravity "
            cmd += "southwest -annotate +5+5 '%s = %.2f' "%(param, pts[i])
            cmd += "%s %s"%(m_file, m_file)
            os.system(cmd)
            cmd = "convert -geometry 640x480! %s %s"%(m_file, m_file)
            os.system(cmd)
            cmd = "convert %s %s"%(m_file, m_file[:-3]+'jpg')
            os.system(cmd)
            ratemaps_file.write(m_file[:-3]+'jpg' + '\n')
            
            # Create activation pattern, annotate it
            p_file = os.path.abspath(
                os.path.join(frame_dir, 'activation_%03d.png'%i))
            array_to_image(self.results['activation'][i], p_file,
                norm=False, cmin=0, cmax=1.15, cmap=cm.gray)
            cmd = "convert -fill white -undercolor '#00000080' -gravity "
            cmd += "southwest -annotate +5+5 '%s = %.2f' "%(param, pts[i])
            cmd += "%s %s"%(p_file, p_file)
            os.system(cmd)
            cmd = "convert %s %s"%(p_file, p_file[:-3]+'jpg')
            os.system(cmd)
            activation_file.write(p_file[:-3]+'jpg' + '\n')
            
            self.out('Saved and annotated:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s'%
                (s_file, n_file, a_file, m_file, p_file))

        smap_list_file.close()
        norms_list_file.close()
        ratemaps_file.close()
        autocorr_file.close()
        activation_file.close()
        
        # Create 12fps movies from the saved frame images
        #
        # This mencoder call string is based on the example found here:
        # http://objectmix.com/perl/253395-creating-movie-file-series-\
        # canvas-snapshots.html
        #
        if not no_encode:
            self.out('Encoding movies from frame image files...')
            cwd = os.getcwd()
            os.chdir(self.datadir)
            cmd = "mencoder mf://@frames/%s.txt -msglevel all=1 -mf type=jpeg:fps=12 "
            cmd += "-vf harddup -o %s.mpg -ofps 24 -of lavf -lavfopts "
            cmd += "format=mpg:i_certify_that_my_video_stream_does_not_use_b_frames "
            cmd += "-ovc lavc -lavcopts vcodec=mpeg1video:vmax_b_frames=0:trell:"
            cmd += "vqscale=2:vhq:mv0:cbp -nosound"
            os.system(cmd%(('stage_repr',)*2))
            os.system(cmd%(('norms',)*2))
            os.system(cmd%(('autocorr',)*2))
            os.system(cmd%(('ratemaps',)*2))
            os.system(cmd%(('activation',)*2))
        
            # Remove all JPEG images, leaving the PNG masters
            cmd = "rm -rf frames/*.jpg"
            os.system(cmd)
            os.chdir(cwd)
        
        # Clear movie data from results memory
        del self.results['stage_repr_map']
        del self.results['norms']
        del self.results['autocorr']
        del self.results['ratemaps']
        del self.results['activation']
        
        # Good-bye
        self.out('All done!')
