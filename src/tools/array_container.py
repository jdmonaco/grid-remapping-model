#encoding: utf-8
"""
array_container -- Basic functionality for classes with large array attributes

Primarily provides an efficient save/load mechanism for python objects containing
large numerical arrays.

Created by Joe Monaco on 2008-06-25.
Copyright (c) 2008 Columbia University. All rights reserved. 
"""

# Library imports
from os import path, unlink, getcwd, chdir, rmdir, makedirs
from sys import stderr
import cPickle, tarfile
import numpy as N

# Traits type
try:
    from enthought.traits.api import HasTraits
except ImportError:
    TRAITS_AVAILABLE = False
else:
    TRAITS_AVAILABLE = True

# Constants for array container save-file format
ARRAYS_FN = 'arrays.npz'
PICKLE_FN = 'object.pickle'


class ArrayContainer(object):
    
    """
    Base class provides some useful methods for subclasses that contain
    large numerical arrays as data attributes.
    
    Available methods:
    tofile -- instance methods that saves an array container object and
        its numerical arrays to an efficient archive data file
    fromfile -- class method that loads a previously saved array 
        container object from an array container archive
    """
    
    def tofile(self, filename):
        """Store array container object in a compressed archive file specified 
        by the filename argument (relative to current directory).
        """
        # Validate archive filename
        if not filename.endswith('.tar.gz'):
            if filename[-1] != '.':
                filename += '.'
            filename += 'tar.gz'
        if path.exists(filename):
            stderr.write('Error: Save file %s already exists'%filename)
            return
            
        # Create a temp directory for component files
        file_path = path.abspath(filename)
        data_path = path.split(file_path)[0]
        tmpdir = path.join(data_path, '_tmp_%06d'%N.random.randint(1000000))
        while path.exists(tmpdir):
            tmpdir = path.join(data_path, '_tmp_%06d'%N.random.randint(1000000))
        
        # Change to new temp directory
        cwd = getcwd()
        makedirs(tmpdir)
        chdir(tmpdir)
        
        # Get attributes list for this object
        if TRAITS_AVAILABLE and isinstance(self, HasTraits):
            attr_list = filter(
                lambda x: not (x.startswith('trait') or x.startswith('_')), 
                    self.traits().keys())
        else:
            attr_list = filter(lambda x: not x.startswith('_'), dir(self))

        # Create dict of numpy arrays in this ratemap object
        map_arrays = {}
        for attr in attr_list:
            value = getattr(self, attr)
            if type(value) is N.ndarray:
                map_arrays[attr] = value
                setattr(self, attr, N.array([])) # empty array for pickling
                
        # Save all numpy arrays to NPZ file
        try:
            N.savez(ARRAYS_FN, **map_arrays)
        except IOError, e:
            raise IOError, 'failed to write arrays file'
        
        # Pickle the rest of this object
        try:
            fd = file(PICKLE_FN, 'w')
            cPickle.dump(self, fd)
        except IOError:
            raise IOError, 'failed to write pickle file'
        finally:
            fd.close()
        
        # Archive the array and pickle files
        tar = tarfile.open(name=file_path, mode='w:gz')
        tar.add(ARRAYS_FN)
        tar.add(PICKLE_FN)
        tar.close()
        unlink(ARRAYS_FN)
        unlink(PICKLE_FN)
        
        # Repopulate the object with its array data
        for arr in map_arrays:
            setattr(self, arr, map_arrays[arr])        
        
        # Remove the temp directory
        chdir(cwd)
        try:
            rmdir(tmpdir)
        except IOError:
            raise IOError, 'could not remove tempdir:\n%s'%tmpdir

    @classmethod
    def fromfile(cls, filename):
        """Retrieve a stored array container object from the data stored in the 
        file specified by filename.
        """        
        # Validate filename
        if not path.isfile(filename):
            raise TypeError, 'bad file name for saved data'
            
        # Get full paths and temp dir path
        file_path = path.abspath(filename)
        data_path = path.split(file_path)[0]
        tmpdir = path.join(data_path, '_tmp_%06d'%N.random.randint(1000000))
        while path.exists(tmpdir):
            tmpdir = path.join(data_path, '_tmp_%06d'%N.random.randint(1000000))
        
        # Create and change to temp dir
        cwd = getcwd()
        makedirs(tmpdir)
        chdir(tmpdir)
        
        # Extract the tar/gzip archive to the temp dir
        try:
            tar = tarfile.open(name=file_path, mode='r:*')
        except:
            raise ValueError, 'failed to open archive file for reading'
        else:
            tar.extractall(path=tmpdir)
            tar.close()
        
        # Verify that the right files exist
        if not (path.exists(ARRAYS_FN) and path.exists(PICKLE_FN)):
            raise ValueError, 'invalid or corrupted archive file'
        
        # Load the array and object data
        array_data = N.load(ARRAYS_FN)
        map_object = N.load(PICKLE_FN)
        unlink(ARRAYS_FN)
        unlink(PICKLE_FN)
        
        # Repopulate the object with its array data
        for arr in array_data.files:
            setattr(map_object, arr, array_data[arr])

        # Clean up and return loaded object
        chdir(cwd)
        try:
            rmdir(tmpdir)
        except IOError:
            raise IOError, 'could not remove tempdir:\n%s'%tmpdir
        return map_object


if TRAITS_AVAILABLE:

    class TraitedArrayContainer(HasTraits, ArrayContainer):
        pass
    

def archive_load(filename):
    """Return the archived ArrayContainer object from the specified archive"""
    return ArrayContainer.fromfile(filename)
