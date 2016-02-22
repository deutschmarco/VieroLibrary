import os
import numpy as np
from VieroLibrary.readcol import readcol
from VieroLibrary.thumbnail_stack import stack_radec



def GetThumbnailRaDec(catalog_in):
    ''' Read in an RA/DEC-only catalog and return the list of coordinates '''
    list_name = catalog_in
    
    #PUT DATA INTO CUBE
    n_sources_max=90000l
    cube = np.zeros([n_sources_max, 2]) # initialize the thumbnail cube
    if os.path.getsize(list_name) > 0:
        ra_to_stack, dec_to_stack = readcol(list_name,fsep=',',twod=False)
        nsources_list=len(ra_to_stack)
        if nsources_list > n_sources_max: 
            raise ValueError('too many sources in catalog: set a larger n_sources_max.')
                
    radec_positions = [ra_to_stack,dec_to_stack]

    print 'Thumbnail-stacked %d sources.' % len(ra_to_stack)
    
    return radec_positions



def syn_kern(map_in, hd_in, catalog_in, side_length):
    ''' Return an empirical, peak-normalized beam profile by stacking on the bright sources of a given map '''
    
    radec_positions = GetThumbnailRaDec(catalog_in)
    cutout = stack_radec(map_in, hd_in, radec_positions, side=side_length)[0]
    cutout_norm = cutout / max(cutout.flatten())
    
    return cutout_norm