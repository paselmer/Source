# Define data structures that can change

import numpy as np

from immutable_data_structs import * #data structs that don't change


def define_MSC_structure(nchans,nbins):
    """ This function will define the final MCS structure once the number
        of channels and number of bins are determined.
        
        INPUTS:
        nchans -> the number of channels
        nbins -> the number of bins	    
    """
    
    MCS_struct = np.dtype([ ('meta',MCS_meta_struct), 
                            ('counts',np.uint16,(nchans,nbins)) ])
    
    return MCS_struct


 
def define_CLS_structure(nchans,nbins,CLS_meta_struct):
    """ This function will define the CLS structure once the number of
        channels and number of bins are determined.
	
        INPUTS:
        nchans -> the number of channels
        nbins -> the number of bins
	
	OUTPUT:
	A Python Numpy data structure for raw CPL "CLS" data.
    """
    
    CLS_struct = np.dtype([ ('meta',CLS_meta_struct), 
                            ('counts',np.uint16,(nchans,nbins)),
			    ('Reserved',np.uint8,(44,))      ])

    return CLS_struct
    
    
def define_CLS_decoded_structure(nchans,nbins):
    """ This function will define the CLS structure once the number of
        channels and number of bins are determined.
	
        INPUTS:
        nchans -> the number of channels
        nbins -> the number of bins
	
	OUTPUT:
	A Python Numpy data structure for CPL "CLS" data, which character
	strings converted into usable data types.
    """
    
    CLS_struct = np.dtype([ ('meta',CLS_decoded_meta_struct), 
                            ('counts',np.uint16,(nchans,nbins)),
			    ('Reserved',np.uint8,(44,))      ])

    return CLS_struct
    
    

def define_CLS_meta_struct(nbytes):
    """ This function creates the meta struct for a raw CLS file once
        CLS_raw_nav_struct is defined.
    """ 

    CLS_raw_nav_struct = np.dtype([ ('Full_rec', np.uint8, (nbytes,) ) ])

    CLS_meta_struct = np.dtype([ ('Header',CLS_raw_header_struct), ('Engineering',CLS_raw_engineering_struct),
                             ('Nav',CLS_raw_nav_struct)                                             ])    
    
    return CLS_meta_struct
    
