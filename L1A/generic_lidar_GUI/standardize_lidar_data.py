# This code will take given lidar data and translate it into a
# standardized structure, standard_lidar_object.

from generic_lidar_GUI_initializations import *

class standard_lidar_object():
    """ This object represents an binnned, elastic backscatter lidar """
    
    # Attribute definitions...
    # nr --> The number of records/profiles (int)
    # nb --> The number of range bins (int)
    # nc --> The number of channels (int)
    # dz --> Vertical bin size in meters (float)
    # dx --> Horizontal distance between profiles (float)
    # nr0 --> The raw resolution number of records (int)
    # nshots --> The number of shots accumulated per raw record (int)
    # pdir --> Is lidar pointing "Up" or "Down?" (string)
    # rec_arr --> Numpy array of record number (np.float32 or int)
    # time_arr --> Numpy array of time (datetime object)
    # bin_arr --> Numpy array of range bin bumber (np.float32 or int)
    # alt_arr --> Numpy array of bin altitude (np.float32)
    # platform_alt_arr --> Numpy array of platform (ie aircraft) altitude
    # counts --> Numpy array of raw signal/counts (np.float32, nr x nb x nc)
    
    
    def __init__(self, nr, nb, nc, dz, dx, nr0, nshots, pdir,
                       rec_arr, time_arr, bin_arr, alt_arr, platform_alt_arr,
                       counts, energy,
                       ingested,
                       ONA_arr, ONA_bw,
                )
        """ Initialize attributes of this class """
        self.nr = nr
        self.nb = nb
        self.nc = nc
        self.dz = dz
        self.dx = dx
        self.nr0 = nr0
        self.nshots = nshots
        self.pdir = pdir
        self.rec_arr = rec_arr
        self.time_arr = time_arr
        self.bin_arr = bin_arr
        self.alt_arr = alt_arr
        self.platform_alt_arr = platform_alt_arr
        self.counts = counts
        self.energy = energy
        self.ingested = ingested
                       
    
