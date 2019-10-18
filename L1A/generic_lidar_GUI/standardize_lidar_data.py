# This code will take given lidar data and translate it into a
# standardized structure, standard_lidar_object.

import numpy as np
import pdb
import matplotlib.pyplot as plt
import datetime as DT
# Custom libraries
from read_routines import *
from generic_lidar_GUI_initializations import *

class standard_lidar_object():
    """ This object represents an binnned, elastic backscatter lidar """
    
    # NOTES:
    #
    # All the ONA stuff was added with CAMAL in mind, but it should be 
    # able to run off of any array of off-nadir angles.
    
    # Attribute definitions...
    # nr --> The number of records/profiles (int)
    # nb --> The number of range bins (int)
    # nc --> The number of channels (int)
    # dz --> Vertical bin size in meters (float)
    # dx --> Horizontal distance between profiles (float)
    # nr0 --> The raw resolution number of records (int)
    # nshots --> The number of shots accumulated per raw record (int)
    # nhori --> The number of profiles to average
    # pdir --> Is lidar pointing "Up" or "Down?" (string)
    # rec_arr --> Numpy array of record number (np.float32 or int)
    # time_arr --> Numpy array of time (datetime object)
    # bin_arr --> Numpy array of range bin bumber (np.float32 or int)
    # alt_arr --> Numpy array of bin altitude (np.float32)
    # platform_alt_arr --> Numpy array of platform (ie aircraft) altitude
    # counts --> Numpy array of raw signal/counts (np.float32, nr x nc x nb)
    # energy --> Numpy array of energy in milli-Joules (x nr)
    # ingested --> Bool that indicates whether data have been loaded
    # avgd_mask --> List or numpy array (x nc) that indicates which channles have been averaged
    # samp_chan_map --> Maps the samp_chan array to the unaveraged data
    # ONA_arr --> Numpy array of off-nadir angles
    # ONA_bw --> Bin width for histogram of ONA_arr
    # ONA_uplim --> Set upper limit for ONA histogram binning 
    # ONA_uplim --> Set lower limit for ONA histogram binning
    # ONA_bin_centers --> Exactly as it sounds
    # ONA_internal_bw --> Spacing between bin centers. Used internally.
    # samp_chan --> Convenient way for GUI to store and maipulate chosen data
    
    def __init__(self, nr, nb, nc, dz, dx, nr0, nshots, nhori, pdir,
                       rec_arr, time_arr, bin_arr, alt_arr, platform_alt_arr,
                       counts, energy,
                       ingested, avgd_mask, samp_chan_map,
                       ONA_arr, ONA_bw, ONA_uplim, ONA_lowlim, 
                       ONA_bin_centers, ONA_internal_bw,
                       samp_chan
                ):
        """ Initialize attributes of this class """
        self.nr = nr
        self.nb = nb
        self.nc = nc
        self.dz = dz
        self.dx = dx
        self.nr0 = nr0
        self.nshots = nshots
        self.nhori = nhori
        self.pdir = pdir
        self.rec_arr = rec_arr
        self.time_arr = time_arr
        self.bin_arr = bin_arr
        self.alt_arr = alt_arr
        self.platform_alt_arr = platform_alt_arr
        self.counts = counts
        self.energy = energy
        self.ingested = ingested
        self.avgd_mask = avgd_mask
        self.samp_chan_map = samp_chan_map
        self.ONA_arr = ONA_arr
        self.ONA_bw = ONA_bw
        self.ONA_uplim = ONA_uplim
        self.ONA_lowlim = ONA_lowlim
        self.ONA_bin_centers = ONA_bin_centers
        self.ONA_internal_bw = ONA_internal_bw
        self.samp_chan = samp_chan
        
    def update_avgd_mask(self,sel_chan):
        """ Update mask that tells GUI which channels have been averaged """
        self.avgd_mask[sel_chan-1] = True
        
    def bin_ONAs(self):
        """ Put the off-nadir angles into histogram bins """
        
        nhbins = int( (self.ONA_uplim - self.ONA_lowlim)/self.ONA_bw )
        dens, edges = np.histogram( self.ONA_arr, bins=nhbins,
                      range=(self.ONA_lowlim, self.ONA_uplim) )
        delta = float(self.ONA_uplim - self.ONA_lowlim) / float(nhbins)
        center = edges[:len(dens)] + delta
        mask = dens != 0
        print("This many bins have a feq of at least 1:", dens[mask].shape)
        print("These are the bins: ",center[mask])
        self.ONA_bin_centers = center[mask]
        self.ONA_internal_bw = delta
        
    def apply_ONA_mask(self,bin_indx,opt):
        """ Method to apply mask to data based on off-nadir angle (ONA). 
            Do not call this method before you load data.
            If off-nadir angles are not in your lidar data, set 
            self.ONA_bw to None.
        """
        
        # This maps the full-length dataset to samp_chan
        samp_chan_map = np.arange(0,self.nr,dtype=np.uint32)
        self.samp_chan_map = samp_chan_map # could be overwritten below
        
        # Indicator that lidar data doesn't have ONAs
        if (self.ONA_bw is None): return
        # Go no farther in no angle is selected
        if (opt == 'no mask'): return
        # If you make it to this line, a mask is going to be applied
        print("An angle mask is being applied.")
        
        cs = self.hist_bin_centers[bin_indx]
        lim1 = cs - self.ONA_internal_bw
        lim2 = cs + self.ONA_internal_bw
        print( 'limit1 = '+str(lim1).strip()+' lim2 = '+str(lim2).strip() )
        
        # Make an array mask to mask "where" angle falls within bin.
        # First reduce array based on averaging, ang_mask then will be
        # in "averaged space." 
        ONA_reduced = self.ONA_arr[::self.nhori]
        ONA_reduced = ONA_reduced[0:self.nr]
        ang_mask = ( (ONA_reduced >= lim1) & 
             (ONA_reduced <= lim2)  )
             
        if opt == 'nogaps':
            
            # Here we make a sample (samp_chan) with no gaps.
            # Remember, ang_mask will be in averaged space.
            self.samp_chan = self.samp_chan[:,:][samp_chan_map[ang_mask]]
            print('The shape of samp_chan is ',self.samp_chan.shape)
            # Finalize the samp_chan_map. The numbers within this array 
            # should correspond to full-ez dataset indicies.
            self.samp_chan_map = samp_chan_map[ang_mask]
            
        elif opt == 'gaps':
            
            # Here we make a sample (samp chan) with gaps
            self.samp_chan[:,:][samp_chan_map[(ang_mask==False)]] = 0
            
        # Wrap up this method, overwrite self.nr
        samp_chan_shape = self.samp_chan.shape
        print('The shape of samp_chan is ', samp_chan_shape)
        self.nr = samp_chan_shape[0]
        print('Length of samp_chan: ',samp_chan_shape[0])
        
        
def select_and_standardize_data(selected_files):
    """ Select which lidar data to use, standardize it, and return it
        in form of a standard_lidar_object.
        
        Add new code to this function for each new instrument.
        
        PLEASE NOTE: This function relies on an initializations import!
    """
    
    # INPUTS:
    # selected_files --> A list of numbers indicating which files to read.
    # --> All other inputs come from initializations at top of this module.
    
    # Instantiate object
    SLO = standard_lidar_object(None,None,None,None,None,None,None,None,
       None,None,None,None,None,None,None,None,None,None,None,None,None,
       None,None,None,None,None,)
    
    if instrument_name == 'Roscoe':
        
        file_list = 'Roscoe_files_for_GUI.txt'
        create_a_file_list(file_list, search_str, raw_dir)
        MCS_data = read_selected_mcs_files_r(file_list,rep_rate,file_len_secs,selected_files)
        
        # The number of seconds needed to convert from the instrument's Unix time
        # to UTC. Typically either 5 hours (cold season) or 4 hours (warm season).
        secs_btwn_instr_UnixT_and_UTC = 18000
        
        SLO.nr = MCS_data['counts'].shape[0]
        SLO.nb = MCS_data['counts'].shape[2]
        SLO.nc = MCS_data['counts'].shape[1]
        SLO.dz = vrZ # or (MCS_data['meta']['binwid'][0] * c) / 2.0
        SLO.dx = dx
        SLO.nr0 = SLO.nr
        SLO.nshots = nshots
        SLO.nhori = nhori
        SLO.pdir = pointing_dir
        SLO.rec_arr = np.arange(0,SLO.nr)
        SLO.time_arr = np.zeros(SLO.nr,dtype=DT.datetime)
        for i in range(0,SLO.nr): 
            SLO.time_arr[i] = DT.datetime.fromtimestamp(MCS_data['meta']['CCSDS']['UnixT'][i]) + \
                              DT.timedelta(seconds=secs_btwn_instr_UnixT_and_UTC)
        SLO.bin_arr = np.arange(0,SLO.nb)
        if SLO.pdir == 'Down':
            SLO.alt_arr = np.arange(SLO.nb,0,-1) * SLO.dz - 2000.0
        elif SLO.pdir == 'Up':
            SLO.alt_arr = SLO.bin_arr * SLO.dz + 0
        else:
            print('pdir can either be "Up" or "Down"')
            pdb.set_trace()
        SLO.platform_alt_arr = np.zeros(SLO.nr) + 20000.0
        SLO.counts = MCS_data['counts'] # should be shape (nr, nc, nb)
        SLO.energy = None#np.zeros(SLO.nr)
        SLO.ingested = True
        SLO.avgd_mask = [False for x in range(0,SLO.nc)]
        SLO.samp_chan_map = None
        SLO.ONA_arr = None
        SLO.ONA_bw = None
        SLO.ONA_uplim = None
        SLO.ONA_lowlim = None
        SLO.ONA_bin_centers = None
        
    elif instrument_name == 'CPL':
        
        bad_cls_nav_time_value = DT.datetime(1970,1,1,0,0,0)
        CLS_data = read_entire_cls_dataset(file_len_recs,raw_dir,nbins,flt_date,bad_cls_nav_time_value,Fcontrol=None)
        SLO.nr = CLS_data['counts'].shape[0]
        SLO.nb = CLS_data['counts'].shape[2]
        SLO.nc = CLS_data['counts'].shape[1]
        SLO.dz = vrZ
        SLO.dx = dx
        SLO.nr0 = nr
        SLO.nshots = nshots
        SLO.nhori = nhori
        SLO.pdir = pointing_dir
        SLO.rec_arr = np.arange(0,SLO.nr)
        SLO.timearr = CLS_data['meta']['Nav']['UTC_Time']
        SLO.bin_arr = np.arange(0,SLO.nb)
        z0 = np.mean(CLS_data['meta']['Nav']['GPS_alt'])
        if SLO.pdir == "Up":
            SLO.alt_arr = z0 + np.arange(0, SLO.nb) * SLO.dz
        elif SLO.pdir == "Down":
            SLO.alt_arr = z0 - np.arange(0, SLO.nb) * SLO.dz
        else:
            SLO.alt_arr = z0 - np.arange(0, SLO.nb) * SLO.dz   
        SLO.platform_alt_arr = z0
        SLO.counts = CLS_data['counts']
        SLO.energy = None
        SLO.ingested = True
        SLO.avgd_mask = [False for x in range(0,SLO.nc)]
        SLO.samp_chan_map = None
        SLO.ONA_arr = None
        SLO.ONA_bw = None
        SLO.ONA_uplim = None
        SLO.ONA_lowlim = None
        SLO.ONA_bin_centers = None        
        
        
    return SLO
        
        
