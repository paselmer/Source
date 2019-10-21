# This code will take given lidar data and translate it into a
# standardized structure, standard_lidar_object.

# Routines from ../lidar.py copied into this module because ../lidar.py
# relies on L1A_initializations.py

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
        
        
def select_and_standardize_data(file_ctrl):
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
        MCS_data = read_selected_mcs_files_r(file_list,rep_rate,file_len_secs,file_ctrl.sel_file_list)
        
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
        CLS_data = read_entire_cls_dataset(file_len_recs,raw_dir,nbins,flt_date,bad_cls_nav_time_value,file_ctrl)
        SLO.nr = CLS_data['counts'].shape[0]
        SLO.nb = CLS_data['counts'].shape[2]
        SLO.nc = CLS_data['counts'].shape[1]
        SLO.dz = vrZ
        SLO.dx = dx
        SLO.nr0 = SLO.nr
        SLO.nshots = nshots
        SLO.nhori = nhori
        SLO.pdir = pointing_dir
        SLO.rec_arr = np.arange(0,SLO.nr)
        SLO.timearr = CLS_data['meta']['Nav']['UTC_Time']
        SLO.bin_arr = np.arange(0,SLO.nb)
        z0 = np.mean(CLS_data['meta']['Nav']['GPS_Altitude'])
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
        

def compute_solar_background(rc,bg_st_bin,bg_ed_bin):
    """This code will compute the solar background of raw counts
       given a 2D counts array and a background region. It will
       return the background as a 2D array, where values are constant
       in the altitude dimension. The altitude dimension should be
       the second dimension (index 1).
    """
    
    # NOTE: This function assumes array of profs x bins [X x Y]
    bg_reg_subset = rc[:,bg_st_bin:bg_ed_bin]
    bg = np.mean(bg_reg_subset, axis=1)
    return bg

    
def correct_raw_counts(rc,E,r,nr,nb,bg_st_bin,bg_ed_bin,opt):
    """This function will produce NRB from raw photon counts"""
    
    # INPUTS:
    #         - rc        = raw counts      [nr,nb]
    #         - E         = energy (Joules) [nr]
    #         - r         = range (meters)  [nb]
    #         - nr        = number of records (scalar-int)
    #         - nb        = number of bins    (scalar-int)
    #         - bg_st_bin = first bin of solar background region (scalar-int)
    #         - bg_ed_bin = last bin of solar background region (scalar-int)
    #         - opt       = string = 'bg_sub' or 'NRB' or 'NRB_no_range' 
    #
    # OUTPUS: 
    #         - NRB = NRB, corrected counts [nr,nb]

    # NOTES:
    # [5/8/18] 'NRB_no_range' opt added. Select if counts are already 
    # background-subtracted and range-corrected.
    
    # Compute the solar background
    
    bg = compute_solar_background(rc,bg_st_bin,bg_ed_bin)  
    
    # Compute NRB

    # Broadcast 1D arrays to 2D so computation is vectorized (read: fast)
    if opt is 'NRB': r = np.broadcast_to(r,(nr,nb))
    E = np.broadcast_to(E,(nb,nr)).transpose()
    bg = np.broadcast_to(bg,(nb,nr)).transpose()
    
    if opt == 'NRB':
        NRB =  ((rc - bg) * r**2) / E 
    elif opt == 'bg_sub':
        NRB = rc - bg
    elif opt == 'NRB_no_range':  
        NRB = rc / E
    else:
        print('Invalid option entered. stopping execution in lidar.py')
        pdb.set_trace()
    
    return NRB
    


def convert_raw_energy_monitor_values(raw_EMs,nwl,instr='CAMAL',e_flg=None):
    """As of 10/18/17, this routine is only specific to CAMAL"""

    # INPUTS:
    # raw_EMs -> Raw, untouched, energy monitor values. In array form.
    #            Dimension should be something like [nchans x nwl] but
    #            this could depend on instrument. I plan on put multiple
    #            instruments' conversions in this code.
    # nwl -> The number of wavelengths.
    # instr -> String identifying the instrument.
    # e_flg -> A carry-over from CPL-IDL world which tells this code which
    #          equations should be used on CPL conversions.
    #
    # OUTPUT:
    # EMall -> array containing converted energy, in micro-joules, for all
    #          wavelengths.
    
    if instr == 'CAMAL':
    
        sz = raw_EMs.shape
        nr = sz[0]
    
        EMall = np.zeros((nr,nwl),dtype=np.int32) #nwl is in initializations.py
    
        EMall[:,0] = raw_EMs[:,5].astype(np.int32)*16**0 + \
            raw_EMs[:,6].astype(np.int32)*16**2 +     \
            raw_EMs[:,7].astype(np.int32)*16**4 +     \
            raw_EMs[:,8].astype(np.int32)*16**6 +     \
            raw_EMs[:,9].astype(np.int32)*16**8
        EMall[:,1] = raw_EMs[:,0].astype(np.int32)*16**0 + \
            raw_EMs[:,1].astype(np.int32)*16**2 +     \
            raw_EMs[:,2].astype(np.int32)*16**4 +     \
            raw_EMs[:,3].astype(np.int32)*16**6 +     \
            raw_EMs[:,4].astype(np.int32)*16**8 
        EMall[:,2] = raw_EMs[:,10].astype(np.int32)*16**0 + \
            raw_EMs[:,11].astype(np.int32)*16**2 +     \
            raw_EMs[:,12].astype(np.int32)*16**4 +     \
            raw_EMs[:,13].astype(np.int32)*16**6 +     \
            raw_EMs[:,14].astype(np.int32)*16**8           
        
        # After this line EMall is now of dtype float
        EMall = EMall.astype(np.float32) 
        # As of 1/26/18, the following 3 lines of calculations are obsolete
        #EMall[:,0] = 8.63101e-4 + 2.88557e-10*EMall[:,0] - 1.76057e-18*EMall[:,0]**2 + 0.0*EMall[:,0]**3
        #EMall[:,1] = 0.003932537 + 4.11985e-11*EMall[:,1] - 8.82687e-21*EMall[:,1]**2 + 0.0*EMall[:,1]**3
        #EMall[:,2] = 0.012758984 + 2.45103e-10*EMall[:,2] - 1.17876e-19*EMall[:,2]**2 + 0.0*EMall[:,2]**3
        # As of 1/26/18, Andrew Kupchock says convert the energy monitors as follows
        EMall[:,0] = 1.29089e-10*EMall[:,0] + 0.024267116
        EMall[:,1] = 2.71388E-11*EMall[:,1] + 0.005726336
        EMall[:,2] = 1.36981E-10*EMall[:,2] + 0.00206153
        # At this point all EM values are in milli-joules. Convert to micro-joules
        # in return statement.
        EMall = EMall*1e3
    
    elif instr == 'CPL':

        # NOTES:
        # This CPL Energy Monitor conversion section was copied & translated from 
        # Mr. Dennis Hlavka's "engin_corr_v3.pro" code. I even included the "e_flg"
        # and kept the meanings the same.
        #
        # CPL's EMs -> 0=355, 1=532, 2=1064

        sz = raw_EMs.shape
        nr = sz[0]   

        if e_flg == 6:
            # TC4-07 settings . Initially valid 6/19/2006
            emon_t = raw_EMs.astype(np.float64)/500.0
            emon_c = np.copy(emon_t)
            m0_355= 0.0
            m1_355= 0.05
            m2_355= 0.0
            m3_355= 0.0
            emon_c[:,0]= m0_355 + m1_355*emon_t[:,0] + m2_355*emon_t[:,0]**2 + m3_355*emon_t[:,0]**3
            m0_532= -0.25721
            m1_532= 0.038868
            m2_532= -1.533e-6
            m3_532= 0.0
            emon_c[:,1]= m0_532 + m1_532*emon_t[:,1] + m2_532*emon_t[:,1]**2 + m3_532*emon_t[:,1]**3
            m0_1064= -26.806
            m1_1064= 0.19851
            m2_1064= -1.2727e-4
            m3_1064= 0.0
            emon_c[:,2]= m0_1064 + m1_1064*emon_t[:,2] + m2_1064*emon_t[:,2]**2 + m3_1064*emon_t[:,2]**3
        elif e_flg == 7:
            # GloPac10 settings (UAV-CPL) (initially valid 02/16/2010)
            emon_t = raw_EMs.astype(np.float64)/500.0
            emon_c = np.copy(emon_t)
            m0_355= -2.4349
            m1_355= 0.10332
            m2_355= 2.5793e-04
            m3_355= 0.0
            emon_c[:,0]= m0_355 + m1_355*emon_t[:,0] + m2_355*emon_t[:,0]**2 + m3_355*emon_t[:,0]**3
            m0_532= 0.36922
            m1_532= 0.013169
            m2_532= 2.4631e-07
            m3_532= 0.0000
            emon_c[:,1]= m0_532 + m1_532*emon_t[:,1] + m2_532*emon_t[:,1]**2 + m3_532*emon_t[:,1]**3
            m0_1064= 5.2746
            m1_1064= 0.047994
            m2_1064= -7.1374e-06
            m3_1064= 0.000
            emon_c[:,2]= m0_1064 + m1_1064*emon_t[:,2] + m2_1064*emon_t[:,2]**2 + m3_1064*emon_t[:,2]**3
        elif e_flg == 9:
            # new laser box settings (ER2-CPL) 26Feb17	  
            emon_t = raw_EMs.astype(np.float64)/500.0
            emon_c = np.copy(emon_t)
            m0_355= -3.618073187e0
            m1_355= 0.334612636e0
            m2_355= -0.001398674e0
            m3_355= 0.0
            emon_c[:,0]= m0_355 + m1_355*emon_t[:,0] + m2_355*emon_t[:,0]**2 + m3_355*emon_t[:,0]**3
            lowlocs = np.where(emon_c[:,0] < 0.001)
            emon_c[lowlocs[0],0] = 0.001
            highlocs = np.where(emon_c[:,0] > 20.0)
            emon_c[highlocs[0],0] = 20.0
            m0_532= -0.08109323e0
            m1_532= 0.074648283e0
            m2_532= -0.0000079241e0
            m3_532= 0.0
            emon_c[:,1]= m0_532 + m1_532*emon_t[:,1] + m2_532*emon_t[:,1]**2 + m3_532*emon_t[:,1]**3
            m0_1064= 2.172313321e0
            m1_1064= 0.063209038e0
            m2_1064= -0.0000237511e0
            m3_1064= 0.0
            emon_c[:,2]= m0_1064 + m1_1064*emon_t[:,2] + m2_1064*emon_t[:,2]**2 + m3_1064*emon_t[:,2]**3
        elif e_flg == 10:
            # ACT-A, C-130 energy conversion
            # (updated 12/09/2016)  EMON_C in units of micro-Joules
            emon_t = raw_EMs.astype(np.float64)/500.0
            emon_c = np.copy(emon_t)
            m0_355= -1.491632131e0
            m1_355= 1.05116511e-1
            m2_355= 0.0
            m3_355= 0.0
            emon_c[:,0]= m0_355 + m1_355*emon_t[:,0] + m2_355*emon_t[:,0]**2 + m3_355*emon_t[:,0]**3
            #if (emon_c(0) lt .001D) then emon_c[:,0]= 0.001 <- IDL equiv. for ref
            emon_c[np.argwhere(emon_c[:,0] < 0.001),0] = 0.001
            #if (emon_c(0) gt 20.00D) then emon_c(0)= 20.00D <- IDL equiv. for ref
            emon_c[np.argwhere(emon_c[:,0] > 20.0),0] = 20.0
            m0_532= 1.45462874e-1
            m1_532= 7.142194e-3
            m2_532= 3.1427e-8
            m3_532= 0.0
            emon_c[:,1]= m0_532 + m1_532*emon_t[:,1] + m2_532*emon_t[:,1]**2 + m3_532*emon_t[:,1]**3
            m0_1064= 1.090521825e0
            m1_1064= 4.1418361e-2
            m2_1064= -4.2348e-7
            m3_1064= 0.0
            emon_c[:,2]= m0_1064 + m1_1064*emon_t[:,2] + m2_1064*emon_t[:,2]**2 + m3_1064*emon_t[:,2]**3 
        elif e_flg == 12:
            # ER2-CPL settings (initially valid 5/1/2018)
            # Transcribed from email sent by Andrew Kupchock on 4/25/18
            # [8/22/18] UPDATE
            # According to Andrew, this EM conversion was changed in June after more
            # tinkering in the lab. The May 1, 2018 version was never used outside of the lab
            # so I am overwriting it with the June version which is currently being used for
            # the first time during the REThinC 2018 campaign. 	    
            emon_t = raw_EMs.astype(np.float64)/250.0
            emon_c = np.copy(emon_t)
            m0_355= -0.000166416385199
            m1_355= 0.161304322
            m2_355= -1.086938077
            m3_355= 0.0
            emon_c[:,0]= m0_355*emon_t[:,0]**2 + m1_355*emon_t[:,0] + m2_355
            m0_532= 0.09809502
            m1_532= -1.017836594
            m2_532= 0.0
            m3_532= 0.0
            emon_c[:,1]=  m0_532*emon_t[:,1] + m1_532
            m0_1064= -2.45355E-06
            m1_1064= 0.001553344
            m2_1064= -0.018753417
            m3_1064= 0.973099563
            emon_c[:,2]= m0_1064*emon_t[:,2]**3 + m1_1064*emon_t[:,2]**2 + m2_1064*emon_t[:,2] + m3_1064                      
        else:
             print(str(e_flg)+' is an invalid e_flg value. Stopping in Energy Monitor conversion funtion.')
             pdb.set_trace()
        EMall = emon_c

    else:
        
        print("Instrument name "+instr+" is not valid.")
        print("Stopping in Energy Monitor conversion function.")
        pdb.set_trace()
    
    return EMall   
        


def calculate_off_nadir_angle(phi0,phi0_azi,Y,P,R, **kwargs):
    """ Calculates the off-nadir angle from given inputs.
        Equations used are derived from matrix rotations.
        Fuselage coordinates are rotated into nadir-fixed coordinates.
        Calculation follows what is shown in {G. Cai et al., Unmanned 
        Rotorcraft Systems, Advances in Industrial Control,
        DOI 10.1007/978-0-85729-635-1_2, © Springer-Verlag London 
        Limited 2011} , {Accuracy of Vertical Air Motions from Nadir-
        Viewing Doppler Airborne Radars, Gerald M. Heymsfield, Journal 
        of Atmospheric and Oceanic Technology, December 1989} , and
        {Mapping of Airborn Doppler Radar Data, Wn-Chau Lee et al, 
        Journal of Atmospheric and Oceanic Technology, 1994}.
        !!! ALL ANGLES ARE RADIANS IN CONTEXT OF THIS FUNCTION !!!
        Y, P, and R angles are assumed to be with respect to the
        instrument coordinate system. The instrument and fuselage 
        coordinate systems are assumed to be perfectly aligned.
    """
    
    # ASSUMPTIONS:
    # - All angles (input and output) are in radians.
    # - Y, P, and R angles are with respect to the instrument coordinate
    #   system.
    # - Instrument coordinate origins: +X points perpendicular to right  
    #   of motion, +Y aligns with motion, and +Z points down.
    # - Y is defined as a rotation about the z-axis.
    # - R is defined as a rotation about the y-axis.
    # - P is defined as a rotation about the x-axis.
    # - Rotation angles (Y,P,R) are positive if they move towards a
    #   positive axis origin.
    
    # INPUTS:
    # phi0 -> The off-nadir angle built into the instrument. When the
    #         instrument coordinate system is perfectly level, this is 
    #         off-nadir angle at which the beam points. phi0 should always
    #         be positive. It's not advised that phi0 be greater than or
    #         equal to 90 degrees. If it is, the code will likely run fine;
    #         however, the accuracy of the off-nadir angle calculation
    #         cannot be garanteed. Negative values of phi0 will also produce
    #         unexpected results.
    # phi0_azi -> The counterclockwise angle from the +X instrument axis
    #             (+X defined as perpendicular and to the right of instrument
    #             motion). Counterclockwise from +X obviously defines
    #             a positive phi0_azi.
    # Y -> Yaw, in radians
    # P -> Pitch, in radians
    # R -> Roll, in radians
    
    # OUTPUTS:
    # ONA -> off-nadir angle (in radians)
    
    # Assuming an initial nadir-pointing vector with length equal to 1,
    # compute the x, y, and z components of the intial vector after it's
    # been rotated to the fuselage-ONA (phi0).
    
    x0 = 0.0
    y0 = 0.0
    z0 = 1.0
    
    # Default point angle of beam, broken into components of instrument
    # coordinate system.
    hori_xy_disp = z0 * math.sin(phi0)
    x = hori_xy_disp * math.cos(phi0_azi) # x comp of def. pointing ang.
    y = hori_xy_disp * math.sin(phi0_azi) # y comp of def. pointing ang.
    z = z0 * math.cos(phi0)
    
    # Now compute the nadir-fixed coordinate system components, after
    # rotating through the 3 attitude angles, yaw (Y, about the z-ax), 
    # roll (R, about the y-ax), and pitch (P, about the x-ax), in that 
    # order.
    x1 = x*math.cos(Y)*math.cos(R) + y*(math.sin(Y)*math.cos(P)+math.sin(P)*math.sin(R)*math.cos(Y)) + \
         z*(math.sin(P)*math.sin(Y)-math.sin(R)*math.cos(Y)*math.cos(P))
    y1 = x*(-1*math.sin(Y)*math.cos(R)) + y*(math.cos(Y)*math.cos(P)-math.sin(P)*math.sin(Y)*math.sin(R)) \
         + z*(math.sin(P)*math.cos(Y)+math.cos(P)*math.sin(Y)*math.sin(R))
    z1 = x*math.sin(R) + y*(-1*math.sin(P)*math.cos(R)) + z*math.cos(P)*math.cos(R)
    
    # Now use the dot-product of the initial, nadir-pointing vector, and
    # the rotated vector, to compute the off-nadir angle.
    dotproduct = x1*x0 + y1*y0 +z1*z0
    magvect = math.sqrt(x1**2 + y1**2 +z1**2)
    magnull = math.sqrt(x0**2 + y0**2 +z0**2)
    ONA = math.acos(dotproduct/(magvect*magnull))
    #print('input phi0 is ',phi0*(180./const.pi),' degrees.')
    #print('Off-nadir angle is ',ONA*(180./const.pi),' degrees.')
    
    if ('xz_ang' in kwargs):
        # Compute the angle in the x-z plane (CAMAL purposes in mind)
        dotproduct = x1*x0 + z1*z0
        magvect = math.sqrt(x1**2 +z1**2)
        magnull = math.sqrt(x0**2 +z0**2)
        xz_ang = math.acos(dotproduct/(magvect*magnull))
        return [ONA, xz_ang]
        
    if ('xy_ang' in kwargs):
        # Compute the angle in the x-y plane (CAMAL purposes in mind)

        # The next if block prevents a div-by-zero error.
        # This occurs when R == 0.0
        if x1 == 0.0: 
            if P >= 0.0: xy_ang = const.pi/2.0 # points forward
            if P < 0.0: xy_ang = -1.0*(const.pi/2.0) # points behind
        else:
            xy_ang = np.arctan(y1/x1)
            if x1 < 0.0:
                xy_ang = xy_ang + const.pi
            elif (x1 > 0.0) and (y1 < 0.0):
                xy_ang = xy_ang + 2.0*const.pi
            elif (x1 > 0.0) and (y1 > 0.0):
                pass

        return [ONA, xy_ang]
    
    return ONA
        
