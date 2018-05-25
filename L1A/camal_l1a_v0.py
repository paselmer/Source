# CAMAL L1A processing code.
# NOTE: The most current version of the L1A code must always be named
#       "camal_l1a.py." All versions of the code must be saved as
#       "camal_l1a_vX.py" but not actually be source file that's exectued.

# UPDATE LOG:
#
# [3/2/18] v0 created
# This version is first attempt to write out to HDF5 file in chunks.
#
# [3/7/18] dead time tables now incorporated.
#
# [3/8/18] Nav_interp_T_float64
# This variable is a key to decreasing this code's runtime.
# This variable gets defined no matter which Nav source is selected.
# This variable is defined to be the instrument's Unix Time and must have
# units of the number of seconds since the Unix Epoch.
#
# [3/9/18] Nav, time
# I think I've finally sorted out the Nav data. The time offsets between 
# the various clocks are now computed automatically, with the option to 
# overwrite with manual values if needed. Even after the time offset are 
# applied to the Nav data, A "nudge" forward in time seems to be neccessary 
# to line up IWG1 with GPS. This is on the order of ~1 sec.
#
# [3/15/18] DEM
# DEM elevations and surface types at coordinates have been added. They 
# DEM data are ingested through a call to a shared library (DLL on
# Windows, SO on Unix). The values are written out to the NRB file
# for both the nadir and laser spot coordinates.
#
# [3/19/18] Laser spot
# DEM nadir surface elevation is now used as the altitude at which to
# compute the laser spot. Formerly, it was LS_alt, an initializations
# file input.
#
# [3/22/18] Nav
# Laserspot computed with GPS vs. Nav wasn't matching too well.
# Long story short, discovered I should just zero out "yaw" in 
# GPS data for now.


# Import libraries <----------------------------------------------------

# Libraries I did not create
import os
import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import time
import timeit
import datetime as DT
import h5py
# Libraries I did create
from initializations import *
from read_routines import *
from lidar import *
from time_conversions import *
from mutable_data_structs import define_MSC_structure
import matplotlib.dates as mdates
import ctypes


# Define functions <----------------------------------------------------

def compute_time_offsets(MCS_file=None,quick='no'):
    """ Compute the offset between the UTC timestamp in the MCS data and
        the UTC timestamp in the IWG1 data.
    """
    
    # NOTES:
    # Set quick == 'quick' to just grab a manually computed number that's
    # been hard-coded into this routine.
    
    if quick == 'quick':
        # Computed manually using 12/7/17 flight dataset. Saved in init. file.
        # Subtract this number from the MCS GPS time to convert to IWG1 time
        time_offset_IWG1 = def_time_offset_IWG1   # in init. file
        # Subtract this number from the Unix Time to convert it to MCS GPS time
        time_offset_UnixT = def_time_offset_UnixT # in init. file       
    else:
        MCS_data_1file = read_in_raw_data(MCS_file)
        nav_data_all = read_entire_nav_dataset()
        UnixT_epoch = DT.datetime(1970,1,1)
        nav_UnixT_timedelta = nav_data_all['UnixT'] - \
            UnixT_epoch - DT.timedelta(seconds=secs_btwn_instr_UnixT_and_UTC)
        nav_UnixT_float = np.zeros(nav_data_all.shape[0],dtype=np.float64)
        IWG1_offsets = [] # sample of offsets
        j = 0
        MCS_UnixT_DT = np.zeros(MCS_data_1file.shape[0],dtype=DT.datetime)
        UnixT_offsets = np.zeros(MCS_data_1file.shape[0],dtype=np.float64)
        for i in range(0,MCS_data_1file.shape[0]):
            MCS_UnixT_DT[i] = DT.datetime.fromtimestamp(MCS_data_1file['meta']['CCSDS']['UnixT'][i]) + \
                DT.timedelta(seconds=secs_btwn_instr_UnixT_and_UTC)
            MCS_UTC = weeksecondstoutc(MCS_data_1file['meta']['GpsWeek'][i]*1.0,
                                MCS_data_1file['meta']['Gps_msec'][i]/1e3,leapsecs )
            UnixT_offsets[i] = (MCS_UnixT_DT[i]- MCS_UTC).total_seconds()
            match_loc = np.where(nav_data_all['UnixT'] == MCS_UnixT_DT[i])
            if (match_loc[0].shape[0] > 0):
                IWG1_offsets.append((nav_data_all['UTC'][match_loc[0][0]] - MCS_UTC).total_seconds())
        if len(IWG1_offsets) == 0:
            print("Time offsets are None. Try a different MCS file...")
            time_offset_UnixT = None
            time_offset_IWG1 = None
        else:
            time_offset_UnixT = DT.timedelta(seconds=np.median(UnixT_offsets))
            time_offset_IWG1 = DT.timedelta(seconds=np.median(np.asarray(IWG1_offsets)))

        
    return [time_offset_UnixT, time_offset_IWG1]

    
# Start main execution here <-------------------------------------------
print('Starting main L1A execution at: ',DT.datetime.now())

# Create and load file list for MCS data
MCS_file_list = 'processing_file_list.txt'
search_str = 'data*'
create_a_file_list(MCS_file_list,search_str)
with open(MCS_file_list) as MCS_list_fobj:
    all_MCS_files = MCS_list_fobj.readlines()
nMCS_files = len(all_MCS_files)

# Load the shared library of C functions
np_clib = np.ctypeslib.load_library(clib,clib_path) 
# DEM shared C++ library
DEM_lib_file = DEM_lib_path + DEM_Cpp_name
DEM_lib = ctypes.cdll.LoadLibrary(DEM_lib_file)
DEM_file = (DEM_dir + DEM_name).encode('utf8')
DEM_file_ctype = ctypes.c_char_p(DEM_file)
# Defined C++ DEM read function output
class DEM_out(ctypes.Structure):
    _fields_ = [ ("Elev",ctypes.c_float), ("LOCI",ctypes.c_ubyte) ]
DEM_lib.get_DEM_point.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_double]
DEM_lib.get_DEM_point.restype = DEM_out    
# Compute the time offset between MCS and IWG1
for MCS_file in all_MCS_files:
    time_offset_UnixT, time_offset_IWG1 = compute_time_offsets(MCS_file.strip(),offset_choice)
    if time_offset_UnixT is not None: break
if time_offset_UnixT is None:
    print("There are not enough matches to compute IWG1 time offset.")
    print("You could try setting the initialization file to grab \
           hard-coded values from the initialization file.")
    print("Time offsets are None! Code will stop!")
print("Time offsets are: ", time_offset_UnixT, time_offset_IWG1)

# Load nav data for entire flight. Select from 1 of several possible
# nav data streams. All these read-in all flight data at once. This
# decision was made to simplify the code and perhaps boost processing speed.

if Nav_source == 'gps':

    # Load the entire gps dataset into memory
    gps_data_all, gps_UTC = read_entire_gps_dataset('scrub')
    
    # Interpolate data records to match CAMAL's data rate.
    # 10 Hz as of 2/2/18. Data rate set in initializations.
    
    del_t = np.zeros(gps_data_all.shape[0],dtype=DT.datetime)
    del_t[1:] = gps_UTC[1:] - gps_UTC[:-1]
    del_t[0] = del_t[1]*0.0
    
    # Compute # of elapsed seconds since first record (float)
    
    del_secs = np.zeros(gps_data_all.shape[0],dtype=np.float64)
    for k in range(0,gps_data_all.shape[0]): del_secs[k] = del_t[k].total_seconds()
    
    elapsed_secs = np.cumsum(del_secs)
    tot_elapsed = np.max(elapsed_secs)
    m = 1.0 / MCS_hz # basically, the time rez in secs you'd like to get to
    ix = np.arange(0,tot_elapsed,m) # x coordinates of the interpolated values
    n_interp = ix.shape[0]
    gps_interp = np.zeros(n_interp,dtype=gps_struct)
    gps_UTC_interp = np.zeros(n_interp,dtype=DT.datetime) # note separate array
    UnixT_epoch = DT.datetime(1970,1,1)
    Nav_interp_T_float64 = np.zeros(n_interp,dtype=np.float64)
    for field in gps_data_all.dtype.names:
        gps_interp[field] = np.interp(ix, elapsed_secs, gps_data_all[field])
    # Now handle population of interp time field
    for k in range(0,n_interp): 
        gps_UTC_interp[k] = gps_UTC[0] + DT.timedelta(seconds=ix[k])
        Nav_interp_T_float64[k] = (gps_UTC_interp[k] - UnixT_epoch + 
            time_offset_UnixT).total_seconds()

    # NOW, set the array that will be used in processing, gps_data_all equal
    # to the interpolated array that was just created, gps_interp.
    gps_data_all = gps_interp
    # Remember, the UTC you'll need later is contained in separate array,
    # gps_UTC_interp.
    
    # This variable will be the same no matter what the Nav_source
    Nav_hz = gps_hz

elif Nav_source == 'nav':
    
    # NOTE: "Nav" variables will apply generically to whichever of the
    #       several Nav_source's is selected. "nav" variables will 

    # Load the entire nav dataset into memory
    nav_data_all = read_entire_nav_dataset()
    
    # Interpolate data records to match CAMAL's data rate.
    # 10 Hz as of 2/2/18. Data rate set in initializations.
    
    del_t = np.zeros(nav_data_all.shape[0],dtype=DT.datetime)
    del_t[1:] = nav_data_all['UTC'][1:] - nav_data_all['UTC'][:-1]
    del_t[0] = del_t[1]*0.0
    
    # Compute # of elapsed seconds since first record (float)
    
    del_secs = np.zeros(nav_data_all.shape[0],dtype=np.float64)
    for k in range(0,nav_data_all.shape[0]): del_secs[k] = del_t[k].total_seconds()
    
    elapsed_secs = np.cumsum(del_secs)
    tot_elapsed = np.max(elapsed_secs)
    m = 1.0 / MCS_hz
    ix = np.arange(0,tot_elapsed,m) # x coordinates of the interpolated values
    n_interp = ix.shape[0]
    nav_interp = np.zeros(n_interp,dtype=nav_struct)
    Nav_interp_T_float64 = np.zeros(n_interp,dtype=np.float64)
    UnixT_epoch = DT.datetime(1970,1,1)
    for field in nav_data_all.dtype.names:
        if (field == 'UTC') or (field == 'UnixT'): continue
        nav_interp[field] = np.interp(ix, elapsed_secs, nav_data_all[field])
    # Now handle population of interp time field
    for k in range(0,n_interp): 
        nav_interp['UnixT'][k] = nav_data_all['UnixT'][0] + DT.timedelta(seconds=ix[k])
        Nav_interp_T_float64[k] = (nav_interp['UnixT'][k] - UnixT_epoch).total_seconds() + nudge
        # offsets in following line convert to MCS UTC
        nav_interp['UTC'][k] = ( nav_data_all['UTC'][0] + DT.timedelta(seconds=ix[k]) 
            - time_offset_IWG1 + DT.timedelta(seconds=nudge) )
    
    # NOW, set the array that will be used in processing, nav_data_all equal
    # to the interpolated array that was just created, nav_interp.
    
    nav_data_all = nav_interp
    
    # This variable will be the same no matter what the Nav_source
    Nav_hz = nav_hz    
    
elif Nav_source == 'iwg1':

    # Read in entire IWG1 dataset. Typically contained in single file.
    IWG1_file = raw_dir + IWG1_file
    IWG1_data = read_in_IWG1_data(IWG1_file,est_IWG1_recs)
    
    # Interpolate data records to match CAMAL's data rate.
    # 10 Hz as of 2/2/18. Data rate set in initializations.
    
    del_t = np.zeros(IWG1_data.shape[0],dtype=DT.datetime)
    del_t[1:] = IWG1_data['UTC'][1:] - IWG1_data['UTC'][:-1]
    del_t[0] = del_t[1]*0.0
    
    # Compute # of elapsed seconds since first record (float)
    
    del_secs = np.zeros(IWG1_data.shape[0],dtype=np.float64)
    for k in range(0,IWG1_data.shape[0]): del_secs[k] = del_t[k].total_seconds()
    
    elapsed_secs = np.cumsum(del_secs)
    tot_elapsed = np.max(elapsed_secs)
    m = 1.0 / MCS_hz
    ix = np.arange(0,tot_elapsed,m) # x coordinates of the interpolated values
    n_interp = ix.shape[0]
    IWG1_interp = np.zeros(n_interp,dtype=IWG1_struct)
    UnixT_epoch = DT.datetime(1970,1,1)
    Nav_interp_T_float64 = np.zeros(n_interp,dtype=np.float64)
    for field in IWG1_data.dtype.names:
        if field == 'UTC': continue
        IWG1_interp[field] = np.interp(ix, elapsed_secs, IWG1_data[field])
    # Now handle population of interp time field
    for k in range(0,n_interp):
        # offsets in following line convert to MCS UTC 
        IWG1_interp['UTC'][k] = ( IWG1_data['UTC'][0] + DT.timedelta(seconds=ix[k]) 
            - time_offset_IWG1 + DT.timedelta(seconds=nudge) )
        Nav_interp_T_float64[k] = (IWG1_interp['UTC'][k] - UnixT_epoch + 
            time_offset_UnixT).total_seconds() + nudge

    # NOW, set the array that will be used in processing, IWG1_data equal
    # to the interpolated array that was just created, IWG1_interp.
    
    IWG1_data = IWG1_interp
    
    # This variable will be the same no matter what the Nav_source
    Nav_hz = IWG1_hz    
    
else:
    
    print("You've entered an invalid Nav_source choice. Halting code.")
    pdb.set_trace()


# Load all the dead time tables

DTT = np.zeros((len(DTT_files),max_counts),dtype=np.float32)
cc = 0
for DTT_file in DTT_files: 
    one_DTT = read_in_dead_time_table(dtt_dir+DTT_file)
    DTT[cc,:] = one_DTT[:max_counts]
    cc += 1   


print('Starting core processing...') # <--------------------------------
first_read = True
for f in range(10,12):#range(0,nMCS_files):
    
    MCS_file = all_MCS_files[f]
    MCS_file = MCS_file.rstrip()
    MCS_data_1file = read_in_raw_data(MCS_file)
    if MCS_data_1file is None:
        print('\n******* Bad file! Skipping file!!! *******')
        print(MCS_file)
        continue  
    nr_1file = MCS_data_1file.shape[0]                    # # of records in current file
    
    # Interpolate MCS times to data rate
    xp = (MCS_data_1file['meta']['Gps_msec'] - MCS_data_1file['meta']['Gps_msec'][0])/1000.0
    fp = MCS_data_1file['meta']['Gps_msec'].astype(np.float64)
    ix = np.arange(np.min(xp),np.max(xp),1.0/MCS_hz)
    interp_Gps_msec = np.interp(ix,xp,fp)
    xp = (MCS_data_1file['meta']['CCSDS']['UnixT'] - MCS_data_1file['meta']['CCSDS']['UnixT'][0])/1.0
    fp = MCS_data_1file['meta']['CCSDS']['UnixT'].astype(np.float64)
    ix = np.arange(np.min(xp),np.max(xp),1.0/MCS_hz)        
    interp_UnixT = np.interp(ix,xp,fp)
    interp2orig_indicies = np.zeros(nr_1file,dtype=np.uint32)
    MCS_UnixT_float64_new = np.zeros(nr_1file,dtype=np.float64)
    # ---> START process of calling a C function <---
    #      This C function populates an array 
    #      mapped to the MCS Unix Time values
    #      with the appropriate interpolated values.
    #      This function also populates an array
    #      that maps the Nav_data indicies to
    #      the MCS data. This drastically speeds
    #      up the profile loop in python code.
    MCS_UnixT_float64_new = np.require(MCS_UnixT_float64_new,float,['ALIGNED','C_CONTIGUOUS'])
    MCS_UnixT_float64_orig = np.require(fp,float,['ALIGNED','C_CONTIGUOUS'])
    interp_UnixT = np.require(interp_UnixT,float,['ALIGNED','C_CONTIGUOUS'])
    Nav_interp_T_float64 = np.require(Nav_interp_T_float64,float,['ALIGNED','C_CONTIGUOUS'])
    interp2orig_indicies = np.require(interp2orig_indicies,int,['ALIGNED','C_CONTIGUOUS'])
    np_clib.map_interp_times_to_orig_frame.restype = None
    np_clib.map_interp_times_to_orig_frame.argtypes = [np.ctypeslib.ndpointer(float,
        ndim=1, flags='aligned'), np.ctypeslib.ndpointer(float, ndim=1,
        flags='aligned'), np.ctypeslib.ndpointer(float, ndim=1,flags='aligned'),
        np.ctypeslib.ndpointer(float, ndim=1,flags='aligned'),
        np.ctypeslib.ndpointer(int, ndim=1,flags='aligned'), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp) ]
    np_clib.map_interp_times_to_orig_frame(MCS_UnixT_float64_new, MCS_UnixT_float64_orig, interp_UnixT,
        Nav_interp_T_float64, interp2orig_indicies, ctypes.c_double(1.0/Nav_hz),
        MCS_UnixT_float64_new.ctypes.strides, MCS_UnixT_float64_new.ctypes.shape, 
        MCS_UnixT_float64_orig.ctypes.strides, MCS_UnixT_float64_orig.ctypes.shape,
        interp_UnixT.ctypes.strides, interp_UnixT.ctypes.shape,
        Nav_interp_T_float64.ctypes.strides, Nav_interp_T_float64.ctypes.shape,  
        interp2orig_indicies.ctypes.strides, interp2orig_indicies.ctypes.shape) 
    # ---> END process of calling a C function <---                       

    # Save/create a few key data parameters...
    if first_read:
        nb =  MCS_data_1file['meta']['nbins'][0]              # # of bins
        vrZ = (MCS_data_1file['meta']['binwid'][0] * c) / 2.0 # vert. rez in m
        nc = MCS_data_1file['meta']['nchans'][0]              # # of channels
        nshots = MCS_data_1file['meta']['nshots'][0]      
        vrT = MCS_data_1file['meta']['binwid'][0]
        phi0_azi = 0.0
        Y = 0.0 * (pi/180.0)                                  # Initialize YPR
        P = 0.0 * (pi/180.0)
        R = 0.0 * (pi/180.0)
        i = 0 # counts total records of entire dataset (flight)
        g = 0
        # Create array to help compute heading if using GPS data...
        m_atan = np.zeros((3,4),dtype=np.uint16)
        m_atan[2,3] = 0
        m_atan[2,1] = 1
        m_atan[0,1] = 2
        m_atan[0,3] = 3
        # Now create a fixed altitude frame onto which to fit the counts
        ffrme = set_fixed_bin_alt_frame(ff_bot_alt,ff_top_alt,vrZ_ff,nb,pointing_dir)
        ffrme = np.require(ffrme,float,['ALIGNED','C_CONTIGUOUS'])
        nb_ff = ffrme.shape[0] # number of bins in the fixed frame
        # Load the overlap table
        overlap_vector = read_in_overlap_table(olap_dir+overlap_file)
        # Overlap is one long array, where first nbins are
        # 1064, next nbins are 532, and next (final) nbins are 355.
        overlaps = overlap_vector.reshape((nwl,nb))
        # Open the hdf5 file and create the datasets
        hdf5_fname = L1_dir+'NRB_'+proj_name+'_'+flt_date+'_'+Nav_source+'.hdf5'
        hdf5_file = h5py.File(hdf5_fname, 'a')         
        try:            
            meta_dset = hdf5_file.create_dataset("meta", (nr_1file,), maxshape=(None,), dtype=MCS_meta_struct)
            PGain_dset = hdf5_file.create_dataset("PGain", (len(PGain),), np.float32)
            nav_dset = hdf5_file.create_dataset("nav", (nr_1file,), maxshape=(None,), dtype=nav_save_struct)
            laserspot_dset = hdf5_file.create_dataset("laserspot", (nr_1file,2) , maxshape=(None,2), dtype=np.float32)
            ONA_dset = hdf5_file.create_dataset("ONA", (nr_1file,), maxshape=(None,), dtype=np.float32)
            bin_alt_dset = hdf5_file.create_dataset("bin_alt_array", ffrme.shape, ffrme.dtype)
            nb_ff_dset = hdf5_file.create_dataset("num_ff_bins", (1,), np.uint32)
            num_recs_dset = hdf5_file.create_dataset("num_recs", (1,), np.uint32)
            DEM_nadir_dset = hdf5_file.create_dataset("DEM_nadir", (nr_1file,), maxshape=(None,), dtype=np.float32)
            DEM_nadir_surftype_dset = hdf5_file.create_dataset("DEM_nadir_surftype", (nr_1file,), maxshape=(None,), dtype=np.int8)            
            DEM_laserspot_dset = hdf5_file.create_dataset("DEM_laserspot", (nr_1file,), maxshape=(None,), dtype=np.float32)
            DEM_laserspot_surftype_dset = hdf5_file.create_dataset("DEM_laserspot_surftype", (nr_1file,), maxshape=(None,), dtype=np.int8)
            EM_dset = hdf5_file.create_dataset("EM", (nr_1file,nwl) , maxshape=(None,nwl), dtype=np.float32)
            NRB_dset = hdf5_file.create_dataset("nrb", (nc,nr_1file,nb_ff), maxshape=(nc,None,nb_ff), dtype=np.float32)
        except RuntimeError:
            print("HDF5 file for this dataset already exists, overwriting old file...")
            hdf5_file.close() #close, delete, reopen...
            delete_file(hdf5_fname)
            hdf5_file = h5py.File(hdf5_fname, 'a')
            meta_dset = hdf5_file.create_dataset("meta", (nr_1file,), maxshape=(None,), dtype=MCS_meta_struct)
            PGain_dset = hdf5_file.create_dataset("PGain", (len(PGain),), np.float32)
            nav_dset = hdf5_file.create_dataset("nav", (nr_1file,), maxshape=(None,), dtype=nav_save_struct)
            laserspot_dset = hdf5_file.create_dataset("laserspot", (nr_1file,2) , maxshape=(None,2), dtype=np.float32)
            ONA_dset = hdf5_file.create_dataset("ONA", (nr_1file,), maxshape=(None,), dtype=np.float32)
            bin_alt_dset = hdf5_file.create_dataset("bin_alt_array", ffrme.shape, ffrme.dtype)
            nb_ff_dset = hdf5_file.create_dataset("num_ff_bins", (1,), np.uint32)
            num_recs_dset = hdf5_file.create_dataset("num_recs", (1,), np.uint32)
            DEM_nadir_dset = hdf5_file.create_dataset("DEM_nadir", (nr_1file,), maxshape=(None,), dtype=np.float32)
            DEM_nadir_surftype_dset = hdf5_file.create_dataset("DEM_nadir_surftype", (nr_1file,), maxshape=(None,), dtype=np.int8)            
            DEM_laserspot_dset = hdf5_file.create_dataset("DEM_laserspot", (nr_1file,), maxshape=(None,), dtype=np.float32)
            DEM_laserspot_surftype_dset = hdf5_file.create_dataset("DEM_laserspot_surftype", (nr_1file,), maxshape=(None,), dtype=np.int8)       
            EM_dset = hdf5_file.create_dataset("EM", (nr_1file,nwl) , maxshape=(None,nwl), dtype=np.float32)
            NRB_dset = hdf5_file.create_dataset("nrb", (nc,nr_1file,nb_ff), maxshape=(nc,None,nb_ff), dtype=np.float32)
        except:
            print("An unanticipated error occurred while trying to create \
                the HDF5 datasets. Stopping execution.")
            pdb.set_trace()
        # Some datasets only need to be written once, right here...
        PGain_dset[:] = np.asarray(PGain) # convert list to array
        bin_alt_dset[:] = ffrme
        nb_ff_dset[:] = nb_ff
    else:
        # Expand dataset sizes to accomodate next input MCS file
        meta_dset.resize(meta_dset.shape[0]+nr_1file, axis=0)
        nav_dset.resize(nav_dset.shape[0]+nr_1file, axis=0)
        laserspot_dset.resize(laserspot_dset.shape[0]+nr_1file, axis=0)
        ONA_dset.resize(ONA_dset.shape[0]+nr_1file, axis=0)
        DEM_nadir_dset.resize(DEM_nadir_dset.shape[0]+nr_1file, axis=0)
        DEM_nadir_surftype_dset.resize(DEM_nadir_surftype_dset.shape[0]+nr_1file, axis=0)
        DEM_laserspot_dset.resize(DEM_laserspot_dset.shape[0]+nr_1file, axis=0)
        DEM_laserspot_surftype_dset.resize(DEM_laserspot_surftype_dset.shape[0]+nr_1file, axis=0)
        EM_dset.resize(EM_dset.shape[0]+nr_1file, axis=0)
        NRB_dset.resize(NRB_dset.shape[1]+nr_1file, axis=1)
		
        
    counts_ff = np.zeros((nc,nr_1file,nb_ff),dtype=np.float32)
    NRB = np.empty_like(counts_ff)
        
    Nav_save = np.zeros(nr_1file,dtype=nav_save_struct) #NOTE THE STRUCTURE TYPE!
    ONA_save = np.zeros(nr_1file,dtype=np.float32)
    laserspot = np.zeros((nr_1file,2),dtype=np.float32) # [lat x lon] 
    DEM_nadir = np.zeros(nr_1file,dtype=np.float32)
    DEM_nadir_surftype = np.zeros(nr_1file,dtype=np.int8)-9
    DEM_laserspot = np.zeros(nr_1file,dtype=np.float32)
    DEM_laserspot_surftype = np.zeros(nr_1file,dtype=np.int8)-9
    
    # Begin main record loop
    
    i1f = 0  # counts the records in one (current) file
    
    while i1f < nr_1file:
    
        if Nav_source == 'gps':
            
            # Use index map to match
            # NOTE: Nav_save is an IWG1-style structure (see immutable_dat_structs.py)
            try:
                Nav_match = gps_data_all[interp2orig_indicies[i1f]]
            except IndexError:
                print("Code stopped via debugger in Nav_source block.")
                print("I am guessing that all Nav times are < min(MCS time)")
                pdb.set_trace()
            Nav_save[i1f]['roll'] = Nav_match['roll']
            Nav_save[i1f]['pitch'] = Nav_match['pitch']
            Nav_save[i1f]['drift'] = 0.0
            Nav_save[i1f]['heading'] = np.arctan2(Nav_match['east'],Nav_match['north'])*(180.0/np.pi)
            Nav_save[i1f]['lon'] = Nav_match['lon']
            Nav_save[i1f]['lat'] = Nav_match['lat'] 
            Nav_save[i1f]['GPS_alt'] = Nav_match['alt']           
            Nav_save[i1f]['UTC'] = np.asarray(list(gps_UTC_interp[interp2orig_indicies[i1f]].strftime(
                "%Y-%m-%dT%H:%M:%S.%f").encode('utf8')),dtype=np.uint8)
            Y = 0.0# Nav_match['yaw'] * (pi/180.0) Woah! This isn't what I think yaw should be [3/22/18]
            P = Nav_match['pitch'] * (pi/180.0)
            R = Nav_match['roll'] * (pi/180.0)             
        
        elif Nav_source == 'nav':
  
            # Use index map to match
            # NOTE: Nav_save is an IWG1-style structure (see immutable_dat_structs.py)
            try:
                Nav_match = nav_data_all[interp2orig_indicies[i1f]]                                                     
            except IndexError:
                print("Code stopped via debugger in Nav_source block.")
                print("I am guessing that all Nav times are < min(MCS time)")
                pdb.set_trace()
            Nav_save[i1f]['roll'] = Nav_match['roll']
            Nav_save[i1f]['pitch'] = Nav_match['pitch']
            Nav_save[i1f]['drift'] = Nav_match['drift']
            Nav_save[i1f]['heading'] = Nav_match['heading']
            Nav_save[i1f]['lon'] = Nav_match['lon']
            Nav_save[i1f]['lat'] = Nav_match['lat']
            Nav_save[i1f]['UTC'] = np.asarray(list(Nav_match['UTC'].strftime(
                "%Y-%m-%dT%H:%M:%S.%f").encode('utf8')),dtype=np.uint8)
            Nav_save[i1f]['GPS_alt'] = Nav_match['GPS_alt']
            Y = Nav_match['drift'] * (pi/180.0) 
            P = Nav_match['pitch'] * (pi/180.0)
            R = Nav_match['roll'] * (pi/180.0)  
            
        elif Nav_source == 'iwg1':
            
            # Find matching IWG1 data point. Use index map to match.
            # NOTE: Nav_save is an IWG1-style structure (see immutable_dat_structs.py)
            Nav_match = IWG1_data[interp2orig_indicies[i1f]]
            Nav_save[i1f]['roll'] = Nav_match['roll']
            Nav_save[i1f]['pitch'] = Nav_match['pitch']
            Nav_save[i1f]['drift'] = Nav_match['drift']
            Nav_save[i1f]['heading'] = Nav_match['heading']
            Nav_save[i1f]['lon'] = Nav_match['lon']
            Nav_save[i1f]['lat'] = Nav_match['lat']
            # The following line is necessary to get a clean, interpretable,
            # time to write out to the HDF5 file. It will essetially be a
            # byte array for IDL.
            Nav_save[i1f]['UTC'] = np.asarray(list(Nav_match['UTC'].strftime(
                "%Y-%m-%dT%H:%M:%S.%f").encode('utf8')),dtype=np.uint8)
            Nav_save[i1f]['GPS_alt'] = Nav_match['GPS_alt']
            Y = 0.0 #Nav_match['drift'] * (pi/180.0) 
            P = Nav_match['pitch'] * (pi/180.0)
            R = Nav_match['roll'] * (pi/180.0)             
    
        # Proceed no farther in processing if altitude is not high enough
        if Nav_save[i1f]['GPS_alt'] < alt_cutoff:
            #print('Bad altitude: ',Nav_save[i1f]['GPS_alt'])
            #print('Will not process until good alt is read.\n')
            i1f = i1f + 1
            continue
            
        # Apply dead time and overlap corrections
        cc = 0 # count the channels
        for DTT_file in DTT_files:
            MCS_data_1file['counts'][i1f][cc,:] = DTT[ cc, np.rint(MCS_data_1file['counts'][i1f][cc,:]).astype(np.uint32) ]
            MCS_data_1file['counts'][i1f][cc,:] = overlaps[wl_map[cc],:] * MCS_data_1file['counts'][i1f][cc,:]
            cc += 1
    
        # Calculate the off nadir angle
        # For CAMAL, positive scan angles are to the right
        phi0 = (MCS_data_1file['meta']['scan_pos'][i1f] + angle_offset) * (pi/180.0)
        if phi0 >= 0.0: 
            phi0_azi = 0.0
        else:
            phi0_azi = pi
        [ONA, xy_ang] = calculate_off_nadir_angle(abs(phi0),phi0_azi,Y,P,R,xy_ang=True)   # accepted way as of [1/31/18]
        #ONA = calculate_off_nadir_angle(0.0,0.0,Y,P,R+(-1.0*phi0)) # accepted way as of [1/25/18]
        ONA_save[i1f] = ONA
        
        # Get the DEM value at aircraft nadir coordinate
        DEM_ret = DEM_lib.get_DEM_point(DEM_file,Nav_save[i1f]['lat'],Nav_save[i1f]['lon'])
        DEM_nadir[i1f] = DEM_ret.Elev
        if DEM_ret.Elev == -9999.0: DEM_nadir[i1f] = 0.0 # ocean surface
        DEM_nadir_surftype[i1f] = DEM_ret.LOCI
        
        # Calculate the laser spot coordinates
        [latLS,lonLS] = laser_spot_aircraft(Nav_save[i1f]['GPS_alt'],Nav_save[i1f]['lat']*(np.pi/180.0),
            Nav_save[i1f]['lon']*(np.pi/180.0),xy_ang,DEM_nadir[i1f],ONA,Nav_save[i1f]['heading']*(np.pi/180.0))
        laserspot[i1f,0] = latLS
        laserspot[i1f,1] = lonLS
        
        # Get the DEM value at the laser spot coordinate
        DEM_ret = DEM_lib.get_DEM_point(DEM_file, laserspot[i1f,0]*(180.0/np.pi),
                  laserspot[i1f,1]*(180.0/np.pi) )
        DEM_laserspot[i1f] = DEM_ret.Elev
        if DEM_ret.Elev == -9999.0: DEM_laserspot[i1f] = 0.0 # ocean surface
        DEM_laserspot_surftype[i1f] = DEM_ret.LOCI
    
        # Compute the altitudes of the raw data bins
        bin_alts = compute_raw_bin_altitudes(nb, pointing_dir,Nav_save[i1f]['GPS_alt'],
                   vrZ, ONA)

        length_bin_alts = bin_alts.shape[0]
        if length_bin_alts > nb:
            print('Length of bin_alts vectors is > than nbins. Stopping code...')
            pdb.set_trace()
    
        # Put the bins in the fixed altitude frame
        counts_ff[:,i1f,:] = put_bins_onto_fixed_frame_C(np_clib,ffrme,
            bin_alts,MCS_data_1file['counts'][i1f],nc) 
    
        i = i + 1    # increment record counter by 1
        i1f = i1f + 1  # increment record counter for current file by 1
        
    
    first_read = False 
    if i == 0: continue
    
    # Compute NRB
    EMs = convert_raw_energy_monitor_values(MCS_data_1file['meta']['EM'],nwl)
    for chan in range(0,nc):
        EMsubmit = EMs[:,wl_map[chan]]     
        NRB[chan,:,:] = correct_raw_counts(counts_ff[chan,:,:],EMsubmit,np.arange(0,nb_ff*vrZ_ff,vrZ_ff),
            i1f,nb_ff,ff_bg_st_bin,ff_bg_ed_bin,'NRB')
            
    # Write out one file's worth of data to the hdf5 file
    meta_dset[i-nr_1file:i] = MCS_data_1file['meta']
    nav_dset[i-nr_1file:i] = Nav_save
    laserspot_dset[i-nr_1file:i,:] = laserspot
    ONA_dset[i-nr_1file:i] = ONA_save
    DEM_nadir_dset[i-nr_1file:i] = DEM_nadir
    DEM_nadir_surftype_dset[i-nr_1file:i] = DEM_nadir_surftype
    DEM_laserspot_dset[i-nr_1file:i] = DEM_laserspot
    DEM_laserspot_surftype_dset[i-nr_1file:i] = DEM_laserspot_surftype    
    EM_dset[i-nr_1file:i,:] = EMs
    NRB_dset[:,i-nr_1file:i,:] = NRB
        
    print('\n**********************************')
    print('\nDone with file: '+MCS_file+'\n')
    print('**********************************\n')

      
print('Main L1A execution finished at: ',DT.datetime.now())

# Write out any final parameters to the HDF5 file that needed to wait

num_recs_dset[:] = i
hdf5_file.close()
print("NRB has been written to the HDF5 file:"+hdf5_fname)

# The final message

print('Total raw profiles processed: '+str(i))
print("camal_l1a.py has finished normally.")
#pdb.set_trace()

######################################################################### BELOW THIS LINE SHALL NOT BE PART OF L1A PROCESS

tit = '532 nm NRB'
xlimits = [0,ONA_save.shape[0]]
ylimits = [1200,200]#[900,500]
samp_chan = NRB[2,:,:] + NRB[3,:,:]
curtain_plot(samp_chan.transpose(), nb_ff, vrZ_ff, ffrme, 0, 3e9, hori_cap, pointing_dir,figW, figL, CPpad, 'records', 'altitude(m)', tit, 'alt',[ylimits[0],ylimits[1]], 'recs',[xlimits[0],xlimits[1]], scale_alt_OofM, 1, out_dir)

#make_custom_plot(samp_chan.transpose(), nb_ff, vrZ_ff, ffrme, 0, 3e9, hori_cap, pointing_dir,
                      #figW, figL, CPpad, 'records', 'altitude(m)', tit, 'alt',  
                      #[ylimits[0],ylimits[1]], 'recs', [xlimits[0],xlimits[1]], scale_alt_OofM, 1, out_dir, str(f).strip()+'.png',
                      #np.arange(xlimits[0],xlimits[1]),MCS_data_1file['meta'][xlimits[0]:xlimits[1]]['scan_pos']+angle_offset,
                      #np.arange(xlimits[0],xlimits[1]),MCS_data_1file['meta'][xlimits[0]:xlimits[1]]['scan_state'],
                      #np.arange(xlimits[0],xlimits[1]),Nav_save['roll'][xlimits[0]:xlimits[1]])
