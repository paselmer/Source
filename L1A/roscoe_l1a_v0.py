# Roscoe L1A processing code.
# NOTE: The most current version of the L1A code must always be named
#       "roscoe_l1a.py." All versions of the code must be saved as
#       "roscoe_l1a_vX.py" but not actually be source file that's exectued.
#
# UPDATE LOG:
# 
# [12/16/19] v0 created
#
# [1/10/20] background subtraction corrected for up-looking processing
# A few more lines of code were added because the background zone
# selection code that worked for downward-looking did not work for
# upward-looking.

# Import libraries <----------------------------------------------------

# Libraries I did not create
import os
import pdb
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import time
import timeit
import datetime as DT
import h5py
import matplotlib.dates as mdates
import ctypes
# Libraries I did create
from initializations import *
from read_routines import *
from lidar import *
from time_conversions import *



# Define functions <----------------------------------------------------

def patch_outliers(data_array, m):
    """ Function to overwrite outliers, in a numpy array, with reasonable
        values.
        Originally written to path over bad Roscoe energy monitor values.
    """
    
    # INPUTS: data_array -> self explanatory
    #         m          -> The number of standard deviations
    
    # RETURNS: A patched version of data_array with the same shape and 
    #          number of elements.
    
    keep_indxs = abs(data_array - np.mean(data_array)) < m * np.std(data_array)
    single_reasonable_value = np.mean( data_array[keep_indxs] )
    data_array[np.invert(keep_indxs)] = single_reasonable_value
    return data_array
    
    
def compute_time_offsets(MCS_file=None,quick='no'):
    """ Compute the offset between the timestamp (UnixT which is converted
        to UTC within this function) in the MCS data and the UTC timestamp
        in the IWG1 data.
    """
    
    # NOTES:
    # Set quick == 'quick' to just grab a manually computed number that's
    # been hard-coded into this routine.
    
    if quick == 'quick':
        # Computed manually using 10/23/19 flight dataset. Saved in init. file.
        # Subtract this number from the Unix Time to convert it to UTC time
        time_offset_UnixT = def_time_offset_UnixT # in init. file       
    else:
        MCS_data_1file = read_in_mcs_data_r(MCS_file)
        nav_data_all = read_entire_nav_dataset('*nav*',raw_dir,est_nav_recs_1file,secs_btwn_instr_UnixT_and_UTC)
        UnixT_offsets = []
        for i in range(0,MCS_data_1file.shape[0]):
            try:
                MCS_UnixT_DT = DT.datetime.fromtimestamp(MCS_data_1file['meta']['CCSDS']['UnixT'][i]) + \
                    DT.timedelta(seconds=secs_btwn_instr_UnixT_and_UTC)
            except:
                print('The UnixT in the MCS data is out-of-bounds!')
                continue
            match_loc = np.where(nav_data_all['UnixT'] == MCS_UnixT_DT)
            if (match_loc[0].shape[0] > 0):
                nav_match_UTC = nav_data_all['UTC'][match_loc[0][0]]
                UnixT_offsets.append( (MCS_UnixT_DT - nav_match_UTC).total_seconds() )
        if len(UnixT_offsets) == 0:
            print("Time offsets are None. Try a different MCS file...")
            time_offset_UnixT = None
        else:
            time_offset_UnixT = DT.timedelta(seconds=np.median(UnixT_offsets))

    return time_offset_UnixT

    
# Start main execution here <-------------------------------------------
print('Starting main L1A execution at: ',DT.datetime.now())

# Create and load file list for raw ('data*') data
MCS_file_list = 'processing_file_list.txt'
create_a_file_list(MCS_file_list,search_str,raw_dir)
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

# Compute the time offset between MCS and IWG1 (which comes in form or "IWG1" or "nav" files)
for MCS_file in all_MCS_files:
    time_offset_UnixT = compute_time_offsets(MCS_file.strip(),offset_choice)
    if time_offset_UnixT is not None: break  

# Load nav data for entire flight. Select from 1 of several possible
# nav data streams. All these read-in all flight data at once. This
# decision was made to simplify the code and perhaps boost processing speed.
# NOTE that for this Roscoe version of L1A, you are able to use files 
# generated by CAMAL/ACATS and IWG1; however, as of 4/5/18, this hasn't
# been tested and there may be time offsets to figure out first.

if Nav_source == 'gps':

    # Load the entire gps dataset into memory
    gps_data_all, gps_UTC = read_entire_gps_dataset('scrub')
    
    # Interpolate data records to match CAMAL's data rate.
    # 10 Hz as of 2/2/18. MCS rate set in initializations.
    
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
    #       refer to CAMAL/ACTS "nav*.data" files. 

    # Load the entire nav dataset into memory
    nav_data_all = read_entire_nav_dataset('*nav*',raw_dir,est_nav_recs_1file,secs_btwn_instr_UnixT_and_UTC)
    
    # Interpolate data records to match Roscoe's data rate.
    # 10 Hz as of 2/2/18. MCS rate set in initializations.
    
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
    for field in nav_data_all.dtype.names:
        if (field == 'UTC') or (field == 'UnixT'): continue
        nav_interp[field] = np.interp(ix, elapsed_secs, nav_data_all[field])
    # Now handle population of interp time field
    for k in range(0,n_interp): 
        #nav_interp['UnixT'][k] = nav_data_all['UnixT'][0] + DT.timedelta(seconds=ix[k])
        # offsets in following line convert to MCS UTC
        nav_interp['UTC'][k] = ( nav_data_all['UTC'][0] + DT.timedelta(seconds=ix[k]) 
            + DT.timedelta(seconds=nudge) )
        Nav_interp_T_float64[k] = (nav_interp['UTC'][k] - UnixT_epoch + time_offset_UnixT).total_seconds()
    
    # NOW, set the array that will be used in processing, nav_data_all equal
    # to the interpolated array that was just created, nav_interp.

    nav_data_all = nav_interp
    
    # This variable will be the same no matter what the Nav_source
    Nav_hz = nav_hz   
    
elif Nav_source == 'iwg1':

    # "Dummy out" these variables, which are only need valid values for MCS-Nav processing
    cls_meta_data_all = None
    nav2cls_indx = np.zeros((1000,2))
    FL_trunc = None 
    usable_file_range = [0,1000]

    # Compute offset between IWG1 time and MCS Instrument Time ("ExactTime")
    iwg1_cls_t_offset = compute_time_offset_IWG1_MCS(all_MCS_files)

    # Read in entire IWG1 dataset. Typically contained in single file.
    IWG1_file = raw_dir + IWG1_file
    IWG1_data = read_in_IWG1_data(IWG1_file,est_IWG1_recs)
    
    # Interpolate data records to match CAMAL's data rate.
    # 10 Hz as of 2/2/18. MCS rate set in initializations.
    
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
    Nav_interp_T_float64 = np.zeros(n_interp,dtype=np.float64)
    for field in IWG1_data.dtype.names:
        if field == 'UTC': continue
        IWG1_interp[field] = np.interp(ix, elapsed_secs, IWG1_data[field])
    # Now handle population of interp time field
    for k in range(0,n_interp):
        # offsets in following line convert to MCS UTC 
        IWG1_interp['UTC'][k] = ( IWG1_data['UTC'][0] + DT.timedelta(seconds=ix[k]) )
        Nav_interp_T_float64[k] = (IWG1_interp['UTC'][k] - UnixT_epoch).total_seconds()

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
hdf5_created = False
trans_bin = [0,0]
trans_total = 0
last_file = False 
for f in range(0,nMCS_files):
    
    # Read in MCS data of one file and initialize a few things
    MCS_file = all_MCS_files[f]
    MCS_file = MCS_file.rstrip()
    # Read in file with EM data
    EM_file = MCS_file.replace('dataup',EM_file_tag)
    MCS_data_1file = read_in_mcs_data_r(EM_file)
    EM_data = np.copy(MCS_data_1file['meta']['EM'])
    MCS_data_1file = read_in_mcs_data_r(MCS_file)
    if EM_data.shape[0] != MCS_data_1file.shape[0]:
        print(attention_bar)
        print('It appears up & down MCS files do not have same shape')
        print('# recs in EM file:  ',EM_data.shape)
        print('# recs in dat file: ', MCS_data_1file.shape)
        print('EM file: ',EM_file)
        print('MCS file: ',MCS_file)
        print(attention_bar)
        if (EM_data.shape[0] < MCS_data_1file.shape[0]):
            MCS_data_1file['meta']['EM'][:EM_data.shape[0]] = EM_data
        else:
            MCS_data_1file['meta']['EM'] = EM_data[:MCS_data_1file.shape[0]]
    else:
        MCS_data_1file['meta']['EM'] = EM_data
    if MCS_data_1file is None:
        print('\n******* Bad file! Skipping file!!! *******')
        print(MCS_file)
        continue  
    nr_1file = MCS_data_1file.shape[0]                          # # of records in current file
    good_rec_bool = np.ones(MCS_data_1file.shape[0],dtype=bool) # is 'True' for good records, 'False' for bad
    satur_flag = False # Will turn True if any bin in current file is potentially saturated.
    
    # Interpolate MCS times to data rate
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
    if interp_UnixT.shape[0] < 1:
        print(attention_bar)
        print("Interpolated UnixT has length of zero.")
        print("Something funky going on with the time.")
        print("Perhaps all UnixT's are the same value.")
        print("Due to this issue, skipping file:")
        print(MCS_file)
        print(attention_bar)
        continue
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

    # Deem as "bad," records outside user-specified processing time range
    # Use nav data time for this filtering
    if Nav_source == 'gps': nav_UTC = gps_UTC_interp[interp2orig_indicies]
    if Nav_source == 'nav': nav_UTC = nav_data_all[interp2orig_indicies]['UTC']
    if Nav_source == 'iwg1': nav_UTC = nav_data_all[interp2orig_indicies]['UTC']
    if ( (nav_UTC.min() < process_start) or (nav_UTC.max() >  process_end) ):
        time_range_mask_low = nav_UTC > process_start
        time_range_mask_high = nav_UTC < process_end
        time_range_mask = time_range_mask_low * time_range_mask_high
        MCS_data_1file = np.extract(time_range_mask,MCS_data_1file)
        interp2orig_indicies = np.extract(time_range_mask,interp2orig_indicies)
        MCS_UnixT_float64_new = np.extract(time_range_mask,MCS_UnixT_float64_new)
        good_rec_bool = np.ones(MCS_data_1file.shape[0],dtype=bool)
        print(attention_bar)
        print('In file '+MCS_file+'...')
        print(str(nr_1file - MCS_data_1file.shape[0]).strip()+' records are out of time range.')
        print(attention_bar)
        nr_1file = MCS_data_1file.shape[0]
        if (nr_1file == 0):
            print(attention_bar)
            print("File number "+str(f).strip()+" has no usable data!")
            print("Skipping processing of this file.")
            print(attention_bar)
            continue   
            
    # Create histogram bin boundaries. If first_read, created farther down in first_read block.
    if not first_read: time_bins = np.arange(trans_bin[0],MCS_UnixT_float64_orig[-1]+3.0*secs2avg,secs2avg)

    # Save/create a few key data parameters...
    if ((first_read) & (not hdf5_created)):
        nb =  MCS_data_1file['meta']['nbins'][0]              # # of bins
        vrZ = (MCS_data_1file['meta']['binwid'][0] * c) / 2.0 # vert. rez in m
        nc = MCS_data_1file['meta']['nchans'][0]              # # of channels
        nshots = MCS_data_1file['meta']['nshots'][0]      
        vrT = MCS_data_1file['meta']['binwid'][0]
        # Range vector, broadcast to nc x nb for vectorized operations later
        bins_range=np.broadcast_to(np.arange(0,nb*vrZ,vrZ,dtype=np.float32),(nc,nb))
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
        # Create histogram bin boundaries.
        time_bins = np.arange(MCS_UnixT_float64_orig[0],MCS_UnixT_float64_orig[-1]+3.0*secs2avg,secs2avg)
        # Load the overlap table
        overlap_vector = read_in_overlap_table(olap_dir+overlap_file)
        # Overlap is one long array, where first nbins (833) are
        # 355, next nbins are 532, and next (final) nbins are 1064.
        overlaps = overlap_vector.reshape((nwl,nb))
        # Make the overlap_vector into an array where the 1st dim. lines up sequentially
        # with the channels
        overlaps_chan_seq = np.ones((nc,nb),dtype=overlaps.dtype)
        for chan in range(0,nc):
            overlaps_chan_seq[chan,:] = overlaps[wl_map[chan],:]
        # Open the hdf5 file and create the datasets
        hdf5_fname = L1_dir+'Roscoe_NRB_'+proj_name+'_'+flt_date+'_'+Nav_source+'.hdf5'
        hdf5_file = h5py.File(hdf5_fname, 'a')
        # Open the lat/lon GEOS met data sampling CSV file
        GEOS_samp_fname = L1_dir+'lon_lat_UTC_'+proj_name+'_'+flt_date+'_'+Nav_source+'.txt'
        GEOS_fobj = open(GEOS_samp_fname, 'wt')
        GEOS_fobj.write('lon,lat,time\n') # file header     
        try:
            PGain_dset = hdf5_file.create_dataset("PGain", (len(PGain),), np.float32)
            nav_dset = hdf5_file.create_dataset("nav", (1,), maxshape=(None,), dtype=nav_save_struct)
            laserspot_dset = hdf5_file.create_dataset("laserspot", (1,2) , maxshape=(None,2), dtype=np.float32)
            ONA_dset = hdf5_file.create_dataset("ONA", (1,), maxshape=(None,), dtype=np.float32)
            bin_alt_dset = hdf5_file.create_dataset("bin_alt_array", ffrme.shape, ffrme.dtype)
            nb_ff_dset = hdf5_file.create_dataset("num_ff_bins", (1,), np.uint32)
            num_recs_dset = hdf5_file.create_dataset("num_recs", (1,), np.uint32)
            DEM_nadir_dset = hdf5_file.create_dataset("DEM_nadir", (1,), maxshape=(None,), dtype=np.float32)
            DEM_nadir_surftype_dset = hdf5_file.create_dataset("DEM_nadir_surftype", (1,), maxshape=(None,), dtype=np.int8)            
            DEM_laserspot_dset = hdf5_file.create_dataset("DEM_laserspot", (1,), maxshape=(None,), dtype=np.float32)
            DEM_laserspot_surftype_dset = hdf5_file.create_dataset("DEM_laserspot_surftype", (1,), maxshape=(None,), dtype=np.int8)
            EM_dset = hdf5_file.create_dataset("EM", (1,nwl) , maxshape=(None,nwl), dtype=np.float32)
            NRB_dset = hdf5_file.create_dataset("nrb", (nc,1,nb_ff), maxshape=(nc,None,nb_ff), dtype=np.float32)
            bg_dset = hdf5_file.create_dataset("bg", (nc,1), maxshape=(nc,None), dtype=np.float32)
            saturate_ht_dset = hdf5_file.create_dataset("saturate_ht", (nc,1), maxshape=(nc,None), dtype=np.float32)
        except RuntimeError:
            print("HDF5 file for this dataset already exists, overwriting old file...")
            hdf5_file.close() #close, delete, reopen...
            delete_file(hdf5_fname)
            hdf5_file = h5py.File(hdf5_fname, 'a')
            PGain_dset = hdf5_file.create_dataset("PGain", (len(PGain),), np.float32)
            nav_dset = hdf5_file.create_dataset("nav", (1,), maxshape=(None,), dtype=nav_save_struct)
            laserspot_dset = hdf5_file.create_dataset("laserspot", (1,2) , maxshape=(None,2), dtype=np.float32)
            ONA_dset = hdf5_file.create_dataset("ONA", (1,), maxshape=(None,), dtype=np.float32)
            bin_alt_dset = hdf5_file.create_dataset("bin_alt_array", ffrme.shape, ffrme.dtype)
            nb_ff_dset = hdf5_file.create_dataset("num_ff_bins", (1,), np.uint32)
            num_recs_dset = hdf5_file.create_dataset("num_recs", (1,), np.uint32)
            DEM_nadir_dset = hdf5_file.create_dataset("DEM_nadir", (1,), maxshape=(None,), dtype=np.float32)
            DEM_nadir_surftype_dset = hdf5_file.create_dataset("DEM_nadir_surftype", (1,), maxshape=(None,), dtype=np.int8)            
            DEM_laserspot_dset = hdf5_file.create_dataset("DEM_laserspot", (1,), maxshape=(None,), dtype=np.float32)
            DEM_laserspot_surftype_dset = hdf5_file.create_dataset("DEM_laserspot_surftype", (1,), maxshape=(None,), dtype=np.int8)       
            EM_dset = hdf5_file.create_dataset("EM", (1,nwl) , maxshape=(None,nwl), dtype=np.float32)
            NRB_dset = hdf5_file.create_dataset("nrb", (nc,1,nb_ff), maxshape=(nc,None,nb_ff), dtype=np.float32)
            bg_dset = hdf5_file.create_dataset("bg", (nc,1), maxshape=(nc,None), dtype=np.float32)
            saturate_ht_dset = hdf5_file.create_dataset("saturate_ht", (nc,1), maxshape=(nc,None), dtype=np.float32)
        except:
            print("An unanticipated error occurred while trying to create the HDF5 datasets. Stopping execution.")
            pdb.set_trace()
        # Some datasets only need to be written once, right here...
        PGain_dset[:] = np.asarray(PGain) # convert list to array
        bin_alt_dset[:] = ffrme
        nb_ff_dset[:] = nb_ff
        hdf5_created = True
        
    counts_ff = np.zeros((nc,nr_1file,nb_ff),dtype=np.float32)
    NRB = np.empty_like(counts_ff)
    bg_save = np.zeros((nc,nr_1file),dtype=np.float32) 
        
    Nav_save = np.zeros(nr_1file,dtype=nav_save_struct) #NOTE THE STRUCTURE TYPE!
    ONA_save = np.zeros(nr_1file,dtype=np.float32)
    laserspot = np.zeros((nr_1file,2),dtype=np.float32) # [lat x lon] 
    DEM_nadir = np.zeros(nr_1file,dtype=np.float32)
    DEM_nadir_surftype = np.zeros(nr_1file,dtype=np.int8)-9
    DEM_laserspot = np.zeros(nr_1file,dtype=np.float32)
    DEM_laserspot_surftype = np.zeros(nr_1file,dtype=np.int8)-9
    saturate_ht = np.zeros((nc,nr_1file),dtype=np.float32)-99999.9

    # Right here, right now, test to see if any bin in current file is potentially saturated

    # Quick, broad test...
    satur_locs = []
    satur_locs_indx = 0
    if (MCS_data_1file['counts'].max() > min(saturation_values)): satur_flag = True
    if satur_flag:
        print('\n'+all_MCS_files[f]+' has saturated bin values in it!\n')
        # satur_locs will have [recs x chan x bin] as each row. recs dim should increase
        # towards higher indexes.
        satur_locs = np.argwhere(MCS_data_1file['counts'] > min(saturation_values))
                
    
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
            time_str = gps_UTC_interp[interp2orig_indicies[i1f]].strftime("%Y-%m-%dT%H:%M:%S.%f")           
            Nav_save[i1f]['UTC'] = np.asarray(list(time_str.encode('utf8')),dtype=np.uint8)
            Y = 0.0# Nav_match['yaw'] * (pi/180.0) Woah! This isn't what I think yaw should be [3/22/18]
            P = Nav_match['pitch'] * (pi/180.0)
            R = Nav_match['roll'] * (pi/180.0) 
            # Write 'lat,lon,UTC' out to file
            GEOS_out_str = '{0:.4f},{1:.4f},{2}\n'.format(Nav_match['lon'],Nav_match['lat'],time_str[:19])
            GEOS_fobj.write(GEOS_out_str)            
        
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
            time_str = Nav_match['UTC'].strftime("%Y-%m-%dT%H:%M:%S.%f")
            Nav_save[i1f]['UTC'] = np.asarray(list(time_str.encode('utf8')),dtype=np.uint8)
            # NOTE on next line: Invalids are -999.9, so max should choose valid if exists
            Nav_save[i1f]['GPS_alt'] = max(Nav_match['GPS_alt_msl'],Nav_match['GPS_alt'])
            Y = Nav_match['drift'] * (pi/180.0) 
            P = Nav_match['pitch'] * (pi/180.0)
            R = Nav_match['roll'] * (pi/180.0)  
            # Write 'lat,lon,UTC' out to file
            GEOS_out_str = '{0:.4f},{1:.4f},{2}\n'.format(Nav_match['lon'],Nav_match['lat'],time_str[:19])
            GEOS_fobj.write(GEOS_out_str)	    
            
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
            time_str = Nav_match['UTC'].strftime("%Y-%m-%dT%H:%M:%S.%f")    
            Nav_save[i1f]['UTC'] = np.asarray(list(time_str.encode('utf8')),dtype=np.uint8)
            Nav_save[i1f]['GPS_alt'] = Nav_match['GPS_alt']
            Y = 0.0 #Nav_match['drift'] * (pi/180.0) 
            P = Nav_match['pitch'] * (pi/180.0)
            R = Nav_match['roll'] * (pi/180.0)
            # Write 'lat,lon,UTC' out to file
            GEOS_out_str = '{0:.4f},{1:.4f},{2}\n'.format(Nav_match['lon'],Nav_match['lat'],time_str[:19])
            GEOS_fobj.write(GEOS_out_str)
    
        # Proceed no farther in processing if altitude is not high enough
        if Nav_save[i1f]['GPS_alt'] < alt_cutoff:
            #print('Bad altitude: ',Nav_save[i1f]['GPS_alt'])
            #print('Will not process until good alt is read.\n')
            good_rec_bool[i1f] = False
            i = i + 1
            i1f = i1f + 1
            continue
            
        # Patch over bins with phantom signal [12/18/19]
        for cc in range(0,nc):
            MCS_data_1file['counts'][i1f,cc,phb1:phb2+1] = ( np.linspace( MCS_data_1file['counts'][i1f,cc,phb1-1],
                MCS_data_1file['counts'][i1f,cc,phb2+1], (phb2+2) - (phb1-1), dtype=np.uint16 )[1:-1] )

        # Convert to float for subsequent manipulation within this loop
        counts_float32 = MCS_data_1file['counts'][i1f].astype(np.float32)
            
        # Apply dead time correction
        cc = 0 # count the channels
        for DTT_file in DTT_files:
            try:
                counts_float32[cc,:] = DTT[ cc, MCS_data_1file['counts'][i1f,cc,:] ]
            except IndexError:
                # Probably a noise spike
                spike_locs = np.where( counts_float32[cc,:] > DTT.shape[1] )
                counts_float32[cc,spike_locs[0][0]] = 0.0
            cc += 1
    
        # Calculate the off nadir angle
        # For CAMAL, positive scan angles are to the right
        # For Roscoe, assume phi0 and phi0_azi are zero
        [ONA, xy_ang] = calculate_off_nadir_angle(0.0,0.0,Y,P,R,xy_ang=True)   # accepted way as of [1/31/18]
        #ONA = calculate_off_nadir_angle(0.0,0.0,Y,P,R+(-1.0*phi0)) # accepted way as of [1/25/18]
        ONA_save[i1f] = ONA
        if ONA > ONA_cutoff:	
            good_rec_bool[i1f] = False
            i = i + 1
            i1f = i1f + 1
            continue
        
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

        # Populate saturate_ht with bin alts if any bins in current profile are saturated
        # NOTE: The order of cond's matters in if blocks in this code paragraph
        if ((satur_flag) and (i1f in satur_locs[:,0])):
            while ((satur_locs_indx < satur_locs.shape[0]) and (satur_locs[satur_locs_indx,0] == i1f)):
                indx = satur_locs[satur_locs_indx,:] #Will always have 3 dims in same order
                if (MCS_data_1file['counts'][indx[0],indx[1],indx[2]] > saturation_values[indx[1]]):
                    saturate_ht[indx[1],i1f] = bin_alts[indx[2]]
                satur_locs_indx += 1

        length_bin_alts = bin_alts.shape[0]
        if length_bin_alts > nb:
            print('Length of bin_alts vectors is > than nbins. Stopping code...')
            pdb.set_trace()

        # Subtract background and correct raw counts for range (in scope of actual bin-alt frame)
        # It's difficult to apply range correction once counts are in fixed frame.
        try:
            if pointing_dir == 'Down':
              bg_loc1 = np.argwhere(bin_alts <= bg_st_alt)[0][0]
              bg_loc2 = np.argwhere(bin_alts >= bg_ed_alt)[-1][0]
            else:
              bg_loc1 = np.argwhere(bin_alts <= bg_st_alt)[-1][0]
              bg_loc2 = np.argwhere(bin_alts >= bg_ed_alt)[0][0]
            bg_1D = np.mean(counts_float32[:,bg_loc1:bg_loc2],axis=1) 
            bg_save[:,i1f] = bg_1D
            bg = np.broadcast_to(bg_1D, (nb,nc)).transpose()
        except:
            print(attention_bar)
            print("!!!!!!! WARNING !!!!!!!")
            print("No data within defined background region! Using background of zero.")
            print(Nav_save[i1f]['UTC'],' ONA: ',ONA*(180.0/pi))
            print("!!!!!!! WARNING !!!!!!!")
            print(attention_bar)
            bg = np.zeros((nc,nb))
        
        range_cor_af_counts_float32 = ( ( counts_float32 - bg )
                                          * bins_range**2  ) 
        
        # Apply overlap correction here. I would have applied it at the same time as the
        # deadtime, except Roscoe's overlaps can be funky in far ranges, which interferes with
        # background sub.
        range_cor_af_counts_float32 = range_cor_af_counts_float32 / overlaps_chan_seq
    
        # Put the bins in the fixed altitude frame
        #counts_ff[:,i1f,:] = put_bins_onto_fixed_frame_C(np_clib,ffrme,
        #    bin_alts,range_cor_af_counts_float32,nc) 
        counts_ff[:,i1f,:] = put_bins_onto_fixed_frame(ffrme,bin_alts,range_cor_af_counts_float32)
    
        i = i + 1    # increment record counter by 1
        i1f = i1f + 1  # increment record counter for current file by 1
         
    if i == 0: continue
    
    # Apply polarization gain ratio to 1064 perpendicular channel
    print(attention_bar)
    print("Applying pgain factor to hard-coded channel of indicies 1,3")
    print("This should represent the 1064 nm channel")
    print(attention_bar)
    counts_ff[1,:,:] = counts_ff[1,:,:] * PGain[1] # 1064 nm perp
    counts_ff[3,:,:] = counts_ff[3,:,:] * PGain[0] # 355 nm perp
    
    # Compute NRB & populate saturate_ht array
    EMs = convert_raw_energy_monitor_values(MCS_data_1file['meta']['EM'],nwl,'Roscoe')
    
    ff_bg_st_bin = np.argwhere(ffrme <= bg_st_alt)[0][0]
    ff_bg_ed_bin = np.argwhere(ffrme >= bg_ed_alt)[-1][0]
    for chan in range(0,nc): 
        EMsubmit = EMs[:,wl_map[chan]]
        # Patch over bad EM readings
        EMsubmit = patch_outliers(EMsubmit,Estds)
        NRB[chan,:,:] = correct_raw_counts(counts_ff[chan,:,:],EMsubmit,None,
            i1f,nb_ff,ff_bg_st_bin,ff_bg_ed_bin,'NRB_no_range')          

    # Delete bad records from output arrays
    tot_good_recs = good_rec_bool.sum()
    if tot_good_recs < nr_1file:
        MCS_UnixT_float64_new = MCS_UnixT_float64_new[good_rec_bool]
        Nav_save = Nav_save[good_rec_bool]
        laserspot = laserspot[good_rec_bool,:]
        ONA_save = ONA_save[good_rec_bool]
        DEM_nadir = DEM_nadir[good_rec_bool]
        DEM_nadir_surftype = DEM_nadir_surftype[good_rec_bool]
        DEM_laserspot = DEM_laserspot[good_rec_bool]
        DEM_laserspot_surftype = DEM_laserspot_surftype[good_rec_bool]
        EMs = EMs[good_rec_bool,:]
        NRB = NRB[:,good_rec_bool,:]
        bg_save = bg_save[:,good_rec_bool]
        saturate_ht = saturate_ht[:,good_rec_bool]
        print('\nXXXXXXXXXXXXXXXXXXXXXXX')
        deleted_recs = nr_1file - tot_good_recs
        print(str(deleted_recs).strip()+' records deleted!')
        print('From file: '+MCS_file)
        print('\nXXXXXXXXXXXXXXXXXXXXXXX')
        i = i - (nr_1file - tot_good_recs)
        nr_1file = tot_good_recs
    if tot_good_recs == 0:
        print(attention_bar)
        print("No records in MCS_file: " + MCS_file)
        print("meet the cutoff criteria.")
        print(attention_bar)
        continue           

    # Average the data
    # Remainder data at end of each file is carried-over to next file
    #
    # Here's the logic behind this averging section:
    # The data are averaged one file at a time, with any leftover records
    # carried-over until the next file is processed. The carried-over records
    # are then averaged with the first 'x' records of the next file, where 'x'
    # represents the number of records required to meet the time-average criterion.
    # The logic is a bit tricky when processing the first and last files, and these
    # details are dealt with using the "first_read" and "cutbegin" variables. Also,
    # n_expand (how much to expand datasets) and expanded_length variables are
    # computed differently for the edge versus middle files. The reason averaging
    # is done one file at a time is to limit usage of RAM on the computer.
    if secs2avg > 0.0: # if < 0, don't average

        bin_numbers = np.digitize(MCS_UnixT_float64_new,time_bins)
        u, ui, ncounts = np.unique(bin_numbers,return_index=True,return_counts=True)
        Nav_save_avg = np.zeros(u.shape[0],dtype=Nav_save.dtype)
        laserspot_avg = np.zeros((u.shape[0],2),dtype=laserspot.dtype)
        ONA_save_avg = np.zeros(u.shape[0],dtype=ONA_save.dtype)
        DEM_nadir_avg = np.zeros(u.shape[0],dtype=DEM_nadir.dtype)
        DEM_nadir_surftype_avg = np.zeros(u.shape[0],dtype=DEM_nadir_surftype.dtype)
        DEM_laserspot_avg = np.zeros(u.shape[0],dtype=DEM_laserspot.dtype)
        DEM_laserspot_surftype_avg = np.zeros(u.shape[0],dtype=DEM_laserspot_surftype.dtype)
        EMs_avg = np.zeros((u.shape[0],nwl),dtype=EMs.dtype)
        NRB_avg = np.zeros((nc,u.shape[0],nb_ff),dtype=NRB.dtype)
        bg_save_avg = np.zeros((nc,u.shape[0]),dtype=bg_save.dtype)
        saturate_ht_max = np.zeros((nc,u.shape[0]),dtype=saturate_ht.dtype)
        rr = 0 # raw record number
        perfectly_aligned = False # True if trans_total == 0
    
        ei = ui.shape[0]-1
        if last_file: ei = ui.shape[0]
        if first_read: 
            si = 0
        # Gotta do an elif cuz cur file might not have any vals within last bin of prev file
        elif time_bins[bin_numbers[ui[0]]-1] == trans_bin[0]:
            si = 1 # start index
            trans_total = float(trans_ncounts+ncounts[0])
            for field in Nav_save_avg.dtype.names:
                if (field == 'UTC'): continue
                Nav_save_avg[field][0] = ( (Nav_save_sum[field] + 
                    Nav_save[field][rr:rr+ncounts[0]].sum())/trans_total )
            Nav_save_avg['UTC'] = Nav_UTC_carryover
            laserspot_avg[0,:] = (laserspot_sum + laserspot[rr:rr+ncounts[0],:].sum(axis=0))/trans_total
            ONA_save_avg[0] = (ONA_save_sum + ONA_save[rr:rr+ncounts[0]].sum())/trans_total
            DEM_nadir_avg[0] = (DEM_nadir_sum + DEM_nadir[rr:rr+ncounts[0]].sum())/trans_total
            DEM_nadir_surftype_avg[0] = stats.mode( np.concatenate((DEM_nadir_surftype_carryover, DEM_nadir_surftype[rr:rr+ncounts[0]])) )[0][0]
            DEM_laserspot_avg[0] = (DEM_laserspot_sum + DEM_laserspot[rr:rr+ncounts[0]].sum())/trans_total
            DEM_laserspot_surftype_avg[0] = stats.mode( np.concatenate((DEM_laserspot_surftype_carryover, DEM_laserspot_surftype[rr:rr+ncounts[0]])) )[0][0]
            EMs_avg[0,:] = (EMs_sum + EMs[rr:rr+ncounts[0],:].sum(axis=0))/trans_total
            NRB_avg[:,0,:] = (NRB_sum + NRB[:,rr:rr+ncounts[0],:].sum(axis=1))/trans_total
            bg_save_avg[:,0] = (bg_save_sum + bg_save[:,rr:rr+ncounts[0]].sum(axis=1))/trans_total
            saturate_ht_max[:,0] = np.asarray( (saturate_ht_carryover.max(axis=1), saturate_ht[:,rr:rr+ncounts[0]].max(axis=1)) ).max(axis=0)
            print('trans_ncounts = ',trans_ncounts)
            rr += ncounts[0]
        else:
            si = 0
            ei = ui.shape[0]
            perfectly_aligned = True
            print(attention_bar)
            print("I guess the time_bins lined up perfectly with edge of previous file")
            print("because there are no values in the previous file's last time bin.")
            print(trans_bin)
            print(attention_bar)

        for tb in range(si,ei):
            for field in Nav_save_avg.dtype.names:
                if (field == 'UTC'): continue
                Nav_save_avg[field][tb] = np.mean(Nav_save[field][rr:rr+ncounts[tb]])
            Nav_save_avg['UTC'][tb] = Nav_save['UTC'][rr+int(ncounts[tb]/2)] # midpoint
            laserspot_avg[tb,:] = np.mean(laserspot[rr:rr+ncounts[tb],:],axis=0)
            ONA_save_avg[tb] = np.mean(ONA_save[rr:rr+ncounts[tb]])
            DEM_nadir_avg[tb] = np.mean(DEM_nadir[rr:rr+ncounts[tb]])
            DEM_nadir_surftype_avg[tb] = stats.mode(DEM_nadir_surftype[rr:rr+ncounts[tb]])[0][0]
            DEM_laserspot_avg[tb] = np.mean(DEM_laserspot[rr:rr+ncounts[tb]])
            DEM_laserspot_surftype_avg[tb] = stats.mode(DEM_laserspot_surftype[rr:rr+ncounts[tb]])[0][0]
            EMs_avg[tb,:] = np.mean(EMs[rr:rr+ncounts[tb],:],axis=0)
            NRB_avg[:,tb,:] = np.mean(NRB[:,rr:rr+ncounts[tb],:],axis=1)	    
            bg_save_avg[:,tb] = np.mean(bg_save[:,rr:rr+ncounts[tb]],axis=1)
            saturate_ht_max[:,tb] = np.max(saturate_ht[:,rr:rr+ncounts[tb]],axis=1)
            #print('iter: ',u[tb],MCS_UnixT_float64_new[rr])
            rr += ncounts[tb]

        # Save the sum and ncounts of last populated time bin to string averaging
        # together between multiple files.
        if ((not last_file) and (ei != ui.shape[0])):
            Nav_save_sum = np.zeros(1,dtype=Nav_save.dtype)
            for field in Nav_save_avg.dtype.names:
                if (field == 'UTC'): continue
                Nav_save_sum[field][0] = Nav_save[field][rr:rr+ncounts[-1]].sum()
            Nav_UTC_carryover = Nav_save['UTC'][rr+ncounts[-1]-1]
            laserspot_sum = laserspot[rr:rr+ncounts[-1],:].sum(axis=0)
            ONA_save_sum = ONA_save[rr:rr+ncounts[-1]].sum()
            DEM_nadir_sum = DEM_nadir[rr:rr+ncounts[-1]].sum()
            DEM_nadir_surftype_carryover = DEM_nadir_surftype[rr:rr+ncounts[-1]]
            DEM_laserspot_sum = DEM_laserspot[rr:rr+ncounts[-1]].sum()
            DEM_laserspot_surftype_carryover = DEM_laserspot_surftype[rr:rr+ncounts[-1]]
            EMs_sum = EMs[rr:rr+ncounts[-1],:].sum(axis=0)
            NRB_sum = NRB[:,rr:rr+ncounts[-1],:].sum(axis=1)
            bg_save_sum = bg_save[:,rr:rr+ncounts[-1]].sum(axis=1)
            saturate_ht_carryover = saturate_ht[:,rr:rr+ncounts[-1]]
            # following line defines the "trasfer" bin; trans_bin
            trans_bin = time_bins[ bin_numbers[ ui[-1] ]-1 : bin_numbers[ ui[-1] ]+1 ]
            trans_ncounts = ncounts[-1]
            
        # Eliminate those avg'd profiles that contain less than min_avg_profs raw
        # profiles right here, exclusively in this single paragraph.
        big_enough_mask = ncounts >= min_avg_profs
        big_enough_mask[0] = True
        big_enough_mask[-1] = True
        if ((not first_read) and (not perfectly_aligned) and (trans_total < min_avg_profs)): big_enough_mask[0] = False
        if big_enough_mask.sum() < ncounts.shape[0]:
            Nav_save_avg = Nav_save_avg[big_enough_mask]
            laserspot_avg = laserspot_avg[big_enough_mask,:]
            ONA_save_avg = ONA_save_avg[big_enough_mask]
            DEM_nadir_avg = DEM_nadir_avg[big_enough_mask]
            DEM_nadir_surftype_avg = DEM_nadir_surftype_avg[big_enough_mask]
            DEM_laserspot_avg = DEM_laserspot_avg[big_enough_mask]
            DEM_laserspot_surftype_avg = DEM_laserspot_surftype_avg[big_enough_mask]
            EMs_avg = EMs_avg[big_enough_mask,:]
            NRB_avg = NRB_avg[:,big_enough_mask,:]
            bg_save_avg = bg_save_avg[:,big_enough_mask]
            saturate_ht_mask = saturate_ht_max[:,big_enough_mask]
            ncounts = ncounts[big_enough_mask]
            ui = ui[big_enough_mask]
            u = u[big_enough_mask]
            print("\nAvg'd profiles eliminated due to min_avg_profs constraint.")
            print(np.argwhere(big_enough_mask == False))
            print(big_enough_mask.shape," reduced to ", u.shape, "\n")

        # Expand dataset sizes to accomodate next input MCS file
        cutbegin = 0
        if ((first_read) and (ncounts[0] < min_avg_profs)): cutbegin = 1
        n_expand = u.shape[0]-1-cutbegin
        if ( (last_file) and (ncounts[-1] > min_avg_profs) ): n_expand = u.shape[0]
        expanded_length = nav_dset.shape[0]+n_expand
        # Now, if it's the first_read, nav_dset has an initialized length of 1; therefore,
        # if you use the expanded_length in the previous line the first_read, you'll 
        # get a dataset size that is too long by 1. The following line fixes this.
        if first_read: expanded_length = n_expand
        nav_dset.resize(expanded_length, axis=0)
        laserspot_dset.resize(expanded_length, axis=0)
        ONA_dset.resize(expanded_length, axis=0)
        DEM_nadir_dset.resize(expanded_length, axis=0)
        DEM_nadir_surftype_dset.resize(expanded_length, axis=0)
        DEM_laserspot_dset.resize(expanded_length, axis=0)
        DEM_laserspot_surftype_dset.resize(expanded_length, axis=0)
        EM_dset.resize(expanded_length, axis=0)
        NRB_dset.resize(expanded_length, axis=1)
        bg_dset.resize(expanded_length, axis=1)
        saturate_ht_dset.resize(expanded_length, axis=1)	
                        
        # Write out one file's worth of data to the hdf5 file
        nav_dset[expanded_length-n_expand:expanded_length] = Nav_save_avg[cutbegin:n_expand+cutbegin]
        laserspot_dset[expanded_length-n_expand:expanded_length,:] = laserspot_avg[cutbegin:n_expand+cutbegin,:]
        ONA_dset[expanded_length-n_expand:expanded_length] = ONA_save_avg[cutbegin:n_expand+cutbegin]
        DEM_nadir_dset[expanded_length-n_expand:expanded_length] = DEM_nadir_avg[cutbegin:n_expand+cutbegin]
        DEM_nadir_surftype_dset[expanded_length-n_expand:expanded_length] = DEM_nadir_surftype_avg[cutbegin:n_expand+cutbegin]
        DEM_laserspot_dset[expanded_length-n_expand:expanded_length] = DEM_laserspot_avg[cutbegin:n_expand+cutbegin]
        DEM_laserspot_surftype_dset[expanded_length-n_expand:expanded_length] = DEM_laserspot_surftype_avg[cutbegin:n_expand+cutbegin]    
        EM_dset[expanded_length-n_expand:expanded_length,:] = EMs_avg[cutbegin:n_expand+cutbegin,:]
        NRB_dset[:,expanded_length-n_expand:expanded_length,:] = NRB_avg[:,cutbegin:n_expand+cutbegin,:]
        bg_dset[:,expanded_length-n_expand:expanded_length] = bg_save_avg[:,cutbegin:n_expand+cutbegin]
        saturate_ht_dset[:,expanded_length-n_expand:expanded_length] = saturate_ht_max[:,cutbegin:n_expand+cutbegin]
        nrecs = expanded_length

    else: # No averaging

        # Expand dataset sizes to accomodate next input MCS file
        nav_dset.resize(i, axis=0)
        laserspot_dset.resize(i, axis=0)
        ONA_dset.resize(i, axis=0)
        DEM_nadir_dset.resize(i, axis=0)
        DEM_nadir_surftype_dset.resize(i, axis=0)
        DEM_laserspot_dset.resize(i, axis=0)
        DEM_laserspot_surftype_dset.resize(i, axis=0)
        EM_dset.resize(i, axis=0)
        NRB_dset.resize(i, axis=1)
        bg_dset.resize(i, axis=1)
        saturate_ht_dset.resize(i, axis=1)
                        
        # Write out one file's worth of data to the hdf5 file
        nav_dset[i-nr_1file:i] = Nav_save
        laserspot_dset[i-nr_1file:i,:] = laserspot
        ONA_dset[i-nr_1file:i] = ONA_save
        DEM_nadir_dset[i-nr_1file:i] = DEM_nadir
        DEM_nadir_surftype_dset[i-nr_1file:i] = DEM_nadir_surftype
        DEM_laserspot_dset[i-nr_1file:i] = DEM_laserspot
        DEM_laserspot_surftype_dset[i-nr_1file:i] = DEM_laserspot_surftype    
        EM_dset[i-nr_1file:i,:] = EMs
        NRB_dset[:,i-nr_1file:i,:] = NRB
        bg_dset[:,i-nr_1file:i] = bg_save
        saturate_ht_dset[:,i-nr_1file:i] = saturate_ht
        nrecs = i

    first_read = False
        
    print('\n**********************************')
    print('\nDone with file: '+MCS_file+'\n')
    print('**********************************\n')

      
print('Main L1A execution finished at: ',DT.datetime.now())

# Write out any final parameters to the HDF5 file that needed to wait

num_recs_dset[:] = nrecs
hdf5_file.close()
GEOS_fobj.close()
print("NRB has been written to the HDF5 file:"+hdf5_fname)

# The final message

print('Total raw profiles processed: '+str(i))
print("roscoe_l1a.py has finished normally.")
