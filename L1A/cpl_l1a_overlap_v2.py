# DESCRIPTION:
#
# This code, derived from "cpl_l1a_v1.py," creates the XDR overlap table
# required for CPL. Set "raob_file" before running.
#
# v1 interpolates a raob profile to each lidar bin for each lidar profile
#
# v2 uses one of our standard ACLS atmospheric profiles and matches each
# lidar bin to the closest ACLS for each lidar profile.

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
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
# Libraries I did create
from initializations import *
from read_routines import *
from lidar import *
from time_conversions import *

# -----------------> SET RAOB FILE HERE <-------------------------------
std_atm_dir = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\Config\\'
std_atm_file = std_atm_dir + 'arctic.winter.acls'
stop_raob_ht = 30.0 # km?
sigma = 2.0 # smoothing stddevs for AMB
sratm = 8.0*np.pi/3.0 # molecular lidar ratio
dz = vrZ

# -----------------> SET RAW OL OUTPUT FILE <---------------------------
out_dir = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\algorithm_dev\\overlap\\'
raw_ol_file = out_dir + 'raw_OL_values_' + proj_name + '_' + flt_date + '.txt'

######### READ IN RAYLEIGH & COMPUTE AMB #########

with open(std_atm_file,'r') as f_obj:
    Ray_data = f_obj.readlines()

Bray_alt = []
Bray355 = []
Bray532 = []
Bray1064 = []
for line in Ray_data:
    linesplit = line.split()
    if len(linesplit) > 7: continue
    Bray_alt.append(float(linesplit[0]))
    Bray355.append(float(linesplit[4]))
    Bray532.append(float(linesplit[5]))
    Bray1064.append(float(linesplit[6]))

Bray_alt = np.asarray(Bray_alt)
Bray355 = np.asarray(Bray355)
Bray532 = np.asarray(Bray532)
Bray1064 = np.asarray(Bray1064)
mtsq_slant355 = np.zeros(Bray355.shape[0])
mtsq_slant532 = np.zeros(Bray532.shape[0])
mtsq_slant1064 = np.zeros(Bray1064.shape[0])

odm355 = - 0.5 * np.log(1.0)
odm532 = - 0.5 * np.log(1.0)
odm1064 = - 0.5 * np.log(1.0)
for j in range(0,Bray355.shape[0]):
    odm355 = odm355 + Bray355[j] * sratm * dz/1e3
    odm532 = odm532 + Bray532[j] * sratm * dz/1e3
    odm1064 = odm1064 + Bray1064[j] * sratm * dz/1e3
    mtsq_slant355[j] = (np.exp(-2.0*odm355))#**(1.0/np.cos(ONA_samp))
    mtsq_slant532[j] = (np.exp(-2.0*odm532))#**(1.0/np.cos(ONA_samp))
    mtsq_slant1064[j] = (np.exp(-2.0*odm1064))#**(1.0/np.cos(ONA_samp))
AMB355 = Bray355 * mtsq_slant355 
AMB532 = Bray532 * mtsq_slant532
AMB1064 = Bray1064 * mtsq_slant1064


# Define functions <----------------------------------------------------

def patch_outliers(data_array, m):
    """ Function to overwrite outliers, in a numpy array, with reasonable
        values.
        Originally written to path over bad CPL energy monitor values.
    """
    
    # INPUTS: data_array -> self explanatory
    #         m          -> The number of standard deviations
    
    # RETURNS: A patched version of data_array with the same shape and 
    #          number of elements.
    
    keep_indxs = abs(data_array - np.mean(data_array)) < m * np.std(data_array)
    single_reasonable_value = np.mean( data_array[keep_indxs] )
    data_array[np.invert(keep_indxs)] = single_reasonable_value
    return data_array

def compute_time_offset_IWG1_CLS(cls_files):
    """ This function computes the time offset between the IWG1 Nav time
        and the CLS Instrument time ("ExactTime").

        Takes the string, full-path names of a list of CLS files as its 
        only argument. Returns the time offset in seconds.

        Add this offset to the Instrument time.
    """

    fc = 0
    time_offsets = []

    while fc < 2:

        cls_data = read_in_cls_data(cls_files[fc].strip())

        diff = cls_data['meta']['Nav']['UTC_Time'] - cls_data['meta']['Header']['ExactTime']
        diff_float_secs = np.zeros(diff.shape[0],dtype=np.float32)
        for i in range(0,diff.shape[0]):
            diff_float_secs[i] = diff[i].total_seconds()

        time_offsets.append(diff_float_secs.mean())

        fc += 1

    print("Computed time offsets: ",time_offsets)

    if abs(time_offsets[0] - time_offsets[1]) > 2.0:
        print("The time offsets between 2 sample files differs by too much.")
        print("Enter 0 for first time offset, 1 for second.")
        print("Time offsets: ",time_offsets)
        choice = int( input("Enter your choice.") )
        time_offset = time_offsets[choice]
    else:
        time_offset = np.mean(time_offsets)

    if abs(time_offset) > 80.0:
        print(attention_bar)
        print("!!!!!!! ATTENTION !!!!!!!")
        print("The computer time offset is "+str(time_offset).strip()+" seconds.")
        print("!!!!!!! ATTENTION !!!!!!!")
        print(attention_bar)
        #pdb.set_trace()
    
    return time_offset

def create_cls_interp_unixt_array(single_cls_file,cls_meta_data_all,
    firstlast_flag,Navi,firstlast_trunc=None):
    """ This function provides the "INTERP" and "ORIG"
        arrarys required for the "map_interp_times_to_orig_frame"
        function in the C library.
        The way in which it does this varies depending on whether 
        your Nav_source is 'cls' or something else.
    """

    # INPUTS:
    #
    # single_cls_file -> String. Name of one CLS file.
    # cls_meta_data_all -> Structure array containing the "straightened" 
    #                      Nav records. If Nav_source != 'cls' set to None. 
    # firstlast_flag  -> -1:first file of data set (flight),
    #                    1:last file of data set (flight),
    #                    0: neither the first nor last file
    # Navi -> List, [straightened Nav start index, end index]
    #         Set to anything if Nav_source != 'cls'
    # firstlast_trunc -> A list, [number of records to delete
    #                    from the beginning of the first file,
    #                    number of records to delete from the end
    #                    of the last file]. Only require for 
    #                    firstlast_flag = -1 or 1. Again, only applies
    #                    to Nav_source == 'cls.' Set to None otherwise.
    #
    # OUTPUTS:
    # 
    # interp_UnixT -> The CLS time, converted to UnixT (standard), interp'd
    #                 to data rate (set in init. file).
    # CLS_UnixT_float64_orig -> The original CLS Nav times.
    # cls_data_1file -> The CLS data from one file, corrected for Nav
    #                   time fluctuations ("straightened").
    #


    # Load the data...

    cls_data_1file = read_in_cls_data(single_cls_file)
    nr0 = cls_data_1file.shape[0]

    if Nav_source == 'cls':

        if firstlast_flag == -1: 
            cls_data_1file = cls_data_1file[firstlast_trunc[0]:]
        if firstlast_flag == 1:
            cls_data_1file = cls_data_1file[:nr0-firstlast_trunc[1]]
        cls_data_1file['meta']['Nav'] = cls_meta_data_all['Nav'][Navi[0]:Navi[1]]
        cls_nav_data_1file = cls_data_1file['meta']['Nav']
        del_t = np.zeros(cls_nav_data_1file.shape[0],dtype=DT.datetime)
        del_t[1:] = cls_nav_data_1file['UTC_Time'][1:] - cls_nav_data_1file['UTC_Time'][:-1]
        del_t[0] = del_t[1]*0.0
    
        # Compute # of elapsed seconds since first record (float)
    
        del_secs = np.zeros(cls_nav_data_1file.shape[0],dtype=np.float64)
        for k in range(0,cls_nav_data_1file.shape[0]): del_secs[k] = del_t[k].total_seconds()
    
        elapsed_secs = np.cumsum(del_secs)
        tot_elapsed = np.max(elapsed_secs)
        UnixT_epoch = DT.datetime(1970,1,1)

        # Also, you'll need "Unix Times" corresponding to the original sampling rate

        CLS_UnixT_float64_orig = (cls_nav_data_1file['UTC_Time'][0]-UnixT_epoch).total_seconds() + elapsed_secs

        interp_UnixT = np.copy(CLS_UnixT_float64_orig)

    else:

        CLS_UnixT_float64_orig = np.zeros(cls_data_1file.shape[0],dtype=np.float64)
        UnixT_epoch = DT.datetime(1970,1,1)

        for k in range(0,cls_data_1file.shape[0]):
            CLS_UnixT_float64_orig[k] = (cls_data_1file['meta']['Header']['ExactTime'][k]-UnixT_epoch).total_seconds()

        xp = CLS_UnixT_float64_orig - CLS_UnixT_float64_orig[0]
        fp = CLS_UnixT_float64_orig
        ix = np.arange(np.min(xp),np.max(xp)+(1.0/CLS_hz),1.0/CLS_hz)        
        interp_UnixT = np.interp(ix,xp,fp)

    return [interp_UnixT,CLS_UnixT_float64_orig,cls_data_1file]
	

    
# Start main execution here <-------------------------------------------
print('Starting main L1A execution at: ',DT.datetime.now())

# Create and load file list for CLS data
CLS_file_list = 'processing_file_list.txt'
search_str = '*.cls'
create_a_file_list(CLS_file_list,search_str)
with open(CLS_file_list) as CLS_list_fobj:
    all_CLS_files = CLS_list_fobj.readlines()
nCLS_files = len(all_CLS_files)

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

# Load nav data for entire flight. Select from 1 of several possible
# nav data streams. All these read-in all flight data at once. This
# decision was made to simplify the code and perhaps boost processing speed.
# NOTE that for this CPL version of L1A, you are able to use files 
# generated by CAMAL/ACATS and IWG1; however, as of 4/5/18, this hasn't
# been tested and there may be time offsets to figure out first.

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
    m = 1.0 / CLS_hz # basically, the time rez in secs you'd like to get to
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
    #       refer to CAMAL/ACTS "nav*.data" files.
    
    # "Dummy out" these variables, which are only need valid values for CLS-Nav processing
    cls_meta_data_all = None
    nav2cls_indx = np.zeros((1000,2))
    FL_trunc = None 
    usable_file_range = [0,1000]

    # Compute offset between IWG1 time and CLS Instrument Time ("ExactTime")
    # This data is the IWG1 data.
    iwg1_cls_t_offset = compute_time_offset_IWG1_CLS(all_CLS_files)    

    # Load the entire nav dataset into memory
    nav_data_all = read_entire_nav_dataset('*.nav*')
    
    # Interpolate data records to match CPL's data rate.
    # 10 Hz as of 2/2/18. Data rate set in initializations.
    
    del_t = np.zeros(nav_data_all.shape[0],dtype=DT.datetime)
    del_t[1:] = nav_data_all['UTC'][1:] - nav_data_all['UTC'][:-1]
    del_t[0] = del_t[1]*0.0
    pdb.set_trace()
    
    # Compute # of elapsed seconds since first record (float)
    
    del_secs = np.zeros(nav_data_all.shape[0],dtype=np.float64)
    for k in range(0,nav_data_all.shape[0]): del_secs[k] = del_t[k].total_seconds()
    
    elapsed_secs = np.cumsum(del_secs)
    tot_elapsed = np.max(elapsed_secs)
    m = 1.0 / CLS_hz
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
        #nav_interp['UnixT'][k] = nav_data_all['UnixT'][0] + DT.timedelta(seconds=ix[k])
        # offsets in following line convert to CLS UTC
        nav_interp['UTC'][k] = ( nav_data_all['UTC'][0] + DT.timedelta(seconds=ix[k]) 
            + DT.timedelta(seconds=nudge) )
        Nav_interp_T_float64[k] = (nav_interp['UTC'][k] - UnixT_epoch).total_seconds()
    
    # NOW, set the array that will be used in processing, nav_data_all equal
    # to the interpolated array that was just created, nav_interp.

    nav_data_all = nav_interp
    
    # This variable will be the same no matter what the Nav_source
    Nav_hz = nav_hz   
    
elif Nav_source == 'iwg1':

    # "Dummy out" these variables, which are only need valid values for CLS-Nav processing
    cls_meta_data_all = None
    nav2cls_indx = np.zeros((1000,2))
    FL_trunc = None 
    usable_file_range = [0,1000]

    # Compute offset between IWG1 time and CLS Instrument Time ("ExactTime")
    iwg1_cls_t_offset = compute_time_offset_IWG1_CLS(all_CLS_files)

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
    m = 1.0 / CLS_hz
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
        # offsets in following line convert to CLS UTC 
        IWG1_interp['UTC'][k] = ( IWG1_data['UTC'][0] + DT.timedelta(seconds=ix[k]) )
        Nav_interp_T_float64[k] = (IWG1_interp['UTC'][k] - UnixT_epoch).total_seconds()

    # NOW, set the array that will be used in processing, IWG1_data equal
    # to the interpolated array that was just created, IWG1_interp.

    IWG1_data = IWG1_interp
    
    # This variable will be the same no matter what the Nav_source
    Nav_hz = IWG1_hz
    
elif Nav_source == 'cls':

    #  This block reads in Nav records from CLS files, then performs
    # the tasks as the other blocks above. This pre-processing is much more
    # involved than the other Nav sources because of the inconsistent sampling/
    # alignment of the Nav data in the CLS files.

    # Make this offset zero, since we're not using IWG1
    iwg1_cls_t_offset = 0.0

    # Load the entire nav dataset into memory. This data will be used in both
    # the 'cls' Nav block and later on in the main file loop.
    cls_meta_data_all, nav2cls_indx, FL_trunc, usable_file_range = read_entire_cls_meta_dataset()
    #meta_save = np.copy(cls_meta_data_all) #uncomment to compare processed to original
    cls_nav_data_all = np.copy(cls_meta_data_all['Nav'])

    UnixT_epoch = DT.datetime(1970,1,1)
    m = 1.0 / CLS_hz
    
    # Identify all the unique UTC_Times and their starting indexes
    u, ui, ncounts = np.unique(cls_nav_data_all['UTC_Time'],return_index=True,return_counts=True)
    # Bad time values will sink to the front of this "unique" array
    if u[0] == bad_cls_nav_time_value:
        print(attention_bar)
        print("\nSome(A) UTC_Time(s) in the Nav records were(was) invalid.")
        print("Attempt is being made to take care of this.") 
        print("Setting first invalid CLS-Nav UTC_Time record equal to")
        print("the first valid CLS-Nav UTC_Time record.\n")
        print(attention_bar)
        cls_nav_data_all['UTC_Time'][0] = u[1]
        u = u[1:]
        ui = ui[1:]
        ncounts = ncounts[1:]

    # I found at least one case where inst. clk. was correct, but data system appeared to be
    # grabbing an old Nav record for a bit. This block of code determines if this is the case
    # and overwrites those values with the first valid values.
    while (u[1]-u[0]).total_seconds() > 3:
        print("\nSome(A) UTC_Time(s) in the Nav records were(was) from stale packets.")
        print("Attempt is being made to take care of this.") 
        print("Setting first stale CLS-Nav UTC_Time record equal to")
        print("the first non-stale CLS-Nav UTC_Time record.\n")
        cls_nav_data_all['UTC_Time'][0] = u[1]
        u = u[1:]
        ui = ui[1:]
        ncounts = ncounts[1:]
    n_u = u.shape[0]      
    
    # Convert the unique UTC_Times to Unix Times
    Nav_UnixT = np.zeros(n_u,dtype=np.float64)
    Nav_unique = np.zeros(n_u,dtype=CLS_decoded_nav_struct)
    for i in range(0,n_u): 
        Nav_UnixT[i] = (cls_nav_data_all['UTC_Time'][ui[i]] - UnixT_epoch).total_seconds()
        Nav_unique[i] = cls_nav_data_all[ui[i]] 

    # Compute # of elapsed seconds since first record (float)    
    u, ui, ncounts = np.unique(cls_meta_data_all['Header']['ExactTime'],return_index=True,return_counts=True)
    if u[0] == bad_cls_nav_time_value:
        u = u[1:]
        ui = ui[1:]
        ncounts = ncounts[1:]
    tot_elapsed = (u.max() - u.min()).total_seconds()
    if tot_elapsed > (max_flt_hours*3600.0):
        tot_hours = tot_elapsed / 3600.0
        print(attention_bar)
        print("HEY! The data span "+str(tot_hours).strip()+" hours!")
        print("Consider trimming processing time range in init. file.")
        print('Enter "c" to contine processing anyway.')
        print('Enter "q" to safely quit.')
        print(attention_bar)
        pdb.set_trace() 

    # Look at the time deltas between unique, you'll use this in a little bit
    ExactTime_ncounts = np.copy(ncounts)
    ExactTime_ui = np.copy(ui)
    inst_clk_deltas = delta_datetime_vector(u)

    # Estimate the data rate from the instr. clock
    computed_CLS_hz = cls_meta_data_all.shape[0] / tot_elapsed
    if abs(computed_CLS_hz - CLS_hz) > CLS_hz_tol:
        print(attention_bar)
        print("!!!!!!! WARNING !!!!!!!")
        print("The estimated (init. file) and computed data rates differ by too much.")
        print("estimated: ",CLS_hz," computed: ",computed_CLS_hz)
        print("!!!!!!! WARNING !!!!!!!")
        print(attention_bar)
        #pdb.set_trace()
    else:
        print('Init. file CLS data rate is sufficiently close to computed data rate.')
        print('CLS_hz: ',CLS_hz, ' computed_CLS_hz: ',computed_CLS_hz)
        CLS_hz = computed_CLS_hz
        m = 1.0 / CLS_hz

    # Create x-ord array of Unix Times at which to compute interpolated values
    ix = np.arange(Nav_UnixT[0],Nav_UnixT[0]+tot_elapsed+1,m) # x coordinates of the interpolated values

    # Create the x-coordinates for the interp data, which will have units of Unix Time
    print("Making the reasonably dangerous assumption that earliest time is first rec")
    print("and latest time is last rec.")
    rough_expected_recs_sec = int(CLS_hz)
    interp_start_time_offset = 0.0
    if ncounts[0] < CLS_hz:
        interp_start_time_offset = (rough_expected_recs_sec - ncounts[0])/CLS_hz
        
    n_interp = ix.shape[0]
    nav_interp = np.zeros(n_interp,dtype=CLS_decoded_nav_struct)
    Nav_interp_T_float64 = ix
    for field in cls_nav_data_all.dtype.names:
        if (field == 'UTC_Time'): continue
        nav_interp[field] = np.interp(ix, Nav_UnixT, Nav_unique[field])
    # Now handle population of interp time field
    for k in range(0,n_interp):
        nav_interp['UTC_Time'][k] = ( cls_nav_data_all['UTC_Time'][0] + DT.timedelta(seconds=(m*k)) )
 
    # Prepare nav_interp 'UTC_Time' array, where time gaps have been
    # masked out, then overwrite original Nav time.
    offset_indx=np.abs((ix - ix[0]) - interp_start_time_offset).argmin()
    if inst_clk_deltas.max() > inst_clk_rez:
        
        print(attention_bar+"Time gaps in instrument clock detected. Handling it!"+attention_bar)
        UTC_ovrwrt = np.copy(nav_interp['UTC_Time'][offset_indx:])
        no_delta = True
        for k in range(0,inst_clk_deltas.shape[0]):
            if inst_clk_deltas[k] > inst_clk_rez:
                if no_delta: 
                    delta_indx = ExactTime_ui[k-1]+ExactTime_ncounts[k-1]
                    no_delta = False
                skip_recs = int(inst_clk_deltas[k]*CLS_hz)
                UTC_ovrwrt[delta_indx:delta_indx+skip_recs] = bad_cls_nav_time_value
                #print(delta_indx,ExactTime_ui[k],k)
                #print('A')
                #pdb.set_trace()
                delta_indx += skip_recs
            elif ((not no_delta) and (inst_clk_deltas[k] <= inst_clk_rez)):
                delta_indx += ExactTime_ncounts[k] #int(CLS_hz)
                #print(delta_indx,ExactTime_ui[k],k)
                #print('B')
                #pdb.set_trace()
        UTC_Time_use = np.extract(UTC_ovrwrt != bad_cls_nav_time_value, UTC_ovrwrt)
        try:
            cls_meta_data_all['Nav']['UTC_Time'] = UTC_Time_use[:cls_meta_data_all.shape[0]]
        except ValueError:
            num_short = str(cls_meta_data_all.shape[0] - UTC_Time_use.shape[0]).strip()
            print(attention_bar+'!!!!!!!!!! WARNING !!!!!!!!!!\n')
            print('UTC_Time_use is shorter than meta array by '+num_short+' elements.')
            print('\n!!!!!!!!!! WARNING !!!!!!!!!!'+attention_bar)
            cls_meta_data_all['Nav']['UTC_Time'][:UTC_Time_use.shape[0]] = UTC_Time_use
            cls_meta_data_all['Nav']['UTC_Time'][UTC_Time_use.shape[0]:] = UTC_Time_use[-1]  
    
    else:
        
        print('\nYou are lucky. No time gaps detected in the middle of the flight! :-)\n')
        cls_meta_data_all['Nav']['UTC_Time'] = nav_interp['UTC_Time'][offset_indx:offset_indx+cls_meta_data_all.shape[0]]

    #plt.plot_date(cls_meta_data_all['Nav']['UTC_Time'],cls_meta_data_all['Nav']['RollAngle'],marker='x')
    #plt.plot_date(nav_interp['UTC_Time'],nav_interp['RollAngle'],marker='o')
    #plt.show()
    #pdb.set_trace() 
    
    # NOW, set the array that will be used in processing, cls_nav_data_all equal
    # to the interpolated array that was just created, nav_interp.
    cls_nav_data_all = nav_interp
    
    # This variable will be the same no matter what the Nav_source
    Nav_hz = CLS_hz    
    
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
    
# Create an array to hold all the overlap for the whole duration    

OL = np.zeros((nCLS_files*file_len_recs,nbins,nwl),dtype=np.float64) - 999.9

print('Starting core processing...') # <--------------------------------
first_read = True
usable_file_indicies = range(usable_file_range[0],usable_file_range[1])
trans_bin = [0,0]
last_file = False 
for f in range(0,nCLS_files):

    if f not in usable_file_indicies:
        print(attention_bar)
        print("File number "+str(f).strip()+" has no usable data!")
        print("Skipping processing of this file.")
        print(attention_bar)
        continue
    
    # Read in CLS data of one file and initialize a few things
    CLS_file = all_CLS_files[f]
    CLS_file = CLS_file.rstrip()
    FL_flag = 0
    if f == usable_file_range[0]: FL_flag=-1
    if ((f == usable_file_range[1]-1) or (f == nCLS_files-1)): 
        FL_flag=1
        last_file = True
    [interp_UnixT,CLS_UnixT_float64_orig,CLS_data_1file] = create_cls_interp_unixt_array(
        CLS_file,cls_meta_data_all,FL_flag,nav2cls_indx[f,:],FL_trunc)
    CLS_UnixT_float64_orig = CLS_UnixT_float64_orig - nudge + iwg1_cls_t_offset
    interp_UnixT = interp_UnixT - nudge + iwg1_cls_t_offset
    if CLS_data_1file is None:
        print('\n******* Bad file! Skipping file!!! *******')
        print(CLS_file)
        continue  
    nr_1file = CLS_data_1file.shape[0]                          # # of records in current file
    good_rec_bool = np.ones(CLS_data_1file.shape[0],dtype=bool) # is 'True' for good records, 'False' for bad
    satur_flag = False # Will turn True if any bin in current file is potentially saturated.

    # Create histogram bin boundaries. If first_read, created farther down in first_read block.
    if not first_read: time_bins = np.arange(trans_bin[0],CLS_UnixT_float64_orig[-1]+3.0*secs2avg,secs2avg)

    # Deem as "bad," records outside user-specified processing time range
    if ( (CLS_data_1file['meta']['Nav']['UTC_Time'].min() < process_start) 
      or (CLS_data_1file['meta']['Nav']['UTC_Time'].max() > process_end) ):
        time_range_mask_low = CLS_data_1file['meta']['Nav']['UTC_Time'] > process_start
        time_range_mask_high = CLS_data_1file['meta']['Nav']['UTC_Time'] < process_end
        time_range_mask = time_range_mask_low * time_range_mask_high
        CLS_data_1file = np.extract(time_range_mask,CLS_data_1file)
        CLS_UnixT_float64_orig = np.extract(time_range_mask,CLS_UnixT_float64_orig)
        good_rec_bool = np.ones(CLS_data_1file.shape[0],dtype=bool) # is 'True' for good records, 'False' for bad
        print(attention_bar)
        print('In file '+CLS_file+'...')
        print(str(nr_1file - CLS_data_1file.shape[0]).strip()+' records are out of time range.')
        print(attention_bar)
        nr_1file = CLS_data_1file.shape[0]
        if (nr_1file == 0): 
            print(attention_bar)
            print("File number "+str(f).strip()+" has no usable data!")
            print("Skipping processing of this file.")
            print(attention_bar)
            continue
    
    # ---> START process of calling a C function <---
    #      This C function populates an array 
    #      mapped to the CLS Unix Time values
    #      with the appropriate interpolated values.
    #      This function also populates an array
    #      that maps the Nav_data indicies to
    #      the CLS data. This drastically speeds
    #      up the profile loop in python code.
    # "interp_UnixT" and "CLS_UnixT_float64_orig" should have been computed
    # earlier in the program, before this file loop.
    interp2orig_indicies = np.zeros(nr_1file,dtype=np.uint32)
    CLS_UnixT_float64_new = np.zeros(nr_1file,dtype=np.float64)
    CLS_UnixT_float64_new = np.require(CLS_UnixT_float64_new,float,['ALIGNED','C_CONTIGUOUS'])
    CLS_UnixT_float64_orig = np.require(CLS_UnixT_float64_orig,float,['ALIGNED','C_CONTIGUOUS'])
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
    np_clib.map_interp_times_to_orig_frame(CLS_UnixT_float64_new, CLS_UnixT_float64_orig, interp_UnixT,
        Nav_interp_T_float64, interp2orig_indicies, ctypes.c_double(1.0/Nav_hz),
        CLS_UnixT_float64_new.ctypes.strides, CLS_UnixT_float64_new.ctypes.shape, 
        CLS_UnixT_float64_orig.ctypes.strides, CLS_UnixT_float64_orig.ctypes.shape,
        interp_UnixT.ctypes.strides, interp_UnixT.ctypes.shape,
        Nav_interp_T_float64.ctypes.strides, Nav_interp_T_float64.ctypes.shape,  
        interp2orig_indicies.ctypes.strides, interp2orig_indicies.ctypes.shape) 
    # ---> END process of calling a C function <---                     

    # Save/create a few key data parameters...
    if first_read:
        nb =  nbins              # # of bins
        #vrZ set in config file. Not in raw data (CLS). Dennis confirmed. [4/6/18]
        nc = CLS_data_1file['meta']['Header']['NumChannels'][0] # # of channels
        #nshots set in config file. not in raw data (CLS) [4/6/18]
        vrT = 1e-7 #CLS_data_1file['meta']['binwid'][0]
        # Range vector, broadcast to nc x nb for vectorized operations later
        bins_range=np.broadcast_to(np.arange(vrZ/2.0,nb*vrZ+vrZ/2.0,vrZ,dtype=np.float32),(nc,nb)) # changed for overlap [7/17/19]
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
        time_bins = np.arange(CLS_UnixT_float64_orig[0],CLS_UnixT_float64_orig[-1]+3.0*secs2avg,secs2avg)
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
        # Open the lat/lon GEOS met data sampling CSV file
        GEOS_samp_fname = L1_dir+'lon_lat_UTC_'+proj_name+'_'+flt_date+'_'+Nav_source+'.txt'
        GEOS_fobj = open(GEOS_samp_fname, 'wt')
        GEOS_fobj.write('lon,lat,time\n') # file header     
        
    counts_ff = np.zeros((nc,nr_1file,nb),dtype=np.float32) # nb instead of nb_ff for overlap code
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
    AMB_save = np.zeros((nr_1file,nb,nwl),dtype=np.float32)-999.9

    # Right here, right now, test to see if any bin in current file is potentially saturated

    # Quick, broad test...
    satur_locs = []
    satur_locs_indx = 0
    if (CLS_data_1file['counts'].max() > min(saturation_values)): satur_flag = True
    if satur_flag:
        print('\n'+all_CLS_files[f]+' has saturated bin values in it!\n')
        # satur_locs will have [recs x chan x bin] as each row. recs dim should increase
        # towards higher indexes.
        satur_locs = np.argwhere(CLS_data_1file['counts'] > min(saturation_values))
                
    
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
                print("I am guessing that all Nav times are < min(CLS time)")
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
                print("I am guessing that all Nav times are < min(CLS time)")
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

        elif Nav_source == 'cls':

            # Find matching Nav data point. Use index map to match.
            # NOTE: Nav_save is an IWG1-style structure (see immutable_dat_structs.py)
            Nav_match = cls_nav_data_all[interp2orig_indicies[i1f]]
            Nav_save[i1f]['roll'] = Nav_match['RollAngle']
            Nav_save[i1f]['pitch'] = Nav_match['PitchAngle']
            #Nav_save[i1f]['drift'] = Nav_match['drift']
            Nav_save[i1f]['heading'] = Nav_match['TrueHeading']
            Nav_save[i1f]['lon'] = Nav_match['GPS_Longitude']
            Nav_save[i1f]['lat'] = Nav_match['GPS_Latitude']
            # The following line is necessary to get a clean, interpretable,
            # time to write out to the HDF5 file. It will essetially be a
            # byte array for IDL.
            time_str = Nav_match['UTC_Time'].strftime("%Y-%m-%dT%H:%M:%S.%f")
            Nav_save[i1f]['UTC'] = np.asarray(list(time_str.encode('utf8')),dtype=np.uint8)
            Nav_save[i1f]['GPS_alt'] = Nav_match['GPS_Altitude']
            Y = 0.0 #Nav_match['drift'] * (pi/180.0) 
            P = Nav_match['PitchAngle'] * (pi/180.0)
            R = Nav_match['RollAngle'] * (pi/180.0)                        
            # Write 'lat,lon,UTC' out to file
            GEOS_out_str = '{0:.4f},{1:.4f},{2}\n'.format(Nav_match['GPS_Longitude'],Nav_match['GPS_Latitude'],time_str[:19])
            GEOS_fobj.write(GEOS_out_str)
    
        # Proceed no farther in processing if altitude is not high enough
        if Nav_save[i1f]['GPS_alt'] < alt_cutoff:
            #print('Bad altitude: ',Nav_save[i1f]['GPS_alt'])
            #print('Will not process until good alt is read.\n')
            good_rec_bool[i1f] = False
            i = i + 1
            i1f = i1f + 1
            continue

        # Convert to float for subsequent manipulation within this loop
        counts_float32 = CLS_data_1file['counts'][i1f].astype(np.float32)
            
        # Apply dead time correction
        cc = 0 # count the channels
        for DTT_file in DTT_files:
            try:
                counts_float32[cc,:] = DTT[ cc, CLS_data_1file['counts'][i1f,cc,:] ]
            except IndexError:
                # Probably a noise spike
                spike_locs = np.where( counts_float32[cc,:] > DTT.shape[1] )
                counts_float32[cc,spike_locs[0][0]] = 0.0
            cc += 1
    
        # Calculate the off nadir angle
        # For CAMAL, positive scan angles are to the right
        # For CPL, assume phi0 and phi0_azi are zero
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
                   
        # Match top bin to closest standard atmosphere profile alt.
        # Then, assume small ONA so that you can assume each bin on down matches.
        TopBinMatch = np.argmin(np.abs(Bray_alt - bin_alts[0]))
        extent = Bray_alt.shape[0] - TopBinMatch # fixed CPL frame is 900 vs. 833 for actual
        if extent > nb: extent = nb
        CPLframe_extent = TopBinMatch + nb
        if CPLframe_extent > Bray_alt.shape[0]: CPLframe_extent = Bray_alt.shape[0]
        AMB_save[i1f,:extent,0] = AMB355[TopBinMatch:CPLframe_extent]
        AMB_save[i1f,:extent,1] = AMB532[TopBinMatch:CPLframe_extent]
        AMB_save[i1f,:extent,2] = AMB1064[TopBinMatch:CPLframe_extent]

        # Populate saturate_ht with bin alts if any bins in current profile are saturated
        # NOTE: The order of cond's matters in if blocks in this code paragraph
        if ((satur_flag) and (i1f in satur_locs[:,0])):
            while ((satur_locs_indx < satur_locs.shape[0]) and (satur_locs[satur_locs_indx,0] == i1f)):
                indx = satur_locs[satur_locs_indx,:] #Will always have 3 dims in same order
                if (CLS_data_1file['counts'][indx[0],indx[1],indx[2]] > saturation_values[indx[1]]):
                    saturate_ht[indx[1],i1f] = bin_alts[indx[2]]
                satur_locs_indx += 1

        length_bin_alts = bin_alts.shape[0]
        if length_bin_alts > nb:
            print('Length of bin_alts vectors is > than nbins. Stopping code...')
            pdb.set_trace()

        # Subtract background and correct raw counts for range (in scope of actual bin-alt frame)
        # It's difficult to apply range correction once counts are in fixed frame.
        try:
            bg_loc1 = np.argwhere(bin_alts <= bg_st_alt)[0][0]
            bg_loc2 = np.argwhere(bin_alts >= bg_ed_alt)[-1][0]	   
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
        
        # No overlap correction applied, cuz, duh.
        
    
        # Store the raw counts in the counts_ff variable for this overlap code.
        # NOTE: In the real L1A code, counts_ff represents the rebinned counts,
        # which doesn't happen in this overlap code.
        counts_ff[:,i1f,:] = range_cor_af_counts_float32
    
        i = i + 1    # increment record counter by 1
        i1f = i1f + 1  # increment record counter for current file by 1
         
    if i == 0: continue
    
    # Apply polarization gain ratio to 1064 perpendicular channel
    print(attention_bar)
    print("Applying pgain factor to hard-coded channel of index 3")
    print("This should represent the 1064 nm channel")
    print(attention_bar)
    counts_ff[3,:,:] = counts_ff[3,:,:] * PGain[1]
    
    # Patch over bad energy monitor values.
    # Different from classic "AM.pro" code ~line 954.
    BadEs = np.argwhere(CLS_data_1file['meta']['Engineering']['LaserEnergyMonitors'] <= 0)
    if BadEs.size > 0:
        medianEs = np.median(CLS_data_1file['meta']['Engineering']['LaserEnergyMonitors'],axis=0)
        for BadE in BadEs:
            # Just replace for all wavelengths, even if only one is bad.
            CLS_data_1file['meta']['Engineering']['LaserEnergyMonitors'][BadE[0],:] = medianEs
    
    # Compute NRB & populate saturate_ht array
    EMs = convert_raw_energy_monitor_values(CLS_data_1file['meta']['Engineering']['LaserEnergyMonitors'],nwl,'CPL',e_flg)
    ff_bg_st_bin = np.argwhere(ffrme <= bg_st_alt)[0][0]
    ff_bg_ed_bin = np.argwhere(ffrme >= bg_ed_alt)[-1][0]
    for chan in range(0,nc): 
        EMsubmit = EMs[:,wl_map[chan]]
        # Patch over bad EM readings
        EMsubmit = patch_outliers(EMsubmit,Estds)
        E = np.broadcast_to(EMsubmit,(nb,nr_1file)).transpose()
        # Calculate NRB right here. Don't pass to function.
        NRB[chan,:,:] = counts_ff[chan,:,:] / E

    # Delete bad records from output arrays
    tot_good_recs = good_rec_bool.sum()
    if tot_good_recs < nr_1file:
        CLS_UnixT_float64_new = CLS_UnixT_float64_new[good_rec_bool]
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
        AMB_save = AMB_save[good_rec_bool,:,:]
        print('\nXXXXXXXXXXXXXXXXXXXXXXX')
        deleted_recs = nr_1file - tot_good_recs
        print(str(deleted_recs).strip()+' records deleted!')
        print('From file: '+CLS_file)
        print('\nXXXXXXXXXXXXXXXXXXXXXXX')
        i = i - (nr_1file - tot_good_recs)
        nr_1file = tot_good_recs
    
    # Calculate the overlap
    OL[i-nr_1file:i,:,0] = AMB_save[:,:,0] / NRB[0,:,:]
    OL[i-nr_1file:i,:,1] = AMB_save[:,:,1] / NRB[1,:,:]
    OL[i-nr_1file:i,:,2] = AMB_save[:,:,2] / (NRB[2,:,:]+NRB[3,:,:])

    # This is the overlap code, no averaging or NRB output

    first_read = False
        
    print('\n**********************************')
    print('\nDone with file: '+CLS_file+'\n')
    print('**********************************\n')

# /////////////// FINAL OVERLAP BIN ARRAY COMPUTATION \\\\\\\\\\\\\\\\\\
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

# Trim off the edge of the OL array & filter out bad OL values

OL = OL[:i,:,:]
OLend = 285 # set an end for the overlap region
SDKeep = 2
# Remove any records where overlap is inf (or just not finite) in overlap region
OL355colmax = OL[:,:,0].max(axis=1)
OL532colmax = OL[:,:,1].max(axis=1)
OL1064colmax = OL[:,:,2].max(axis=1)
good355 = OL355colmax != np.inf
good532 = OL532colmax != np.inf
good1064 = OL1064colmax != np.inf
OL355_2D = OL[good355,:,0]
OL532_2D = OL[good532,:,1]
OL1064_2D = OL[good1064,:,2]
# Only retain SDKeep stddevs to retain per bin
for j in range(0,nbins):    
    OL355_2D[:,j] = patch_outliers(OL355_2D[:,j], SDKeep)
    OL532_2D[:,j] = patch_outliers(OL532_2D[:,j], SDKeep)
    OL1064_2D[:,j] = patch_outliers(OL1064_2D[:,j], SDKeep)

# Finesse and finalize the overlap calculation
OL355 = OL355_2D.mean(axis=0)
OL532 = OL532_2D.mean(axis=0)
OL1064 = OL1064_2D.mean(axis=0)
o5 = (1/OL532) / 4e+14 # Manually determined norm factor
o3 = (1/OL355) / 7e+12
o1 = (1/OL1064) / 1e13
# Write raw overlap values to a file
with open(raw_ol_file,'w') as f_obj:
    f_obj.write('OL 355 nm, OL 532 nm, OL 1064 nm\n')
    for a,b,c in zip(o3,o5,o1):
        f_obj.write(str(a)+','+str(b)+','+str(c)+'\n')
pdb.set_trace()
def fit_func(x,a,b,c,d):
    return a * np.arctan(b*x + c) + d
popt, pcov = curve_fit(fit_func, np.arange(0,OLend), o5[:OLend])
plt.plot( np.arange(0,OLend), fit_func(np.arange(0,OLend),*popt) )
plt.plot( np.arange(0,OLend), o5[:OLend] )
plt.show()
pdb.set_trace()
# ////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////
      
print('Main L1A execution finished at: ',DT.datetime.now())
pdb.set_trace()

# Close output files

GEOS_fobj.close()
print("NRB has been written to the HDF5 file:"+hdf5_fname)

# The final message

print('Total raw profiles processed: '+str(i))
print("camal_l1a.py has finished normally.")

######################################################################### BELOW THIS LINE SHALL NOT BE PART OF L1A PROCESS

tit = '355 nm NRB'
xlimits = [0,ONA_save.shape[0]]
ylimits = [810,0]#[900,500]
samp_chan = NRB[0,:,:]# + NRB[3,:,:]
curtain_plot(samp_chan.transpose(), nb_ff, vrZ_ff, ffrme, 0, 5e8, hori_cap, pointing_dir,
                      figW, figL, CPpad, 'records', 'alt (m)', tit, 'alt', ylimits, 'recs', 
                      xlimits, scale_alt_OofM, 1, out_dir)
pdb.set_trace()

make_custom_plot(samp_chan.transpose(), nb_ff, vrZ_ff, ffrme, 0, 1e8, hori_cap, pointing_dir,
                      figW, figL, CPpad, 'records', 'altitude(m)', tit, 'alt',  
                      [ylimits[0],ylimits[1]], 'recs', [xlimits[0],xlimits[1]], scale_alt_OofM, 1, out_dir, str(f).strip()+'.png',
                      np.arange(xlimits[0],xlimits[1]),Nav_save['pitch'][xlimits[0]:xlimits[1]],
                      np.arange(xlimits[0],xlimits[1]),ONA_save[xlimits[0]:xlimits[1]]*(180.0/np.pi),
                      np.arange(xlimits[0],xlimits[1]),Nav_save['roll'][xlimits[0]:xlimits[1]])
pdb.set_trace()
