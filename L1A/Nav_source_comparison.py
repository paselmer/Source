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
            print("There are not enough matches to compute IWG1 time offset.")
            print("You could try setting the initialization file to grab \
                   hard-coded values from the initialization file.")
            print("Time offsets are None! Code will stop!")
            time_offset_UnixT = None
            time_offset_IWG1 - None
            exit()
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
# Compute the time offset between MCS and IWG1
time_offset_UnixT, time_offset_IWG1 = compute_time_offsets(all_MCS_files[0].strip(),offset_choice)
print("Time offsets are: ", time_offset_UnixT, time_offset_IWG1)

# Load nav data for entire flight. Select from 1 of several possible
# nav data streams. All these read-in all flight data at once. This
# decision was made to simplify the code and perhaps boost processing speed.

if True:

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
    gps_interp_T64 = Nav_interp_T_float64
    gps_data_all = gps_interp
    # Remember, the UTC you'll need later is contained in separate array,
    # gps_UTC_interp.
    
    # This variable will be the same no matter what the Nav_source
    Nav_hz = gps_hz


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
        nav_interp['UTC'][k] = nav_data_all['UTC'][0] + DT.timedelta(seconds=ix[k]) - time_offset_IWG1 + DT.timedelta(seconds=nudge)
    
    # NOW, set the array that will be used in processing, nav_data_all equal
    # to the interpolated array that was just created, nav_interp.
    
    nav_interp_T64 = Nav_interp_T_float64
    nav_data_all = nav_interp
    
    # This variable will be the same no matter what the Nav_source
    Nav_hz = nav_hz    
    
anchor_loc1 = np.where(gps_interp_T64 == nav_interp_T64[0])
anchor_loc2 = np.where(gps_interp_T64 == nav_interp_T64[-1])
gps_UTC_interp = gps_UTC_interp[anchor_loc1[0][0]:anchor_loc2[0][0]+1]
gps_data_all = gps_data_all[anchor_loc1[0][0]:anchor_loc2[0][0]+1]


roll_diff = nav_data_all['roll'] - gps_data_all['roll']
pitch_diff = nav_data_all['pitch'] - gps_data_all['pitch']
alt_diff = nav_data_all['GPS_alt'] - gps_data_all['alt']
lat_diff = nav_data_all['lat'] - gps_data_all['lat']
lon_diff = nav_data_all['lon'] - gps_data_all['lon']

ax1 = plt.subplot(2,1,1)
plt.plot_date(nav_data_all['UTC'],nav_data_all['roll'],marker='o')
plt.plot_date(gps_UTC_interp,gps_data_all['roll'],marker='x')
ax1.xaxis.set_visible(False)
plt.ylabel('roll (deg)')
plt.title('IWG1 roll (blue) and GPS roll (orange): '+flt_date)
plt.subplot(2,1,2,sharex=ax1)
plt.plot_date(gps_UTC_interp,roll_diff)
plt.xlabel('UTC')
plt.ylabel('IWG1 roll minus GPS roll (deg)')
plt.title('Difference between IWG1 roll and GPS roll: '+flt_date)
plt.show()

ax1 = plt.subplot(2,1,1)
plt.plot_date(nav_data_all['UTC'],nav_data_all['pitch'],marker='o')
plt.plot_date(gps_UTC_interp,gps_data_all['pitch'],marker='x')
ax1.xaxis.set_visible(False)
plt.ylabel('pitch (deg)')
plt.title('IWG1 pitch (blue) and GPS pitch (orange): '+flt_date)
plt.subplot(2,1,2,sharex=ax1)
plt.plot_date(gps_UTC_interp,pitch_diff)
plt.xlabel('UTC')
plt.ylabel('IWG1 pitch minus GPS pitch (deg)')
plt.title('Difference between IWG1 pitch and GPS pitch: '+flt_date)
plt.show()

ax1 = plt.subplot(2,1,1)
plt.plot_date(nav_data_all['UTC'],nav_data_all['GPS_alt'],marker='o')
plt.plot_date(gps_UTC_interp,gps_data_all['alt'],marker='x')
ax1.xaxis.set_visible(False)
plt.ylabel('alt (m)')
plt.title('IWG1 alt (blue) and GPS alt (orange): '+flt_date)
plt.subplot(2,1,2,sharex=ax1)
plt.plot_date(gps_UTC_interp,alt_diff)
plt.xlabel('UTC')
plt.ylabel('IWG1 alt minus GPS alt (m)')
plt.title('Difference between IWG1 alt and GPS alt: '+flt_date)
plt.show()

ax1 = plt.subplot(2,1,1)
plt.plot_date(nav_data_all['UTC'],nav_data_all['lat'],marker='o')
plt.plot_date(gps_UTC_interp,gps_data_all['lat'],marker='x')
ax1.xaxis.set_visible(False)
plt.ylabel('lat (deg)')
plt.title('IWG1 lat (blue) and GPS lat (orange): '+flt_date)
plt.subplot(2,1,2,sharex=ax1)
plt.plot_date(gps_UTC_interp,lat_diff)
plt.xlabel('UTC')
plt.ylabel('IWG1 lat minus GPS lat (deg)')
plt.title('Difference between IWG1 lat and GPS lat: '+flt_date)
plt.show()

ax1 = plt.subplot(2,1,1)
plt.plot_date(nav_data_all['UTC'],nav_data_all['lon'],marker='o')
plt.plot_date(gps_UTC_interp,gps_data_all['lon'],marker='x')
ax1.xaxis.set_visible(False)
plt.ylabel('lon (deg)')
plt.title('IWG1 lon (blue) and GPS lon (orange): '+flt_date)
plt.subplot(2,1,2,sharex=ax1)
plt.plot_date(gps_UTC_interp,lon_diff)
plt.xlabel('UTC')
plt.ylabel('IWG1 lon minus GPS lon (deg)')
plt.title('Difference between IWG1 lon and GPS lon: '+flt_date)
plt.show()



pdb.set_trace()
