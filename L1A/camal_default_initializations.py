from scipy import constants as const
import os
from subprocess import check_output
import numpy as np
import datetime as DT

# NOTE: Directories are at bottom because structure different
#       for Unix vs. Windows

# ******** SET THE NAME OF THE PROJECT & FLIGHT DATE ********
proj_name = 'CAMAL_17'
flt_date = '20171205'

# ******** DEFINE CONSTANTS & DATA PARAMS ********
# Speed of light in meters per second
c = const.c
pi =const.pi
# The number of leap seconds since 1980-01-06 00:00:00
leapsecs = 16
# Each MCS data file should contain 5 minutes of data
file_len_secs = 5 * 60
# The laser rep rate in Hertz
rep_rate = 5000.0
# Start and end bin of solar background region
bg_st_bin = 850
bg_ed_bin = 950
# Stand and end bin of solar background region in fixed frame
ff_bg_st_bin = 770
ff_bg_ed_bin = 820
# The bin resolution in the fixed frame (m)
vrZ_ff = 30.0
# List containing top and bottom altitudes (m) of the fixed frame
ff_bot_alt,ff_top_alt = [-15e3,22.005e3]
# The number of wavelengths
nwl = 3
# Which channels #'s are which wavelengths? (0=355,1=532,2=1064)
wl_map = [2, 2, 1, 1, 0]
# Create a list where index # of WL corresponds to string name of WL
wl_str = ['355','532','1064']
# The bin width for the histogram used to determine scan angles
scan_pos_bw = 1.75
# The max scan angle, for use with histogram code
scan_pos_uplim = 17.0
# The min scan angle, for use with histogram code
scan_pos_lowlim = -17.0
# Scan angle offset, degrees. This number gets added to scan_pos.
angle_offset = 2.1
# Minimum GPS week value. All < than are filtered out.
minweek = 1000
# Housekeeping record size in bytes
hk_rec_size = 269
# GPS data rate in GPS files (records per second). Please make a float!
gps_hz = 2.0
# An estimate of the number of records in 1 gps file
est_gps_recs_1file = 5*60*10 + 30
# A maximum # of IWG1 records. Data are 1 Hz and flight can be 8 hours.
est_IWG1_recs = 30000
# An estimate of the number of nav records in 1 file
est_nav_recs_1file = 5*60*10 + 30
# MCS data records per second. Please make it a float!
MCS_hz = 10.0
# IWG1 records per second. Please make it a float!
IWG1_hz = 1.0
# nav (refering to file type) records per second. (Please make it a float!) 
nav_hz = 1.0
# Set the Polarization Gain ratios for the wavelenghts [532nm, 1064nm]
PGain = [0.00,0.00]
# Set this to the maximum possible count rate. Don't set > recs in DTT file!
max_counts = 16000
# Dead time tables [list]. List in channel order OR ELSE!!!
DTT_files = ['dttable_camal_chan1_27238-022318.xdr',
    'dttable_camal_chan2_27243-022318.xdr', 'dttable_camal_chan3_27239-022318.xdr',
    'dttable_camal_chan4_27242-022318.xdr', 'dttable_355_9999-022410.xdr']
# The overlap file to use
overlap_file = 'OLOnes_CAMAL.xdr' #'olaptable_cpl-ccviceax_comb_iceland12.xdr'      
# The number of seconds needed to convert from the instrument's Unix time
# to UTC. Typically either 5 hours (cold season) or 4 hours (warm season).
secs_btwn_instr_UnixT_and_UTC = 18000
# Offset that might be required if Unix Time between nav files and MSC
# files doesn't provide an exact match to link to other clocks (seconds).
nudge = 1.0
# Set this to 'quick' to just grab the hard-coded time offsets below
offset_choice = 'no'
def_time_offset_UnixT = DT.timedelta(seconds=1645.798)
def_time_offset_IWG1 = DT.timedelta(seconds=33.0)
# Roll and pitch offsets for GPS (degrees). Subtract these from read-in vals.
gps_roll_offset = 0.033
gps_pitch_offset = 0.088


# ******** CONFIGURATION PARAMETERS ********
# horizontal averaging (# of raw profiles)
nhori = 1
# default scale of color bar
CBar_max = 50.0
CBar_max_NRB = 5e13
CBar_min = 0.0
# The order of magnitude of the alt scale (help make round ticks)
scale_alt_OofM = 1e3 #order of mag.
# Default vertical range of bins to be plotted in curtain
minbin=1000
maxbin=0
# curtain plot width and height (inches)
figW = 18
figL = 10
# profile plot width and height (inches)
pp_figW = 6.5
pp_figL = 7
# "Up" or "Down?" Which direction is lidar pointed?
pointing_dir = "Down"
# default axes limits for profile plot [xmin,xmax]/[ymin,ymax]
pp_xax_bounds_raw_counts = [0,100]
pp_xax_bounds_bgsub_counts = [-15,45]
pp_xax_bounds_NRB = [-7e12,5e13]
pp_yax_bounds_bins = [1000,0]
pp_yax_bounds_alt = [-10e3,20e3]
# Y-axis bounds of the energy monitor plot
EMp_yax_bounds = [0,150]
# Channel # entries by user cannot be outside this range
min_chan = 1
max_chan = 5
# Cap on the number of profiles that can be used to generate curtain.
hori_cap = 20000
# The padding around the curtain plot in inches
CPpad = 0.1
# The format specification string for dealing with the angles dropbox list
angle_list_format_str = '{:9.5f}'
# Set default housekeeping axis values
hk_y_max = 100
hk_y_min = -100
# For actual data processing (not the GUI), nav data source
Nav_source = 'nav' #'nav' 'gps' or 'iwg1'
# IWG1 data file
IWG1_file = "IWG1.08Dec2017-0031.txt"
# Don't process any data when below this alt (m). Doesn't apply to GUI.
alt_cutoff = 10000

if (os.name != 'nt'): # IF UNIX-BASED MACHINE, DEFINE DIRECTORIES HERE

    # ******** DEFINE ALL DIRECTORIES ********
    # Set the directory of the raw data
    raw_dir = '/cpl3/CAMAL/Data/'+proj_name+'/'+flt_date+'/raw/'
    # Set the directory of the L0 data
    L0_dir = '/cpl3/CAMAL/Data/'+proj_name+'/'+flt_date+'/L0/'
    # Set the directory of the L1 data
    L1_dir = '/cpl3/CAMAL/Data/'+proj_name+'/'+flt_date+'/L1/'
    # Set the directory of the L2 data
    L2_dir = '/cpl3/CAMAL/Data/'+proj_name+'/'+flt_date+'/L2/'
    # Set the directory that contains configuration files
    config_dir = '/cpl3/CAMAL/Config/'
    # The source directory
    source_dir = '/cpl3/CAMAL/Source/L1A/'
    # Directory to put output
    out_dir = '/cpl3/CAMAL/Analysis/'+proj_name+'/'+flt_date+'/'
    # Directory and name of library containing C codes
    clib_path = source_dir + 'C_lib_unix/'
    clib = 'CAMAL_C_lib_v0.so'
    dtt_dir = config_dir + 'dttables/'
    # Directory containing overlap files
    olap_dir = '/home/selmer/compute_OL/' 
    # Directory and name of the DEM
    DEM_dir = config_dir + 'DEM/'
    DEM_name = 'JPL-CloudSat.dem'
    # Directory and name of C++ library with functions to read in DEM
    DEM_lib_path = DEM_dir + 'JPL_DEM_CPP_SO/'
    DEM_Cpp_name = 'JPL_DEM_CPP_FUNCTIONS.so'
    

else:                 # IF WINDOWS-BASED MACHINE, DEFINE DIRECTORIES HERE

	    # ******** DEFINE ALL DIRECTORIES ********
    # Set the directory of the raw data
    #raw_dir = 'F:\\CAMAL\\from_datakey\\'+flt_date+'\\'
    raw_dir = 'C:\\Users\\pselmer\\Documents\\CAMAL\\camal\\data\\'+proj_name+'\\'+flt_date+'\\raw\\'
    # Set the directory of the L0 data
    L0_dir = 'C:\\Users\\pselmer\\Documents\\CAMAL\\camal\\data\\'+proj_name+'\\'+flt_date+'\\L0\\'
    # Set the directory of the L1 data
    L1_dir = 'C:\\Users\\pselmer\\Documents\\CAMAL\\camal\\data\\'+proj_name+'\\'+flt_date+'\\L1\\'
    # Set the directory of the L2 data
    L2_dir = 'C:\\Users\\pselmer\\Documents\\CAMAL\\camal\\data\\'+proj_name+'\\'+flt_date+'\\L2\\'
    # Set the directory that contains configuration files
    config_dir = 'C:\\Users\\pselmer\\Documents\\CAMAL\\camal\\config_source\\'
    # The source directory
    source_dir = 'C:\\Users\\pselmer\\Documents\\CAMAL\\camal\\source\\'
    # Directory to put output
    out_dir = 'C:\\Users\\pselmer\\Documents\\CAMAL\\camal\\analysis\\'+proj_name+'\\'+flt_date+'\\'
    # Directory and name of library containing C codes
    clib_path = 'C:\\Users\\pselmer\\Documents\\CAMAL\\camal\\source\\C_lib_Win64\\CAMAL_C_lib_Win64\\x64\\Release\\'
    clib = 'CAMAL_C_lib_Win64.dll'
    # Directory of dead time tables
    dtt_dir = config_dir + 'dttables\\'
    # Directory containing overlap files
    olap_dir = config_dir
    # Directory and name of the DEM
    DEM_dir = config_dir + '\\DEM\\'
    DEM_name = 'JPL-CloudSat.dem'
    # Directory and name of C++ library with functions to read in DEM
    DEM_lib_path = 'C:\\Users\\pselmer\\Documents\\CAMAL\\camal\\config_source\\DEM\\JPL_DEM_CPP_DLL\\x64\\Release\\'
    DEM_Cpp_name = 'JPL_DEM_CPP_DLL.dll'

