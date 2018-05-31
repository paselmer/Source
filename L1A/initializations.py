from scipy import constants as const
import os
from subprocess import check_output
import numpy as np
import datetime as DT

# NOTE: Directories are at bottom because structure different
#       for Unix vs. Windows

# ******** SET THE NAME OF THE PROJECT & FLIGHT DATE ********
proj_name = 'Wallops_12'
flt_date = '20sep12' # in "CPL" form, not "CAMAL" form
sortie = '12-922'

# ******** SET THE TIME RANGE BOUNDARIES FOR DATA PROCESSING ********
process_start = DT.datetime(2000,9,1,0,0,0) #yr,mon,dy,hr,min,sec
process_end   = DT.datetime(2020,10,13,0,0,0)

# ******** DEFINE CONSTANTS & DATA PARAMS ********
# Speed of light in meters per second
c = const.c
pi =const.pi
# The number of leap seconds since 1980-01-06 00:00:00
leapsecs = 15
# Each CLS data file should contain 9000 records
file_len_recs = 9000
# Amazingly, the number of vertical bins is not reported in the
# CLS data don't not include this. So, define it here.
nbins = 833
# Number of shots accumulated per record
nshots = 500
# CLS data also don't include bin size. So define here (m).
vrZ = 29.98
# The laser rep rate in Hertz
rep_rate = 5000.0
# Start and end bin of solar background region (GUI only)
bg_st_bin = 730
bg_ed_bin = 800
# Start and end altitudes of solar background region (meters)
bg_st_alt = -500.0
bg_ed_alt = -1500.0
# The bin resolution in the fixed frame (m)
vrZ_ff = 30.0
# List containing top and bottom altitudes (m) of the fixed frame
ff_bot_alt,ff_top_alt = [-10e3,22.5e3]
# This flag tells the code which equations to use to convert energy
e_flg = 6
# The number of wavelengths
nwl = 3
# Which channels #'s are which wavelengths? (0=355,1=532,2=1064)
wl_map = [0, 1, 2, 2]
# Create a list where index # of WL corresponds to string name of WL
wl_str = ['355','532','1064']
# Minimum GPS week value. All < than are filtered out.
minweek = 1000
# GPS data rate in GPS files (records per second). Please make a float!
gps_hz = 2.0
# An estimate of the number of records in 1 gps file
est_gps_recs_1file = 5*60*10 + 30
# A maximum # of IWG1 records. Data are 1 Hz and flight can be 8 hours.
est_IWG1_recs = 30000
# An estimate of the number of nav records in 1 file
est_nav_recs_1file = 5*60*10 + 30
# CLS data records per second. Please make it a float!
CLS_hz = float(rep_rate) / float(nshots)
# The resolution of CPL's instrument clock in seconds
inst_clk_rez = 1.0
# Computed and init. file difference can't exceed this tolerance.
CLS_hz_tol = 0.05
# IWG1 records per second. Please make it a float!
IWG1_hz = 1.0
# nav (refering to file type) records per second. (Please make it a float!) 
nav_hz = 1.0
# Set the Polarization Gain ratios for the wavelenghts [532nm, 1064nm]
PGain = [0.00,0.766]
# Set this to the maximum possible count rate. Don't set > recs in DTT file!
max_counts = 16000
# Dead time tables [list]. List in channel order OR ELSE!!!
DTT_files = ['dttable_355_0135-072102.xdr','dttable_532_11296-021009_30m.xdr',
    'dttable_1064par_8351-021009_30m.xdr','dttable_1064per_8346-021009_30m.xdr']
# The overlap file to use
overlap_file = 'olaptable_cpl-ccviceax_comb_iceland12.xdr'
# The number of seconds needed to convert from the instrument's Unix time
# to UTC. Typically either 5 hours (cold season) or 4 hours (warm season).
secs_btwn_instr_UnixT_and_UTC = 14400
# Roll and pitch offsets for GPS (degrees). Subtract these from read-in vals.
gps_roll_offset = 0.0
gps_pitch_offset = 0.0
# Subtract this number of seconds from CLS data to manually fudge a better
# match to the Nav data.
nudge = -1.0


# ******** CONFIGURATION PARAMETERS ********
# Average the data to this time. In seconds. Set to -99.9 for no averaging.
secs2avg = 1.0
# Minimum number of raw profiles than can be used in an average profile
min_avg_profs = 4
# horizontal averaging (# of raw profiles)
nhori = 1
# default scale of color bar
CBar_max = 50.0
CBar_max_NRB = 5e13
CBar_min = 0.0
# The order of magnitude of the alt scale (help make round ticks)
scale_alt_OofM = 1e3 #order of mag.
# Default vertical range of bins to be plotted in curtain
minbin=833
maxbin=0
# curtain plot width and height (inches)
figW = 16
figL = 9
# profile plot width and height (inches)
pp_figW = 6.5
pp_figL = 7
# "Up" or "Down?" Which direction is lidar pointed?
pointing_dir = "Down"
# default axes limits for profile plot [xmin,xmax]/[ymin,ymax]
pp_xax_bounds_raw_counts = [0,100]
pp_xax_bounds_bgsub_counts = [-15,45]
pp_xax_bounds_NRB = [-7e12,5e10]
pp_yax_bounds_bins = [1000,0]
pp_yax_bounds_alt = [-10e3,20e3]
# Y-axis bounds of the energy monitor plot
EMp_yax_bounds = [0,150]
# Channel # entries by user cannot be outside this range
min_chan = 1
max_chan = 4
# Cap on the number of profiles that can be used to generate curtain.
hori_cap = 20000
# The padding around the curtain plot in inches
CPpad = 0.1
# For actual data processing (not the GUI), nav data source
Nav_source = 'cls' #'nav' 'gps' 'iwg1' or 'cls'
# IWG1 data file
IWG1_file = "IWG1.20Jan2013-2149"
# Don't process any data when below this alt (m). Doesn't apply to GUI.
alt_cutoff = 500
# Don't process any profiles where off-nadir angle exceeds this many radians.
ONA_cutoff = 30.0 * (pi/180.0)
# Invalid/bad Nav values get overwritten with this value - datetime.datetime object.
bad_cls_nav_time_value = DT.datetime(1970,1,1,0,0,0)
# The attention bar
attention_bar = '\n******************************\n'
# The maximum number of hours a dataset is expected to be
max_flt_hours = 30

if (os.name != 'nt'): # IF UNIX-BASED MACHINE, DEFINE DIRECTORIES HERE

    # ******** DEFINE ALL DIRECTORIES ********
    # Set the directory of the raw data
    raw_dir = '/cpl3/'+proj_name+'/Raw_data/'+sortie+'/'
    # Set the directory of the L0 data
    L0_dir = '/cpl3/'+proj_name+'/L0/'
    # Set the directory of the L1 data
    L1_dir = '/cpl3/'+proj_name+'/L1/'
    # Set the directory of the L2 data
    L2_dir = '/cpl3/'+proj_name+'/L2/'
    # Set the directory that contains configuration files
    config_dir = '/cpl/dhlavka/Cpl/Config/'
    # The source directory
    source_dir = '/cpl3/CAMAL/Source/L1A/'
    # Directory to put output
    out_dir = '/cpl3/'+proj_name+'/analysis/'
    # Directory and name of library containing C codes
    clib_path = source_dir + 'C_lib_unix/'
    clib = 'CAMAL_C_lib_v0.so'
    # Directory containing dead time files
    dtt_dir = '/cpl/dhlavka/Cpl/Config/'
    # Directory containing overlap files
    olap_dir = '/cpl/dhlavka/Cpl/Config/'
    # Directory and name of the DEM
    DEM_dir = '/cpl3/CAMAL/Config/DEM/'
    DEM_name = 'JPL-CloudSat.dem'
    # Directory and name of C++ library with functions to read in DEM
    DEM_lib_path =  '/cpl3/CAMAL/Source/DEM_reader_code/JPL_DEM_CPP_SO/'
    DEM_Cpp_name = 'JPL_DEM_CPP_FUNCTIONS.so'
    

else:                 # IF WINDOWS-BASED MACHINE, DEFINE DIRECTORIES HERE

    # ******** DEFINE ALL DIRECTORIES ********
    # Set the directory of the raw data
    raw_dir = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\ACT-A\\'+proj_name+'\\'+flt_date+'\\'
    # Set the directory of the L0 data
    L0_dir = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\ACT-A\\'+proj_name+'\\'+flt_date+'\\L0\\'
    # Set the directory of the L1 data
    L1_dir = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\ACT-A\\'+proj_name+'\\'+flt_date+'\\L1\\'
    # Set the directory of the L2 data
    L2_dir = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\ACT-A\\'+proj_name+'\\'+flt_date+'\\L2\\'
    # Set the directory that contains configuration files
    config_dir = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\Config\\'
    # The source directory
    source_dir = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\Source\\L1A\\'
    # Directory to put output
    out_dir = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\ACT-A\\'+proj_name+'\\'+flt_date+'\\analysis\\'
    # Directory and name of library containing C codes
    clib_path = source_dir + 'C_lib_Win64\\CAMAL_C_lib_Win64\\x64\\Release\\'
    clib = 'CAMAL_C_lib_Win64.dll'
    # Directory containing dead time files    
    dtt_dir = config_dir
    # Directory containing overlap files
    olap_dir = config_dir    
    # Directory and name of the DEM
    DEM_dir = 'C:\\Users\\pselmer\\Documents\\CAMAL\\camal\\config_source\\DEM\\'
    DEM_name = 'JPL-CloudSat.dem'
    # Directory and name of C++ library with functions to read in DEM
    DEM_lib_path = DEM_dir + 'JPL_DEM_CPP_DLL\\x64\\Release\\'
    DEM_Cpp_name = 'JPL_DEM_CPP_DLL.dll'

