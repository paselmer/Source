from scipy import constants as const
import os

# NOTE: Directories are at bottom because structure different
#       for Unix vs. Windows

# ******** SET THE NAME OF THE PROJECT & FLIGHT DATE ********
proj_name = 'Project_Name'
flt_date = '20170926'

# ******** DEFINE CONSTANTS & DATA PARAMS ********
# Speed of light in meters per second
c = const.c
pi =const.pi
# The number of leap seconds since 1980-01-06 00:00:00
leapsecs = 16
# Each file should contain 5 minutes of data
file_len_secs = 5 * 60
# The laser rep rate in Hertz
rep_rate = 5000.0
# Start and end bin of solar background region
bg_st_bin = 3
bg_ed_bin = 23
# The number of wavelengths
nwl = 3
# Which channels #'s are which wavelengths? (0=355,1=532,2=1064)
wl_map = [2, 2, 1, 1, 0]
# Create a list where index # of WL corresponds to string name of WL
wl_str = ['355','532','1064']
# The bin width for the histogram used to determine scan angles
scan_pos_bw = 0.1
# The max scan angle, for use with histogram code
scan_pos_uplim = 10.0
# The min scan angle, for use with histogram code
scan_pos_lowlim = -10
# Scan angle offset, degrees. This number gets added to scan_pos.
angle_offset = 2.1
# Housekeeping record size in bytes
hk_rec_size = 269

# ******** CONFIGURATION PARAMETERS ********
# horizontal averaging (# of raw profiles)
nhori = 1
# default scale of color bar
CBar_max = 75.0
CBar_max_NRB = 5e6
CBar_min = 0.0
# the value of the min & max y tick marks (in meters)
min_scale_alt = -5e3
max_scale_alt = 20e3
scale_alt_OofM = 1e3 #order of mag.
# curtain plot width and height (inches)
figW = 11
figL = 7
# profile plot width and height (inches)
pp_figW = 4.5
pp_figL = 7
# "Up" or "Down?" Which direction is lidar pointed?
pointing_dir = "Up"
# default axes limits for profile plot [xmin,xmax]/[ymin,ymax]
pp_xax_bounds_raw_counts = [0,100]
pp_xax_bounds_bgsub_counts = [-15,100]
pp_xax_bounds_NRB = [-5e5,5e6]
pp_yax_bounds_bins = [0,1000]
pp_yax_bounds_alt = [0,25e3]
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

if (os.name != 'nt'): # IF UNIX-BASED MACHINE, DEFINE DIRECTORIES HERE

    # ******** DEFINE ALL DIRECTORIES ********
    # Set the directory of the raw data
    raw_dir = '/home/litespar/camal/data/'+proj_name+'/'+flt_date+'/raw/'
    # Set the directory of the L0 data
    L0_dir = '/home/litespar/camal/data/'+proj_name+'/'+flt_date+'/L0/'
    # Set the directory of the L1 data
    L1_dir = '/home/litespar/camal/data/'+proj_name+'/'+flt_date+'/L1/'
    # Set the directory of the L2 data
    L2_dir = '/home/litespar/camal/data/'+proj_name+'/'+flt_date+'/L2/'
    # Set the directory that contains configuration files
    config_dir = '/home/litespar/camal/config_source/'
    # The source directory
    source_dir = '/home/litespar/camal/source/'
    # Directory to put output
    out_dir = '/home/litespar/camal/analysis/'+proj_name+'/'+flt_date+'/'

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

