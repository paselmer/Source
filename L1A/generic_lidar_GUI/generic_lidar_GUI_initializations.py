import os
from subprocess import check_output
import numpy as np
import datetime as DT

# NOTES: 
#
# Directories are at bottom because structure different
# for Unix vs. Windows.
#
# Some variables below are only required for certain instruments. As of 
# 10/8/19, I've noted a few that aren't required for CPL. As this GUI is
# used for more and more instruements, there will inevitably be ample 
# variation in what variables are required. Please use this initializations
# file intelligently and thoughtfully. 
#
# ONLY ADD variables. NEVER subtract. All new variables will be dealt with
# in "translate_to_standard_lidar_structure.py."


# ******** WHAT IS THE NAME OF THE INSTRUMENT? USE EXACT CODE  *********
# CPL, CAMAL, Roscoe, etc.
instrument_name = 'Roscoe' 

# ************* SET THE NAME OF THE PROJECT & FLIGHT DATE **************
proj_name = 'Uplook_19'
flt_date = '04oct19'     # in "CPL" form, not "CAMAL" form
sortie = '20191004'   # Needn't actually be a sortie. Just a folder name.
search_str = '*datadown*.data'     # How code will identify file names of files

# *********** DEFINE CONSTANTS, DATA PARAMS, CONFIG OPTIONS  ***********
# horizontal averaging (# of raw profiles)
nhori = 1
# Some instrument's read routines require this. CAMAL,Roscoe don't.
file_len_recs = 9000
# MCS data files (CAMAL, Roscoe, ACATS) should contain 5 minutes of data
file_len_secs = 5 * 60
# Only required for CPL, for which it's not reported in the raw data.
nbins = 900
# Number of shots accumulated per record - again only required for CPL.
nshots = 500
# CLS data also don't include bin size. So define here (m).
vrZ = 29.98
# The laser rep rate in Hertz
rep_rate = 5000.0
# Horizontal distance between profiles (m)
dx = 200.0
# Start and end bin of solar background region
bg_st_bin = 870
bg_ed_bin = 900
# This flag tells the code which equations to use to convert energy
e_flg = 6
# The number of wavelengths
nwl = 2
# The bin width for the histogram used to determine off-nadir angles
ONA_bw = 3.0
# The max off-nadir angle, for use with histogram code
ONA_uplim = 17.0
# The min off-nadir angle, for use with histogram code
ONA_lowlim = -17.0
# Which channels #'s are which wavelengths? (0=355,1=532,2=1064)
wl_map = [1, 1, 0, 0]
# Create a list where index # of WL corresponds to string name of WL
wl_str = ['355','1064']
# default scale of color bar
CBar_max = 70.0
CBar_max_NRB = 5e13
CBar_min = 0.0
# The order of magnitude of the alt scale (help make round ticks)
scale_alt_OofM = 1e3 #order of mag.
# Default vertical range of bins to be plotted in curtain
minbin=900
maxbin=0
# curtain plot width and height (inches)
figW = 9.8
figL = 7
# profile plot width and height (inches)
pp_figW = 6
pp_figL = 6.5
# "Up" or "Down?" Which direction is lidar pointed?
pointing_dir = "Up"
# default axes limits for profile plot [xmin,xmax]/[ymin,ymax]
pp_xax_bounds_raw_counts = [0,500]
pp_xax_bounds_bgsub_counts = [-15,145]
pp_xax_bounds_NRB = [-7e12,5e13]
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
# The format specification string for dealing with the angles dropbox list
angle_list_format_str = '{:9.5f}'
# The attention bar
attention_bar = '\n******************************\n'

if (os.name != 'nt'): # IF UNIX-BASED MACHINE, DEFINE DIRECTORIES HERE

    # ******** DEFINE ALL DIRECTORIES ********
    # Set the directory of the raw data
    raw_dir = '/cpl3/'+proj_name+'/Raw_data/'+sortie+'/'
    # Directory to put output
    out_dir = '/cpl3/'+proj_name+'/analysis/'   

else:                 # IF WINDOWS-BASED MACHINE, DEFINE DIRECTORIES HERE

    # ******** DEFINE ALL DIRECTORIES ********
    # Set the directory of the raw data
    raw_dir = 'C:\\Users\\pselmer\\Documents\\Roscoe\\'+proj_name+'\\'+sortie+'\\'
    # Directory to put output
    out_dir = 'C:\\Users\\pselmer\\Documents\\Roscoe\\'+proj_name+'\\analysis\\'
