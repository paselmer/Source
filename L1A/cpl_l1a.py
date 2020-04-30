# CPL L1A processing code.
# NOTE: The most current version of the L1A code must always be named
#       "cpl_l1a.py." All versions of the code must be saved as
#       "cpl_l1a_vX.py" but not actually be source file that's exectued.

# DESCRIPTION:
#
# [6/21/18] Version 0 (v0)
#
# This code processes raw CPL ".cls" files into NRB HDF5 files. It also
# produces a text file of lon,lat,&time for the purpose of sampling
# model data.
#
# This code can be thought of as having several main parts:
#
# 1) Nav ingestion, cleanup, and interpolation 
# 2) Nav-to-CLS time syncronization {per file basis}
# 3) Initialization of arrays & HDF5 parameters {per file basis}
# 4) Nav-to-CLS record matching {per raw record basis}
# 5) Compuation/population of various parameters,
#    preparation for NRB computation {per raw record basis}
# 6) PGain application & actual NRB computation {per file basis}
# 7) Deletion of bad records {per file basis}
# 8) Averaging of data & final HDF5 output {per file basis}
#
# Explanation of each part:
#
# 1) Nav ingestion, cleanup, and interpolation
# Oh boy, this part is a doozy. There are a couple functions defined at
# the top. One computes a time offset between IWG1 Nav data and the
# instrument. The other function reads in one CLS file and provides the
# arrays needed for time syncronization. The later function isn't used
# until the file loop (part 2).
# The Nav_source is chosen via the initialization file and the
# corresponding block of code is chosen. In each block, the same
# general things happen - the data are read in, then interpolated to the
# CLS data rate.
# The 'cls' block is more complex due to the inconsistent sampling of
# the IWG1 feed by the CPL data system (a single IWG1 record does not
# always repeat for 10 +/-3 records). In the 'cls' block the data
# rate is computed and compared to the data rate that the user input
# in the initialization file. If the difference is too much, the code
# just prints out a warning message then proceeds to use the user input
# data rate. Since the Nav time has a 1 Hz  rez and the photon
# data have a ~10 Hz rez, an array of interpolated (to ~10 Hz) times 
# replaces the original CLS Nav times. Additional data cleanup
# is done within the 'cls' block due to things like "status records"
# being at the beginnings/ends of the flight.
#
# 2) Nav-to-CLS time syncronization {per file basis}
# At the top of the file loop, data from the first CLS file are read.
# The interpolated instrument time and the interpolated Nav time
# (from whatever Nav_source was chosen) are sent to a C function
# in an external library. This C function uses those arrays to return
# an array of subscripts which map the Nav_source records to the CLS
# records. There should be a subscript for each record of the current
# CLS file. Doing this means that the matching implied by part 4 is 
# essentially already done.
# 
# 3) Initialization of arrays & HDF5 parameters {per file basis}
# During the first iteration of the file loop, various arrays and
# HDF5 dataset objects are initialized. The two output files, HDF5 and
# text, are opened. A small number of parameters are reinitialized each
# iteration of the file loop.
#
# 4) Nav-to-CLS record matching {per raw record basis}
# The work of finding matches is already done by this point, as
# mentioned in part 2's explanation. The Nav record that matches the
# current CLS profile is found by applying the subscript computed
# in part 2. A different 'if' block is executed depending on the Nav_source.
# Each if block does the same things - stores the matched Nav into
# the Nav_save array and writes the lon,lat,UTC out to the GEOS text
# file. Nav_save is a standardized structure. It gets populated here
# no matter the Nav_source.
#
# 5) Compuation/population of various parameters,
#    preparation for NRB computation {per raw record basis}
# Now that you have the matching Nav data, a whole host of 
# computations take place...dead time correction, off-nadir angle, overlap
# geolocation, bin altitude. Background subtraction and range correction
# are applied here because they're hard to do after the counts are
# put into a fixed altitude frame. Immediately  after this, the counts
# are put into a fixed altitude frame.
#
# 6) PGain application & actual NRB computation {per file basis}
# At this point, you've gone through every raw rez profile in a single
# CLS file and are back in the file loop. PGain is applied here, and the
# final NRB is computed. These things are done here to increase speed by
# vectorizing the computation (ie - implicit loops across recs & bins
# instead of explicit loop); you'll notice, there's no bins loop in
# this code.
#
# 7) Deletion of bad records {per file basis}
# A boolean (bool) array, good_recs_bool, is created anew for each
# iteration of the file loop. It is initialized to all True (1's).
# At several points in the code, the current record can be flagged as
# bad, False (0), under certain conditions (altitude too low, off-nadir
# angle too high, etc.). In this relatively small block of code the
# bad records are deleted by logic of keeping only the remaining good 
# records. Arrays are adjusted accordingly.
#
# 8) Averaging of data & final HDF5 output {per file basis}
# Data are averaged here, file by file...
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
# Within the averaging section, the averaged data are dumped into the HDF5 file.
# The code concludes outside the loops by simply recording the total number of
# CLS records and closing the output files.


# UPDATE LOG:
#
# [4/2/18] Beta version started.
# Using camal_l1a_v0.py as a starting template.
#
# [4/6/18] NRB success
# NRB first successfully calculated. Used Podex_13 data.
#
# [5/7/18] Variable time-averaging of data successful
#
# [5/31/18] "Banding" issue fixed
# Fixed issue where banding would appear in NRB profiles
# (very obvious in images) due to integer rounding of
# counts array. Basically fixed by converting counts
# array from an int array to float array.
#
# [6/12/18] Code modified to accomodate wonky overlap
# Noticeable vertical "banding" was present in images of NRB.
# Determined that it was caused by CPL's strange overlaps.
# Fixed the issue by applying overlap after background sub.
# GEOS-lon-lat-UTC output file is now created by this code.
# Prior to this date it had not been.
#
# [6/18/18] Nav_source of "nav" now tested for first time
# Tested using "nav" text files as the source of nav data for processing
# for ACTA-18 27apr18 flight. Worked out a few kinks - Windows being
# dumb with its 'dir' command, slight difference between CAMAL and CPL
# style nav text files. Visually compared images of NRB processed using
# CLS and NRB processed using Nav. They appeared identical. 20-degree 
# rolls throughout.
#
# [6/20/18] All-zeros first record fixed
# Rebecca discovered that the first record of every file I gave her
# had all zeros for NRB in the first record. I discovered that this
# was due to a flaw in logic in how the averaging section of this 
# code handled the HDF5 dataset resizing of the first file. I added
# a couple more lines of logic and the problem seems to be fixed.
#
# [6/21/18] v0 description written
#
# [6/22/18] More "nav" Nav_source capability
# Those nav text files produced by the various instruments (CAMAL,CPLs)
# have different formats depending on the instrument & aircraft.
# Code can now read & process with nav text files produced by ER2-CPL,
# UAV-CPL, and CAMAL - all of which have different formats. Processing
# with all 3 types mentioned here has been tested.
# ...
# Improved some logic in the "nav" Nav_source block within the profile
# loop. Code will now pick the max between GPS_alt and GPS_alt_msl.
# They can be various combinations of invalid. Right now, "nav" should
# be a second or third processing choice, in my opinion. Use IWG1 files
# for best results.
#
# [7/9/18] Overlap correction fixed
# This is embarrassing. I was multiplying by the overlap factor in those
# XDR files instead of dividing by it. This was causing the weird slope
# in the calibration region (especially in 532 nm). This has been fixed
# so that now division is performed.
#  The following is overlap application in AM.pro...
#  "
#  RSqOverLap.RsqOL(0:(NumSamples-1),i)=$
#             RSqOverLap.RsqOL(0:(NumSamples-1),i)/OL(*,IndxOL)
#  "
#
# [7/30/18] Started addition of saturate_ht
# This variable is required by current CPL L2 process.
#
# [7/31/18] Finished adding saturate_ht
# Each iteration of the file loop, a flag is set (satur_flag) based on a 
# broad check of whether any bins might be saturated.
# Within the profile loop saturate_ht is populated with heights of bins
# that exceed the saturation_values. These bins are referenced to the actual
# bin height and not the fixed-frame altitude. This is fine by me.
#
# [9/24/18] Fixed time range boundaries error
# There was a pretty obvious error in the block of code that masks out
# records not within the user specified time range. The error specifically
# the np.extract() call. Also, np.extract() needed to be applied to CLS_UnixT_float64_orig.
# I must have never actually tested this feature of the L1A code. It works
# now.
#
# [10/1/18] saturate_ht issue fixed
# This issue has yet to have caused a crash; however, when looking at
# camal_l1a_v1 and cpl_l1a_v0 side-by-side, I noticed that I neglected
# filter out the bad records for saturate_ht within the CPL code.
# I'm guessing this never caused a crash because all CPL flights
# used to test since saturate_ht was added probably had nice, continuous
# data.
#
# [11/13/18] Nav_save averaging error fixed
# Nav_save_sum was averaging with itself instead of being averaged with
# Nav_save from the current file. This is now fixed.
#
# [11/14/18] Another averaging error discovered and fixed
# rr was not being incremented in the trans_bin block 
#
# [1/22/19] background added as output variable
# The background signal per channel was added as an output variable. This was done
# because it is needed later to compute attenuated total backscatter uncertainty. 
#
# [1/31/19] fixed a bug that affects datasets with very large time gaps
# I fixed a bug with essentially a band-aid. In cases where datasets have very large
# timegaps between contiguous data, such as the Guam 2014 ATTREX deployment, an issue
# arose in the averaging section when the data gap occured between 2 consecutive CLS files.
# The codes essentially thought the time bins lines up perfectly because no times in the 
# current CLS file belonged in the last time bin of the previous file. Therefore it set
# 'ei' (end index) to the full shape (as opposed to -1) which resulted in rr going
# out of bounds by going thru its loop and then into the "not last_file" block right afterwords.
# The bandaid is this: if ((not last_file) and (ei != ui.shape[0]))
# Both can't be true, the way I have things coded, so this should prevent the code from crashing
# in this situation.
#
# [5/3/19] v1 first started
# Created a new function (inside this file), "patch_outliers," that is used
# to patch over wonky energy monitor (EM) values with reasonable values.
# The impetus for creating this was the 28Jan13 PODEX flight.
#
# [5/16/19] min_avg_prof enforcement block added in averaging section.
#
# [7/1/19] trans_total initialized to zero. Error corrected WRT last update which added
# the min_avg_prof enforcement block.
#
# [10/18/19] adjustments for initializations-free read_routines
#
# [4/22/20] ***** Version 2 started, cpl_l1a_v2 *****
# This version was started to incorporate attitude correction parameters
# from Andrew. The option to incorporate multiple overlaps, started in
# code cpl_l1a_v1a.py, was also included in this new version.
# The initializations file has been updated to include attitude correction
# parameters, and parameters in file have been regrouped.
#
# [4/30/20] Issue with setting 'ix' using computed data rate. Improved.
# During processing of PELIcoe_19 21Oct19 data, Andrew Kupchock noticed
# that nudge improved the Nav sync for one part of the flight, which making
# it worse for another portion of the flight.
# After some discussion we think it was the result of creating the interpolated
# Nav oordinate array, "ix," with the improper data rate. I was comparing
# a computed data rate with the initialized one, and using the computed
# one if it was different from the computed one by an arbitrary tolerance.
# With a tolerance (CLS_hz_tol in initializations) set to 0.05 Hz, the
# aforementioned PELIcoe flight's computed data rate of 0.03833 Hz slipped
# through. A difference of 0.03833 Hz will result in a 690 record difference
# if it's really off from the true data rate by that much! This can explain
# the symptom Andrew Kupchock first noticed.
# The following changes were made for now...
# - tot_elapsed seconds now computed using Nav time instead of the CLS clock time.
# - Data rate in initializations file is now always used to create oordinates
#   of interpolated Nav array. User is now only warned if computed data
#   rate differs too much from input data rate.


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

        cls_data = read_in_cls_data(cls_files[fc].strip(),nbins,flt_date,bad_cls_nav_time_value)

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

    cls_data_1file = read_in_cls_data(single_cls_file,nbins,flt_date,bad_cls_nav_time_value)
    nr0 = cls_data_1file.shape[0]

    if Nav_source == 'cls':

        if firstlast_flag == -1: 
            cls_data_1file = cls_data_1file[firstlast_trunc[0]:]
        if firstlast_flag == 1:
            cls_data_1file = cls_data_1file[:nr0-firstlast_trunc[1]]
        good_data_len = Navi[1] - Navi[0]
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
	

# Dictionary of letters and integers
AlphaNum = {}
for i in range(65,91): AlphaNum[chr(i)] = i-65 
	
def load_all_possible_overlaps(overlap_dict,nb,nc):
    """ Function that will load all possible overlap tables into an array
        of [nOLs, nchans, nbins]. That array is returned.
    """
    
    # INPUTS:
    # overlap_dict -> From load_overlap_configuration(). Contains the 
    #    names of all overlap tables.
    # nb -> Number of bins (833 for CPL)
    # nc -> Number of channels (4 for CPL)
    
    # OUTPUT:
    # all_overlaps -> array containing the OL data from all OL table files
    #    [nOLs, nc, nb]
    
    nOLs = len(overlap_dict)
    all_overlaps = np.zeros((nOLs,nc,nb),dtype=np.float32)
    
    for i in range(0,nOLs):
        # Load the overlap table
        overlap_file = overlap_dict[i]
        overlap_vector = read_in_overlap_table(olap_dir+overlap_file)
        # Overlap is one long array, where first nbins (833) are
        # 355, next nbins are 532, and next (final) nbins are 1064.
        overlaps = overlap_vector.reshape((nwl,nb))
        # Make the overlap_vector into an array where the 1st dim. lines up sequentially
        # with the channels
        overlaps_chan_seq = np.ones((nc,nbins),dtype=overlaps.dtype)
        for chan in range(0,nc):
            overlaps_chan_seq[chan,:] = overlaps[wl_map[chan],:]
        all_overlaps[i,:,:] = overlaps_chan_seq
        
    return all_overlaps
    
def load_overlap_configuration(overlaps_config_file):
    """ Function that will look at special (not L1A initializations)
        user-defined config file to determine which overlap tables go
        to which segments (defined by time) of the data.
    """
    
    # INPUT:
    # overlaps_config_file -> specially formatted file that identifies the
    #    overlap tables and the time spans to which to apply them.
    
    # OUTPUTS:
    # overlap_dict -> dictionary that uses integer codes as keys to
    #    identify the string names of overlap table files.
    # overlap_map -> list [[], []] that maps time spans to the OL table
    #    to use, by wavelength.
    #    [ [datetime1, datetime2], [355,532,1064] ]
    
    
    with open(overlaps_config_file) as f_obj:
        line = f_obj.readline()
        Lc = 0 # line counter
        Cc = 0 # codes element counter
        Tc = 0 # times element counter
        OL_codes = []
        OL_times = []
        XXXflag = False
        flt_indx = None
        while line:
            split_line = line.split(',')
            if split_line[0].strip() == 'XXX': XXXflag = True
            if (Lc > 0) & (not XXXflag):
                OL_codes.append(split_line)
                if flt_date in OL_codes[Cc]: flt_indx = Cc
                Cc += 1
            elif (split_line[0] == 'Times for overlaps') | (Lc == 0) | (split_line[0].strip() == 'XXX'):
                pass # Skip this line
            elif (Lc > 1) & (len(OL_codes) > 0):
                OL_times.append(split_line)
                Tc += 1
            else:
                print('Overlaps configuration not formatted as expected...')
                pdb.set_trace()
            line = f_obj.readline()
            Lc += 1
        
        # Figure out the number of overlap tables the user has programmed.
        list_len = len(OL_times)
        Ic = 1 # list Item counter
        for item in OL_times:
            if item[0] == 'Overlaps key':
                nOL = list_len - Ic
            Ic += 1
        OL_keys = OL_times[-1*nOL:]  # Last 
        OL_times = OL_times[:-1*nOL-1]
        print('{0:d} overlap tables detected.'.format(nOL))
        
        overlap_dict = {}
        for item in OL_keys: overlap_dict[AlphaNum[item[0].strip()]] = item[1].strip()
        
        OL_codes = OL_codes[flt_indx]
        OL_times = OL_times[flt_indx]
        
        Ic = 0 # list Item counter
        overlap_map = []
        for codes, times in zip(OL_codes,OL_times):
            if (codes == '') | (codes == '\n'): break
            if (Ic > 0):
                split_codes = codes.split('/')
                split_times = times.split('/')
                overlap_map.append( [
                      [ DT.datetime.strptime(split_times[0],"%Y-%m-%dT%H:%M:%S"), DT.datetime.strptime(split_times[1],"%Y-%m-%dT%H:%M:%S") ], 
                      [ AlphaNum[split_codes[1]], AlphaNum[split_codes[2]], AlphaNum[split_codes[0]] ] 
                    ] )
            Ic += 1
    
    return overlap_dict, overlap_map

    
# Start main execution here <-------------------------------------------
print('Starting main L1A execution at: ',DT.datetime.now())

# File detailing how to map OLs to different portions of the flight
overlaps_config_file = 'overlaps_configuration.csv'

# Create and load file list for CLS data
CLS_file_list = 'processing_file_list.txt'
search_str = '*.cls'
create_a_file_list(CLS_file_list,search_str,raw_dir)
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
    cls_meta_data_all, nav2cls_indx, FL_trunc, usable_file_range = read_entire_cls_meta_dataset(raw_dir,file_len_recs,nbins,flt_date,bad_cls_nav_time_value)
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
    tot_elapsed = Nav_UnixT[-1] - Nav_UnixT[0]
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
    u, ui, ncounts = np.unique(cls_meta_data_all['Header']['ExactTime'],return_index=True,return_counts=True)
    if u[0] == bad_cls_nav_time_value:
        u = u[1:]
        ui = ui[1:]
        ncounts = ncounts[1:]    
    ExactTime_ncounts = np.copy(ncounts)
    ExactTime_ui = np.copy(ui)
    inst_clk_deltas = delta_datetime_vector(u)

    # Estimate the data rate from the Nav time and length of data array.
    # Warn user if computed data doesn't match initialized data rate.
    computed_CLS_hz = cls_meta_data_all.shape[0] / tot_elapsed
    if abs(computed_CLS_hz - CLS_hz) > CLS_hz_tol:
        print(attention_bar)
        print("!!!!!!! WARNING !!!!!!!")
        print("The estimated (init. file) and computed data rates differ by too much.")
        print("estimated: ",CLS_hz," computed: ",computed_CLS_hz)
        print("Make sure CLS_hz from initializations is correct")
        print("!!!!!!! WARNING !!!!!!!")
        print('Enter "c" to contine processing anyway.')
        print('Enter "q" to safely quit.')
        print('If you continue, data rate from initializations will be used.')
        print(attention_bar)
        pdb.set_trace()
    else:
        print('Init. file CLS data rate is sufficiently close to computed data rate.')
        print('CLS_hz: ',CLS_hz, ' computed_CLS_hz: ',computed_CLS_hz)
        print('To be clear, code will use CLS data rate from initializations (CLS_hz)\n')
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
        
        print('If you are seeing this message, tell Patrick Selmer to')
        print('investigate why this block of code is necessary.')
        # Maybe this relates to how NEW is populated in the C code.
        # Wish I would have take better notes about this block.
        pdb.set_trace()
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

print('Starting core processing...') # <--------------------------------
first_read = True
usable_file_indicies = range(usable_file_range[0],usable_file_range[1])
trans_bin = [0,0]
trans_total = 0
last_file = False 
OL_map_indx = 0 # count thru overlap map as required
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
    if ( (CLS_data_1file['meta']['Header']['ExactTime'].min() < process_start) 
      or (CLS_data_1file['meta']['Header']['ExactTime'].max() > process_end) ):
        time_range_mask_low = CLS_data_1file['meta']['Header']['ExactTime'] > process_start
        time_range_mask_high = CLS_data_1file['meta']['Header']['ExactTime'] < process_end
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
        bins_range=np.broadcast_to(np.arange(0,nb*vrZ,vrZ,dtype=np.float32),(nc,nb))
        Y = 0.0 * (np.pi/180.0)                                  # Initialize YPR
        P = 0.0 * (np.pi/180.0)
        R = 0.0 * (np.pi/180.0)
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
        # Make the overlap_vector into an array where the 1st dim. lines up sequentially
        # with the channels
        overlaps_chan_seq = np.ones((nc,nb),dtype=np.float32)
        # Code has option of whether to process using single or multiple overlaps  
        if multi_OLs:
            # Load all possible overlap tables into array [nOLs, nchans, nbins]
            # First user-defined defined configuration for this flight
            overlap_dict, overlap_map = load_overlap_configuration(overlaps_config_file) 
            # Now load all the OL tables into an array
            all_overlaps = load_all_possible_overlaps(overlap_dict,nb,nc)        
        else: # intializations say only using 1 overlap table
            # Load the overlap table
            overlap_vector = read_in_overlap_table(olap_dir+overlap_file)
            # Overlap is one long array, where first nbins (833) are
            # 355, next nbins are 532, and next (final) nbins are 1064.
            overlaps = overlap_vector.reshape((nwl,nb))
            for chan in range(0,nc):
                overlaps_chan_seq[chan,:] = overlaps[wl_map[chan],:]
        # Open the hdf5 file and create the datasets
        hdf5_fname = L1_dir+'NRB_'+proj_name+'_'+flt_date+'_'+Nav_source+'.hdf5'
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
            current_Nav_UTC = Nav_match['UTC_Time'] # Nav_match fields dependant on dataset
            Y = 0.0# Nav_match['yaw'] * (pi/180.0) Woah! This isn't what I think yaw should be [3/22/18]
            P = Nav_match['pitch'] * (np.pi/180.0)
            R = Nav_match['roll'] * (np.pi/180.0) 
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
            current_Nav_UTC = Nav_match['UTC_Time'] # Nav_match fields dependant on dataset
            # NOTE on next line: Invalids are -999.9, so max should choose valid if exists
            Nav_save[i1f]['GPS_alt'] = max(Nav_match['GPS_alt_msl'],Nav_match['GPS_alt'])
            Y = Nav_match['drift'] * (np.pi/180.0) 
            P = Nav_match['pitch'] * (np.pi/180.0)
            R = Nav_match['roll'] * (np.pi/180.0)  
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
            current_Nav_UTC = Nav_match['UTC_Time'] # Nav_match fields dependant on dataset
            Nav_save[i1f]['GPS_alt'] = Nav_match['GPS_alt']
            Y = 0.0 #Nav_match['drift'] * (np.pi/180.0) 
            P = Nav_match['pitch'] * (np.pi/180.0)
            R = Nav_match['roll'] * (np.pi/180.0)
            # Write 'lat,lon,UTC' out to file
            GEOS_out_str = '{0:.4f},{1:.4f},{2}\n'.format(Nav_match['lon'],Nav_match['lat'],time_str[:19])
            GEOS_fobj.write(GEOS_out_str)

        elif Nav_source == 'cls':

            # Find matching Nav data point. Use index map to match.
            # NOTE: Nav_save is an IWG1-style structure (see immutable_dat_structs.py)
            Nav_match = cls_nav_data_all[interp2orig_indicies[i1f]]
            Nav_save[i1f]['roll'] = (Nav_match['RollAngle'] + cpl_roll_offset)
            Nav_save[i1f]['pitch'] = ((Nav_match['PitchAngle'] + cpl_pitch_offset) + ((Nav_match['RollAngle'] + cpl_roll_offset) * cpl_pitch_roll_factor))
            Nav_save[i1f]['drift'] = (Nav_match['TrackAngleTrue'] - Nav_match['TrueHeading'])
            # Need to make sure drift is in the correct units after Track/Heading Difference Calculated
            if Nav_save[i1f]['drift'] > 180.0 :
                Nav_save[i1f]['drift'] = (Nav_save[i1f]['drift'] - 360.0)
            
            Nav_save[i1f]['heading'] = Nav_match['TrueHeading']
            Nav_save[i1f]['lon'] = Nav_match['GPS_Longitude']
            Nav_save[i1f]['lat'] = Nav_match['GPS_Latitude']
            # The following line is necessary to get a clean, interpretable,
            # time to write out to the HDF5 file. It will essetially be a
            # byte array for IDL.
            time_str = Nav_match['UTC_Time'].strftime("%Y-%m-%dT%H:%M:%S.%f")
            Nav_save[i1f]['UTC'] = np.asarray(list(time_str.encode('utf8')),dtype=np.uint8)
            current_Nav_UTC = Nav_match['UTC_Time'] # Nav_match fields dependant on dataset
            Nav_save[i1f]['GPS_alt'] = Nav_match['GPS_Altitude']
            Y = Nav_save[i1f]['drift'] * (np.pi/180.0) 											     #Yaw now incorporated!! 3/30/2020
            P = ((Nav_match['PitchAngle'] + cpl_pitch_offset) + ((Nav_match['RollAngle'] + cpl_roll_offset) * cpl_pitch_roll_factor))  * (np.pi/180.0)  #Offsets added 3/30/2020
            R = (Nav_match['RollAngle'] + cpl_roll_offset) * (np.pi/180.0)
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
            print(Nav_save[i1f]['UTC'],' ONA: ',ONA*(180.0/np.pi))
            print("!!!!!!! WARNING !!!!!!!")
            print(attention_bar)
            bg = np.zeros((nc,nb))
        
        range_cor_af_counts_float32 = ( ( counts_float32 - bg )
                                          * bins_range**2  ) 
        
        # Apply overlap correction here. I would have applied it at the same time as the
        # deadtime, except CPL's overlaps can be funky in far ranges, which interferes with
        # background sub.
        if multi_OLs: # Only if using multiple overlap tables
            if (current_Nav_UTC > overlap_map[OL_map_indx][0][1]) & (OL_map_indx < len(overlap_map)-1):
                print('UTC = '+current_Nav_UTC.strftime("%Y-%m-%dT%H:%M:%S")+'. Now using the following overlap tables...')
                print('355 nm: '+overlap_dict[overlap_map[OL_map_indx][1][0]])
                print('532 nm: '+overlap_dict[overlap_map[OL_map_indx][1][1]])
                print('1064 nm: '+overlap_dict[overlap_map[OL_map_indx][1][2]]+'\n')
                OL_map_indx += 1
                # Just do the channels w/o a loop for now [11/20/19]
                overlaps_chan_seq[0,:] = all_overlaps[overlap_map[OL_map_indx][1][0],0,:]
                overlaps_chan_seq[1,:] = all_overlaps[overlap_map[OL_map_indx][1][1],1,:]
                overlaps_chan_seq[2,:] = all_overlaps[overlap_map[OL_map_indx][1][2],2,:]
                overlaps_chan_seq[3,:] = all_overlaps[overlap_map[OL_map_indx][1][2],3,:]
            else:
                overlaps_chan_seq[0,:] = all_overlaps[overlap_map[OL_map_indx][1][0],0,:]
                overlaps_chan_seq[1,:] = all_overlaps[overlap_map[OL_map_indx][1][1],1,:]
                overlaps_chan_seq[2,:] = all_overlaps[overlap_map[OL_map_indx][1][2],2,:]
                overlaps_chan_seq[3,:] = all_overlaps[overlap_map[OL_map_indx][1][2],3,:]
        range_cor_af_counts_float32 = range_cor_af_counts_float32 / overlaps_chan_seq
    
        # Put the bins in the fixed altitude frame
        # Updated to use split_bins_onto_fixed_frame, as of v2 [4/22/20]
        counts_ff[:,i1f,:] = split_bins_onto_fixed_frame(ffrme,bin_alts,range_cor_af_counts_float32)
    
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
        NRB[chan,:,:] = correct_raw_counts(counts_ff[chan,:,:],EMsubmit,None,
            i1f,nb_ff,ff_bg_st_bin,ff_bg_ed_bin,'NRB_no_range')

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
        print('\nXXXXXXXXXXXXXXXXXXXXXXX')
        deleted_recs = nr_1file - tot_good_recs
        print(str(deleted_recs).strip()+' records deleted!')
        print('From file: '+CLS_file)
        print('\nXXXXXXXXXXXXXXXXXXXXXXX')
        i = i - (nr_1file - tot_good_recs)
        nr_1file = tot_good_recs

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

        bin_numbers = np.digitize(CLS_UnixT_float64_new,time_bins)
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
            #print('iter: ',u[tb],CLS_UnixT_float64_new[rr])
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

        # Expand dataset sizes to accomodate next input CLS file
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

        # Expand dataset sizes to accomodate next input CLS file
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
    print('\nDone with file: '+CLS_file+'\n')
    print('**********************************\n')

      
print('Main L1A execution finished at: ',DT.datetime.now())

# Write out any final parameters to the HDF5 file that needed to wait

num_recs_dset[:] = nrecs
hdf5_file.close()
GEOS_fobj.close()
print("NRB has been written to the HDF5 file:"+hdf5_fname)

# The final message

print('Total raw profiles processed: '+str(i))
print("cpl_l1a.py has finished normally.")
