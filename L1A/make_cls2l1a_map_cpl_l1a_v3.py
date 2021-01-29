# CPL L1A processing code.
# NOTE: The most current version of the L1A code must always be named
#       "cpl_l1a.py." All versions of the code must be saved as
#       "cpl_l1a_vX.py" but not actually be source file that's executed.

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
# 2) Nav-to-CLS time synchronization {per file basis}
# 3) Initialization of arrays & HDF5 parameters {per file basis}
# 4) Nav-to-CLS record matching {per raw record basis}
# 5) Computation/population of various parameters,
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
# arrays needed for time synchronization. The later function isn't used
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
# 2) Nav-to-CLS time synchronization {per file basis}
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
# HDF5 data set objects are initialized. The two output files, HDF5 and
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
# 5) Computation/population of various parameters,
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
# Here's the logic behind this averaging section:
# The data are averaged one file at a time, with any leftover records
# carried-over until the next file is processed. The carried-over records
# are then averaged with the first 'x' records of the next file, where 'x'
# represents the number of records required to meet the time-average criterion.
# The logic is a bit tricky when processing the first and last files, and these
# details are dealt with using the "first_read" and "cutbegin" variables. Also,
# n_expand (how much to expand data sets) and expanded_length variables are
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
# [6/12/18] Code modified to accommodate wonky overlap
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
# code handled the HDF5 data set resizing of the first file. I added
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
# arose in the averaging section when the data gap occurred between 2 consecutive CLS files.
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
# Nav coordinate array, "ix," with the improper data rate. I was comparing
# a computed data rate with the initialized one, and using the computed
# one if it was different from the computed one by an arbitrary tolerance.
# With a tolerance (CLS_hz_tol in initializations) set to 0.05 Hz, the
# aforementioned PELIcoe flight's computed data rate of 0.03833 Hz slipped
# through. A difference of 0.03833 Hz will result in a 690 record difference
# if it's really off from the true data rate by that much! This can explain
# the symptom Andrew Kupchock first noticed.
# The following changes were made for now...
# - tot_elapsed seconds now computed using Nav time instead of the CLS clock time.
# - Data rate in initializations file is now always used to create coordinates
#   of interpolated Nav array. User is now only warned if computed data
#   rate differs too much from input data rate.
#
# [5/12/20] *v3 start*
# NaN GPS altitudes experienced in 22 Oct 19 data. Numpy isfinite() function
# is now used in alt_cutoff check block.
#
# [6/18/2020]
# Added code to automate calculation of average frequency in the entire dataset
# Also added variable widths of time between records based on the calculated frequency
# M_array created for scenario of data missing and not missing.
# nav_interp Time array definition was redefined. SEAC4RS campaign flights are
# good candidates to test functionality of this code.
#
# [7/8/20]
# Fixed logic in 'if' block of code that handles invalid times in the
# Nav_source=='cls' pre-processing block. Whatever I did before didn't make
# any sense. It seems like I was assuming that invalid times could only
# occur at the beginning of flights. Now, a more sophisticated method of
# patching over the invalid times takes place.
# Also in this update, minor tweak to 'make_frequency_array' code.
#
# [7/22/20] Flaw in averaging code logic
# I fixed a flaw in the averaging code logic today.
# This flaw was caused by the band-aid created on [1/31/19].
# In those situations, I actually want it to go into the "not last_file"
# IF block because it is unknown whether the data in the last averaging bin
# (aka time bin) needs to go into the first averaging bin of the next file.
# The flaw caused the code to lock out of using a transfer bin (trans_bin)
# after the first instance that an averaging bin boundary was exactly 
# synced with a file boundary.
# I changed the logic so that a transfer bin is used EVERY file, except
# the last file (not last_file condition block). The first averaging bin,
# corresponding to element zero in the averaged-data arrays, is ALWAYS
# reserved for data from the last file's bin. Now, when the time boundary
# and file boundary line up, the array's get expanded to accomodate the
# transfer bin. So it's now possible for a transfer bin to contain a full
# secs2avg worth of raw profiles.
# ALSO, posgap_Nav now used in  CPL Nav pre-processing block.
#
# [1/29/21]
# New cpl_l1a version? 
# I am writing this note from "make_cls2l1a_map_cpl_l1a_v3.py."
# This code is cpl_l1a_v3.py, but with additional arrays that track 
# exactly which CLS records go to which L1A records.
# I fixed an error with this code that occurred when the time bin boundaries
# line up with a file boundary. In this case, trans_bin and first bin become
# separate L1A records instead of merging. This created havoc tracking CLS-to-L1A.
# Even more special care around masking out min_avg_profs records had to be 
# taken with the tracking arrays (L1Arec_nums_1file!).
# This code will form the underpinnings of cpl_l1a_v4.py because I want this
# record mapping functionality to be part of the standard L1A process.


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
    single_reasonable_value = np.mean(data_array[keep_indxs])
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

        cls_data = read_in_cls_data(cls_files[fc].strip(), nbins, flt_date, bad_cls_nav_time_value)

        diff = cls_data['meta']['Nav']['UTC_Time'] - cls_data['meta']['Header']['ExactTime']
        diff_float_secs = np.zeros(diff.shape[0], dtype=np.float32)
        for i in range(0, diff.shape[0]):
            diff_float_secs[i] = diff[i].total_seconds()

        time_offsets.append(diff_float_secs.mean())

        fc += 1

    print("Computed time offsets: ", time_offsets)

    if abs(time_offsets[0] - time_offsets[1]) > 2.0:
        print("The time offsets between 2 sample files differs by too much.")
        print("Enter 0 for first time offset, 1 for second.")
        print("Time offsets: ", time_offsets)
        choice = int(input("Enter your choice."))
        time_offset = time_offsets[choice]
    else:
        time_offset = np.mean(time_offsets)

    if abs(time_offset) > 80.0:
        print(attention_bar)
        print("!!!!!!! ATTENTION !!!!!!!")
        print("The computer time offset is " + str(time_offset).strip() + " seconds.")
        print("!!!!!!! ATTENTION !!!!!!!")
        print(attention_bar)
        # pdb.set_trace()

    return time_offset


# New Function added to fix time offset in Instrument Time
def Convert_CLS_Time_for_IWG1(cls_files):
    """ This function computes the time offset between the IWG1 Nav time
      and the CLS Instrument time ("ExactTime").

      Takes the string, full-path names of a list of CLS files as its
      only argument. Returns the time offset in seconds.

      Add this offset to the Instrument time.
  """

    fc = 0
    time_offsets = []

    while fc < 2:

        cls_data = read_in_cls_data(cls_files[fc].strip(), nbins, flt_date, bad_cls_nav_time_value)

        diff = cls_data['meta']['Nav']['UTC_Time'] - cls_data['meta']['Header']['ExactTime']
        diff_float_secs = np.zeros(diff.shape[0], dtype=np.float32)
        for i in range(0, diff.shape[0]):
            diff_float_secs[i] = diff[i].total_seconds()

        time_offsets.append(diff_float_secs.mean())

        fc += 1

    print("Computed time offsets: ", time_offsets)

    if abs(time_offsets[0] - time_offsets[1]) > 2.0:
        print("The time offsets between 2 sample files differs by too much.")
        print("Enter 0 for first time offset, 1 for second.")
        print("Time offsets: ", time_offsets)
        choice = int(input("Enter your choice."))
        time_offset = time_offsets[choice]
    else:
        time_offset = np.mean(time_offsets)

    if abs(time_offset) > 80.0:
        print(attention_bar)
        print("!!!!!!! ATTENTION !!!!!!!")
        print("The computer time offset is " + str(time_offset).strip() + " seconds.")
        print("!!!!!!! ATTENTION !!!!!!!")
        print(attention_bar)
        # pdb.set_trace()

    return time_offset


def create_cls_interp_unixt_array(single_cls_file, cls_meta_data_all,
                                  firstlast_flag, Navi, firstlast_trunc=None):
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

    cls_data_1file = read_in_cls_data(single_cls_file, nbins, flt_date, bad_cls_nav_time_value)
    nr0 = cls_data_1file.shape[0]

    if Nav_source == 'cls':

        if firstlast_flag == -1:
            cls_data_1file = cls_data_1file[firstlast_trunc[0]:]
        if firstlast_flag == 1:
            cls_data_1file = cls_data_1file[:nr0 - firstlast_trunc[1]]
        good_data_len = Navi[1] - Navi[0]
        cls_data_1file['meta']['Nav'] = cls_meta_data_all['Nav'][Navi[0]:Navi[1]]
        cls_nav_data_1file = cls_data_1file['meta']['Nav']
        del_t = np.zeros(cls_nav_data_1file.shape[0], dtype=DT.datetime)
        del_t[1:] = cls_nav_data_1file['UTC_Time'][1:] - cls_nav_data_1file['UTC_Time'][:-1]
        del_t[0] = del_t[1] * 0.0

        # Compute # of elapsed seconds since first record (float)

        del_secs = np.zeros(cls_nav_data_1file.shape[0], dtype=np.float64)
        for k in range(0, cls_nav_data_1file.shape[0]): del_secs[k] = del_t[k].total_seconds()

        elapsed_secs = np.cumsum(del_secs)
        tot_elapsed = np.max(elapsed_secs)
        UnixT_epoch = DT.datetime(1970, 1, 1)

        # Also, you'll need "Unix Times" corresponding to the original sampling rate

        CLS_UnixT_float64_orig = (cls_nav_data_1file['UTC_Time'][0] - UnixT_epoch).total_seconds() + elapsed_secs

        interp_UnixT = np.copy(CLS_UnixT_float64_orig)

    else:

        CLS_UnixT_float64_orig = np.zeros(cls_data_1file.shape[0], dtype=np.float64)
        UnixT_epoch = DT.datetime(1970, 1, 1)

        for k in range(0, cls_data_1file.shape[0]):
            CLS_UnixT_float64_orig[k] = (cls_data_1file['meta']['Header']['ExactTime'][k] - UnixT_epoch).total_seconds()

        xp = CLS_UnixT_float64_orig - CLS_UnixT_float64_orig[0]
        fp = CLS_UnixT_float64_orig
        ix = np.arange(np.min(xp), np.max(xp) + (1.0 / CLS_hz), 1.0 / CLS_hz)
        interp_UnixT = np.interp(ix, xp, fp)

    return [interp_UnixT, CLS_UnixT_float64_orig, cls_data_1file]


# Dictionary of letters and integers
AlphaNum = {}
for i in range(65, 91): AlphaNum[chr(i)] = i - 65


def load_all_possible_overlaps(overlap_dict, nb, nc):
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
    all_overlaps = np.zeros((nOLs, nc, nb), dtype=np.float32)

    for i in range(0, nOLs):
        # Load the overlap table
        overlap_file = overlap_dict[i]
        overlap_vector = read_in_overlap_table(olap_dir + overlap_file)
        # Overlap is one long array, where first nbins (833) are
        # 355, next nbins are 532, and next (final) nbins are 1064.
        overlaps = overlap_vector.reshape((nwl, nb))
        # Make the overlap_vector into an array where the 1st dim. lines up sequentially
        # with the channels
        overlaps_chan_seq = np.ones((nc, nbins), dtype=overlaps.dtype)
        for chan in range(0, nc):
            overlaps_chan_seq[chan, :] = overlaps[wl_map[chan], :]
        all_overlaps[i, :, :] = overlaps_chan_seq

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
        Lc = 0  # line counter
        Cc = 0  # codes element counter
        Tc = 0  # times element counter
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
                pass  # Skip this line
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
        Ic = 1  # list Item counter
        for item in OL_times:
            if item[0] == 'Overlaps key':
                nOL = list_len - Ic
            Ic += 1
        OL_keys = OL_times[-1 * nOL:]  # Last
        OL_times = OL_times[:-1 * nOL - 1]
        print('{0:d} overlap tables detected.'.format(nOL))

        overlap_dict = {}
        for item in OL_keys: overlap_dict[AlphaNum[item[0].strip()]] = item[1].strip()

        OL_codes = OL_codes[flt_indx]
        OL_times = OL_times[flt_indx]

        Ic = 0  # list Item counter
        overlap_map = []
        for codes, times in zip(OL_codes, OL_times):
            if (codes == '') | (codes == '\n'): break
            if Ic > 0:
                split_codes = codes.split('/')
                split_times = times.split('/')
                overlap_map.append([
                    [DT.datetime.strptime(split_times[0], "%Y-%m-%dT%H:%M:%S"),
                     DT.datetime.strptime(split_times[1], "%Y-%m-%dT%H:%M:%S")],
                    [AlphaNum[split_codes[1]], AlphaNum[split_codes[2]], AlphaNum[split_codes[0]]]])
            Ic += 1

    return overlap_dict, overlap_map


# Start main execution here <-------------------------------------------
print('Starting main L1A execution at: ', DT.datetime.now())

# File detailing how to map OLs to different portions of the flight
overlaps_config_file = 'overlaps_configuration.csv'
if 'make_frequency_array' not in locals():
    make_frequency_array = 1
# Create and load file list for CLS data
CLS_file_list = 'processing_file_list.txt'
search_str = '*.cls'
create_a_file_list(CLS_file_list, search_str, raw_dir)
with open(CLS_file_list) as CLS_list_fobj:
    all_CLS_files = CLS_list_fobj.readlines()
nCLS_files = len(all_CLS_files)

# Load the shared library of C functions
np_clib = np.ctypeslib.load_library(clib, clib_path)
# DEM shared C++ library
DEM_lib_file = DEM_lib_path + DEM_Cpp_name
DEM_lib = ctypes.cdll.LoadLibrary(DEM_lib_file)
DEM_file = (DEM_dir + DEM_name).encode('utf8')
DEM_file_ctype = ctypes.c_char_p(DEM_file)


# Defined C++ DEM read function output


class DEM_out(ctypes.Structure):
    _fields_ = [("Elev", ctypes.c_float), ("LOCI", ctypes.c_ubyte)]


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

    del_t = np.zeros(gps_data_all.shape[0], dtype=DT.datetime)
    del_t[1:] = gps_UTC[1:] - gps_UTC[:-1]
    del_t[0] = del_t[1] * 0.0

    # Compute # of elapsed seconds since first record (float)

    del_secs = np.zeros(gps_data_all.shape[0], dtype=np.float64)
    for k in range(0, gps_data_all.shape[0]): del_secs[k] = del_t[k].total_seconds()

    elapsed_secs = np.cumsum(del_secs)
    tot_elapsed = np.max(elapsed_secs)
    m = 1.0 / CLS_hz  # basically, the time rez in secs you'd like to get to
    ix = np.arange(0, tot_elapsed, m)  # x coordinates of the interpolated values
    n_interp = ix.shape[0]
    gps_interp = np.zeros(n_interp, dtype=gps_struct)
    gps_UTC_interp = np.zeros(n_interp, dtype=DT.datetime)  # note separate array
    UnixT_epoch = DT.datetime(1970, 1, 1)
    Nav_interp_T_float64 = np.zeros(n_interp, dtype=np.float64)
    for field in gps_data_all.dtype.names:
        gps_interp[field] = np.interp(ix, elapsed_secs, gps_data_all[field])
    # Now handle population of interp time field
    for k in range(0, n_interp):
        gps_UTC_interp[k] = gps_UTC[0] + DT.timedelta(seconds=ix[k])
        Nav_interp_T_float64[k] = (gps_UTC_interp[k] - UnixT_epoch + time_offset_UnixT).total_seconds()

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
    nav2cls_indx = np.zeros((1000, 2))
    FL_trunc = None
    usable_file_range = [0, 1000]

    # Compute offset between IWG1 time and CLS Instrument Time ("ExactTime")
    # This data is the IWG1 data.
    iwg1_cls_t_offset = compute_time_offset_IWG1_CLS(all_CLS_files)

    # Load the entire nav dataset into memory
    nav_data_all = read_entire_nav_dataset('*.nav*')

    # Interpolate data records to match CPL's data rate.
    # 10 Hz as of 2/2/18. Data rate set in initializations.

    del_t = np.zeros(nav_data_all.shape[0], dtype=DT.datetime)
    del_t[1:] = nav_data_all['UTC'][1:] - nav_data_all['UTC'][:-1]
    del_t[0] = del_t[1] * 0.0
    pdb.set_trace()

    # Compute # of elapsed seconds since first record (float)

    del_secs = np.zeros(nav_data_all.shape[0], dtype=np.float64)
    for k in range(0, nav_data_all.shape[0]): del_secs[k] = del_t[k].total_seconds()

    elapsed_secs = np.cumsum(del_secs)
    tot_elapsed = np.max(elapsed_secs)
    m = 1.0 / CLS_hz
    ix = np.arange(0, tot_elapsed, m)  # x coordinates of the interpolated values
    n_interp = ix.shape[0]
    nav_interp = np.zeros(n_interp, dtype=nav_struct)
    Nav_interp_T_float64 = np.zeros(n_interp, dtype=np.float64)
    UnixT_epoch = DT.datetime(1970, 1, 1)
    for field in nav_data_all.dtype.names:
        if (field == 'UTC') or (field == 'UnixT'): continue
        nav_interp[field] = np.interp(ix, elapsed_secs, nav_data_all[field])
    # Now handle population of interp time field
    for k in range(0, n_interp):
        # nav_interp['UnixT'][k] = nav_data_all['UnixT'][0] + DT.timedelta(seconds=ix[k])
        # offsets in following line convert to CLS UTC
        nav_interp['UTC'][k] = (nav_data_all['UTC'][0] + DT.timedelta(seconds=ix[k]) + DT.timedelta(seconds=nudge))
        Nav_interp_T_float64[k] = (nav_interp['UTC'][k] - UnixT_epoch).total_seconds()

    # NOW, set the array that will be used in processing, nav_data_all equal
    # to the interpolated array that was just created, nav_interp.

    nav_data_all = nav_interp

    # This variable will be the same no matter what the Nav_source
    Nav_hz = nav_hz

elif Nav_source == 'iwg1':

    # "Dummy out" these variables, which are only need valid values for CLS-Nav processing
    cls_meta_data_all = None
    nav2cls_indx = np.zeros((1000, 2))
    FL_trunc = None
    usable_file_range = [0, 1000]

    # Compute offset between IWG1 time and CLS Instrument Time ("ExactTime")
    iwg1_cls_t_offset = compute_time_offset_IWG1_CLS(all_CLS_files)

    # Read in entire IWG1 dataset. Typically contained in single file.
    IWG1_file = raw_dir + IWG1_file
    IWG1_data = read_in_IWG1_data(IWG1_file, est_IWG1_recs)

    # Interpolate data records to match CAMAL's data rate.
    # 10 Hz as of 2/2/18. Data rate set in initializations.

    del_t = np.zeros(IWG1_data.shape[0], dtype=DT.datetime)
    del_t[1:] = IWG1_data['UTC'][1:] - IWG1_data['UTC'][:-1]
    del_t[0] = del_t[1] * 0.0

    # Compute # of elapsed seconds since first record (float)

    del_secs = np.zeros(IWG1_data.shape[0], dtype=np.float64)
    for k in range(0, IWG1_data.shape[0]): del_secs[k] = del_t[k].total_seconds()

    elapsed_secs = np.cumsum(del_secs)
    tot_elapsed = np.max(elapsed_secs)
    m = 1.0 / CLS_hz
    ix = np.arange(0, tot_elapsed, m)  # x coordinates of the interpolated values
    n_interp = ix.shape[0]
    IWG1_interp = np.zeros(n_interp, dtype=IWG1_struct)
    UnixT_epoch = DT.datetime(1970, 1, 1)
    Nav_interp_T_float64 = np.zeros(n_interp, dtype=np.float64)
    for field in IWG1_data.dtype.names:
        if field == 'UTC': continue
        IWG1_interp[field] = np.interp(ix, elapsed_secs, IWG1_data[field])
    # Now handle population of interp time field
    for k in range(0, n_interp):
        # offsets in following line convert to CLS UTC
        IWG1_interp['UTC'][k] = (IWG1_data['UTC'][0] + DT.timedelta(seconds=ix[k]))
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
    cls_meta_data_all, nav2cls_indx, FL_trunc, usable_file_range = read_entire_cls_meta_dataset(raw_dir, file_len_recs,
                                                                                                nbins, flt_date,
                                                                                                bad_cls_nav_time_value)
    # meta_save = np.copy(cls_meta_data_all) #uncomment to compare processed to original
    cls_nav_data_all = np.copy(cls_meta_data_all['Nav'])

    UnixT_epoch = DT.datetime(1970, 1, 1)

    # Identify all the unique UTC_Times and their starting indexes
    u, ui, ncounts = np.unique(cls_nav_data_all['UTC_Time'], return_index=True, return_counts=True)
    # Bad time values will sink to the front of this "unique" array

    # The follow code deals with there being bad Nav clock times in the data
    # One example of where this code will execute is ACEPOL 18-604 25Oct17
    # Another example is SEAC4RS 13-955 19aug13
    # As of 7/8/2020, this code block is a sophisticated band-aid.
    # It takes the nav record before a stretch of invalid times and simply
    # repeats the nav record across that stretch.
    if u[0] == bad_cls_nav_time_value:
        print(attention_bar)
        print("\nA(Some) UTC_Time(s) in the Nav record(s) was(were) invalid.")
        print("Attempt is being made to take care of this.")
        print("Setting first invalid CLS-Nav UTC_Time record equal to")
        print("the first valid CLS-Nav UTC_Time record.\n")
        print(attention_bar)
        # Hopefully there are very few bad cls times.
        # Patch across the bad_cls_nav_time_value regions in the array
        good_time_indexes = np.where(cls_nav_data_all['UTC_Time'] != bad_cls_nav_time_value)[0]
        index_diffs = np.arange(good_time_indexes.shape[0], dtype=np.uint32)
        index_diffs[0] = 0
        index_diffs[1:] = np.diff(good_time_indexes)
        missing_time_index = np.where(index_diffs > 1)[0]
        for mti in missing_time_index:
            pregap_Nav = cls_nav_data_all[good_time_indexes][mti - 1]
            posgap_Nav = cls_nav_data_all[good_time_indexes][mti]
            for i in range(1, index_diffs[mti]):
                ii = i + good_time_indexes[mti - 1]
                cls_nav_data_all[ii] = posgap_Nav
        # Once the bad times are seen again in the individual-file CLS data in
        # the file loop ('core processing'), then the process_start initializations
        # variable should take care of things.
        # Rerun line from before the if statements to have the correct indices if there were problem records
        u, ui, ncounts = np.unique(cls_nav_data_all['UTC_Time'], return_index=True, return_counts=True)

    # I found at least one case where inst. clk. was correct, but data system appeared to be
    # grabbing an old Nav record for a bit. This block of code determines if this is the case
    # and overwrites those values with the first valid values.
    if (u[1] - u[0]).total_seconds() > 3:
        print("\nSome(A) UTC_Time(s) in the Nav records were(was) from stale packets.")
        print("Attempt is being made to take care of this.")
        print("Setting first stale CLS-Nav UTC_Time record equal to")
        print("the first non-stale CLS-Nav UTC_Time record.\n")
        print('IF the code stops here you can proceed by pressing c and then Enter but notify Andrew/Patrick')
        pdb.set_trace()
        cls_nav_data_all['UTC_Time'][0] = u[1]
        u = u[1:]
        ui = ui[1:]
        ncounts = ncounts[1:]

    n_u = u.shape[0]
    # Convert the unique UTC_Times to Unix Times
    Nav_UnixT = np.zeros(n_u, dtype=np.float64)
    Nav_unique = np.zeros(n_u, dtype=CLS_decoded_nav_struct)
    for i in range(0, n_u):
        Nav_UnixT[i] = (cls_nav_data_all['UTC_Time'][ui[i]] - UnixT_epoch).total_seconds()
        Nav_unique[i] = cls_nav_data_all[ui[i]]

        # Compute # of elapsed seconds since first record (float)
    tot_elapsed = Nav_UnixT[-1] - Nav_UnixT[0]
    if tot_elapsed > (max_flt_hours * 3600.0):
        tot_hours = tot_elapsed / 3600.0
        print(attention_bar)
        print("HEY! The data span " + str(tot_hours).strip() + " hours!")
        print("Consider trimming processing time range in init. file.")
        print('Enter "c" to continue processing anyway.')
        print('Enter "q" to safely quit.')
        print(attention_bar)
        pdb.set_trace()

        # Look at the time deltas between unique, you'll use this in a little bit
    # On 5/4/2020 updated to separate Instrument time from Nav Time u ui and ncounts arrays
    u_inst, ui_inst, ncounts_inst = np.unique(cls_meta_data_all['Header']['ExactTime'], return_index=True, return_counts=True)
    if u_inst[0] == bad_cls_nav_time_value:
        u_inst = u_inst[1:]
        ui_inst = ui_inst[1:]
        ncounts_inst = ncounts_inst[1:]
    ExactTime_ncounts = np.copy(ncounts_inst)
    ExactTime_ui = np.copy(ui_inst)
    inst_clk_deltas = delta_datetime_vector(u_inst)
    # Estimate the data rate from the Nav time and length of data array.
    # Warn user if computed data doesn't match initialized data rate.
    # computed_CLS_hz = cls_meta_data_all.shape[0] / tot_elapsed
    # new computed hz does not have as much of an issue with gaps in data from the laser being turned off 6/14/2020
    # starts at 1 to remove partial second at the beginning of data collection - USING FROM NAV! not Instrument clock
    computed_CLS_hz = np.mean(np.diff(ui[1:]))
    if abs(computed_CLS_hz - CLS_hz) > CLS_hz_tol:
        print(attention_bar)
        print("!!!!!!! WARNING !!!!!!!")
        print("The estimated (init. file) and computed data rates differ by too much.")
        print("estimated: ", CLS_hz, " computed: ", computed_CLS_hz)
        print("Make sure CLS_hz from initializations is correct")
        print("!!!!!!! WARNING !!!!!!!")
        print('Enter "c" to continue processing anyway.')
        print('Enter "q" to safely quit.')
        print('If you continue, data rate from initializations will be used.')
        print(attention_bar)
        pdb.set_trace()
    else:
        print('Init. file CLS data rate is sufficiently close to computed data rate.')
        print('CLS_hz: ', CLS_hz, ' computed_CLS_hz: ', computed_CLS_hz)
        print('To be clear, code will use CLS data rate from initializations (CLS_hz)\n')
    m = 1.0 / CLS_hz
    # interp_start_time_offset calculation moved up on 5/4/2020 for availability in ix generation
    # Create the x-coordinates for the interp data, which will have units of Unix Time
    print("Making the reasonably dangerous assumption that earliest time is first rec")
    print("and latest time is last rec.")
    rough_expected_recs_sec = int(CLS_hz)
    interp_start_time_offset = 0.0
    # if ncounts[0] < CLS_hz:      # 5/5/2020 Removed for instances that ncounts is larger than expected 1 second
    interp_start_time_offset = (rough_expected_recs_sec - ncounts[0]) / CLS_hz  # Now calculated using the ncounts from NAV DATA

    # Added section below so that we can correctly determine the end of time in the data
    if ncounts[-1] < CLS_hz:  # Now calculated using the ncounts from NAV DATA
        interp_end_time_offset = (ncounts[-1] - 1) / CLS_hz
    '''
    This block of code finds the correct average frequency of the instrument for the flight based on the records 
    that were generated if there are no data gaps (and even if there are, provided there is enough time for a measurement)
    '''
    instrument_clk_threshold = inst_clk_rez
    # Check for gaps in instrument clock from threshold set with instrument_clk_threshold
    if np.max(np.diff(u_inst[1:])).total_seconds() <= instrument_clk_threshold:
        # Case where everything is perfect and there are no sizable gaps
        ratio = 2
        # Calculate mean frequency within tolerance by while loop
        while abs(ratio - 1) > freq_tol:  # Limitation set to 1e-10 - could be lower for more rigidity
            total_time = (Nav_UnixT[-1] - Nav_UnixT[0]) + ((ncounts[-1] - 1) - (rough_expected_recs_sec - ncounts[0])) * m
            new_m = total_time / np.sum(ncounts)
            ratio = new_m / m
            m = new_m
            print('M equals: ', m)
            # Update time offsets for future use
            interp_start_time_offset = (rough_expected_recs_sec - ncounts[0]) * m
            interp_end_time_offset = (ncounts[-1] - 1) * m
        # ************************************************************** #
        if make_frequency_array:
            # Making array of time deltas based on local average data collection frequency
            # Removing the first sub-second with calculation starting at seconds unique time
            records_between_seconds = np.diff(ui[1:])
            time_between_seconds = np.diff(Nav_UnixT[1:])
            frequency = records_between_seconds / time_between_seconds
            # Take out the high frequency elements of the timing data - Averages to the best of it's ability at least
            # 'gap' number of time-updates (typically seconds)
            temp_freq_orig = np.zeros(len(frequency), dtype=float)
            for step in np.arange(0, len(frequency)):
                if step < gap:
                    temp_freq_orig[step] = np.mean(frequency[0:step + gap])
                else:
                    if step <= (len(frequency) - gap):
                        temp_freq_orig[step] = np.mean(frequency[step - gap:step + gap])
                    else:
                        temp_freq_orig[step] = np.mean(frequency[step - gap:-1])

            # from the now semi-cleaned up data. Take the average records per second to define a rate for spans of time
            # of the flight
            # step size is typically 100
            timespans = np.arange(0, len(temp_freq_orig), step_size)
            temp_freq = np.zeros(len(timespans), dtype=float)
            # Assign Averaged Frequency for all but the end
            for span in np.arange(0, len(timespans) - 1): temp_freq[span] = np.mean(
                temp_freq_orig[timespans[span]:timespans[span + 1] - 1])
            # Assign Average Frequency for the last bit less than step_size
            if span != len(temp_freq) + 1:
                if len(temp_freq_orig[timespans[-1]:-1]) > 0:
                    temp_freq[-1] = np.mean(temp_freq_orig[timespans[-1]:-1])
                else:
                    temp_freq[-1] = temp_freq[-2]
            # Add initial sub-second to array of m
            m_array = np.full((1, ui[1]), 1 / temp_freq[0])
            # Add calculated frequency averages by step to array of m values
            for section_index in np.arange(0, len(temp_freq) - 1):
                m_array = np.concatenate((m_array, np.full(
                    (1, ui[timespans[section_index + 1] + 1] - ui[timespans[section_index] + 1]),
                    1 / temp_freq[section_index])), axis=1)
            # Add last section of values to array of m values separated since width is different
            m_array = np.concatenate((m_array,
                                      np.full((1, (cls_nav_data_all['UTC_Time'].shape[0] - ui[timespans[-1] + 1])),
                                              1 / temp_freq[-1])), axis=1)
            # Adjust frequencies so that average frequency matches known average frequency
            m_array = np.squeeze(m_array)
            frequency_error = (m - np.mean(m_array))
            print("After averaging time and resetting frequency subsets the error was found to be ", frequency_error,
                  " seconds per record")
            m_array = m_array + frequency_error
            ix_variable = np.zeros(len(m_array))
            ix_variable[0] = (Nav_UnixT[0] + interp_start_time_offset)
            for element in np.arange(0, len(ix_variable) - 1):
                ix_variable[element + 1] = ix_variable[element] + m_array[element]
    else:
        # case where there is at least 1 gap of duration longer than instrument_clk_threshold
        # Build continuous instrument time array to pull gaps out efficiently (array of datetimes hard to deal with)
        Inst_Time = np.zeros(len(u_inst), dtype=np.float64)
        for i in np.arange(0, len(u_inst)):
            Inst_Time[i] = (cls_meta_data_all['Header']['ExactTime'][ui_inst[i]] - UnixT_epoch).total_seconds()
        # Use instrument clock to find separations - index in separation is before time jump
        inst_jumps = np.where(np.diff(Inst_Time) > instrument_clk_threshold)
        inst_jumps = inst_jumps[0]
        separations = ui_inst[inst_jumps]
        # steps is the number of records until the instrument time changes drastically
        # Location of clock deltas should be placed at record numbers (separations + steps_to_the_edge + 1)
        steps_to_the_edge = ncounts_inst[inst_jumps] - 1
        # Transfer into Navigation Clock Time frame
        before_jump_times = cls_nav_data_all['UTC_Time'][separations + steps_to_the_edge]
        after_jump_times = cls_nav_data_all['UTC_Time'][separations + steps_to_the_edge + 1]
        # Determine subseconds of before and after jump by finding ncounts of before_jump_times and after_jump_times
        nav_time_indices = np.zeros(len(separations), dtype=int)
        for j in np.arange(0, len(separations)):
            nav_time_indices[j] = np.where(ui <= (separations[j] + steps_to_the_edge[j]))[0][-1]
        # subseconds before or after could be as many as 10, subtract 1 for 0 case
        subseconds_before = ncounts[nav_time_indices] - 1
        subseconds_after = ncounts[nav_time_indices + 1]
        # Makes it so there is no subsecond subtraction if there is a jump right after another
        subseconds_after[np.concatenate(
            ((nav_time_indices[0:-1] == nav_time_indices[1:] - 1), np.array([False])))] = rough_expected_recs_sec
        separation_deltas = []
        for i in np.arange(0, len(separations)):
            separation_deltas.append((after_jump_times[i] - before_jump_times[i]).total_seconds() +
                                     (rough_expected_recs_sec - subseconds_after[i]) * m - (subseconds_before[i] * m))
        nav_jump_indices = np.zeros(len(separations), dtype=np.uint64)
        for j in np.arange(0, len(separations)):
            nav_jump_indices[j] = np.where(ui < separations[j])[0][-1]
        inst_clk_offsets = [0]
        if len(nav_jump_indices) > 1:
            # Removes jumps that immediately follow other jumps in the data
            repeat_removal = 1
            while repeat_removal <= (len(nav_jump_indices) - 1):
                if abs(nav_jump_indices[repeat_removal] - nav_jump_indices[repeat_removal - 1]) < 2:
                    inst_clk_offsets[repeat_removal - 1] = inst_clk_offsets[repeat_removal - 1] + \
                                                           (inst_jumps[repeat_removal] - inst_jumps[repeat_removal - 1])
                    inst_jumps = np.delete(inst_jumps, repeat_removal, axis=0)
                    nav_jump_indices = np.delete(nav_jump_indices, repeat_removal, axis=0)
                else:
                    repeat_removal += 1
                    inst_clk_offsets.append(0)
        # Add beginning and end to determine locations for analysis - note one at beginning to remove sub-second
        end_idxs = np.concatenate((nav_jump_indices, np.asarray([len(ui) - 1], dtype=np.uint64)))
        # Jumping over gap: 1 for gap 1 for difference between instrument and nav clock
        start_idxs = np.asarray(nav_jump_indices + inst_clk_offsets + 2, dtype=int)
        # skipping partial second from index 0 to 1;
        start_idxs = np.concatenate((np.asarray([1], dtype=int), start_idxs))
        total_times = np.zeros(len(end_idxs), dtype=float)
        # Minimum number of seconds for a timespan to use the calculated m value
        # min_jump typically 400
        # Calculate mean frequency in the time spans without instrument time gaps
        for k in np.arange(0, len(end_idxs)):
            total_times[k] = (cls_nav_data_all['UTC_Time'][ui[end_idxs[k]]] -
                              cls_nav_data_all['UTC_Time'][ui[start_idxs[k]]]).total_seconds()
        rec_spans = (ui[end_idxs] - ui[start_idxs])
        new_ms = total_times / rec_spans
        # Removes m from calculation if it is in a region where the time span is less than min_jump*inst_clk_freq
        span_weights = (rec_spans * (rec_spans >= (min_jump * rough_expected_recs_sec))) / \
                       np.sum(rec_spans * (rec_spans >= (min_jump * rough_expected_recs_sec)))
        new_m = np.sum(new_ms * span_weights)
        if new_m != 0:
            m = new_m
            ratio = 2
            # Calculate mean frequency within tolerance by while loop
            while abs(ratio - 1) > freq_tol:  # Limitation set to 1e-10 - could be lower for more rigidity
                total_time = (Nav_UnixT[-1] - Nav_UnixT[0]) - np.sum(separation_deltas) + \
                             (ncounts[-1] - (rough_expected_recs_sec - ncounts[0])) * m
                new_m = total_time / np.sum(ncounts)
                ratio = new_m / m
                m = new_m
                print('M equals: ', m)
                # Update time offsets for future use
                interp_start_time_offset = (rough_expected_recs_sec - ncounts[0]) * m
                interp_end_time_offset = (ncounts[-1] - 1) * m
            print('Calculated with Gaps - M average equals: ', m, ' using {0:.4f}% of the data '.
                  format(
                np.sum(rec_spans * (rec_spans >= (min_jump * rough_expected_recs_sec))) / len(cls_meta_data_all)))
        else:
            print("The GAPS in data are too frequent to update the frequency automatically!!!")
            print("Manual effort may be required for correction")
        # ************************************************************** #
        # Develop an m-array with time gaps
        if make_frequency_array:
            # Making array of time deltas based on local average data collection frequency with data gap issue
            records_between_seconds = np.zeros(1, dtype=np.uint64)
            time_between_seconds = np.zeros(1, dtype=float)
            for idx in np.arange(0, len(end_idxs)):
                records_between_seconds = np.concatenate((records_between_seconds,
                                                          np.diff(ui[start_idxs[idx]:int(end_idxs[idx] + 1)])))
                time_between_seconds = np.concatenate((time_between_seconds,
                                                       np.diff(Nav_UnixT[start_idxs[idx]:int(end_idxs[idx] + 1)])))
            # Remove first 0 from when the array was initialized
            records_between_seconds = records_between_seconds[1:]
            time_between_seconds = time_between_seconds[1:]
            frequency = records_between_seconds / time_between_seconds
            # Take out the high frequency elements of the timing data - Averages to the best of it's ability at least
            # 'gap' number of time-updates (typically seconds) - may be better if not spanning missing data
            temp_freq_orig = np.zeros(len(frequency), dtype=float)
            for step in np.arange(0, len(frequency)):
                if step < gap:
                    temp_freq_orig[step] = np.mean(frequency[0:step + gap])
                else:
                    if step <= (len(frequency) - gap):
                        temp_freq_orig[step] = np.mean(frequency[step - gap:step + gap])
                    else:
                        temp_freq_orig[step] = np.mean(frequency[step - gap:-1])
            # from the now semi-cleaned up data. Take the average records per second to define a rate for spans of time
            # of the flight
            timespans = np.arange(0, len(temp_freq_orig), step_size)
            temp_freq = np.zeros(len(timespans), dtype=float)
            # Assign Averaged Frequency for all but the end
            for span in np.arange(0, len(timespans) - 1): temp_freq[span] = np.mean(
                temp_freq_orig[timespans[span]:timespans[span + 1] - 1])
            # Assign Average Frequency for the last bit less than step_size
            if span != len(temp_freq) + 1:
                if len(temp_freq_orig[timespans[-1]:-1]) > 0:
                    temp_freq[-1] = np.mean(temp_freq_orig[timespans[-1]:-1])
                else:
                    temp_freq[-1] = temp_freq[-2]
            # Add initial sub-second to array of m
            m_array = np.full((1, ui[1]), 1 / temp_freq[0])
            # Add calculated frequency averages by step to array of m values
            # Skipped subseconds are assigned the same frequency of the surrounding data concurrently
            for section_index in np.arange(0, len(temp_freq) - 1):
                m_array = np.concatenate((m_array, np.full(
                    (1, ui[timespans[section_index + 1] + 1] - ui[timespans[section_index] + 1]),
                    1 / temp_freq[section_index])), axis=1)
            # Add last section of values to array of m values separated since width is different
            m_array = np.concatenate((m_array,
                                      np.full((1, (cls_nav_data_all['UTC_Time'].shape[0] - ui[timespans[-1] + 1])),
                                              1 / temp_freq[-1])), axis=1)
            # Adjust frequencies so that average frequency matches known average frequency
            m_array = np.squeeze(m_array)
            frequency_error = (m - np.mean(m_array))
            print("After averaging time and resetting frequency subsets the error was found to be ", frequency_error,
                  " seconds per record")
            m_array = m_array + frequency_error
            ix_variable = np.zeros(len(m_array))
            ix_variable[0] = (Nav_UnixT[0] + interp_start_time_offset)
            for element in np.arange(0, len(ix_variable) - 1):
                ix_variable[element + 1] = ix_variable[element] + m_array[element]
            # Recreating array of every single jump location -
            # removed repeat items previously to generate clean array for calculating steps of m-array
            inst_jumps = np.where(np.diff(Inst_Time) > instrument_clk_threshold)
            inst_jumps = inst_jumps[0]
            # Adds in expected time difference from the jumps in data collection one at a time for simplicity
            for jump_location_idx in np.arange(0, len(separation_deltas)):
                ix_variable[ui_inst[inst_jumps[jump_location_idx] + ncounts_inst[inst_jumps[jump_location_idx]]]:] += \
                    separation_deltas[jump_location_idx]
    '''
    This terminates the block of code where the instrument frequency is automatically calculated and an 
    m_array is created if desired with the make_frequency_array selector
    '''
    # 5/4/2020 ix array calculation updated to include start_time_offset and end_time_offset
    # Create x-ord array of Unix Times at which to compute interpolated values
    print('If frequency is manual then try ', computed_CLS_hz * float(nshots), 'Hz if no data is missing')

    ix = np.arange((Nav_UnixT[0] + interp_start_time_offset), ((Nav_UnixT[-1] + interp_end_time_offset) + m),
                   m)  # x coordinates of the interpolated values
    if make_frequency_array:
        # Use variable time-width profile array for interpolation
        ix = ix_variable
        # Reassign frequency to calculated average frequency which matches array
        CLS_hz = 1.0 / m
    n_interp = ix.shape[0]
    nav_interp = np.zeros(n_interp, dtype=CLS_decoded_nav_struct)
    Nav_interp_T_float64 = ix
    for field in cls_nav_data_all.dtype.names:
        if field == 'UTC_Time': continue
        nav_interp[field] = np.interp(ix, Nav_UnixT, Nav_unique[field])
    # Now handle population of interp time field
    # 5/4/2020 added timedelta of interp_start_time_offset to nav_interp UTC Time values
    # for k in range(0, n_interp):
    #     nav_interp['UTC_Time'][k] = (cls_nav_data_all['UTC_Time'][0] + DT.timedelta(seconds=interp_start_time_offset)
    #                                  + DT.timedelta(seconds=(m*k)))

    # Time interpolation changed on 6/17/2020
    for k in range(0, n_interp):
        nav_interp['UTC_Time'][k] = (UnixT_epoch + DT.timedelta(seconds=ix[k]))

    # Added 6/30/2020 to remove Data from around time jumps. Utilizes the altitude cutoff to assign an altitude below
    # the threshold - removing the data from being utilized
    if make_frequency_array and ('inst_jumps' in locals()):
        nudge_offset = int(nudge / m)
        for jump_location_idx in np.arange(0, len(inst_jumps)):
            nav_interp['GPS_Altitude'][(ui_inst[inst_jumps[jump_location_idx] + 
            ncounts_inst[inst_jumps[jump_location_idx]]] + 
            nudge_offset - 9 * rough_expected_recs_sec):(ui_inst[inst_jumps[jump_location_idx] + 
            ncounts_inst[inst_jumps[jump_location_idx]]] + nudge_offset + 1 * rough_expected_recs_sec)] = (alt_cutoff - 100.0)
    # Prepare nav_interp 'UTC_Time' array, where time gaps have been
    # masked out, then overwrite original Nav time.
    # 5/4/2020 replaced calculation offset_indx = np.abs((ix - ix[0]) - interp_start_time_offset).argmin()
    offset_indx = (rough_expected_recs_sec - ncounts[0])  # Uses Nav Time for ncounts
    # Added or to conditional for make frequency array for updates on 6/19/2020
    # DESCRIPTION OF FOLLOWING CODE...
    # The purpose of the code in this next paragraph is create the proper mapping
    # of the CLS data to the interpolated Nav data, which comes from the CLS data
    # itself when Nav_source=='cls', and there are gaps in the data.
    # The UTC from the interpolated Nav data is pulled into a new array, UTC_ovrwrt.
    # At this point, cls_meta_data_all has exactly the # of total valid (non-watchdog)
    # records in the entire flight. But before this paragraph was written, it only did the
    # thing in the else block - and this allows the CLS records to map to the Nav data inside
    # itself. cls_meta_data_all is thrown into the create_cls_interp_unixt_array function
    # in order allow the code to match the CLS data to itself (ie the interpolated Nav data). The
    # ORIG and INTERP Unix Times that go into the C mapping function come from cls_meta_data_all['Nav']
    # implanted into the record-space of a single CLS file via the map variable, "Navi."
    # If there are data gaps, the interpolated data (pre-version 3) is interpolated across the gaps,
    # resulting in a much longer (many more profiles) array than cls_meta_data_all. This means
    # that interpolated Nav data now residing in cls_meta_data_all will map incorrectly, the error
    # increasing after each data gap. These are large errors. By deleting the interpolated
    # gap-spanning data, cls_meta_data_all will implant each CLS file with correctly
    # mapped interpolated navigation (to match to itself). Simple, right?
    # Test this code out on SEAC4RS data.
    # NOTE: As of version 3, this code will not execute if make_frequency_array
    #       set to 1. Version 3 computes an 'm_array' that is exactly defined
    #       for every CLS record and does not span data gaps.
    if (inst_clk_deltas.max() > inst_clk_rez) and not make_frequency_array:

        print('If you are seeing this message, tell Patrick Selmer to')
        print('investigate why this block of code is necessary.')
        # Maybe this relates to how NEW is populated in the C code.
        # Wish I would have take better notes about this block.
        # pdb.set_trace()
        print(attention_bar + "Time gaps in instrument clock detected. Handling it!" + attention_bar)
        UTC_ovrwrt = np.copy(nav_interp['UTC_Time'][offset_indx:])
        no_delta = True
        for k in range(0, inst_clk_deltas.shape[0]):
            if inst_clk_deltas[k] > inst_clk_rez:
                if no_delta:
                    delta_indx = ExactTime_ui[k - 1] + ExactTime_ncounts[k - 1]
                    no_delta = False
                skip_recs = int(inst_clk_deltas[k] * CLS_hz)
                UTC_ovrwrt[delta_indx:delta_indx + skip_recs] = bad_cls_nav_time_value
                # print(delta_indx,ExactTime_ui[k],k)
                # print('A')
                # pdb.set_trace()
                delta_indx += skip_recs
            elif (not no_delta) and (inst_clk_deltas[k] <= inst_clk_rez):
                delta_indx += ExactTime_ncounts[k]  # int(CLS_hz)
                # print(delta_indx,ExactTime_ui[k],k)
                # print('B')
                # pdb.set_trace()
        UTC_Time_use = np.extract(UTC_ovrwrt != bad_cls_nav_time_value, UTC_ovrwrt)
        try:
            print(cls_meta_data_all['Nav']['UTC_Time'].shape[0])
            print(UTC_Time_use.shape[0])
            cls_meta_data_all['Nav']['UTC_Time'] = UTC_Time_use[:cls_meta_data_all.shape[0]]
        except ValueError:
            num_short = str(cls_meta_data_all.shape[0] - UTC_Time_use.shape[0]).strip()
            print(attention_bar + '!!!!!!!!!! WARNING !!!!!!!!!!\n')
            print('UTC_Time_use is shorter than meta array by ' + num_short + ' elements.')
            print('\n!!!!!!!!!! WARNING !!!!!!!!!!' + attention_bar)
            cls_meta_data_all['Nav']['UTC_Time'][:UTC_Time_use.shape[0]] = UTC_Time_use
            cls_meta_data_all['Nav']['UTC_Time'][UTC_Time_use.shape[0]:] = UTC_Time_use[-1]

    else:

        print('\nYou are lucky. No time gaps detected in the middle of the flight! :-)\n')
        # 5/4/2020 offset_indx removed from assignment below and will need to be checked for lines above
        print(cls_meta_data_all['Nav']['UTC_Time'].shape[0])
        print(nav_interp['UTC_Time'].shape[0])
        cls_meta_data_all['Nav']['UTC_Time'] = nav_interp['UTC_Time'][0:cls_meta_data_all.shape[0]]

    # plt.plot_date(cls_meta_data_all['Nav']['UTC_Time'],cls_meta_data_all['Nav']['RollAngle'],marker='x')
    # plt.plot_date(nav_interp['UTC_Time'],nav_interp['RollAngle'],marker='o')
    # plt.show()
    # pdb.set_trace()

    # NOW, set the array that will be used in processing, cls_nav_data_all equal
    # to the interpolated array that was just created, nav_interp.
    cls_nav_data_all = nav_interp

    # This variable will be the same no matter what the Nav_source
    Nav_hz = CLS_hz

else:

    print("You've entered an invalid Nav_source choice. Halting code.")
    pdb.set_trace()

# Load all the dead time tables

DTT = np.zeros((len(DTT_files), max_counts), dtype=np.float32)
cc = 0
for DTT_file in DTT_files:
    one_DTT = read_in_dead_time_table(dtt_dir + DTT_file)
    DTT[cc, :] = one_DTT[:max_counts]
    cc += 1

print('Starting core processing...')  # <--------------------------------
first_read = True
usable_file_indicies = range(usable_file_range[0], usable_file_range[1])
trans_bin = [0, 0]
trans_total = 0
last_file = False
OL_map_indx = 0  # count through overlap map as required
######## CLS 2 L1A MAP-SAVING CODE ########
true_undoctored_num_CLS_recs = 0      # not a typical L1A variable
last_true_undoctored_num_CLS_recs = 0 # not a typical L1A variable
n_L1A = 0 # counts the number of L1A-resolution records
n_expand_runsum = 0
for f in range(0, nCLS_files):

    if f not in usable_file_indicies:
        print(attention_bar)
        print("File number " + str(f).strip() + " has no usable data!")
        print("Skipping processing of this file.")
        print(attention_bar)
        continue

    # Read in CLS data of one file and initialize a few things
    CLS_file = all_CLS_files[f]
    CLS_file = CLS_file.rstrip()
    FL_flag = 0
    if f == usable_file_range[0]: FL_flag = -1
    if (f == usable_file_range[1] - 1) or (f == nCLS_files - 1):
        FL_flag = 1
        last_file = True
    [interp_UnixT, CLS_UnixT_float64_orig, CLS_data_1file] = create_cls_interp_unixt_array(
        CLS_file, cls_meta_data_all, FL_flag, nav2cls_indx[f, :], FL_trunc)
    CLS_UnixT_float64_orig = CLS_UnixT_float64_orig - nudge + iwg1_cls_t_offset
    interp_UnixT = interp_UnixT - nudge + iwg1_cls_t_offset
    if CLS_data_1file is None:
        print('\n******* Bad file! Skipping file!!! *******')
        print(CLS_file)
        continue
    nr_1file = CLS_data_1file.shape[0]  # # of records in current file
    good_rec_bool = np.ones(CLS_data_1file.shape[0], dtype=bool)  # is 'True' for good records, 'False' for bad
    satur_flag = False  # Will turn True if any bin in current file is potentially saturated.
    ################## MASTER CLS RECORD TRACKER FOR SINGLE FILE ###############################
    if (FL_flag == -1) and (Nav_source != 'iwg1'):  # first file of flight
        true_undoctored_num_CLS_recs = nr_1file + FL_trunc[0] + last_true_undoctored_num_CLS_recs
        rectrack_1file = np.arange(last_true_undoctored_num_CLS_recs, true_undoctored_num_CLS_recs, dtype=np.uint32)
        rectrack_1file = rectrack_1file[FL_trunc[0]:] 
    elif (FL_flag == 1) and (Nav_source != 'iwg1'): # last file of flight
        true_undoctored_num_CLS_recs = nr_1file + FL_trunc[1] + last_true_undoctored_num_CLS_recs
        rectrack_1file = np.arange(last_true_undoctored_num_CLS_recs, true_undoctored_num_CLS_recs, dtype=np.uint32)
        if FL_trunc[1] != 0: rectrack_1file = rectrack_1file[:int(-1*FL_trunc[1])]
    else:              # flight in the middle
        true_undoctored_num_CLS_recs = nr_1file + last_true_undoctored_num_CLS_recs
        rectrack_1file = np.arange(last_true_undoctored_num_CLS_recs, true_undoctored_num_CLS_recs, dtype=np.uint32)
    #print(true_undoctored_num_CLS_recs,last_true_undoctored_num_CLS_recs)
    #print(rectrack_master.min(), rectrack_master.max(), rectrack_master.shape)
    #pdb.set_trace()
    last_true_undoctored_num_CLS_recs = true_undoctored_num_CLS_recs
    ############################################################################################    

    # Create histogram bin boundaries. If first_read, created farther down in first_read block.
    if not first_read: time_bins = np.arange(trans_bin[0], CLS_UnixT_float64_orig[-1] + 3.0 * secs2avg, secs2avg)

    # Deem as "bad," records outside user-specified processing time range
    if ((CLS_data_1file['meta']['Header']['ExactTime'].min() < process_start)
            or (CLS_data_1file['meta']['Header']['ExactTime'].max() > process_end)):
        time_range_mask_low = CLS_data_1file['meta']['Header']['ExactTime'] > process_start
        time_range_mask_high = CLS_data_1file['meta']['Header']['ExactTime'] < process_end
        time_range_mask = time_range_mask_low * time_range_mask_high
        CLS_data_1file = np.extract(time_range_mask, CLS_data_1file)
        rectrack_1file = np.extract(time_range_mask,rectrack_1file)
        CLS_UnixT_float64_orig = np.extract(time_range_mask, CLS_UnixT_float64_orig)
        good_rec_bool = np.ones(CLS_data_1file.shape[0], dtype=bool)  # is 'True' for good records, 'False' for bad
        print(attention_bar)
        print('In file ' + CLS_file + '...')
        print(str(nr_1file - CLS_data_1file.shape[0]).strip() + ' records are out of time range.')
        print(attention_bar)
        nr_1file = CLS_data_1file.shape[0]
        if nr_1file == 0:
            print(attention_bar)
            print("File number " + str(f).strip() + " has no usable data!")
            print("Skipping processing of this file.")
            print(attention_bar)
            continue

    # ---> START process of calling a C function <---
    #      This C function populates an array
    #      mapped to the CLS Unix Time values
    #      with the appropriate interpolated values.
    #      This function also populates an array
    #      that maps the Nav_data indices to
    #      the CLS data. This drastically speeds
    #      up the profile loop in python code.
    # "interp_UnixT" and "CLS_UnixT_float64_orig" should have been computed
    # earlier in the program, before this file loop.
    interp2orig_indicies = np.zeros(nr_1file, dtype=np.uint32)
    CLS_UnixT_float64_new = np.zeros(nr_1file, dtype=np.float64)
    CLS_UnixT_float64_new = np.require(CLS_UnixT_float64_new, float, ['ALIGNED', 'C_CONTIGUOUS'])
    CLS_UnixT_float64_orig = np.require(CLS_UnixT_float64_orig, float, ['ALIGNED', 'C_CONTIGUOUS'])
    interp_UnixT = np.require(interp_UnixT, float, ['ALIGNED', 'C_CONTIGUOUS'])
    Nav_interp_T_float64 = np.require(Nav_interp_T_float64, float, ['ALIGNED', 'C_CONTIGUOUS'])
    interp2orig_indicies = np.require(interp2orig_indicies, int, ['ALIGNED', 'C_CONTIGUOUS'])
    np_clib.map_interp_times_to_orig_frame.restype = None
    np_clib.map_interp_times_to_orig_frame.argtypes = [np.ctypeslib.ndpointer(float, ndim=1, flags='aligned'),
                                                       np.ctypeslib.ndpointer(float, ndim=1, flags='aligned'),
                                                       np.ctypeslib.ndpointer(float, ndim=1, flags='aligned'),
                                                       np.ctypeslib.ndpointer(float, ndim=1, flags='aligned'),
                                                       np.ctypeslib.ndpointer(int, ndim=1, flags='aligned'),
                                                       ctypes.POINTER(ctypes.c_double),
                                                       ctypes.POINTER(np.ctypeslib.c_intp),
                                                       ctypes.POINTER(np.ctypeslib.c_intp),
                                                       ctypes.POINTER(np.ctypeslib.c_intp),
                                                       ctypes.POINTER(np.ctypeslib.c_intp),
                                                       ctypes.POINTER(np.ctypeslib.c_intp),
                                                       ctypes.POINTER(np.ctypeslib.c_intp),
                                                       ctypes.POINTER(np.ctypeslib.c_intp),
                                                       ctypes.POINTER(np.ctypeslib.c_intp),
                                                       ctypes.POINTER(np.ctypeslib.c_intp),
                                                       ctypes.POINTER(np.ctypeslib.c_intp)]
    np_clib.map_interp_times_to_orig_frame(CLS_UnixT_float64_new, CLS_UnixT_float64_orig, interp_UnixT,
                                           Nav_interp_T_float64, interp2orig_indicies, ctypes.c_double(1.0 / Nav_hz),
                                           CLS_UnixT_float64_new.ctypes.strides, CLS_UnixT_float64_new.ctypes.shape,
                                           CLS_UnixT_float64_orig.ctypes.strides, CLS_UnixT_float64_orig.ctypes.shape,
                                           interp_UnixT.ctypes.strides, interp_UnixT.ctypes.shape,
                                           Nav_interp_T_float64.ctypes.strides, Nav_interp_T_float64.ctypes.shape,
                                           interp2orig_indicies.ctypes.strides, interp2orig_indicies.ctypes.shape)
    # ---> END process of calling a C function <---

    # Save/create a few key data parameters...
    if first_read:
        nb = nbins  # # of bins
        # vrZ set in config file. Not in raw data (CLS). Dennis confirmed. [4/6/18]
        nc = CLS_data_1file['meta']['Header']['NumChannels'][0]  # # of channels
        # nshots set in config file. not in raw data (CLS) [4/6/18]
        vrT = 1e-7  # CLS_data_1file['meta']['binwid'][0]
        # Range vector, broadcast to nc x nb for vectorized operations later
        bins_range = np.broadcast_to(np.arange(0, nb * vrZ, vrZ, dtype=np.float32), (nc, nb))
        Y = 0.0 * (np.pi / 180.0)  # Initialize YPR
        P = 0.0 * (np.pi / 180.0)
        R = 0.0 * (np.pi / 180.0)
        i = 0  # counts total records of entire dataset (flight)
        g = 0
        # Create array to help compute heading if using GPS data...
        m_atan = np.zeros((3, 4), dtype=np.uint16)
        m_atan[2, 3] = 0
        m_atan[2, 1] = 1
        m_atan[0, 1] = 2
        m_atan[0, 3] = 3
        # Now create a fixed altitude frame onto which to fit the counts
        ffrme = set_fixed_bin_alt_frame(ff_bot_alt, ff_top_alt, vrZ_ff, nb, pointing_dir)
        ffrme = np.require(ffrme, float, ['ALIGNED', 'C_CONTIGUOUS'])
        nb_ff = ffrme.shape[0]  # number of bins in the fixed frame
        # Create histogram bin boundaries.
        time_bins = np.arange(CLS_UnixT_float64_orig[0], CLS_UnixT_float64_orig[-1] + 3.0 * secs2avg, secs2avg)
        # Make the overlap_vector into an array where the 1st dim. lines up sequentially
        # with the channels
        overlaps_chan_seq = np.ones((nc, nb), dtype=np.float32)
        # Code has option of whether to process using single or multiple overlaps
        if multi_OLs:
            # Load all possible overlap tables into array [nOLs, nchans, nbins]
            # First user-defined defined configuration for this flight
            overlap_dict, overlap_map = load_overlap_configuration(overlaps_config_file)
            # Now load all the OL tables into an array
            all_overlaps = load_all_possible_overlaps(overlap_dict, nb, nc)
        else:  # initializations say only using 1 overlap table
            # Load the overlap table
            overlap_vector = read_in_overlap_table(olap_dir + overlap_file)
            # Overlap is one long array, where first nbins (833) are
            # 355, next nbins are 532, and next (final) nbins are 1064.
            overlaps = overlap_vector.reshape((nwl, nb))
            for chan in range(0, nc):
                overlaps_chan_seq[chan, :] = overlaps[wl_map[chan], :]
        # Open the hdf5 file and create the data sets
        hdf5_fname = L1_dir + 'NRB_' + proj_name + '_' + flt_date + '_' + Nav_source + '.hdf5'
        hdf5_file = h5py.File(hdf5_fname, 'a')
        # Open the lat/lon GEOS met data sampling CSV file
        GEOS_samp_fname = L1_dir + 'lon_lat_UTC_' + proj_name + '_' + flt_date + '_' + Nav_source + '.txt'
        GEOS_fobj = open(GEOS_samp_fname, 'wt')
        GEOS_fobj.write('lon,lat,time\n')  # file header
        try:
            PGain_dset = hdf5_file.create_dataset("PGain", (len(PGain),), np.float32)
            nav_dset = hdf5_file.create_dataset("nav", (1,), maxshape=(None,), dtype=nav_save_struct)
            laserspot_dset = hdf5_file.create_dataset("laserspot", (1, 2), maxshape=(None, 2), dtype=np.float32)
            ONA_dset = hdf5_file.create_dataset("ONA", (1,), maxshape=(None,), dtype=np.float32)
            bin_alt_dset = hdf5_file.create_dataset("bin_alt_array", ffrme.shape, ffrme.dtype)
            nb_ff_dset = hdf5_file.create_dataset("num_ff_bins", (1,), np.uint32)
            num_recs_dset = hdf5_file.create_dataset("num_recs", (1,), np.uint32)
            DEM_nadir_dset = hdf5_file.create_dataset("DEM_nadir", (1,), maxshape=(None,), dtype=np.float32)
            DEM_nadir_surftype_dset = hdf5_file.create_dataset("DEM_nadir_surftype", (1,), maxshape=(None,),dtype=np.int8)
            DEM_laserspot_dset = hdf5_file.create_dataset("DEM_laserspot", (1,), maxshape=(None,), dtype=np.float32)
            DEM_laserspot_surftype_dset = hdf5_file.create_dataset("DEM_laserspot_surftype", (1,), maxshape=(None,), dtype=np.int8)
            EM_dset = hdf5_file.create_dataset("EM", (1, nwl), maxshape=(None, nwl), dtype=np.float32)
            NRB_dset = hdf5_file.create_dataset("nrb", (nc, 1, nb_ff), maxshape=(nc, None, nb_ff), dtype=np.float32)
            bg_dset = hdf5_file.create_dataset("bg", (nc, 1), maxshape=(nc, None), dtype=np.float32)
            saturate_ht_dset = hdf5_file.create_dataset("saturate_ht", (nc, 1), maxshape=(nc, None), dtype=np.float32)
        except RuntimeError:
            print("HDF5 file for this dataset already exists, overwriting old file...")
            hdf5_file.close()  # close, delete, reopen...
            delete_file(hdf5_fname)
            hdf5_file = h5py.File(hdf5_fname, 'a')
            PGain_dset = hdf5_file.create_dataset("PGain", (len(PGain),), np.float32)
            nav_dset = hdf5_file.create_dataset("nav", (1,), maxshape=(None,), dtype=nav_save_struct)
            laserspot_dset = hdf5_file.create_dataset("laserspot", (1, 2), maxshape=(None, 2), dtype=np.float32)
            ONA_dset = hdf5_file.create_dataset("ONA", (1,), maxshape=(None,), dtype=np.float32)
            bin_alt_dset = hdf5_file.create_dataset("bin_alt_array", ffrme.shape, ffrme.dtype)
            nb_ff_dset = hdf5_file.create_dataset("num_ff_bins", (1,), np.uint32)
            num_recs_dset = hdf5_file.create_dataset("num_recs", (1,), np.uint32)
            DEM_nadir_dset = hdf5_file.create_dataset("DEM_nadir", (1,), maxshape=(None,), dtype=np.float32)
            DEM_nadir_surftype_dset = hdf5_file.create_dataset("DEM_nadir_surftype", (1,), maxshape=(None,), dtype=np.int8)
            DEM_laserspot_dset = hdf5_file.create_dataset("DEM_laserspot", (1,), maxshape=(None,), dtype=np.float32)
            DEM_laserspot_surftype_dset = hdf5_file.create_dataset("DEM_laserspot_surftype", (1,), maxshape=(None,), dtype=np.int8)
            EM_dset = hdf5_file.create_dataset("EM", (1, nwl), maxshape=(None, nwl), dtype=np.float32)
            NRB_dset = hdf5_file.create_dataset("nrb", (nc, 1, nb_ff), maxshape=(nc, None, nb_ff), dtype=np.float32)
            bg_dset = hdf5_file.create_dataset("bg", (nc, 1), maxshape=(nc, None), dtype=np.float32)
            saturate_ht_dset = hdf5_file.create_dataset("saturate_ht", (nc, 1), maxshape=(nc, None), dtype=np.float32)
        except:
            print("An unanticipated error occurred while trying to create the HDF5 data sets. Stopping execution.")
            pdb.set_trace()
        # Some data sets only need to be written once, right here...
        PGain_dset[:] = np.asarray(PGain)  # convert list to array
        bin_alt_dset[:] = ffrme
        nb_ff_dset[:] = nb_ff

    counts_ff = np.zeros((nc, nr_1file, nb_ff), dtype=np.float32)
    NRB = np.empty_like(counts_ff)
    bg_save = np.zeros((nc, nr_1file), dtype=np.float32)

    Nav_save = np.zeros(nr_1file, dtype=nav_save_struct)  # NOTE THE STRUCTURE TYPE!
    ONA_save = np.zeros(nr_1file, dtype=np.float32)
    laserspot = np.zeros((nr_1file, 2), dtype=np.float32)  # [lat x lon]
    DEM_nadir = np.zeros(nr_1file, dtype=np.float32)
    DEM_nadir_surftype = np.zeros(nr_1file, dtype=np.int8) - 9
    DEM_laserspot = np.zeros(nr_1file, dtype=np.float32)
    DEM_laserspot_surftype = np.zeros(nr_1file, dtype=np.int8) - 9
    saturate_ht = np.zeros((nc, nr_1file), dtype=np.float32) - 99999.9

    # Right here, right now, test to see if any bin in current file is potentially saturated

    # Quick, broad test...
    satur_locs = []
    satur_locs_indx = 0
    if CLS_data_1file['counts'].max() > min(saturation_values): satur_flag = True
    if satur_flag:
        print('\n' + all_CLS_files[f] + ' has saturated bin values in it!\n')
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
            Nav_save[i1f]['heading'] = np.arctan2(Nav_match['east'], Nav_match['north']) * (180.0 / np.pi)
            Nav_save[i1f]['lon'] = Nav_match['lon']
            Nav_save[i1f]['lat'] = Nav_match['lat']
            Nav_save[i1f]['GPS_alt'] = Nav_match['alt']
            time_str = gps_UTC_interp[interp2orig_indicies[i1f]].strftime("%Y-%m-%dT%H:%M:%S.%f")
            Nav_save[i1f]['UTC'] = np.asarray(list(time_str.encode('utf8')), dtype=np.uint8)
            current_Nav_UTC = Nav_match['UTC_Time']  # Nav_match fields dependant on data set
            Y = 0.0  # Nav_match['yaw'] * (pi/180.0) Woah! This isn't what I think yaw should be [3/22/18]
            P = Nav_match['pitch'] * (np.pi / 180.0)
            R = Nav_match['roll'] * (np.pi / 180.0)
            # Write 'lat,lon,UTC' out to file
            GEOS_out_str = '{0:.4f},{1:.4f},{2}\n'.format(Nav_match['lon'], Nav_match['lat'], time_str[:19])
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
            Nav_save[i1f]['UTC'] = np.asarray(list(time_str.encode('utf8')), dtype=np.uint8)
            current_Nav_UTC = Nav_match['UTC_Time']  # Nav_match fields dependant on data set
            # NOTE on next line: Invalids are -999.9, so max should choose valid if exists
            Nav_save[i1f]['GPS_alt'] = max(Nav_match['GPS_alt_msl'], Nav_match['GPS_alt'])
            Y = Nav_match['drift'] * (np.pi / 180.0)
            P = Nav_match['pitch'] * (np.pi / 180.0)
            R = Nav_match['roll'] * (np.pi / 180.0)
            # Write 'lat,lon,UTC' out to file
            GEOS_out_str = '{0:.4f},{1:.4f},{2}\n'.format(Nav_match['lon'], Nav_match['lat'], time_str[:19])
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
            # time to write out to the HDF5 file. It will essentially be a
            # byte array for IDL.
            time_str = Nav_match['UTC'].strftime("%Y-%m-%dT%H:%M:%S.%f")
            Nav_save[i1f]['UTC_Time'] = np.asarray(list(time_str.encode('utf8')), dtype=np.uint8)
            # Nav_match fields dependant on data set    Switch to UTC for IWG1 and UTC_Time for cls
            current_Nav_UTC = Nav_match['UTC']
            Nav_save[i1f]['GPS_alt'] = Nav_match['GPS_alt']
            Y = 0.0  # Nav_match['drift'] * (np.pi/180.0)
            P = Nav_match['pitch'] * (np.pi / 180.0)
            R = Nav_match['roll'] * (np.pi / 180.0)
            # Write 'lat,lon,UTC' out to file
            GEOS_out_str = '{0:.4f},{1:.4f},{2}\n'.format(Nav_match['lon'], Nav_match['lat'], time_str[:19])
            GEOS_fobj.write(GEOS_out_str)

        elif Nav_source == 'cls':

            # Find matching Nav data point. Use index map to match.
            # NOTE: Nav_save is an IWG1-style structure (see immutable_dat_structs.py)
            Nav_match = cls_nav_data_all[interp2orig_indicies[i1f]]
            Nav_save[i1f]['roll'] = (Nav_match['RollAngle'] + cpl_roll_offset)
            Nav_save[i1f]['pitch'] = ((Nav_match['PitchAngle'] + cpl_pitch_offset) + ((Nav_match['RollAngle'] + cpl_roll_offset) * cpl_pitch_roll_factor))
            Nav_save[i1f]['drift'] = (Nav_match['TrackAngleTrue'] - Nav_match['TrueHeading'])
            # Need to make sure drift is in the correct units after Track/Heading Difference Calculated
            if Nav_save[i1f]['drift'] > 180.0:
                Nav_save[i1f]['drift'] = (Nav_save[i1f]['drift'] - 360.0)

            Nav_save[i1f]['heading'] = Nav_match['TrueHeading']
            Nav_save[i1f]['lon'] = Nav_match['GPS_Longitude']
            Nav_save[i1f]['lat'] = Nav_match['GPS_Latitude']
            # The following line is necessary to get a clean, interpretable,
            # time to write out to the HDF5 file. It will essentially be a
            # byte array for IDL.
            time_str = Nav_match['UTC_Time'].strftime("%Y-%m-%dT%H:%M:%S.%f")
            Nav_save[i1f]['UTC'] = np.asarray(list(time_str.encode('utf8')), dtype=np.uint8)
            current_Nav_UTC = Nav_match['UTC_Time']  # Nav_match fields dependant on data set
            Nav_save[i1f]['GPS_alt'] = Nav_match['GPS_Altitude']
            Y = Nav_save[i1f]['drift'] * (np.pi / 180.0)  # Yaw now incorporated!! 3/30/2020
            P = ((Nav_match['PitchAngle'] + cpl_pitch_offset) + ((Nav_match['RollAngle'] + cpl_roll_offset) * cpl_pitch_roll_factor)) * (np.pi / 180.0)  # Offsets added 3/30/2020
            R = (Nav_match['RollAngle'] + cpl_roll_offset) * (np.pi / 180.0)
            # Write 'lat,lon,UTC' out to file
            GEOS_out_str = '{0:.4f},{1:.4f},{2}\n'.format(Nav_match['GPS_Longitude'], Nav_match['GPS_Latitude'],
                                                          time_str[:19])
            GEOS_fobj.write(GEOS_out_str)

        # Proceed no farther in processing if altitude is not high enough
        if (Nav_save[i1f]['GPS_alt'] < alt_cutoff) | (not np.isfinite(Nav_save[i1f]['GPS_alt'])):
            print('Bad altitude: ', Nav_save[i1f]['GPS_alt'])
            # print('Will not process until good alt is read.\n')
            good_rec_bool[i1f] = False
            i = i + 1
            i1f = i1f + 1
            continue

        # Convert to float for subsequent manipulation within this loop
        counts_float32 = CLS_data_1file['counts'][i1f].astype(np.float32)

        # Apply dead time correction
        cc = 0  # count the channels
        for DTT_file in DTT_files:
            try:
                counts_float32[cc, :] = DTT[cc, CLS_data_1file['counts'][i1f, cc, :]]
            except IndexError:
                # Probably a noise spike
                spike_locs = np.where(counts_float32[cc, :] > DTT.shape[1])
                counts_float32[cc, spike_locs[0][0]] = 0.0
            cc += 1

        # Calculate the off nadir angle
        # For CAMAL, positive scan angles are to the right
        # For CPL, assume phi0 and phi0_azi are zero
        [ONA, xy_ang] = calculate_off_nadir_angle(0.0, 0.0, Y, P, R, xy_ang=True)  # accepted way as of [1/31/18]
        # ONA = calculate_off_nadir_angle(0.0,0.0,Y,P,R+(-1.0*phi0)) # accepted way as of [1/25/18]
        ONA_save[i1f] = ONA
        if ONA > ONA_cutoff:
            good_rec_bool[i1f] = False
            i = i + 1
            i1f = i1f + 1
            continue

        # Get the DEM value at aircraft nadir coordinate
        DEM_ret = DEM_lib.get_DEM_point(DEM_file, Nav_save[i1f]['lat'], Nav_save[i1f]['lon'])
        DEM_nadir[i1f] = DEM_ret.Elev
        if DEM_ret.Elev == -9999.0: DEM_nadir[i1f] = 0.0  # ocean surface
        DEM_nadir_surftype[i1f] = DEM_ret.LOCI

        # Calculate the laser spot coordinates
        [latLS, lonLS] = laser_spot_aircraft(Nav_save[i1f]['GPS_alt'], Nav_save[i1f]['lat'] * (np.pi / 180.0),
                                             Nav_save[i1f]['lon'] * (np.pi / 180.0), xy_ang, DEM_nadir[i1f], ONA,
                                             Nav_save[i1f]['heading'] * (np.pi / 180.0))
        laserspot[i1f, 0] = latLS
        laserspot[i1f, 1] = lonLS

        # Get the DEM value at the laser spot coordinate
        DEM_ret = DEM_lib.get_DEM_point(DEM_file, laserspot[i1f, 0] * (180.0 / np.pi),
                                        laserspot[i1f, 1] * (180.0 / np.pi))
        DEM_laserspot[i1f] = DEM_ret.Elev
        if DEM_ret.Elev == -9999.0: DEM_laserspot[i1f] = 0.0  # ocean surface
        DEM_laserspot_surftype[i1f] = DEM_ret.LOCI

        # Compute the altitudes of the raw data bins
        bin_alts = compute_raw_bin_altitudes(nb, pointing_dir, Nav_save[i1f]['GPS_alt'],
                                             vrZ, ONA)

        # Populate saturate_ht with bin alts if any bins in current profile are saturated
        # NOTE: The order of cond's matters in if blocks in this code paragraph
        if satur_flag and (i1f in satur_locs[:, 0]):
            while (satur_locs_indx < satur_locs.shape[0]) and (satur_locs[satur_locs_indx, 0] == i1f):
                indx = satur_locs[satur_locs_indx, :]  # Will always have 3 dims in same order
                if CLS_data_1file['counts'][indx[0], indx[1], indx[2]] > saturation_values[indx[1]]:
                    saturate_ht[indx[1], i1f] = bin_alts[indx[2]]
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
            bg_1D = np.mean(counts_float32[:, bg_loc1:bg_loc2], axis=1)
            bg_save[:, i1f] = bg_1D
            bg = np.broadcast_to(bg_1D, (nb, nc)).transpose()
        except:
            print(attention_bar)
            print("!!!!!!! WARNING !!!!!!!")
            print("No data within defined background region! Using background of zero.")
            print(Nav_save[i1f]['UTC'], ' ONA: ', ONA * (180.0 / np.pi))
            print('End Bin Altitudes ,', bin_alts[800:833])
            print('Reported Altitude ,', Nav_save[i1f]['GPS_alt'])
            print("!!!!!!! WARNING !!!!!!!")
            print(attention_bar)
            bg = np.zeros((nc, nb))

        range_cor_af_counts_float32 = ((counts_float32 - bg) * bins_range ** 2)

        # Apply overlap correction here. I would have applied it at the same time as the
        # deadtime, except CPL's overlaps can be funky in far ranges, which interferes with
        # background sub.
        if multi_OLs:  # Only if using multiple overlap tables
            if (current_Nav_UTC > overlap_map[OL_map_indx][0][1]) & (OL_map_indx < len(overlap_map) - 1):
                print('UTC = ' + current_Nav_UTC.strftime("%Y-%m-%dT%H:%M:%S") +
                      '. Now using the following overlap tables...')
                print('355 nm: ' + overlap_dict[overlap_map[OL_map_indx][1][0]])
                print('532 nm: ' + overlap_dict[overlap_map[OL_map_indx][1][1]])
                print('1064 nm: ' + overlap_dict[overlap_map[OL_map_indx][1][2]] + '\n')
                OL_map_indx += 1
                # Just do the channels w/o a loop for now [11/20/19]
                overlaps_chan_seq[0, :] = all_overlaps[overlap_map[OL_map_indx][1][0], 0, :]
                overlaps_chan_seq[1, :] = all_overlaps[overlap_map[OL_map_indx][1][1], 1, :]
                overlaps_chan_seq[2, :] = all_overlaps[overlap_map[OL_map_indx][1][2], 2, :]
                overlaps_chan_seq[3, :] = all_overlaps[overlap_map[OL_map_indx][1][2], 3, :]
            else:
                overlaps_chan_seq[0, :] = all_overlaps[overlap_map[OL_map_indx][1][0], 0, :]
                overlaps_chan_seq[1, :] = all_overlaps[overlap_map[OL_map_indx][1][1], 1, :]
                overlaps_chan_seq[2, :] = all_overlaps[overlap_map[OL_map_indx][1][2], 2, :]
                overlaps_chan_seq[3, :] = all_overlaps[overlap_map[OL_map_indx][1][2], 3, :]
        range_cor_af_counts_float32 = range_cor_af_counts_float32 / overlaps_chan_seq

        # Put the bins in the fixed altitude frame
        # Updated to use split_bins_onto_fixed_frame, as of v2 [4/22/20]
        counts_ff[:, i1f, :] = split_bins_onto_fixed_frame(ffrme, bin_alts, range_cor_af_counts_float32)

        i = i + 1  # increment record counter by 1
        i1f = i1f + 1  # increment record counter for current file by 1

    if i == 0: continue

    # Apply polarization gain ratio to 1064 perpendicular channel
    print(attention_bar)
    print("Applying pgain factor to hard-coded channel of index 3")
    print("This should represent the 1064 nm channel")
    print(attention_bar)
    counts_ff[3, :, :] = counts_ff[3, :, :] * PGain[1]

    # Patch over bad energy monitor values.
    # Different from classic "AM.pro" code ~line 954.
    BadEs = np.argwhere(CLS_data_1file['meta']['Engineering']['LaserEnergyMonitors'] <= 0)
    if BadEs.size > 0:
        medianEs = np.median(CLS_data_1file['meta']['Engineering']['LaserEnergyMonitors'], axis=0)
        for BadE in BadEs:
            # Just replace for all wavelengths, even if only one is bad.
            CLS_data_1file['meta']['Engineering']['LaserEnergyMonitors'][BadE[0], :] = medianEs

    # Compute NRB & populate saturate_ht array
    EMs = convert_raw_energy_monitor_values(CLS_data_1file['meta']['Engineering']['LaserEnergyMonitors'], nwl, 'CPL',e_flg)
    ff_bg_st_bin = np.argwhere(ffrme <= bg_st_alt)[0][0]
    ff_bg_ed_bin = np.argwhere(ffrme >= bg_ed_alt)[-1][0]
    for chan in range(0, nc):
        EMsubmit = EMs[:, wl_map[chan]]
        # Patch over bad EM readings
        EMsubmit = patch_outliers(EMsubmit, Estds)
        NRB[chan, :, :] = correct_raw_counts(counts_ff[chan, :, :], EMsubmit, None,
                            i1f, nb_ff, ff_bg_st_bin, ff_bg_ed_bin, 'NRB_no_range')

    # Delete bad records from output arrays
    tot_good_recs = good_rec_bool.sum()
    if tot_good_recs < nr_1file:
        CLS_UnixT_float64_new = CLS_UnixT_float64_new[good_rec_bool]
        rectrack_1file = rectrack_1file[good_rec_bool]
        Nav_save = Nav_save[good_rec_bool]
        laserspot = laserspot[good_rec_bool, :]
        ONA_save = ONA_save[good_rec_bool]
        DEM_nadir = DEM_nadir[good_rec_bool]
        DEM_nadir_surftype = DEM_nadir_surftype[good_rec_bool]
        DEM_laserspot = DEM_laserspot[good_rec_bool]
        DEM_laserspot_surftype = DEM_laserspot_surftype[good_rec_bool]
        EMs = EMs[good_rec_bool, :]
        NRB = NRB[:, good_rec_bool, :]
        bg_save = bg_save[:, good_rec_bool]
        saturate_ht = saturate_ht[:, good_rec_bool]
        print('\nXXXXXXXXXXXXXXXXXXXXXXX')
        deleted_recs = nr_1file - tot_good_recs
        print(str(deleted_recs).strip() + ' records deleted!')
        print('From file: ' + CLS_file)
        print('\nXXXXXXXXXXXXXXXXXXXXXXX')
        i = i - (nr_1file - tot_good_recs)
        nr_1file = tot_good_recs

    # Average the data
    # Remainder data at end of each file is carried-over to next file
    #
    # Here's the logic behind this averaging section:
    # The data are averaged one file at a time, with any leftover records
    # carried-over until the next file is processed. The carried-over records
    # are then averaged with the first 'x' records of the next file, where 'x'
    # represents the number of records required to meet the time-average criterion.
    # The logic is a bit tricky when processing the first and last files, and these
    # details are dealt with using the "first_read" and "cutbegin" variables. Also,
    # n_expand (how much to expand data sets) and expanded_length variables are
    # computed differently for the edge versus middle files. The reason averaging
    # is done one file at a time is to limit usage of RAM on the computer.
    if secs2avg > 0.0:  # if < 0, don't average
        
        ######## CLS 2 L1A MAP-SAVING CODE ########
        L1Arec_nums_1file = np.empty_like(rectrack_1file)
        ONA_1file = np.zeros(rectrack_1file.shape[0], dtype=ONA_save.dtype)
        Plane_Alt_1file = np.zeros(rectrack_1file.shape[0], dtype=np.float32)
        E_1file = np.zeros((rectrack_1file.shape[0], 3), dtype=EMs.dtype)        

        bin_numbers = np.digitize(CLS_UnixT_float64_new, time_bins)
        u, ui, ncounts = np.unique(bin_numbers, return_index=True, return_counts=True)
        Nav_save_avg = np.zeros(u.shape[0], dtype=Nav_save.dtype)
        laserspot_avg = np.zeros((u.shape[0], 2), dtype=laserspot.dtype)
        ONA_save_avg = np.zeros(u.shape[0], dtype=ONA_save.dtype)
        DEM_nadir_avg = np.zeros(u.shape[0], dtype=DEM_nadir.dtype)
        DEM_nadir_surftype_avg = np.zeros(u.shape[0], dtype=DEM_nadir_surftype.dtype)
        DEM_laserspot_avg = np.zeros(u.shape[0], dtype=DEM_laserspot.dtype)
        DEM_laserspot_surftype_avg = np.zeros(u.shape[0], dtype=DEM_laserspot_surftype.dtype)
        EMs_avg = np.zeros((u.shape[0], nwl), dtype=EMs.dtype)
        NRB_avg = np.zeros((nc, u.shape[0], nb_ff), dtype=NRB.dtype)
        bg_save_avg = np.zeros((nc, u.shape[0]), dtype=bg_save.dtype)
        saturate_ht_max = np.zeros((nc, u.shape[0]), dtype=saturate_ht.dtype)
        rr = 0  # raw record number
        perfectly_aligned = False 

        ei = ui.shape[0] - 1           # End Index (ei) set only in 
        if last_file: ei = ui.shape[0] # these two lines.
        # Three different types of cases in following chain of IF blocks...
        # 1 - It's the first file.
        # 2 - You're in the middle of processing and a time bin spans across files
        # 3 - You're in the middle of processing and a time bin doesn't span files
        # If it's the last file, 2 or 3 could be executed with additional
        # controls on the end index and trans_bin (via "last_file") outside
        # these IF blocks.
        if first_read:
            si = 0
        # Gotta do an elif cuz cur file might not have any vals within last bin of prev file
        elif time_bins[bin_numbers[ui[0]] - 1] == trans_bin[0]:
            si = 1  # start index
            trans_total = float(trans_ncounts + ncounts[0])
            for field in Nav_save_avg.dtype.names:
                if field == 'UTC': continue
                Nav_save_avg[field][0] = ((Nav_save_sum[field] + Nav_save[field][rr:rr + ncounts[0]].sum()) / trans_total)
            Nav_save_avg['UTC'][0,:] = Nav_UTC_carryover
            laserspot_avg[0, :] = (laserspot_sum + laserspot[rr:rr + ncounts[0], :].sum(axis=0)) / trans_total
            ONA_save_avg[0] = (ONA_save_sum + ONA_save[rr:rr + ncounts[0]].sum()) / trans_total
            DEM_nadir_avg[0] = (DEM_nadir_sum + DEM_nadir[rr:rr + ncounts[0]].sum()) / trans_total
            DEM_nadir_surftype_avg[0] = stats.mode(np.concatenate((DEM_nadir_surftype_carryover,DEM_nadir_surftype[rr:rr + ncounts[0]])))[0][0]
            DEM_laserspot_avg[0] = (DEM_laserspot_sum + DEM_laserspot[rr:rr + ncounts[0]].sum()) / trans_total
            DEM_laserspot_surftype_avg[0] = stats.mode(np.concatenate((DEM_laserspot_surftype_carryover,DEM_laserspot_surftype[rr:rr + ncounts[0]])))[0][0]
            EMs_avg[0, :] = (EMs_sum + EMs[rr:rr + ncounts[0], :].sum(axis=0)) / trans_total
            NRB_avg[:, 0, :] = (NRB_sum + NRB[:, rr:rr + ncounts[0], :].sum(axis=1)) / trans_total
            bg_save_avg[:, 0] = (bg_save_sum + bg_save[:, rr:rr + ncounts[0]].sum(axis=1)) / trans_total
            saturate_ht_max[:, 0] = np.asarray((saturate_ht_carryover.max(axis=1),saturate_ht[:, rr:rr + ncounts[0]].max(axis=1))).max(axis=0)
            ######## CLS 2 L1A MAP-SAVING CODE ########
            # This block of code is where the transfer bin is carried over,
            # so set the current L1A record number equal to the last value
            # from the previous file iteration. In the transfer save block,
            # n_L1A is not incremented, so the following line is good.
            L1Arec_nums_1file[rr:rr+ncounts[0]] = n_L1A 
            # Now increment n_L1A because you're done with the transfer record
            n_L1A += 1
            ONA_1file[rr:rr+ncounts[0]] = ONA_save[rr:rr+ncounts[0]]
            Plane_Alt_1file[rr:rr+ncounts[0]] = Nav_save['GPS_alt'][rr:rr+ncounts[0]]
            E_1file[rr:rr+ncounts[0],:] = EMs[rr:rr+ncounts[0],:]            
            print('trans_ncounts = ', trans_ncounts)
            rr += ncounts[0]
        else:
            # So the 'time bin' of the first record of this file does not
            # match the 'time bin' of the last record of the previous file
            # (aka "trans_bin"). BUT we still have the saved trans_bin data
            # from the previous file. We need to put it into our current avg
            # arrays. First, we'll expand the size of our averaged arrays by one.
            # Then we'll feed the data from the trans_bin (remember, previous file)
            # into the first element (element zero). We'll then set si (start index)
            # to one and code outside this block should proceed normally.
            Nav_save_avg = np.zeros(u.shape[0]+1, dtype=Nav_save.dtype)
            laserspot_avg = np.zeros((u.shape[0]+1, 2), dtype=laserspot.dtype)
            ONA_save_avg = np.zeros(u.shape[0]+1, dtype=ONA_save.dtype)
            DEM_nadir_avg = np.zeros(u.shape[0]+1, dtype=DEM_nadir.dtype)
            DEM_nadir_surftype_avg = np.zeros(u.shape[0]+1, dtype=DEM_nadir_surftype.dtype)
            DEM_laserspot_avg = np.zeros(u.shape[0]+1, dtype=DEM_laserspot.dtype)
            DEM_laserspot_surftype_avg = np.zeros(u.shape[0]+1, dtype=DEM_laserspot_surftype.dtype)
            EMs_avg = np.zeros((u.shape[0]+1, nwl), dtype=EMs.dtype)
            NRB_avg = np.zeros((nc, u.shape[0]+1, nb_ff), dtype=NRB.dtype)
            bg_save_avg = np.zeros((nc, u.shape[0]+1), dtype=bg_save.dtype)
            saturate_ht_max = np.zeros((nc, u.shape[0]+1), dtype=saturate_ht.dtype)
            # Now insert trans_bin data into element zero. NOTE rr remains zero.
            trans_total = float(trans_ncounts)
            for field in Nav_save_avg.dtype.names:
                if field == 'UTC': continue
                Nav_save_avg[field][0] = Nav_save_sum[field] / trans_total
            Nav_save_avg['UTC'][0,:] = Nav_UTC_carryover
            laserspot_avg[0, :] = laserspot_sum / trans_total
            ONA_save_avg[0] = ONA_save_sum / trans_total
            DEM_nadir_avg[0] = DEM_nadir_sum / trans_total
            DEM_nadir_surftype_avg[0] = stats.mode(DEM_nadir_surftype_carryover)[0][0]
            DEM_laserspot_avg[0] = DEM_laserspot_sum / trans_total
            DEM_laserspot_surftype_avg[0] = stats.mode(DEM_laserspot_surftype_carryover)[0][0]
            EMs_avg[0, :] = EMs_sum / trans_total
            NRB_avg[:, 0, :] = NRB_sum / trans_total
            bg_save_avg[:, 0] = bg_save_sum / trans_total
            saturate_ht_max[:, 0] = saturate_ht_carryover.max(axis=1)
            # Now do these critical arrays as well
            u = np.concatenate(([-99],u))        # now elem 0 corr. with trans_bin
            ui = np.concatenate(([-99],ui))
            ncounts = np.concatenate(([trans_ncounts],ncounts))
            ######## CLS 2 L1A MAP-SAVING CODE ########
            # This block of code is where the transfer bin is carried over,
            # so set the current L1A record number equal the last value
            # from the previous file iteration. In the transfer save block,
            # n_L1A is not incremented. But THIS is a SPECIAL CASE. Just above
            # alignment with the file boundary essentially added a new record.
            # Therefore n_L1A needs to be incremented both before AND after
            # assignment into array!!! [1/27/21]
            n_L1A += 1
            #L1Arec_nums_1file[rr:rr+ncounts[0]] = n_L1A    
            # Now increment n_L1A because you're done with the transfer record           
            #ONA_1file[rr:rr+ncounts[0]] = ONA_save[rr:rr+ncounts[0]]
            #Plane_Alt_1file[rr:rr+ncounts[0]] = Nav_save['GPS_alt'][rr:rr+ncounts[0]]
            #E_1file[rr:rr+ncounts[0],:] = EMs[rr:rr+ncounts[0],:]            
            print('trans_ncounts = ', trans_ncounts)                       
            si = 1
            ei = Nav_save_avg.shape[0] - 1 # You expanded the array in this block
            if last_file: ei = Nav_save_avg.shape[0] # loop thru til end if last file
            perfectly_aligned = True
            print(attention_bar)
            print("I guess the time_bins lined up perfectly with edge of previous file")
            print("because there are no values in the previous file's last time bin.")
            print(trans_bin)
            print(attention_bar)

        for tb in range(si, ei):
            for field in Nav_save_avg.dtype.names:
                if field == 'UTC': continue
                Nav_save_avg[field][tb] = np.mean(Nav_save[field][rr:rr + ncounts[tb]])
            Nav_save_avg['UTC'][tb] = Nav_save['UTC'][rr + int(ncounts[tb] / 2)]  # midpoint
            laserspot_avg[tb, :] = np.mean(laserspot[rr:rr + ncounts[tb], :], axis=0)
            ONA_save_avg[tb] = np.mean(ONA_save[rr:rr + ncounts[tb]])
            DEM_nadir_avg[tb] = np.mean(DEM_nadir[rr:rr + ncounts[tb]])
            DEM_nadir_surftype_avg[tb] = stats.mode(DEM_nadir_surftype[rr:rr + ncounts[tb]])[0][0]
            DEM_laserspot_avg[tb] = np.mean(DEM_laserspot[rr:rr + ncounts[tb]])
            DEM_laserspot_surftype_avg[tb] = stats.mode(DEM_laserspot_surftype[rr:rr + ncounts[tb]])[0][0]
            EMs_avg[tb, :] = np.mean(EMs[rr:rr + ncounts[tb], :], axis=0)
            NRB_avg[:, tb, :] = np.mean(NRB[:, rr:rr + ncounts[tb], :], axis=1)
            bg_save_avg[:, tb] = np.mean(bg_save[:, rr:rr + ncounts[tb]], axis=1)
            saturate_ht_max[:, tb] = np.max(saturate_ht[:, rr:rr + ncounts[tb]], axis=1)
            ######## CLS 2 L1A MAP-SAVING CODE ########
            L1Arec_nums_1file[rr:rr+ncounts[tb]] = n_L1A
            n_L1A += 1    
            ONA_1file[rr:rr+ncounts[tb]] = ONA_save[rr:rr+ncounts[tb]]
            Plane_Alt_1file[rr:rr+ncounts[tb]] = Nav_save['GPS_alt'][rr:rr+ncounts[tb]]
            E_1file[rr:rr+ncounts[tb],:] = EMs[rr:rr+ncounts[tb],:]   
            # print('iter: ',u[tb], CLS_UnixT_float64_new[rr])
            rr += ncounts[tb]

        # Save the sum and ncounts of last populated time bin to string averaging
        # together between multiple files.
        if not last_file:
            Nav_save_sum = np.zeros(1, dtype=Nav_save.dtype)
            for field in Nav_save_avg.dtype.names:
                if field == 'UTC': continue
                Nav_save_sum[field][0] = Nav_save[field][rr:rr + ncounts[-1]].sum()
            Nav_UTC_carryover = Nav_save['UTC'][rr + ncounts[-1] - 1]
            laserspot_sum = laserspot[rr:rr + ncounts[-1], :].sum(axis=0)
            ONA_save_sum = ONA_save[rr:rr + ncounts[-1]].sum()
            DEM_nadir_sum = DEM_nadir[rr:rr + ncounts[-1]].sum()
            DEM_nadir_surftype_carryover = DEM_nadir_surftype[rr:rr + ncounts[-1]]
            DEM_laserspot_sum = DEM_laserspot[rr:rr + ncounts[-1]].sum()
            DEM_laserspot_surftype_carryover = DEM_laserspot_surftype[rr:rr + ncounts[-1]]
            EMs_sum = EMs[rr:rr + ncounts[-1], :].sum(axis=0)
            NRB_sum = NRB[:, rr:rr + ncounts[-1], :].sum(axis=1)
            bg_save_sum = bg_save[:, rr:rr + ncounts[-1]].sum(axis=1)
            saturate_ht_carryover = saturate_ht[:, rr:rr + ncounts[-1]]
            ######## CLS 2 L1A MAP-SAVING CODE ########
            L1Arec_nums_1file[rr:rr+ncounts[-1]] = n_L1A   
            print(L1Arec_nums_1file[rr:rr+ncounts[-1]])
            print('n_L1A = ',n_L1A)
            ONA_1file[rr:rr+ncounts[-1]] = ONA_save[rr:rr+ncounts[-1]]
            Plane_Alt_1file[rr:rr+ncounts[-1]] = Nav_save['GPS_alt'][rr:rr+ncounts[-1]]
            E_1file[rr:rr+ncounts[-1],:] = EMs[rr:rr+ncounts[-1],:]
            # following line defines the "transfer" bin; trans_bin
            trans_bin = time_bins[bin_numbers[ui[-1]] - 1: bin_numbers[ui[-1]] + 1]
            trans_ncounts = ncounts[-1]

        # [1/28/21]
        if perfectly_aligned:
            L1Arec_nums_for_current_L1Aarrs = np.concatenate((L1Arec_nums[-1*ncounts[0]:],L1Arec_nums_1file))
            mo = ncounts[0]# mask offset. will apply to mask that filters out min_avg_profs profs.
        else:
            L1Arec_nums_for_current_L1Aarrs = np.copy(L1Arec_nums_1file)
            mo = 0

        # Eliminate those avg'd profiles that contain less than min_avg_profs raw
        # profiles right here, exclusively in this single paragraph.
        big_enough_mask = ncounts >= min_avg_profs
        big_enough_mask[0] = True
        big_enough_mask[-1] = True
        if (not first_read) and (not perfectly_aligned) and (trans_total < min_avg_profs): big_enough_mask[0] = False
        if big_enough_mask.sum() < ncounts.shape[0]:	
            Nav_save_avg = Nav_save_avg[big_enough_mask]
            laserspot_avg = laserspot_avg[big_enough_mask, :]
            ONA_save_avg = ONA_save_avg[big_enough_mask]
            DEM_nadir_avg = DEM_nadir_avg[big_enough_mask]
            DEM_nadir_surftype_avg = DEM_nadir_surftype_avg[big_enough_mask]
            DEM_laserspot_avg = DEM_laserspot_avg[big_enough_mask]
            DEM_laserspot_surftype_avg = DEM_laserspot_surftype_avg[big_enough_mask]
            EMs_avg = EMs_avg[big_enough_mask, :]
            NRB_avg = NRB_avg[:, big_enough_mask, :]
            bg_save_avg = bg_save_avg[:, big_enough_mask]
            saturate_ht_mask = saturate_ht_max[:, big_enough_mask]
            ncounts = ncounts[big_enough_mask]
            ui = ui[big_enough_mask]
            u = u[big_enough_mask]
            ######## CLS 2 L1A MAP-SAVING CODE ########
            # 1) Identify the L1A rec #'s that have been deleted.
            # 2) Find where the L1Arec_nums_1file array == those rec #'s
            # 3) Remove those elements from L1Arec_nums_1file & rectrack_master_1file
            unique_L1Arec_nums = np.unique(L1Arec_nums_for_current_L1Aarrs)
            too_small_mask = np.invert(big_enough_mask)
            deleted_L1A_indicies = unique_L1Arec_nums[too_small_mask]
            for di in deleted_L1A_indicies:
                mask = L1Arec_nums_for_current_L1Aarrs != di # locations not equal to deleted index
                L1Arec_nums_for_current_L1Aarrs = L1Arec_nums_for_current_L1Aarrs[mask]
                L1Arec_nums_1file = L1Arec_nums_1file[mask[mo:]] # running count of # of L1A records
                rectrack_1file = rectrack_1file[mask[mo:]]
                ONA_1file = ONA_1file[mask[mo:]]
                Plane_Alt_1file = Plane_Alt_1file[mask[mo:]]
                E_1file = E_1file[mask[mo:],:]                
                # Some L1A recs were deleted, so decrement appropriately.
                # pushback_mask shouldn't need an offset.
                pushback_mask = L1Arec_nums_1file > di
                L1Arec_nums_1file[pushback_mask] = L1Arec_nums_1file[pushback_mask] - 1
                # Now, you also have to apply diff-length pushback_mask to array that matches L1A-rez variables
                pushback_mask = L1Arec_nums_for_current_L1Aarrs > di
                L1Arec_nums_for_current_L1Aarrs[pushback_mask] = L1Arec_nums_for_current_L1Aarrs[pushback_mask] - 1
                n_L1A -= 1
            print("\nAvg'd profiles eliminated due to min_avg_profs constraint.")
            print(np.argwhere(big_enough_mask == False))
            print(big_enough_mask.shape, " reduced to ", u.shape, "\n")

        # Expand data set sizes to accommodate next input CLS file
        cutbegin = 0
        if (first_read) and (ncounts[0] < min_avg_profs): 
            cutbegin = 1
            L1Arec_nums_1file = L1Arec_nums_1file[ncounts[0]:]
            L1Arec_nums_1file = L1Arec_nums_1file - 1
            n_L1A -= 1
            rectrack_1file = rectrack_1file[ncounts[0]:]
            ONA_1file = ONA_1file[ncounts[0]:]
            Plane_Alt_1file = Plane_Alt_1file[ncounts[0]:]
            E_1file = E_1file[ncounts[0]:,:]              
        n_expand = u.shape[0] - 1 - cutbegin
        if ( (last_file) and (ncounts[-1] > min_avg_profs) ): 
            n_expand = u.shape[0]
        elif ((last_file) and (ncounts[-1] <= min_avg_profs)):
            # added 7/30/2020 cuz PELIcoe 22oct19 flight did not have min_avg_profs in last rec.
            L1Arec_nums_1file = L1Arec_nums_1file[:-1*ncounts[-1]]
            rectrack_1file = rectrack_1file[:-1*ncounts[-1]]
            ONA_1file = ONA_1file[:-1*ncounts[-1]]
            Plane_Alt_1file = Plane_Alt_1file[:-1*ncounts[-1]]
            E_1file = E_1file[:-1*ncounts[-1],:]             
            n_L1A -= 1
        expanded_length = nav_dset.shape[0] + n_expand
        # Now, if it's the first_read, nav_dset has an initialized length of 1; therefore,
        # if you use the expanded_length in the previous line the first_read, you'll
        # get a data set size that is too long by 1. The following line fixes this.
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
        nav_dset[expanded_length - n_expand:expanded_length] = Nav_save_avg[cutbegin:n_expand + cutbegin]
        laserspot_dset[expanded_length - n_expand:expanded_length, :] = laserspot_avg[cutbegin:n_expand + cutbegin, :]
        ONA_dset[expanded_length - n_expand:expanded_length] = ONA_save_avg[cutbegin:n_expand + cutbegin]
        DEM_nadir_dset[expanded_length - n_expand:expanded_length] = DEM_nadir_avg[cutbegin:n_expand + cutbegin]
        DEM_nadir_surftype_dset[expanded_length - n_expand:expanded_length] = DEM_nadir_surftype_avg[cutbegin:n_expand + cutbegin]
        DEM_laserspot_dset[expanded_length - n_expand:expanded_length] = DEM_laserspot_avg[cutbegin:n_expand + cutbegin]
        DEM_laserspot_surftype_dset[expanded_length - n_expand:expanded_length] = DEM_laserspot_surftype_avg[cutbegin:n_expand + cutbegin]
        EM_dset[expanded_length - n_expand:expanded_length, :] = EMs_avg[cutbegin:n_expand + cutbegin, :]
        NRB_dset[:, expanded_length - n_expand:expanded_length, :] = NRB_avg[:, cutbegin:n_expand + cutbegin, :]
        bg_dset[:, expanded_length - n_expand:expanded_length] = bg_save_avg[:, cutbegin:n_expand + cutbegin]
        saturate_ht_dset[:, expanded_length - n_expand:expanded_length] = saturate_ht_max[:,cutbegin:n_expand + cutbegin]
        nrecs = expanded_length
        n_expand_runsum = n_expand + n_expand_runsum
        print('Nav_save_avg shape: ',Nav_save_avg.shape)
        print('cutbegin: ',cutbegin)
        print('n_expand: ',n_expand)
        print('n_expand_runsum: ',n_expand_runsum)
        print('n_L1A: ',n_L1A)
        print(L1Arec_nums_1file.shape)
        print(rectrack_1file.shape)

    else:  # No averaging

        # Expand data set sizes to accommodatee next input CLS file
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
        nav_dset[i - nr_1file:i] = Nav_save
        laserspot_dset[i - nr_1file:i, :] = laserspot
        ONA_dset[i - nr_1file:i] = ONA_save
        DEM_nadir_dset[i - nr_1file:i] = DEM_nadir
        DEM_nadir_surftype_dset[i - nr_1file:i] = DEM_nadir_surftype
        DEM_laserspot_dset[i - nr_1file:i] = DEM_laserspot
        DEM_laserspot_surftype_dset[i - nr_1file:i] = DEM_laserspot_surftype
        EM_dset[i - nr_1file:i, :] = EMs
        NRB_dset[:, i - nr_1file:i, :] = NRB
        bg_dset[:, i - nr_1file:i] = bg_save
        saturate_ht_dset[:, i - nr_1file:i] = saturate_ht
        nrecs = i
        print('** WARNING ** WARNING ** WARNING ** WARNING ** WARNING **')
        print('No code written to deal with no averaging')   
        
    # At this point, no more raw CLS records will be deleted. So now, concatenate
    # rectrack_master with rectrack_1file.
    if 'rectrack_master' in locals(): # Will not be in locals() first iter.
        rectrack_master = np.concatenate((rectrack_master,rectrack_1file))
        L1Arec_nums = np.concatenate((L1Arec_nums, L1Arec_nums_1file))
        L1A_Time = np.concatenate((L1A_Time, Nav_save_avg['UTC'][cutbegin:n_expand+cutbegin]), axis=0)
        Plane_Alt = np.concatenate((Plane_Alt, Plane_Alt_1file))
        ONA_4map = np.concatenate((ONA_4map, ONA_1file))
        E_4map = np.concatenate((E_4map, E_1file),axis=0)
    else:
        rectrack_master = np.copy(rectrack_1file)
        L1Arec_nums = np.copy(L1Arec_nums_1file)
        L1A_Time = np.copy(Nav_save_avg['UTC'][cutbegin:n_expand+cutbegin])
        Plane_Alt = np.copy(Plane_Alt_1file)
        ONA_4map = np.copy(ONA_1file)
        E_4map = np.copy(E_1file)
    
    first_read = False
        
    print('\n**********************************')
    print('\nDone with file: '+CLS_file+'\n')
    print('**********************************\n')

      
print('Main L1A execution finished at: ',DT.datetime.now())

######## CLS 2 L1A MAP-SAVING CODE ########
# Now write out the CLS-to-L1A index maps to a text file
print('Shape of rectrack_master: ',rectrack_master.shape)
print('Shape of L1Arec_nums: ',L1Arec_nums.shape)
print('Shape of L1A_Time: ',L1A_Time.shape)
index_map_file = raw_dir + 'CLS2L1A_map_'+proj_name+'_'+flt_date+'_'+Nav_source+'.csv'

with open(index_map_file, 'w') as map_f_obj:
    oord = np.arange(0, rectrack_master.shape[0])
    for CLSi, L1Ai, x in zip(rectrack_master, L1Arec_nums, oord):
        str_time = str( L1A_Time[L1Ai,:].view('S26') )[3:29]
        str_alt = '{0:11.4f}'.format( Plane_Alt[x] )
        str_ONA = '{0:13.9f}'.format( ONA_4map[x] )
        str_E0 = '{0:21.6f}'.format( E_4map[x,0] )
        str_E1 = '{0:21.6f}'.format( E_4map[x,1] )
        str_E2 = '{0:21.6f}'.format( E_4map[x,2] )
        outlist = [str(CLSi),str(L1Ai),str_time,str_alt,str_ONA,str_E0,str_E1,str_E2]
        string2write = ','.join(outlist) + '\n'
        map_f_obj.write(string2write)


num_recs_dset[:] = nrecs
print('nrecs: ',nrecs)
hdf5_file.close()
GEOS_fobj.close()
print("NRB has been written to the HDF5 file:" + hdf5_fname)

# The final message

print('Total raw profiles processed: ' + str(i))
print("cpl_l1a.py has finished normally.")

print('Nudge Used: ' + str(nudge))
