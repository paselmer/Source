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
# Libraries I did create
from initializations import *
from read_routines import *
from lidar import *
from time_conversions import *
from mutable_data_structs import define_MSC_structure
import matplotlib.dates as mdates
import ctypes


# Define functions <----------------------------------------------------

def compute_MCS_IWG1_time_offset_IWG1():
    """ Compute the offset between the UTC timestamp in the MCS data and
        the UTC timestamp in the IWG1 data.
    """
    
    # Computed manually using 12/7/17 flight dataset
    # Subtract this number from the MCS GPS time
    time_offset_IWG1 = 33.0
    time_offset_IWG1 = DT.timedelta(seconds=33)
    
    return time_offset_IWG1
    
def compute_Unix_time_offset():
    """ This is the time offset between the 2 internal clocks on the
        CAMAL instrument; the Unix Time and the GPS time.
    """
    
    # Computed manually using 12/7/17 flight dataset
    # Subtract this number from the Unix time
    time_offset_UnixT = 1645.798
    time_offset_UnixT = DT.timedelta(seconds=1645.798)
    
    return time_offset_UnixT  
        
    
# Start main execution here <-------------------------------------------
print('Starting main L1A execution at: ',DT.datetime.now())

# Create and load file lists (MCS,GPS,nav)...

MCS_file_list = 'processing_file_list.txt'
search_str = 'data*'
create_a_file_list(MCS_file_list,search_str)
GPS_file_list = 'gps_file_list.txt'
search_str = 'gps*'
create_a_file_list(GPS_file_list,search_str)
nav_file_list = 'nav_file_list.txt'
search_str = 'nav*'
create_a_file_list(nav_file_list,search_str)
# Load the shared library
np_clib = np.ctypeslib.load_library(clib,clib_path) 
# Compute the time offset between MCS and IWG1
time_offset_IWG1 = compute_MCS_IWG1_time_offset_IWG1()
# Compute the time offset between GPS and UnixT within MCS data stream
time_offset_UnixT = compute_Unix_time_offset()

with open(MCS_file_list) as MCS_list_fobj:
    all_MCS_files = MCS_list_fobj.readlines()
nMCS_files = len(all_MCS_files)

# Load nav data for entire flight. Select from 1 of several possible
# nav data streams. All these read-in all flight data at once. This
# decision was made to simplify the code and perhaps boost processing speed.

if Nav_source == 'gps':
    
    with open(GPS_file_list) as GPS_list_fobj:
        all_GPS_files = GPS_list_fobj.readlines()
    nGPS_files = len(all_GPS_files)

elif Nav_source == 'nav':
    
    # NOTE: "Nav" variables will apply generically to whichever of the
    #       several Nav_source's is selected. "nav" variables will 
    
    with open(nav_file_list) as nav_list_fobj:
        all_nav_files = nav_list_fobj.readlines()
    nnav_files = len(all_nav_files)
    
    est_nav_recs = est_nav_recs_1file*nnav_files
    nav_data_all = np.zeros(est_nav_recs,dtype=nav_struct)
    j = 0
    for i in range(0,nnav_files):
        nav_file = all_nav_files[i]
        nav_file = nav_file.strip()
        nav_data_1file = read_in_nav_data(nav_file, est_nav_recs_1file, DT.timedelta(seconds=18000))
        try:
            nav_data_all[j:j+nav_data_1file.shape[0]] = nav_data_1file
        except ValueError:
            pdb.set_trace()
        j += nav_data_1file.shape[0]
    nav_data_all = nav_data_all[0:j] # trim the fat
    
    # Interpolate data records to match CAMAL's data rate.
    # 10 Hz as of 2/2/18. Data rate set initializations.
    
    del_t = np.zeros(nav_data_all.shape[0],dtype=DT.datetime)
    del_t[1:] = nav_data_all['UTC'][1:] - nav_data_all['UTC'][:-1]
    del_t[0] = del_t[1]*0.0
    
    # Compute # of elapsed seconds since first record (float)
    
    del_secs = np.zeros(nav_data_all.shape[0],dtype=np.float64)
    for k in range(0,nav_data_all.shape[0]): del_secs[k] = del_t[k].total_seconds()
    
    elapsed_secs = np.cumsum(del_secs)
    tot_elapsed = np.max(elapsed_secs)
    m = nav_hz / MCS_hz
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
        Nav_interp_T_float64[k] = (nav_interp['UnixT'][k] - UnixT_epoch).total_seconds()
        nav_interp['UTC'][k] = nav_data_all['UTC'][0] + DT.timedelta(seconds=ix[k])
    
    # NOW, set the array that will be used in processing, nav_data_all equal
    # to the interpolated array that was just created, nav_interp.
    
    nav_data_all = nav_interp
    
elif Nav_source == 'iwg1':

    IWG1_file = raw_dir + "IWG1.08Dec2017-0031.txt"
    IWG1_data = read_in_IWG1_data(IWG1_file,est_IWG1_recs)
    
    # Interpolate data records to match CAMAL's data rate.
    # 10 Hz as of 2/2/18. Data rate set initializations.
    
    del_t = np.zeros(IWG1_data.shape[0],dtype=DT.datetime)
    del_t[1:] = IWG1_data['UTC'][1:] - IWG1_data['UTC'][:-1]
    del_t[0] = del_t[1]*0.0
    
    # Compute # of elapsed seconds since first record (float)
    
    del_secs = np.zeros(IWG1_data.shape[0],dtype=np.float64)
    for k in range(0,IWG1_data.shape[0]): del_secs[k] = del_t[k].total_seconds()
    
    elapsed_secs = np.cumsum(del_secs)
    tot_elapsed = np.max(elapsed_secs)
    m = IWG1_hz / MCS_hz
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
        IWG1_interp['UTC'][k] = IWG1_data['UTC'][0] + DT.timedelta(seconds=ix[k])
        Nav_interp_T_float64[k] = (IWG1_interp['UTC'][k] - UnixT_epoch + 
            time_offset_UnixT).total_seconds()

    # NOW, set the array that will be used in processing, IWG1_data equal
    # to the interpolated array that was just created, IWG1_interp.
    
    IWG1_data = IWG1_interp
    
else:
    
    print("You've entered an invalid Nav_source choice. Halting code.")
    pdb.set_trace()


print('Starting core processing...') # <--------------------------------
first_read = False
for f in range(0,nMCS_files):
    
    MCS_file = all_MCS_files[f]
    MCS_file = MCS_file.rstrip()
    MCS_data_1file = read_in_raw_data(MCS_file)
    if MCS_data_1file is None:
        print('\n******* Bad file! Skipping file!!! *******')
        print(MCS_file)
        continue  
    nr_1file = MCS_data_1file.shape[0]                    # # of records in current file
    
    ## Interpolate MCS times to data rate
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
        np.ctypeslib.ndpointer(int, ndim=1,flags='aligned'), 
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp) ]
    np_clib.map_interp_times_to_orig_frame(MCS_UnixT_float64_new, MCS_UnixT_float64_orig, interp_UnixT,
        Nav_interp_T_float64, interp2orig_indicies,
        MCS_UnixT_float64_new.ctypes.strides, MCS_UnixT_float64_new.ctypes.shape, 
        MCS_UnixT_float64_orig.ctypes.strides, MCS_UnixT_float64_orig.ctypes.shape,
        interp_UnixT.ctypes.strides, interp_UnixT.ctypes.shape,
        Nav_interp_T_float64.ctypes.strides, Nav_interp_T_float64.ctypes.shape,  
        interp2orig_indicies.ctypes.strides, interp2orig_indicies.ctypes.shape) 
    # ---> END process of calling a C function <---                          

    # Save a few key data parameters...
    if first_read is False:
        nb =  MCS_data_1file['meta']['nbins'][0]              # # of bins
        vrZ = (MCS_data_1file['meta']['binwid'][0] * c) / 2.0 # vert. rez in m
        nc = MCS_data_1file['meta']['nchans'][0]              # # of channels
        nshots = MCS_data_1file['meta']['nshots'][0]      
        vrT = MCS_data_1file['meta']['binwid'][0]
        tot_est_recs = int(rep_rate/nshots)*file_len_secs*nMCS_files
        nr = tot_est_recs
        phi0_azi = 0.0
        Y = 0.0 * (pi/180.0)                                  # Initialize YPR
        P = 0.0 * (pi/180.0)
        R = 0.0 * (pi/180.0)

        # Now create a fixed altitude frame onto which to fit the counts
        ffrme = set_fixed_bin_alt_frame(ff_bot_alt,ff_top_alt,vrZ_ff,nb,pointing_dir)
        nb_ff = ffrme.shape[0] # number of bins in the fixed frame
        counts_ff = np.zeros((nc,nr,nb_ff),dtype=np.float32)
        NRB = np.empty_like(counts_ff)
        mult = np.zeros(nb_ff,dtype=np.float32)
        
        i = 0
        g = 0
        Nav_save = np.zeros(nr,dtype=nav_save_struct) #NOTE THE STRUCTURE TYPE!
        ONA_save = np.zeros(nr,dtype=np.float32)
        xy_ang_save = np.zeros(nr,dtype=np.float32)
        laserspot = np.zeros((nr,2),dtype=np.float32) # [lat x lon]
        MCS_struct = define_MSC_structure(nc,nb) 
        MCS_metadata_all = np.zeros(nr,dtype = MCS_meta_struct)
        
    i1f = 0 
    j = 0
    
    while i1f < nr_1file:
    
        if Nav_source == 'gps':
            
            # Find mathing GPS data point
            print("Please write some code. Try to follow IWG1/nav examples")
            pdb.set_trace()
            #Y = GPS_data_all['yaw'][Gind[0]] * (pi/180.0)           #*****----> DON'T DELETE THESE LINES, WILL USE EVENTUALLY
            #P = GPS_data_all['pitch'][Gind[0]] * (pi/180.0)
            #R = GPS_data_all['roll'][Gind[0]] * (pi/180.0)            
        
        elif Nav_source == 'nav':
  
            # Use index map to match
            Nav_match = nav_data_all[interp2orig_indicies[i1f]]                                                     
            Nav_save[i]['roll'] = Nav_match['roll']
            Nav_save[i]['pitch'] = Nav_match['pitch']
            Nav_save[i]['drift'] = Nav_match['drift']
            Nav_save[i]['lon'] = Nav_match['lon']
            Nav_save[i]['lat'] = Nav_match['lat']
            Nav_save[i]['UTC'] = Nav_match['UTC'] 
            Nav_save[i]['GPS_alt'] = Nav_match['GPS_alt']
            Y = 0.0 #Nav_match['drift'] * (pi/180.0) 
            P = Nav_match['pitch'] * (pi/180.0)
            R = Nav_match['roll'] * (pi/180.0)  
            
        elif Nav_source == 'iwg1':
            
            # Find matching IWG1 data point. Use index map to match.
            # NOTE: Nav_save is an IWG1-style structure (see immutable_dat_structs.py)
            Nav_match = IWG1_data[interp2orig_indicies[i1f]]
            Nav_save[i]['roll'] = Nav_match['roll']
            Nav_save[i]['pitch'] = Nav_match['pitch']
            Nav_save[i]['drift'] = Nav_match['drift']
            Nav_save[i]['heading'] = Nav_match['heading']
            Nav_save[i]['lon'] = Nav_match['lon']
            Nav_save[i]['lat'] = Nav_match['lat']
            # The following line is necessary to get a clean, interpretable,
            # time to write out to the HDF5 file. It will essetially be a
            # byte array for IDL.
            Nav_save[i]['UTC'] = np.asarray(list(Nav_match['UTC'].strftime("%Y-%m-%dT%H:%M:%S.%f").encode('utf8')),dtype=np.uint8)
            Nav_save[i]['GPS_alt'] = Nav_match['GPS_alt']
            Y = 0.0 #Nav_match['drift'] * (pi/180.0) 
            P = Nav_match['pitch'] * (pi/180.0)
            R = Nav_match['roll'] * (pi/180.0)             
    
        # Proceed no farther in processing if altitude is not high enough
        if Nav_match['GPS_alt'] < alt_cutoff:
            print(Nav_match['GPS_alt'])
            i1f = i1f + 1
            continue
    
        # Calculate the off nadir angle
        # For CAMAL, positive scan angles are to the right
        phi0 = (MCS_data_1file['meta']['scan_pos'][i1f] + angle_offset) * (pi/180.0)
        if phi0 >= 0.0: 
            phi0_azi = 0.0
        else:
            phi0_azi = pi
        [ONA, xy_ang] = calculate_off_nadir_angle(abs(phi0),phi0_azi,Y,P,R,xy_ang=True)   # accepted way as of [1/31/18]
        #ONA = calculate_off_nadir_angle(0.0,0.0,Y,P,R+(-1.0*phi0)) # accepted way as of [1/25/18]
        ONA_save[i] = ONA
        xy_ang_save[i] = xy_ang
        
        # Calculate the laser spot coordinates
        [latLS,lonLS] = laser_spot_aircraft(Nav_match['GPS_alt'],Nav_match['lat']*(np.pi/180.0),
            Nav_match['lon']*(np.pi/180.0),xy_ang,100.0,ONA,Nav_match['heading']*(np.pi/180.0))
        laserspot[i,0] = latLS
        laserspot[i,1] = lonLS
    
        # Compute the altitudes of the raw data bins
        bin_alts = compute_raw_bin_altitudes(nb, pointing_dir,MCS_data_1file['meta']['GpsAlt'][i1f],
                   vrZ, ONA)

        length_bin_alts = bin_alts.shape[0]
        if length_bin_alts > nb:
            print('Length of bin_alts vectors is > than nbins. Stopping code...')
            pdb.set_trace()
    
        # Put the bins in the fixed altitude frame
        #counts_ff[:,i,:] = put_bins_onto_fixed_frame_C(np_clib,ffrme,
        #    bin_alts,MCS_data_1file['counts'][i1f],nc,mult) 
        counts_ff[:,i,:] = put_bins_onto_fixed_frame(ffrme,
            bin_alts,MCS_data_1file['counts'][i1f],nc)  
            
        # Store all meta data from raw files in array   
        MCS_metadata_all[i] = MCS_data_1file['meta'][i1f]
    
        i = i + 1    # increment record counter by 1
        i1f = i1f + 1  # increment record counter for current file by 1
        
    
    if i == 0: continue
    
    # Compute NRB
    EMs = convert_raw_energy_monitor_values(MCS_data_1file['meta']['EM'],nwl)
    for chan in range(0,nc):
        EMsubmit = EMs[:,wl_map[chan]]     
        NRB[chan,i-i1f:i,:] = correct_raw_counts(counts_ff[chan,i-i1f:i,:],EMsubmit,np.arange(0,nb_ff*vrZ_ff,vrZ_ff),
            i1f,nb_ff,ff_bg_st_bin,ff_bg_ed_bin,'NRB')      
        
    first_read = True
    print('\n**********************************')
    print('\nDone with file: '+MCS_file+'\n')
    print('**********************************\n')
    
    
# Trim the excess of arrays, due to over-estimation of size

MCS_metadata_all = MCS_metadata_all[:i]
laserspot = laserspot[:i,:]
counts_ff = counts_ff[:,:i,:]
NRB = NRB[:,:i,:]
Nav_save = Nav_save[:i]
      
print('Main L1A execution finished at: ',DT.datetime.now())

# Write NRB data out to HDF5 file

print("Now writing data out to HDF5 file.")
import h5py

# Create the file    
hdf5_fname = out_dir+'NRB_'+proj_name+'_'+flt_date+'.hdf5'
hdf5_file = h5py.File(hdf5_fname, 'w')
# Create the datasets and write out data
meta_dset = hdf5_file.create_dataset("meta", (i,), MCS_meta_struct)
meta_dset[:] = MCS_metadata_all
PGain = 0.0
PGain_dset = hdf5_file.create_dataset("PGain", (1,), np.float32)
PGain_dset[:] = PGain
nav_dset = hdf5_file.create_dataset("nav", (i,), nav_save_struct)
nav_dset[:] = Nav_save
laserspot_dset = hdf5_file.create_dataset("laserspot", laserspot.shape, laserspot.dtype)
laserspot_dset[:] = laserspot
bin_alt_dset = hdf5_file.create_dataset("bin_alt_array", ffrme.shape, ffrme.dtype)
bin_alt_dset[:] = ffrme
nb_ff_dset = hdf5_file.create_dataset("num_ff_bins", (1,), np.uint32)
nb_ff_dset[:] = nb_ff
num_recs_dset = hdf5_file.create_dataset("num_recs", (1,), np.uint32)
num_recs_dset[:] = i
NRB_dset = hdf5_file.create_dataset("nrb", (nc,i,nb_ff), NRB.dtype) 
NRB_dset[:] = NRB
hdf5_file.close()


#pdb.set_trace()
#tit = '532 nm NRB'
#xlimits = [0,i]
#ylimits = [1100,200]
#samp_chan = NRB[2,:,:] + NRB[3,:,:]
#curtain_plot(samp_chan.transpose(), nb_ff, vrZ_ff, ffrme, 0, 3e9, hori_cap, pointing_dir,
    #figW, figL, CPpad, 'records', 'altitude(m)', tit, 'alt',[ylimits[0],ylimits[1]], 'recs',
    #[xlimits[0],xlimits[1]], scale_alt_OofM, 1, out_dir)
#pdb.set_trace()

## Pick channel(s)
#samp_chan = counts_ff[2,:i,:] + counts_ff[3,:i,:] # 532, in fixed frame

## Background subtract the rebinned counts and then plot
#wl = 1 # numpy index 1 is 532 nm wl
#opt = 'bg_sub'
#EMs = convert_raw_energy_monitor_values(MCS_metadata_all['EM'][:i],nwl)
#EMsubmit = EMs[:,wl] 
#samp_chan = correct_raw_counts(samp_chan,EMsubmit,np.arange(0,nb_ff*vrZ_ff,vrZ_ff),
                               #i,nb_ff,975,1045,opt)

                               
## Code used to make masked scan angle plots ----------------------------
## ----------------------------------------------------------------------                               

## Save samp_chan at raw rez, before any averaging.
#samp_chanX = np.empty_like(samp_chan)
#samp_chanX[:,:] = samp_chan                               
#print('The shape of samp_chan before averaging for total image = ',samp_chan.shape)                                  
##samp_chan = average_lidar_data(samp_chan[15:,:],10,0)                                           
#print('The shape after averaging = ',samp_chan.shape)
#tit = 'CAMAL 532 nm background subtracted counts, all look angles'
#xlimits = [0,i]
##xlimits = [mdates.date2num(weeksecondstoutc(MCS_data_1file['meta']['GpsWeek'][0]*1.0,
           ##MCS_data_1file['meta']['Gps_msec'][0]*1e-3,leapsecs)),
           ##mdates.date2num(weeksecondstoutc(MCS_data_1file['meta']['GpsWeek'][i-1]*1.0,
           ##MCS_data_1file['meta']['Gps_msec'][i-1]*1e-3,leapsecs))]        
#ylimits = [1300,800]
#curtain_plot(samp_chan.transpose(), nb_ff, vrZ_ff, ffrme, 0, 35.0, hori_cap, pointing_dir,
                      #figW, figL, CPpad, 'records', 'altitude(m)', tit, 'alt',  
                      #[ylimits[0],ylimits[1]], 'recs', [xlimits[0],xlimits[1]], scale_alt_OofM, 1, out_dir)
#pdb.set_trace()                      
#make_custom_plot(samp_chan.transpose(), nb_ff, vrZ_ff, ffrme, 0, 60.0, hori_cap, pointing_dir,
                      #figW, figL, CPpad, 'records', 'altitude(m)', tit, 'alt',  
                      #[ylimits[0],ylimits[1]], 'recs', [xlimits[0],xlimits[1]], scale_alt_OofM, 1, out_dir, str(f).strip()+'.png',
                      #np.arange(xlimits[0],xlimits[1]),MCS_metadata_all[xlimits[0]:xlimits[1]]['scan_pos']+angle_offset,
                      #np.arange(xlimits[0],xlimits[1]),MCS_metadata_all[xlimits[0]:xlimits[1]]['scan_state'],
                      #np.arange(xlimits[0],xlimits[1]),Nav_save[xlimits[0]:xlimits[1]]['roll'])
#pdb.set_trace()                      
### Apply a mask for each look angle. Plot the curtain w/gaps.
#angle_arr = MCS_data_1file['meta']['scan_pos'] + angle_offset                               
#[hist_bin_width, hist_bin_centers] = determine_look_angles(
                                     #angle_arr,
                                     #scan_pos_uplim, scan_pos_lowlim, 
                                     #scan_pos_bw )     

#pdb.set_trace()
