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
# Libraries I did create
from initializations import *
from read_routines import *
from lidar import *
from time_conversions import *


# Define functions <----------------------------------------------------



    
# Start main execution here <-------------------------------------------

# Create and load file lists (MCS,GPS)...

MCS_file_list = 'processing_file_list.txt'
search_str = 'data*'
create_a_file_list_in_Windows(MCS_file_list,search_str)
GPS_file_list = 'gps_file_list.txt'
search_str = 'gps*'
create_a_file_list_in_Windows(GPS_file_list,search_str)
# Load the shared library
np_clib = np.ctypeslib.load_library(clib,clib_path) 

with open(MCS_file_list) as MCS_list_fobj:
    all_MCS_files = MCS_list_fobj.readlines()
nMCS_files = len(all_MCS_files)

with open(GPS_file_list) as GPS_list_fobj:
    all_GPS_files = GPS_list_fobj.readlines()
nGPS_files = len(all_GPS_files)


# Load IWG1 data for entire flight

IWG1_file = raw_dir + "IWG1.08Dec2017-0031.txt"
IWG1_data = read_in_IWG1_data(IWG1_file,est_IWG1_recs)


# Finally read in MCS data from first file...

MCS_file = all_MCS_files[25]
MCS_file = MCS_file.rstrip()
MCS_data_1file = read_in_raw_data(MCS_file)
GPS_file = all_GPS_files[40]
GPS_file = GPS_file.rstrip()
GPS_data_1file = read_in_gps_data(GPS_file)

# Save a few key data parameters...

nr = MCS_data_1file.shape[0]                          # # of records
nb =  MCS_data_1file['meta']['nbins'][0]              # # of bins
vrZ = (MCS_data_1file['meta']['binwid'][0] * c) / 2.0 # vert. rez in m
nc = MCS_data_1file['meta']['nchans'][0]              # # of channels
nshots = MCS_data_1file['meta']['nshots'][0]      
vrT = MCS_data_1file['meta']['binwid'][0]
phi0_azi = 0.0
Y = 0.0 * (pi/180.0)                                  # Initialize YPR
P = 0.0 * (pi/180.0)
R = 0.0 * (pi/180.0)

## Create a plot of the raw BGsub data.
 
#EMs = convert_raw_energy_monitor_values(MCS_data_1file['meta']['EM'],nwl)
#EMsubmit = EMs[:,1]
#samp_chan = MCS_data_1file['counts'][:,2,:] + MCS_data_1file['counts'][:,3,:]
#samp_chan = correct_raw_counts(samp_chan,EMsubmit,np.arange(0,nb*vrZ,vrZ),
                               #nr,nb,bg_st_bin,bg_ed_bin,'bg_sub')
#angle_arr = MCS_data_1file['meta']['scan_pos'] + angle_offset                               
#[hist_bin_width, hist_bin_centers] = determine_look_angles(
                                     #angle_arr,
                                     #scan_pos_uplim, scan_pos_lowlim, 
                                     #scan_pos_bw )  

#sti = np.array([57,158,259,360,462])                                     
#edi = np.array([157,258,359,461,561])
#samp_chan_save = np.empty_like(samp_chan)
#samp_chan_save[:,:] = samp_chan                                     
#for i in range(0,5):   
    #samp_chan = np.empty_like(samp_chan_save)
    #samp_chan[:,:] = samp_chan_save     
    #mean_prof = np.mean(samp_chan[sti[i]:edi[i],:],axis=0)
    #pdb.set_trace()
    ##samp_chan = apply_look_angle_mask(samp_chan,angle_arr,hist_bin_centers[i],hist_bin_width,nr,nhori,'gaps') 
    ##samp_chan = samp_chan.transpose()                               
    ##curtain_plot(samp_chan, nb, vrZ, np.arange(0,1000), 0, 100.0, hori_cap, pointing_dir,
                      ##figW, figL, CPpad, 'recs', 'bins', 'test', 'bins',  
                      ##[1000,0], 'recs', [0,nr], scale_alt_OofM, 1, out_dir)
#pdb.set_trace()


# Now create a fixed altitude frame onto which to fit the counts

bot_alt,top_alt = [-15e3,30.005e3]
vrZ_ff = 30.0
ffrme = set_fixed_bin_alt_frame(bot_alt,top_alt,vrZ_ff,nb,pointing_dir)
nb_ff = ffrme.shape[0] # number of bins in the fixed frame
counts_ff = np.zeros((nc,nr,nb_ff))
mult = np.zeros(nb_ff)

print('Starting core processing...')
i = 0
g = 0
match = False
nav_source = 'GPS' # enter 'GPS' or 'IWG1'
IWG1_match = IWG1_data[0]
IWG1_save = np.zeros(nr,dtype=IWG1_struct)
ONA_save = np.zeros(nr,dtype=np.float32)
xzang_save = np.zeros(nr,dtype=np.float32)
while i < nr:
    
    ## Check to see if new GPS is need with this new record [i]
    #if match is True:
        #while GPS_data_1file['Gps_sec'][Gind[0]]*1e3 < MCS_data_1file['meta']['Gps_msec'][i]:
            #Gind[0] += 1
            #if Gind[0] >= GPS_data_1file.shape[0]:
                #match = False
                #break
    
    ## Find GPS record that matches
    #while match is False:
        ## Read in a GPS file
        #GPS_file = all_GPS_files[g]
        #GPS_file = GPS_file.rstrip()
        #GPS_data_1file = read_in_gps_data(GPS_file)
        #if GPS_data_1file is None: # If no data in file, skip
            #print("Skipping file " + all_GPS_files[g])
            #g += 1
            #continue
        #GPS_ind = np.arange(GPS_data_1file.shape[0])+1 # +1 is very important here!
        #condlist = [ GPS_data_1file['Gps_sec']*1e3 > MCS_data_1file['meta']['Gps_msec'][i] ]
        #choicelist = [ GPS_ind ]
        #selected_inds = np.select(condlist,choicelist,default=0)
        #if np.sum(selected_inds) > 0: 
            #Gind = GPS_ind[GPS_ind == selected_inds] - 1
            #match = True
        #g += 1
    ##print(inds)
    
    # Find matching IWG1 data point
    MCS_datetime_rec = weeksecondstoutc(MCS_data_1file['meta']['GpsWeek'][i]*1.0,
                       MCS_data_1file['meta']['Gps_msec'][i]*1e-3,leapsecs)
    if IWG1_match['UTC'] < MCS_datetime_rec:
        IWG1_match_bool = IWG1_data['UTC'] >= MCS_datetime_rec
        IWG1_match = IWG1_data[IWG1_match_bool][0]
    IWG1_save[i]['roll'] = IWG1_match['roll']
    IWG1_save[i]['pitch'] = IWG1_match['pitch']
    IWG1_save[i]['drift'] = IWG1_match['drift']
    
    # Calculate the off nadir angle
    # For CAMAL, positive scan angles are to the right
    phi0 = (MCS_data_1file['meta']['scan_pos'][i] + angle_offset) * (pi/180.0)
    if phi0 >= 0.0: 
        phi0_azi = 0.0
    else:
        phi0_azi = pi
    #Y = GPS_data_1file['yaw'][Gind[0]] * (pi/180.0)           #*****----> DON'T DELETE THESE LINES, WILL USE EVENTUALLY
    #P = GPS_data_1file['pitch'][Gind[0]] * (pi/180.0)
    #R = GPS_data_1file['roll'][Gind[0]] * (pi/180.0)
    Y = 0.0
    P = IWG1_match['pitch'] * (pi/180.0)
    R = IWG1_match['roll'] * (pi/180.0)
    ONA = calculate_off_nadir_angle(-1.0*phi0,phi0_azi,Y,P,R)
    ONA_save[i] = ONA
    xzang_save[i] = calculate_off_nadir_angle(-1.0*phi0,phi0_azi,Y,P,R,xz_ang=True)
    #et = timeit.timeit("calculate_off_nadir_angle(-1.0*phi0,phi0_azi,Y,P,R,xz_ang=True)",
        #setup="from lidar import calculate_off_nadir_angle\nfrom __main__ import phi0, \
        #phi0_azi,Y,P,R",number=1)
    #print("ONA calc time = %g " % et)
    
    # Compute the altitudes of the raw data bins
    bin_alts = compute_raw_bin_altitudes(nb, pointing_dir,MCS_data_1file['meta']['GpsAlt'][i],
                   vrZ, ONA)
    #et = timeit.timeit("compute_raw_bin_altitudes(nb, pointing_dir,MCS_data_1file['meta']['GpsAlt'][i], \
        #vrZ, ONA)",setup="from lidar import compute_raw_bin_altitudes\nfrom __main__ import nb, \
        #pointing_dir,MCS_data_1file,i,vrZ,ONA",number=1)
    #print("Raw bin calc time = %g" % et)
    length_bin_alts = bin_alts.shape[0]
    if length_bin_alts > nb:
        #print('Length of bin_alts vectors is > than nbins. Stopping code...')
        pdb.set_trace()
    
    # Put the bins in the fixed altitude frame
    #counts_ff[:,i,:] = put_bins_onto_fixed_frame(ffrme,bin_alts,MCS_data_1file['counts'][i],nc) 
    #et = timeit.timeit("put_bins_onto_fixed_frame(ffrme,bin_alts,MCS_data_1file['counts'][i], \
        #nc,ONA*(180.0/pi)) ",setup="from lidar import put_bins_onto_fixed_frame\nfrom __main__ \
        #import ffrme,bin_alts,MCS_data_1file,i,nc,ONA,pi",number=1)
    counts_ff[:,i,:] = put_bins_onto_fixed_frame_C(np_clib,ffrme,
        bin_alts,MCS_data_1file['counts'][i],nc,mult) 
    #et = timeit.timeit("put_bins_onto_fixed_frame_C(np_clib,ffrme, \
        #bin_alts,MCS_data_1file['counts'][i],nc,mult)",
        #setup="from lidar import put_bins_onto_fixed_frame_C\nfrom __main__ \
        #import ffrme,bin_alts,MCS_data_1file,i,nc,ONA,pi,mult,np_clib",number=1)        
    #print("Rebinning time = %g " % et)
    
    i += 1 # increment record counter by 1

# Pick channel(s)
samp_chan = counts_ff[2,:,:] + counts_ff[3,:,:] # 532, in fixed frame

# Background subtract the rebinned counts and then plot
wl = 1 # numpy index 1 is 532 nm wl
opt = 'bg_sub'
EMs = convert_raw_energy_monitor_values(MCS_data_1file['meta']['EM'],nwl)
EMsubmit = EMs[:,wl] 
samp_chan = correct_raw_counts(samp_chan,EMsubmit,np.arange(0,nb_ff*vrZ_ff,vrZ_ff),
                               nr,nb_ff,975,1045,opt)
                               
# Code used to make masked scan angle plots ----------------------------
# ----------------------------------------------------------------------                               

# Save samp_chan at raw rez, before any averaging.
samp_chanX = np.empty_like(samp_chan)
samp_chanX[:,:] = samp_chan                               
print('The shape of samp_chan before averaging for total image = ',samp_chan.shape)                                  
samp_chan = average_lidar_data(samp_chan[15:,:],10,0)                                           
print('The shape after averaging = ',samp_chan.shape)
tit = 'CAMAL 532 nm background subtracted counts, all look angles'
pdb.set_trace()
curtain_plot(samp_chan[0:81].transpose(), nb_ff, vrZ_ff, ffrme, 0, 50.0, hori_cap, pointing_dir,
                      figW, figL, CPpad, 'records', 'altitude(m)', tit, 'alt',  
                      [1040,400], 'recs', [0,80], scale_alt_OofM, 1, out_dir)
## Apply a mask for each look angle. Plot the curtain w/gaps.
angle_arr = MCS_data_1file['meta']['scan_pos'] + angle_offset                               
[hist_bin_width, hist_bin_centers] = determine_look_angles(
                                     angle_arr,
                                     scan_pos_uplim, scan_pos_lowlim, 
                                     scan_pos_bw )     

n_angles = hist_bin_centers.shape[0]
print('There are '+str(n_angles).strip()+' angles.')                              
sti = np.array([57,158,259,360,462])                                     
edi = np.array([157,258,359,461,561])                     
mean_profs = np.zeros((n_angles,nb_ff))      
mean_ONA_along_transect = np.zeros(n_angles)               
mean_xzang_along_transect = np.zeros(n_angles)
for i in range(0,n_angles):   
    samp_chan = np.empty_like(samp_chanX)
    samp_chan[:,:] = samp_chanX  
    print('The shape of reloaded samp_chan is ',samp_chan.shape)
    mean_prof = np.mean(samp_chan[sti[i]:edi[i],:],axis=0)
    mean_profs[i,:] = mean_prof
    mean_ONA_along_transect[i] = np.mean(ONA_save[sti[i]:edi[i]])
    mean_xzang_along_transect[i] = np.mean(xzang_save[sti[i]:edi[i]])
    print('Done with profile for scan angle '+str(i))
    #plt.plot(mean_prof,ffrme)
    #plt.show()                                
    #samp_chan = apply_look_angle_mask(samp_chan,angle_arr,hist_bin_centers[i],hist_bin_width,nr,nhori,'gaps')  
    #samp_chan = average_lidar_data(samp_chan[15:],10,0)                              
    #samp_chan = samp_chan[0:81].transpose()
    #print('Then shape of samp_chan right before image is...',samp_chan.shape)
    #pdb.set_trace()                               
    #tit = 'CAMAL 532 nm background subtracted counts, look angle # '+str(i+1).strip()
    #curtain_plot(samp_chan, nb_ff, vrZ_ff, ffrme, 0, 50.0, hori_cap, pointing_dir,
    #                  figW, figL, CPpad, 'records', 'altitude (m)', tit, 'alt',  
    #                  [1040,400], 'recs', [0,80], scale_alt_OofM, 1, out_dir)
# ----------------------------------------------------------------------                  
# ----------------------------------------------------------------------
             
# Open a file in which to save profile arrays
pdb.set_trace()
arr_file = out_dir + 'save_prof_arrays_file26.npz'
# Save necessary arrays to file
np.savez(arr_file,p=mean_profs,z=ffrme,Oang=mean_ONA_along_transect,
         XZang=mean_xzang_along_transect)
print('Profile arrays have been saved. Code has finished!')
             
                      
pdb.set_trace()
