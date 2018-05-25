# Libraries I have not created...
import numpy as np
import pdb
import matplotlib.pyplot as plt
import os
from subprocess import check_output
from mpl_toolkits.axes_grid1 import make_axes_locatable
# Libraries I have created...
from initializations import * #all directories and constants defined here
from read_routines import *
from mutable_data_structs import define_MSC_structure
from GUI_function import *



# Some declarations

MCS_file_list = config_dir + 'MCS_file_list.txt'


# Create file list using system commands.
# Send different sys commands depending on OS.
# At end of this block, a list should contain the full-path file names
# of all MCS data files to loop thru.

if (os.name != 'nt'):                                # Linux/Unix
	cmd = 'touch ' + MCS_file_list
	cmd_feedback = check_output(cmd, shell=True)
	cmd = 'rm ' + MCS_file_list
	cmd_feedback = check_output(cmd, shell=True)
	cmd = './make_file_list_unix ' + raw_dir + ' > ' + MCS_file_list
	cmd_feedback = check_output(cmd, shell=True)
else:                                                # Windows
	print('Nothing yet. Working on it...')
	pdb.set_trace()	

with open(MCS_file_list) as MCS_list_fobj:
    all_MCS_files = MCS_list_fobj.readlines()
nMCS_files = len(all_MCS_files)

# Read in the MCS (science) data

first_read = 1
r=0
for MCS_file in all_MCS_files:
    MCS_file = MCS_file.rstrip()
    MCS_data_1file = read_in_raw_data(MCS_file)
    if first_read:
        first_read = 0	
        # Put the parameters that won't change during 1 flight into variables
        nc = MCS_data_1file['meta']['nchans'][0]
        nb = MCS_data_1file['meta']['nbins'][0]
        nshots = MCS_data_1file['meta']['nshots'][0]
        vrT = MCS_data_1file['meta']['binwid'][0]
        vrZ = (vrT * c) / 2.0
        # declare data structure to hold all data. estimate tot # recs first
        tot_est_recs = int(rep_rate/nshots)*file_len_secs*nMCS_files
        MCS_struct = define_MSC_structure(nc,nb)
        MCS_data = np.zeros(tot_est_recs, dtype=MCS_struct)
    nr_1file = MCS_data_1file.shape[0]
    MCS_data[r:r+nr_1file] = MCS_data_1file
    r += nr_1file              
    #NOTE: You could put conditional break
    #      statement in this loop to read-in
    #      data from time segment only.
    
    
MCS_data = MCS_data[0:r]

# Prepare data for a plot

samp_chan = MCS_data['counts'][:,4,:]
samp_chan = average_lidar_data(samp_chan,10,0)
samp_chan = samp_chan.transpose()
samp_chan = np.flipud(samp_chan) # reverse the order in columns (not rows)
z = np.flipud(np.arange(0,nb*vrZ,vrZ))

# Plot the data

CPlot = make_curtain_plot(samp_chan,nb,vrZ,z)

pdb.set_trace()



