# Libraries I have not created...
from tkinter import *
from tkinter import ttk
from PIL import ImageTk, Image
import numpy as np
import pdb
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import os
from subprocess import check_output
from mpl_toolkits.axes_grid1 import make_axes_locatable

import datetime as DT
# Libraries I have created...
from initializations import * #all directories and constants defined here
from read_routines import *
from mutable_data_structs import define_MSC_structure
from GUI_function import *


def close_program():
    root.quit()
    root.destroy()
    
def save_curtain_plot():
    if (os.name != 'nt'):                                # Linux/Unix
        cmd = 'cp '+ source_dir+'current_curtain.png ' + out_dir + \
               'curtain_chan_' + chan_sel.get() + '_' + str(DT.datetime.now().time())
        cmd_feedback = check_output(cmd, shell=True)
        print(cmd_feedback)
        print('Curtain should have been saved in '+out_dir)
    else:                                                # Windows
        print('Nothing yet. Working on it...')
        pdb.set_trace()

def load_and_plot(*args):
    
    # Check user channel input for error

    selected_chan = int(chan_sel.get())
    if selected_chan not in range(1,6):
        print('Channel entered is outside selected range.')
        print('Please re-enter a channel number.')
        return None
    
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
    print('All MCS data are loaded.')

    # Prepare data for a plot

    samp_chan = MCS_data['counts'][:,selected_chan-1,:]
    if nhori not in range(1,100):
        print('Set a more reasonable averaging # in initialization file.')
        return None
    elif (nhori > 1):
        samp_chan = average_lidar_data(samp_chan,nhori,0)
        print('Finished averaging channel '+chan_sel.get()+' to '+ \
               str(nhori) + ' profiles.')
        
    samp_chan = samp_chan.transpose()
    samp_chan = np.flipud(samp_chan) # reverse the order in columns (not rows)
    z = np.flipud(np.arange(0,nb*vrZ,vrZ))

    # Plot the data

    plt.clf()
    CPlot = make_curtain_plot_update_fig(samp_chan,nb,vrZ,z)
    canvas.show()
    
root = Tk()
root.title("Experiments")

# Frame to hold all button on RHS
RHSframe = ttk.Frame(root, borderwidth=2)
RHSframe.grid(column=0,row=0,stick=(N, W, E, S))
RHSframe.rowconfigure(0,weight=1)

# Image Canvas
fig99 = plt.figure(99,figsize=(figW,figL))
canvas =FigureCanvasTkAgg(fig99,master=root)
canvas.show()
canvas.get_tk_widget().grid(column=1,row=0)

# Buttons in RHS frame
load_b = Button(RHSframe,text='Load & plot',command=load_and_plot)
load_b.grid(column=0, row=0, sticky=W)
savcurt_b = Button(RHSframe,text='Save curtain',command=save_curtain_plot).grid(column=0, row=3, sticky=W)
close_b = Button(RHSframe,text='Exit',command=close_program).grid(column=0, row=4, sticky=W)

# Text entry boxes in RHS frame
chan_sel_l = ttk.Label(RHSframe, text='Selected channel').grid(column=0, row=1, sticky=S)
chan_sel = StringVar()
chan_entry = ttk.Entry(RHSframe, width=2, textvariable=chan_sel)
chan_entry.insert(0,'1') # put in a default value
chan_entry.grid(column=0, row=2, sticky=N)


# The following commands tell the "RHSframe" how to change the shape of its
# columns if the user resizes the window.
for child in RHSframe.winfo_children(): child.grid_configure(padx=5,pady=5)


root.mainloop()
