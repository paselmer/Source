# Libraries I have not created...
from tkinter import *
from tkinter import ttk
#from PIL import ImageTk, Image
import numpy as np
import pdb
import matplotlib.pyplot as plt
import os
from subprocess import check_output
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime as DT
import matplotlib.dates as mdates
# Libraries I have created...
from initializations import * #all directories and constants defined here
from read_routines import *
from immutable_data_structs import hk_struct
from GUI_function import *
from lidar import *
from time_conversions import *

def close_program():
    root.destroy()

def load_data():
    """ This function loads data into an object, and performs conversions """

    # Some declarations ---

    hk_file_list = config_dir + 'hk_file_list.txt'

    # Create file list using system commands. ---
    # Send different sys commands depending on OS.
    
    if (os.name != 'nt'):                                # Linux/Unix
        cmd = 'touch ' + hk_file_list
        cmd_feedback = check_output(cmd, shell=True)
        cmd = 'rm ' + hk_file_list
        cmd_feedback = check_output(cmd, shell=True)
        cmd = './make_hk_file_list_unix ' + raw_dir + ' > ' + hk_file_list
        cmd_feedback = check_output(cmd, shell=True)
    else:                                                # Windows
        cmd = 'type nul > ' + hk_file_list
        cmd_feedback = check_output(cmd, shell=True)
        cmd = 'del ' + hk_file_list
        cmd_feedback = check_output(cmd, shell=True)
        cmd = 'dir ' + raw_dir + 'hk* /B/S/OD > ' + hk_file_list 
        cmd_feedback = check_output(cmd, shell=True)

    rtn = ingest_housekeeping_data(hk_file_list, CAMAL_hk_obj)

    CAMAL_hk_obj.convert_temps_and_press()

def make_standard_plots():
    """ Make the plots you'll always want, no matter what 
        Try not to touch this function too much (I've probably cursed it
        to being modified constantly).
    """
    
    ybounds = [-20,50]       # <---- Change y-bounds
    ybounds_press = [-20,50] # <---- Change y-bounds
    ntemps = 15              # <---- Number of temperatures to be plotted
    tlab = ['Laser housing', 'Laser heat sink by center BackbB', 
            'Laser heat sing at vessel wall', 'Transceiver SS structure by laser',
            'Backbone by electronics housing', 'Aft optics bench', 
            'SPCM 1064perp', 'SPCM 1064 parallel',
            'SPCM 532 perp', 'SPCM 532 parallel',
            '355 PMT', 'Power convertors',
            'Encoder housing', 'Backbone by heater connector',
            'Lower SS structure by scanner' ]# <---- List of thermistor labels
    
    fig1 = plt.figure(1,figsize=(figW,figL))
    ax = plt.gca()
    box = ax.get_position()
    # The follow command set the size of the plotting area (axes) within 
    # the figure.
    ax.set_position([box.x0*0.5, box.y0, box.x1*0.70, box.y1*0.90])
    
    xax_type = xax_drop.get()      # get user-selected axis type
    if xax_type == "time":
        
        # Set the full axis to datetime objects, initially preferrably
        xvals = np.zeros(CAMAL_hk_obj.nrecs0,dtype=DT.datetime)
        for i in range(0,CAMAL_hk_obj.nrecs0):
            xvals[i] = weeksecondstoutc(CAMAL_hk_obj.data['GpsWeek'][i]*1.0,
                              CAMAL_hk_obj.data['Gps_msec'][i]/1e3,leapsecs )
        xtit = 'Time (UTC)' 
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                                     
    else:   
                 
        xvals = np.arange(0,CAMAL_hk_obj.nrecs0)
        xtit = 'Record number'
    
    # Temperatures plot ------------------------------------------------
    
    therm_ax_bounds = [ xvals[0], max(xvals) ] + ybounds
    # The first "ScanTemp" is actually a pressure reading, so start at 1
    var1, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,1], '#AC9E39', label = tlab[0], figure=fig1)
    var2, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,2], '#E99920', label = tlab[1], figure=fig1)
    var3, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,3], '#D38B78', label = tlab[2], figure=fig1)
    var4, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,4], '#5EFB56', label = tlab[3], figure=fig1)
    var5, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,5], '#56F8FB', label = tlab[4], figure=fig1)
    var6, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,6], '#94C1FF', label = tlab[5], figure=fig1)
    var7, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,7], '#1721EE', label = tlab[6], figure=fig1)
    var8, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,8], '#A017EE', label = tlab[7], figure=fig1)
    var9, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,9], '#322C35', label = tlab[8], figure=fig1)
    var10, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,10], '#818183', label = tlab[9], figure=fig1)
    var11, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,11], '#F31B44', label = tlab[10], figure=fig1)
    var12, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,12], '#9A1834', label = tlab[11], figure=fig1)
    var13, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,13], '#25B67C', label = tlab[12], figure=fig1)
    var14, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,14], '#F0EB24', label = tlab[13], figure=fig1)
    var15, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,15], '#218312', label = tlab[14], figure=fig1)
    var1.axes.axis( therm_ax_bounds ) 
    plt.legend( handles = [var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12,
                           var13, var14, var15] )
    # Place legend outside plotting area
    plt.legend(bbox_to_anchor=(1.01,1), loc="upper left")
    plt.xlabel(xtit)
    plt.ylabel('Temperature (deg C)')
    plt.title('CAMAL temperatures')
    fig1.canvas.draw_idle()
    
    # Pressure plot ----------------------------------------------------
    
    press_ax_bounds = [ xvals[0], max(xvals) ] + ybounds_press
    fig2 = plt.figure(2,figsize=(figW,figL))
    ax = plt.gca()
    box = ax.get_position()
    # The follow command set the size of the plotting area (axes) within 
    # the figure.
    ax.set_position([box.x0*0.5, box.y0, box.x1*0.8, box.y1*0.90])
    varp1, = plt.plot(xvals, CAMAL_hk_obj.data['ScanTemp'][:,0], '#1721EE', label = 'pressure', figure=fig2)
    varp1.axes.axis( press_ax_bounds ) 
    plt.legend( handles = [varp1] )
    # Place legend outside plotting area
    plt.legend(bbox_to_anchor=(1.01,1), loc="upper left")
    plt.xlabel(xtit)
    plt.ylabel('Pressure (?)')
    plt.title('CAMAL pressures')
    fig2.canvas.draw_idle()
    plt.show()
    plt.close(fig1)
    plt.close(fig2)
    
def make_custom_plots():
    """ Feel free to change here as needed 
        This is the place you should touch (modify), not the standard 
        plot function above.
    """
    pass
    
def enter_debug_mode():
    """ This will set a trace inside the code so that you can interrogate    
        the data as desired.
    """
    
    print('\n*****************************************************\n',end="")
    print('YOU HAVE ENTERED PYTHON DEBUGGER! \nFrom here you can ',end="")
    print('print variables, and interrogate \nthe data as desired. \n',end="")
    print("Enter 'c' to escape and return to GUI menu.\n",end="")
    print('*****************************************************\n',end="")
    pdb.set_trace()
    
    
# Object initializations go here

CAMAL_hk_obj = HandS_obj(None,hk_rec_size,None) #<-- Instantiation
rtn = load_data() # <--- Data get loaded right here
shp = CAMAL_hk_obj.data.shape
CAMAL_hk_obj.nrecs0 = shp[0]
print("There are " + str(CAMAL_hk_obj.nrecs0) + " records.")

# ---- Set up a very basic GUI using tkinter ----

root = Tk()
root.title("CAMAL HK plot toolbar")

ax_dbox_group = LabelFrame(root)
xax_drop_l = Label(ax_dbox_group,text="Xaxis type")
xax_drop = ttk.Combobox(ax_dbox_group,values=("recs","time"),width=7)
xax_drop.current(0)
xax_drop_l.grid(column=0, row=0)
xax_drop.grid(column=0, row=1)
#
stand_b = Button(root,text='Standard plots',command=make_standard_plots)
custm_b = Button(root,text='Custom plot(s)',command=make_custom_plots)
debug_b = Button(root,text='Debug/Interrogate',command=enter_debug_mode)
close_b = Button(root,text='Exit',command=close_program)

# Set the grid organization for widgets within root

ax_dbox_group.grid(column=0, row=0)
stand_b.grid(column=0, row=1)
custm_b.grid(column=0, row=2)
debug_b.grid(column=0, row=3)
close_b.grid(column=0, row=4)

# The following commands tell the root how to change the shape of its
# columns if the user resizes the window.
for child in root.winfo_children(): child.grid_configure(padx=4,pady=4)

root.mainloop()

# ----------------------------------------------------------------------
#For CAMAL please use the following conversions for the scanner pressure (1) and temperatures (2:16).

#ScannerTemps(:,2:end) = (ScannerTemps(:,2:end).^2*(-1.9724214))+(ScannerTemps(:,2:end).*(120.992655))-443.417057;

#ScannerTemps(:,1) = (ScannerTemps(:,1).*(7.4134))-3.20342

