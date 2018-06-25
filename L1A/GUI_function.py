from tkinter import *
from tkinter import ttk
from PIL import ImageTk, Image
import numpy as np
import pdb
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as tf
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.dates as mdates
from matplotlib.ticker import FormatStrFormatter
import datetime as DT
import math as math
# Stuff that I've created...
from lidar import average_lidar_data, get_a_color_map
from initializations import * #all directories and constants defined here
from read_routines import *
from mutable_data_structs import define_MSC_structure
from time_conversions import *

class file_control():
    """Object for managing the ingestion of various files"""
    
    def __init__(self, sel_file_list):
        self.sel_file_list = sel_file_list

class canvas_control():
    """Object for all modifiable properties related to canvas"""
    
    def __init__(self, img_bounds, pix2prof, prof_window_open, prof_plot_top,
                 prof_plot_canvas, lines2D_obj, ppi_this_screen,
                 var_window_open, var_plot_top, var_plot_canvas):
        self.img_bounds = img_bounds
        self.pix2prof = pix2prof
        self.prof_window_open = prof_window_open
        self.prof_plot_top = prof_plot_top
        self.prof_plot_canvas = prof_plot_canvas
        self.lines2D_obj = lines2D_obj
        self.ppi_this_screen = ppi_this_screen
        self.var_window_open = var_window_open
        self.var_plot_top = var_plot_top
        self.var_plot_canvas = var_plot_canvas
        
    def create_prof_plot_window(self,fig,fig_title):
        self.prof_plot_top = Toplevel()
        self.prof_plot_top.title(fig_title)
        self.prof_plot_canvas = FigureCanvasTkAgg(fig,master=self.prof_plot_top)
        self.prof_plot_canvas.draw()
        self.prof_plot_canvas.get_tk_widget().grid(column=0,row=0)
        
    def create_var_plot_window(self,fig,tit):
        """ Method to create a side plot in its own window. """
        self.var_plot_top = Toplevel()
        self.var_plot_top.title(tit)
        self.var_plot_canvas = FigureCanvasTkAgg(fig,master=self.var_plot_top)
        self.var_plot_canvas.draw()
        self.var_plot_canvas.get_tk_widget().grid(column=0,row=0)
        
    def update_any_canvas_ctrl_plot(self,x,y,fig):
        self.lines2D_obj.set_xdata(x)
        self.lines2D_obj.set_ydata(y)
        fig.canvas.draw_idle()        


class HandS_obj():
    """ This object holds and manipulates health and status data """
    
    def __init__(self, data, rec_size, nrecs0):
        """ Initialize attributes of this class """
        self.data = data
        self.rec_size = rec_size
        self.nrecs0 = nrecs0
        
    def convert_temps_and_press(self):
        """ Convert the raw temperatures and pressures into standard units """
        
        # Converts raw thermistor data to degrees Celcius
        self.data['ScanTemp'][:,1:] = ( (self.data['ScanTemp'][:,1:]**2)*(-1.9724214) ) + \
                                      ( self.data['ScanTemp'][:,1:]*(120.992655) ) - 443.417057 
                                      
        self.data['ScanTemp'][:,0] = (self.data['ScanTemp'][:,0]*(7.4134))-3.20342                                             
                                      

class lidar_data_obj():
    """ This object holds and manipulates lidar data """
    
    def __init__(self, data, std_params, ingested, nprofs, z, avgd_mask,
                 avgd_E, samp_chan, time_ax_arr, hist_bin_centers,
                 hist_bin_width, samp_chan_map, nprofs0):
       """Initialize attributes of this class"""
       self.data = data
       # std_params is list -> [nc,nb,nshots,vrT,vrZ] 
       self.std_params = std_params
       self.ingested = ingested
       self.nprofs = nprofs
       self.z = z
       self.avgd_mask = avgd_mask
       self.avgd_E = avgd_E
       self.samp_chan = samp_chan # current data being used in plots
       self.time_ax_arr = time_ax_arr
       self.hist_bin_centers = hist_bin_centers
       self.hist_bin_width = hist_bin_width
       self.samp_chan_map = samp_chan_map
       self.nprofs0 = nprofs0
        
    def update_avgd_mask(self,sel_chan):
        self.avgd_mask[sel_chan-1] = True
        
    def set_bin_alt_array(self,pdir,yax,ang_bin_indx):
        """pdir has string literal values of "Up" or "Down"
           depending on the direction the lidar is pointed.
        """
        
        # First bin index should correspond to 'all' angles. In this case,
        # apply no angle correction 
        if ang_bin_indx < 0:
            cosfact = 1.0 # same as applying to angle
            print('No cosine factor applied.')
            
        else:
            cosfact = math.cos(self.hist_bin_centers[ang_bin_indx]*(pi/180.0))
            print('------> cosine factor of '+str(self.hist_bin_centers[ang_bin_indx])+ \
                  ' being applied.')
        
        if yax == 'alt':
            if pdir == "Up":
                self.z = np.median(self.data['meta']['GpsAlt']) + np.arange(0,self.std_params[1]*self.std_params[4],
                                   self.std_params[4])
            elif pdir == "Down":
                self.z = np.median(self.data['meta']['GpsAlt']) - np.arange(0,
                                   self.std_params[1]*self.std_params[4],
                                   self.std_params[4]*cosfact)
            else:
                print("You entered a pointing_dir other than Up or Down. Assuming Down...")
                self.z = np.median(self.data['meta']['GpsAlt']) - np.arange(0,
                                   self.std_params[1]*self.std_params[4],
                                   self.std_params[4]*cosfact)
            #self.z = self.z * cosfact # Apply cosine no matter which above
        elif yax == 'bins':
            self.z = np.arange(1,self.std_params[1]+1)
            #if pdir == "Up":
            #    self.z = np.arange(1,self.std_params[1]+1)
            #elif pdir == "Down":
            #    self.z = np.arange(self.std_params[1],0,-1)
            #else:
            #    print("You entered a pointing_dir other than Up or Down. Assuming Down...")   
            #    self.z = np.arange(self.std_params[1],0,-1)      
            #pdb.set_trace()
        else:
            print("Something is seriously wrong with your code, Patrick. Stopping in \
                   GUI_function.py")
            pdb.set_trace()
                               
    def set_time_axis_array(self):
        """ This method will create the time axis of samp_chan. 
            Unless you know what you're doing, only call this immediately
            after initial data read at raw resolution.
        """
        
        # Set the full axis to datetime objects, initially preferrably
        self.time_ax_arr = np.zeros(self.nprofs,dtype=DT.datetime)
        for i in range(0,self.nprofs):
            #pdb.set_trace()
            self.time_ax_arr[i] = weeksecondstoutc(self.data['meta']['GpsWeek'][i]*1.0,
                                self.data['meta']['Gps_msec'][i]/1e3,leapsecs )
                                
    def determine_look_angles(self):
        """ 
        This method will determine all the look angles within data.
        Originally wittten with CAMAL in mind because CAMAL has
        programmable, continuously variable scan angles.
        
        This method will produce a list of the float values of the
        angles. 
        """
        
        # Correct scan angles for fixed, known offset. Provided in 
        # initialization file.
        
        self.data['meta']['scan_pos'] = self.data['meta']['scan_pos'] + angle_offset
        
        # The histogram parameters here come from the initialization file
        nhbins = int( (scan_pos_uplim - scan_pos_lowlim)/scan_pos_bw )
        dens, edges = np.histogram(self.data['meta']['scan_pos'],
                      bins=nhbins,range=(scan_pos_lowlim,scan_pos_uplim))
        delta = float(scan_pos_uplim - scan_pos_lowlim) / float(nhbins)
        print("delta is ",delta)
        center = edges[:len(dens)] + delta
        mask = dens != 0
        #print("The bin edges are :",edges)
        print("This many bins have a freq of at least 1: ",dens[mask].shape)
        print("These are the bins: ",center[mask])
        self.hist_bin_centers = center[mask]
        self.hist_bin_width = delta
        #plt.plot(center,dens)
        #plt.show()
        
    def apply_look_angle_mask(self,bin_indx,opt):
        """ Code to apply mask to science data based on look angle 
            Like other methods in this class, don't call before you load
            data.
        """   
        
        # This maps the full-length dataset to samp_chan
        samp_chan_map = np.arange(0,self.nprofs).astype('uint32')
        self.samp_chan_map = samp_chan_map # could be overwritten below
        
        # Go no farther if no angle is selected
        if (opt == 'no_mask'): return
        # If you make it to this line, a mask is going to be applied
        print("An angle mask is being applied.")
        
        cs = self.hist_bin_centers[bin_indx]
        lim1 = cs - self.hist_bin_width
        lim2 = cs + self.hist_bin_width
        print( 'limit1 = '+str(lim1).strip()+' limit2 = '+str(lim2).strip() )

        # Make an array mask to mask "where" angle fall within bin
        # First reduce array based on averaging, ang_mask then will be
        # in "averaged space."
        scan_pos_reduced = self.data[::nhori]['meta']['scan_pos']
        scan_pos_reduced = scan_pos_reduced[0:self.nprofs]
        ang_mask = ( (scan_pos_reduced >= lim1) & 
             (scan_pos_reduced <= lim2)  )
        #pdb.set_trace()
        
        if opt == 'nogaps':
            
            # Here we make a sample (samp_chan) with no gaps
            # Remember, ang_mask will be in averaged space.
            self.samp_chan = self.samp_chan[:,:][samp_chan_map[ang_mask]]
            print('The shape of samp_chan is ',self.samp_chan.shape)
            # Finalize the samp_chan_map. The numbers within this array
            # should correspond to full-rez dataset indicies.
            self.samp_chan_map = samp_chan_map[ang_mask]
            
        elif opt == 'gaps':
            
            # Here we make a sample (samp chan) with gaps
            self.samp_chan[:,:][samp_chan_map[(ang_mask==False)]] = 0
        
        # Wrap up this method, overwrite self.nprofs
        samp_chan_shape = self.samp_chan.shape
        print('The shape of samp_chan is ',samp_chan_shape) 
        self.nprofs = samp_chan_shape[0]
        print('Length of samp_chan: ',samp_chan_shape[0])    
        
        
        

def ingest_entire_dataset(MCS_file_list, raw_data_object, Fcontrol):
    """Function (module, whatev) that will read in entire
       CAMAL dataset (usually in form of a single flight).
       Will take list file as arg, and return list of 
       data parameters.
    """
    
    with open(MCS_file_list) as MCS_list_fobj:
        all_MCS_files = MCS_list_fobj.readlines()
    nMCS_files = len(all_MCS_files)

    # Read in the MCS (science) data

    first_read = 1
    r=0
    for i in range(0,nMCS_files):
        if i not in Fcontrol.sel_file_list: continue
        print("File # ====== ",i)
        MCS_file = all_MCS_files[i]
        MCS_file = MCS_file.rstrip()
        MCS_data_1file = read_in_raw_data(MCS_file)
        if MCS_data_1file is None:
            print('Skipping this file. Potentially corrupt data')
            continue
        skip_whole_file = filter_out_bad_recs(MCS_data_1file)
        if skip_whole_file == 1:
            print('Skipping this file...') 
            continue
        if first_read:
            first_read = 0	
            # Put the parameters that won't change during 1 flight into variables
            nc = MCS_data_1file['meta']['nchans'][0]
            nb = MCS_data_1file['meta']['nbins'][0]
            nshots = MCS_data_1file['meta']['nshots'][0]
            vrT = MCS_data_1file['meta']['binwid'][0]
            vrZ = (vrT * c) / 2.0
            data_params = [nc,nb,nshots,vrT,vrZ] #store in list
            # declare data structure to hold all data. estimate tot # recs first
            n_files = min(nMCS_files,len(Fcontrol.sel_file_list))
            tot_est_recs = int(rep_rate/nshots)*file_len_secs*n_files
            MCS_struct = define_MSC_structure(nc,nb)
            MCS_data = np.zeros(tot_est_recs, dtype=MCS_struct)
        nr_1file = MCS_data_1file.shape[0]
        MCS_data[r:r+nr_1file] = MCS_data_1file
        r += nr_1file   
        #NOTE: You could put conditional break
        #      statement in this loop to read-in
        #      data from time segment only. 
    
    try:
        MCS_data = MCS_data[0:r]
    except UnboundLocalError:
        print("\n\n******************************************\n")
        print("There is no valid data in selected files. Pick other files.\n")
        print("******************************************\n")
        return False
    
    raw_data_object.data = MCS_data
    raw_data_object.std_params = data_params
    raw_data_object.ingested = True
    raw_data_object.nprofs = r
    raw_data_object.avgd_mask = np.zeros(nc,dtype=bool)
    print('All MCS data are loaded.')    
    return True     


def make_curtain_plot(counts_imgarr, nb, vrZ, z, canvas_ctrl, cb_min, cb_max,
                      xlab, ylab, tit, yax, yax_lims, xax, xax_lims, mpl_flg):
    """ Function that will create a curtain plot """
    # USE THIS FUNCTION FOR IMAGE
    # The variable counts_imgarr is local to this function, and will not
    # mess up anything outside this scope. Therefore it can be sliced and
    # diced here, without crashing anything outside this scope.

    # expand vertical dimension of image by using np.repeat [retired Oct 2017]
    # subsetting array by user-input bins installed [12/6/17]
    counts_imgarr = counts_imgarr[int(min(yax_lims)):int(max(yax_lims)),:]
    img_arr_shape = counts_imgarr.shape
    print('The shape of counts_imgarr is: ',img_arr_shape)
    
    # horizontal dimension capped at certain # of profiles (in init file)
    # doing this to save memory & CPU work when rendering image
    if img_arr_shape[1] > hori_cap:
        nthin = int(img_arr_shape[1] / hori_cap)
        counts_imgarr = counts_imgarr[:,::nthin]
        img_arr_shape = counts_imgarr.shape
        print('The shape of counts_imgarr is: ',img_arr_shape)     

    # you need to manipulate altitude array in order to properly
    # label image plot
    if pointing_dir == "Up":
        alt_imgarr = np.flipud(z)
    else:
        alt_imgarr = z
    newsize=alt_imgarr.size
    alt_ind = np.linspace(0,newsize-1,newsize,dtype='uint32')
    print('The shape of alt_ind is: ',alt_ind.shape)
    #yax_lims = [ alt_ind[newsize-1],alt_ind[0] ]
    #yax_lims = [ y1,y2 ] # y1, y2 determine zoom along vertical axis

    # Actually plot the photon counts
    fig1 = plt.figure(1,figsize=(figW,figL))
    ax = plt.gca() # ax is now the "handle" to the figure
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(tit)
    cm = get_a_color_map()
    im = ax.imshow(counts_imgarr, cmap=cm, clim=(cb_min,cb_max),
                   interpolation='nearest',aspect='auto', 
                   extent=xax_lims+yax_lims)                 
                   
    # Format the x-axis   
    if (xax == 'time'):
        ax.xaxis_date()
        time_format = mdates.DateFormatter('%H:%M:%S')
        ax.xaxis.set_major_formatter(time_format)
        fig1.autofmt_xdate()
                   
    # Format the y-axis
    locs, labels = plt.yticks()
    if (yax == 'alt'): 
        delta_check = vrZ
        y2_ind = int(max(yax_lims))
        if y2_ind >= nb: 
            y2 = z[nb-1]
        else:
            y2 = z[y2_ind]
        y1 = z[int(min(yax_lims))]            
        ys = [y1,y2]
        min_y = math.ceil(min(ys)/scale_alt_OofM)*scale_alt_OofM
        max_y = math.floor(max(ys)/scale_alt_OofM)*scale_alt_OofM
        ndivi = int( (max_y - min_y) / scale_alt_OofM )
        ytick_lab = np.linspace(min_y,max_y, ndivi+1)
    else:
        delta_check = 2
        ndivi = 20
        ytick_lab = np.linspace(yax_lims[0],yax_lims[1],ndivi+1)
    ytick_ind = np.zeros(ytick_lab.size) + 999
    k=0
    for e in ytick_lab:
        m = abs(e - alt_imgarr) < delta_check # where diff is smaller than 60 meters
        if (np.sum(m) == 0):  #all false
            k=k+1
            continue
        else:
            ytick_ind[k] = alt_ind[m][0]
            k=k+1
    actual_ticks_mask = ytick_ind != 999
    ytick_lab = ytick_lab[actual_ticks_mask]
    ytick_ind = ytick_ind[actual_ticks_mask]          	    
    plt.yticks(ytick_ind,ytick_lab)
    # Format the ytick labels
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  
        
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right",size="5%",pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.autoscale
    plt.savefig(source_dir+'current_curtain.png',bbox_inches='tight',pad_inches=CPpad)
    CPpad_pix = CPpad * canvas_ctrl.ppi_this_screen # CPpad in pixels 
    # The following 1 line stores the "dpi" coordinates of the image corners
    im_Bbox = tf.Bbox.get_points(im.get_window_extent())
    im_Bbox = (im_Bbox / fig1.dpi) * canvas_ctrl.ppi_this_screen
    width_im_box_pix = im_Bbox[1,0] - im_Bbox[0,0]
    # Now get the pixel coordinates bounding the y-axis label
    yax = ax.get_yaxis()
    yax_Bbox = yax.get_tightbbox(fig1.canvas.get_renderer())
    yax_bb = tf.Bbox.get_points(yax_Bbox)
    yax_bb = (yax_bb/fig1.dpi) * canvas_ctrl.ppi_this_screen
    # Use the yax_Bbox to set the left edge of the "tight" curtain image
    im_Bbox[0:2,0] = im_Bbox[0:2,0] - yax_bb[0,0] + CPpad_pix
    if mpl_flg == 1: plt.show()
    plt.close(fig1)
    
    return [fig1, width_im_box_pix, im_Bbox]


def ingest_housekeeping_data(hk_file_list, raw_data_object):
    """Function (module, whatev) that will read in entire
       CAMAL dataset (usually in form of a single flight).
       Will take list file as arg, and populate the object
       which is passed to it. Returns True.
    """
    
    with open(hk_file_list) as hk_list_fobj:
        all_hk_files = hk_list_fobj.readlines()
    nhk_files = len(all_hk_files)
    
    # Allocate array size based on total estimated # of records
    
    statinfo = os.stat(all_hk_files[0].rstrip())
    size1file = statinfo.st_size
    sample_record = np.zeros(1,dtype=hk_struct)
    recs_per_file = int( size1file / raw_data_object.rec_size )
    est_recs = recs_per_file * nhk_files
    hk_data = np.zeros(est_recs,dtype=hk_struct)

    # Read in the MCS (science) data

    r=0
    for i in range(0,nhk_files):
        print("hk File # ====== ",i)
        hk_file = all_hk_files[i]
        hk_file = hk_file.rstrip()
        hk_data_1file = read_in_housekeeping_data(hk_file)
        nr_1file = hk_data_1file.shape[0]
        try:
            hk_data[r:r+nr_1file] = hk_data_1file
        except ValueError:
            pdb.set_trace()
        r += nr_1file   
        #NOTE: You could put conditional break
        #      statement in this loop to read-in
        #      data from time segment only. 
          
    hk_data = hk_data[0:r]
    
    raw_data_object.data = hk_data
    print('All housekeeping data are loaded.')    
    return True     


def ingest_gps_data(gps_file_list):
    """Function (module, whatev) that will read in entire
       CAMAL GPS dataset (usually in form of a single flight).
       Will take list file as arg, and return a numpy array
       of the GPS data.
    """
    
    with open(gps_file_list) as gps_list_fobj:
        all_gps_files = gps_list_fobj.readlines()
    ngps_files = len(all_gps_files)
    
    # Allocate array size based on total estimated # of records
    
    est_recs_1file = file_len_secs * (gps_hz+1) # add a little buffer
    est_recs = int(est_recs_1file * ngps_files)
    gps_struct_1file = np.zeros(est_recs, dtype=gps_struct)
    gps_data = np.zeros(est_recs,dtype=gps_struct)

    # Read in the MCS (science) data

    r=0
    for i in range(0,ngps_files):
        print("hk File # ====== ",i)
        gps_file = all_gps_files[i]
        gps_file = gps_file.rstrip()
        gps_data_1file = read_in_gps_data(gps_file)
        nr_1file = gps_data_1file.shape[0]
        try:
            gps_data[r:r+nr_1file] = gps_data_1file
        except ValueError:
            pdb.set_trace()
        r += nr_1file   
        #NOTE: You could put conditional break
        #      statement in this loop to read-in
        #      data from time segment only. 
          
    gps_data = gps_data[0:r]

    print('All GPS data are loaded.')    
    return gps_data  






    
def make_curtain_plot_update_fig(samp_chan,nb,vrZ,z):
    """ Function that will create a curtain plot """
    # USE THIS VERSION FOR FIGURE

    # expand vertical dimension of image by using np.repeat
    counts_imgarr = np.repeat(samp_chan,vstretch,axis=0)
    img_arr_shape = counts_imgarr.shape
    print('The shape of counts_imgarr is: ',img_arr_shape)
    
    # horizontal dimension capped at certain # of profiles (in init file)
    if img_arr_shape[1] > hori_cap:
        nthin = int(img_arr_shape[1] / hori_cap)
        counts_imgarr = counts_imgarr[:,::nthin]
        print('The shape of counts_imgarr is: ',counts_imgarr.shape)

    # you need to manipulate altitude array in order to properly
    # label image plot
    alt_imgarr = np.repeat(z,vstretch)
    print('The shape of alt_imgarr is: ',alt_imgarr.shape)
    newsize=alt_imgarr.size
    alt_ind = np.linspace(0,newsize-1,newsize,dtype='uint32')
    print('The shape of alt_ind is: ',alt_ind.shape)

    # Actually plot the photon counts
    ax = plt.gca() # ax is now the "handle" to the figure
    im = ax.imshow(counts_imgarr, cmap='gist_ncar', clim=(0.0,cbar_max))
    locs, labels = plt.yticks()
    ytick_lab = np.linspace(0,max_scale_alt,11)
    ytick_ind = np.zeros(ytick_lab.size)
    k=0
    for e in ytick_lab:
        m = abs(e - alt_imgarr) < vrZ # where diff is smaller than 60 meters
        ytick_ind[k] = alt_ind[m][0]
        k=k+1	    
    plt.yticks(ytick_ind,ytick_lab)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right",size="5%",pad=0.05)
    plt.colorbar(im, cax=cax)
    #plt.show()
    ax.autoscale
    plt.savefig(source_dir+'current_curtain.png',bbox_inches='tight')
    #plt.close(fig1)
    
    return None  
