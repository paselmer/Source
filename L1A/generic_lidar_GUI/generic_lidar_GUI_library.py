import sys
sys.path.append('..')     # Code is one level down*
sys.path.append('../../') # L1A initializations inside "Source" directory
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
from read_routines import *
from mutable_data_structs import define_MSC_structure
from time_conversions import *
from generic_lidar_GUI_initializations import *

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



def make_curtain_plot(counts_imgarr, nb, vrZ, z, canvas_ctrl, cb_min, cb_max,
                      xlab, ylab, tit, yax, yax_lims, xax, xax_lims, mpl_flg):
    """ Function that will create a curtain plot """
    # USE THIS FUNCTION FOR IMAGE
    # The variable counts_imgarr is local to this function, and will not
    # mess up anything outside this scope. Therefore it can be sliced and
    # diced here, without crashing anything outside this scope.

    # Expand vertical dimension of image by using np.repeat [retired Oct 2017]
    # Subsetting array by user-input bins installed [12/6/17]
    # Bug fixed where formatting yaxis cancelled out appropriate labels for [9/18/18]
    # the case where the user wanted an 'alt' y-axis in the GUI curtain.
    
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
    plt.yticks(ytick_ind)
    # The above plt.yticks command doesn't seem to work in the GUI, so...
    ax.set_yticklabels(['{:.1f}'.format(item) for item in ytick_lab])
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right",size="5%",pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.autoscale
    plt.savefig('current_curtain.png',bbox_inches='tight',pad_inches=CPpad)
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
    plt.savefig('current_curtain.png',bbox_inches='tight')
    #plt.close(fig1)
    
    return None  
