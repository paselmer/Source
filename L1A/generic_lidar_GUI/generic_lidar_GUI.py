# Main program level of GUI to perform basic visualization and analysis of
# lidars crafted in the vein of Goddard Space Flight Center's Cloud Physics
# Lidar (CPL).

import sys
sys.path.append('..')     # Code is one level down*
sys.path.append('../../') # L1A initializations inside "Source" directory
import numpy as np
import matplotlib.pyplot as plt
import pdb
# Custom libraries/modules
from generic_lidar_GUI_library import *
from standardize_lidar_data import *
from read_routines import * # no longer imports L1A initializations [10/9/19]
from generic_lidar_GUI_initializations import * # !!NEEDS TO BE LAST!!

##file_list = [x for x in range(7,10)]
##lidar_obj = select_and_standardize_data(file_list)


def close_program():
    root.destroy()
    
def close_prof_plot_window():
    canvas_ctrl.prof_plot_top.destroy()
    canvas_ctrl.prof_window_open = False 
    plt.clf()    
    
def close_var_plot_window():
    canvas_ctrl.var_plot_top.destroy() 
    canvas_ctrl.var_window_open = False 
    
def parse_file_selection(sel_str):
    """ Take a list of string integers and convert them to actual
        integers. Treat '-' as a range and ',' and meaning individual
        files.
    """
    
    # First, check for a hypen '-' or comma ','
    indx_hy = sel_str.find('-')
    indx_cm = sel_str.find(',')
    
    if (indx_hy != -1) and (indx_cm != -1):
        # Will not allow both a '-' and ','
        print("File selection input cannot contain both a '-' and ','\n")
        print("Selecting all files...")
        return list(range(0,99))
    elif (indx_hy != -1):
        # If '-' detected, split string
        substrs = sel_str.split('-')
        return(list(range(int(substrs[0])-1,int(substrs[1])))) # assume only 2 values entered
    elif (indx_cm != -1):
        substrs = sel_str.split(',')
        for i in range(0,len(substrs)):
            substrs[i] = int(substrs[i])-1 # will be used for subscript            
        return substrs
    else:
        return [int(sel_str)-1]        
    #pdb.set_trace()
    
def get_selected_channels(return_mode):
    """ This function gets the selected channel input. It allows the user
        to enter multiple channels to be summed, as long as they are
        separated by -'s. For example, if the user wants to display a 
        curtain of channels 1 and 2, they would enter "1-2" in the entry
        box. The first channel entered by the user (even if only one is
        entered) is that the channel that is used for label retrieval
    """
    
    while True:
        
        inString = "".join(chan_sel.get().split())
        sel_chans = inString.split('-')
        
        try:
            
            if return_mode == 'first':
                return int(sel_chans[0])
            elif return_mode == 'all':
                for i in range(0,len(sel_chans)): sel_chans[i] = int(sel_chans[i])
                return sel_chans
            else:
                print("return_mode entry is invalid. Stoppping code.") 
                pdb.set_trace()
                
        except ValueError:
            
            print('*****************************************************\n')
            print("Something is wrong with the selected channel input\n")
            print('You entered "' + inString + '"\n')
            print('If entering multiple channels (to sum), use a hypen (-)\n')
            print('For example, 3-4.\n')
            print('*****************************************************\n')
            response = input('Enter your new channel selection right here at the prompt...\n') 
            chan_sel.set(response)
        
                    
    
def save_curtain_plot():
    if (os.name != 'nt'):                                # Linux/Unix
        cmd = 'cp ' + 'current_curtain.png ' + out_dir + \
               'curtain_chan_' + "".join(chan_sel.get().split()) + '_' +  \
               str(DT.datetime.now().time()) + '.png'
        cmd_feedback = check_output(cmd, shell=True)
        print(cmd_feedback)
        print('Curtain should have been saved in '+out_dir)
    else:                                                # Windows
        file_tag = str(DT.datetime.now().time()).split(':')
        sep = '-'
        Win_file_tag = sep.join(file_tag)
        cmd = 'copy ' + 'current_curtain.png ' + out_dir + \
              'curtain_chan_' + "".join(chan_sel.get().split()) + '_' + Win_file_tag + '.png'
        cmd_feedback = check_output(cmd, shell=True)
        print(cmd_feedback)
        print('Curtain should have been saved in '+out_dir)
 
        
def save_EM_plot():
    if (os.name != 'nt'):                                # Linux/Unix
        cmd = 'cp ' + 'current_EM_plot.png ' + out_dir + \
               'EM_plot_' + str(DT.datetime.now().time())
        cmd_feedback = check_output(cmd, shell=True)
        print(cmd_feedback)
        print('EM plot should have been saved in '+out_dir)
    else:                                                # Windows
        file_tag = str(DT.datetime.now().time()).split(':')
        sep = '-'
        Win_file_tag = sep.join(file_tag)
        cmd = 'copy ' + 'current_EM_plot.png ' + out_dir + \
              'EM_plot_' + Win_file_tag + '.png'
        cmd_feedback = check_output(cmd, shell=True)
        print(cmd_feedback)
        print('EM plot should have been saved in '+out_dir)     
        
def save_prof_plot():
    """Saves figure directly as opposed to other similar functions.
       Also note file name is the same for all OS's here.
    """
    
    # fig2 should be the current figure right now, but just in case it
    # isn't, set the figure
    plt.figure(fig2.number)
    
    file_tag = str(DT.datetime.now().time()).split(':')
    sep = '-'
    file_tag = sep.join(file_tag)
    plt.savefig(out_dir+'prof_plot_'+file_tag+'.png',bbox_inches='tight',pad_inches=CPpad)
    print('Profile plot should have been saved in '+out_dir)
        

def save_popup_window_plot():
    """The way I have this GUI coded, only one pop-up window plot
       may be displayebackgroundd at any given time. This function will save
       whichever is currently open. The objects I have built
       contain boolean attributes to tell me which type of window is
       open. Based on the type, different functions specific to window
       names are called to save their respective plots.
    """
    
    if canvas_ctrl.var_window_open is True:
        save_EM_plot()
    elif canvas_ctrl.prof_window_open is True:
        save_prof_plot()
    else:
        print('No window plot is open. Nothing saved.')
        
        
def get_plot_parameters_for_current_data_type(curtain_type, yax_type, selected_chan):
    """ This code gets the parameters that change each time you select
        a new channel or data type (counts, NRB, ect.)
    """
    
    wl = wl_map[selected_chan-1]
    
    opt = 'bg_sub' #default. opt is for lidar data type (counts, NRB...)
    lab_tag = 'raw counts'
    pp_xax_lims = pp_xax_bounds_raw_counts
    if curtain_type == 2: 
        opt='bg_sub'
        lab_tag = 'counts (background subtracted)'
        pp_xax_lims = pp_xax_bounds_bgsub_counts
    if curtain_type == 3: 
        opt='NRB'
        lab_tag = 'NRB'
        pp_xax_lims = pp_xax_bounds_NRB
    
    if yax_type == 'bins':
	    pp_yax_lims = pp_yax_bounds_bins
    elif yax_type == 'alt':
	    pp_yax_lims = pp_yax_bounds_alt
    else:
	    print('Stopping. Something wrong with yax_type.')
	    pdb.set_trace()
    
    pp_ax_lims = pp_xax_lims + pp_yax_lims
    return [opt, lab_tag, wl, pp_ax_lims]       

def load_and_plot(*args):
    
    # Check user channel input for error ---

    selected_chan = get_selected_channels( 'all' )
    if (min(selected_chan) not in range(1,6)) or (max(selected_chan) not in range(1,6)):
        print('Channel entered is outside hard-coded limits.')
        print('Please re-enter a channel number.')
        return None

    # Check user input for selected files ---
    
    if ts_switch.get() == 0:
        new_sfl = parse_file_selection(sel_files.get())
    else:
        # Create a list so long, it's bound to include all files
        new_sfl = list(range(0,99))
    if (file_ctrl.sel_file_list == new_sfl) is False:
        lidar_obj.ingested = False
    print("-----> Ingested = ",lidar_obj.ingested)        
    file_ctrl.sel_file_list = new_sfl        
        
    # Read in entire dataset ---

    reinitialize_avgd_mask = False
    if lidar_obj.ingested is False:
        # Read in the data from selected instrument & files (controlled by file_ctrl)
        lidar_obj_temp = select_and_standardize_data(file_ctrl.sel_file_list)
        lidar_obj.nr = lidar_obj_temp.nr
        lidar_obj.nb = lidar_obj_temp.nb
        lidar_obj.nc = lidar_obj_temp.nc
        lidar_obj.dz = lidar_obj_temp.dz
        lidar_obj.dx = lidar_obj_temp.dx
        lidar_obj.nr0 = lidar_obj_temp.nr0
        lidar_obj.nshots = lidar_obj_temp.nshots
        lidar_obj.nhori = lidar_obj_temp.nhori
        lidar_obj.pdir = lidar_obj_temp.pdir
        lidar_obj.rec_arr = lidar_obj_temp.rec_arr
        lidar_obj.time_arr = lidar_obj_temp.time_arr
        lidar_obj.bin_arr = lidar_obj_temp.bin_arr
        lidar_obj.alt_arr = lidar_obj_temp.alt_arr
        lidar_obj.platform_alt_arr = lidar_obj_temp.platform_alt_arr
        lidar_obj.counts = lidar_obj_temp.counts
        lidar_obj.energy = lidar_obj_temp.energy
        lidar_obj.ingested = lidar_obj_temp.ingested
        lidar_obj.avgd_mask = lidar_obj_temp.avgd_mask
        lidar_obj.samp_chan_map = lidar_obj_temp.samp_chan_map
        lidar_obj.ONA_arr = lidar_obj_temp.ONA_arr
        lidar_obj.ONA_bw = lidar_obj_temp.ONA_bw
        lidar_obj.ONA_uplim = lidar_obj_temp.ONA_uplim
        lidar_obj.ONA_lowlim = lidar_obj_temp.ONA_lowlim
        lidar_obj.ONA_bin_centers = lidar_obj_temp.ONA_bin_centers
        lidar_obj.ONA_internal_bw = lidar_obj_temp.ONA_internal_bw
        lidar_obj.samp_chan = lidar_obj_temp.samp_chan        
        # Now that data have been read-in, set the ingested flag to True
        lidar_obj.ingested = True
        # In addition to data ingestion, compute PPI here (so only does 1x)
        ppi_img = ImageTk.PhotoImage(Image.open('PPI_ref_img.png'))
        # actual computation of PPI
        canvas_ctrl.ppi_this_screen = ppi_img.width()/1.0
        if lidar_obj.ONA_bw is not None: lidar_obj.bin_ONAs()
        # Write the look angles to strings, then insert in GUI dropbox
        cur_angles = ['all']
        # Save the original nprofs. NOTE: This will be overwritten if
        # averaging is set to be performed.
        lidar_obj.nr0 = lidar_obj.counts.shape[0]
        if lidar_obj.ONA_bin_centers is not None:
            for angle in lidar_Obj.ONA_bin_centers:
                # angle_list_format_str is in initialization file
                angle_str = angle_list_format_str.format(angle)
                cur_angles.append(angle_str)
        ang_list.config(values=(cur_angles))
        ang_list.current(0)
        reinitialize_avgd_mask = True	

    # Set the y-axis array and save other essential params ---  
    yax_type = yax_drop.get()      # get user-selected axis type
    xax_type = xax_drop.get()      # get user-selected axis type
    print('ang_list.current() = ',ang_list.current())
    nc = lidar_obj.nc
    nb = lidar_obj.nb
    # Set the averaged mask here. In CAMAL GUI this is done outside this file.
    if reinitialize_avgd_mask: lidar_obj.avgd_mask = np.zeros(nc,dtype=bool)
    nshots = lidar_obj.nshots
    nr = lidar_obj.nr
    if yax_type != 'alt': 
        z = lidar_obj.bin_arr
    else:
        z = lidar_obj.alt_arr
    curtain_type = curt_type.get() # radio button number

    # Prepare data for a plot (average data, if configured to do so) ---

    # Average each of the channels entered by user (could only be 1 chan)
    # If nhori set to 1, nothing should be executed in following loop
    for i in range(0,len(selected_chan)):
        if not lidar_obj.avgd_mask[selected_chan[i]-1]:
            if nhori not in range(1,200):
                print('Set a more reasonable averaging # in initialization file.')
                return None
            elif (nhori > 1):
                samp_chan = lidar_obj.data['counts'][:,selected_chan[i]-1,:]
                samp_chan = average_lidar_data(samp_chan,nhori,0)
                new_shape = samp_chan.shape
                lidar_obj.nr = new_shape[0]
                lidar_obj.update_avgd_mask(selected_chan[i])
                lidar_obj.counts[0:new_shape[0],selected_chan[i]-1,:] = samp_chan
                print('Finished averaging channel '+str(selected_chan[0])+' to '+ \
                     str(nhori) + ' profiles.')
                # Save the nprofs parameter. Masking out the data by look angle will
                # alter nprofs. This will cause a bug in the code IF YOU DON'T SAVE
                # IT HERE!!!!!!!!
                lidar_obj.nr0 = lidar_obj.nr                                             
    # Use numpy method to sum together channels. Store in "samp_chan"
    # This following statement is valid whether or not data were averaged                               
    samp_chan = lidar_obj.counts[0:lidar_obj.nr0,
                selected_chan[0]-1:selected_chan[len(selected_chan)-1],:].sum(1)
    lidar_obj.nr = lidar_obj.nr0 # RESET nprofs here!!!                                  
                   
    # Retreive plot params                                     
    [opt, lab_tag, wl, pp_ax_lims] = get_plot_parameters_for_current_data_type(curtain_type,
                              yax_type,selected_chan[0])      
    # Correct photon counts if user has check box checked.
    # I've decided to only apply this to samp_chan and not entire dataset.
    if curtain_type > 1:
        wl = wl_map[selected_chan[0]-1]
        if lidar_obj.energy is not None:
            EMs = convert_raw_energy_monitor_values(lidar_obj.data['meta']['Engineering']['LaserEnergyMonitors'],nwl,instrument_name,e_flg)
            EMsubmit = EMs[:,wl]
        else:
            print('No energy data. Setting energy array to all ones.')
            EMsubmit = np.ones(lidar_obj.nr,nwl)
        if nhori > 1:
            nelems = EMsubmit.shape
            trunc_nelems = int(nelems[0]/nhori)
            leftover = nelems[0] % nhori
            EMsubmit = EMsubmit[0:nelems[0]-leftover]
            print("The shape of EMsubmit is ",EMsubmit.shape)
            EMsubmit = np.mean( EMsubmit.reshape(-1,nhori), axis=1 )*1e-6
            print("The shape of EMsubmit is ",EMsubmit.shape)
        samp_chan = correct_raw_counts(samp_chan,EMsubmit,np.arange(0,nb*vrZ,vrZ),
                                       lidar_obj.nr,nb,bg_st_bin,bg_ed_bin,
                                       opt)
    # Store samp_chan for possible profile plotting (samp_chan_mask set in here!)                                      
    lidar_obj.samp_chan = samp_chan   
    # Plot only data of current angle
    ang_mask_opt = 'nogaps'
    if ang_list.current() == 0: 
        ang_mask_opt = 'no_mask'
    elif gaps.get() == 1:
        ang_mask_opt = 'gaps'
    else:
        ang_mask_opt = 'nogaps'
    print('------>',ang_list.current(),gaps.get(),ang_mask_opt)                       
    lidar_obj.apply_ONA_mask(ang_list.current()-1,ang_mask_opt)                                                                       
    # Manipulate samp_chan to get correct image orientation
    samp_chan = lidar_obj.samp_chan
    samp_chan = samp_chan.transpose()
    if pointing_dir == "Up":
        samp_chan = np.flipud(samp_chan) # reverse the order in columns (not rows)

    # Plot the data ---

    yax_min = float(ylim_min.get().strip())
    yax_max = float(ylim_max.get().strip()) 
    yax_lims = [yax_min,yax_max]   
    color_bar_min = float(cbar_min.get().strip())
    color_bar_max = float(cbar_max.get().strip())
    if (palt.get() != 'N') & (pointing_dir == 'Down'):
        z = float(palt.get()) - np.arange(0,lidar_obj.nb*lidar_obj.dz,lidar_obj.dz)
    if opt == 'NRB' and color_bar_max < 1e3: 
        print('Setting color bar max to default NRB value.')
        cbar_max.set(str(CBar_max_NRB).strip())
    if opt != 'NRB' and color_bar_max > 5e4:
        print('Setting color bar max to default counts value')
        cbar_max.set(str(CBar_max).strip())
    color_bar_max = float(cbar_max.get().strip())        
    print('color_bar_max is ', color_bar_max)
    ylabel = 'Bin number'
    if yax_type == 'alt': ylabel = 'Altitude (m)'
    title = instrument_name + 'Curtain Plot - channel ' + "".join(chan_sel.get().split()) + \
            ' - ' + wl_str[wl] + ' nm ' + lab_tag
    if xax_type == 'time':
        xax_lims = [ mdates.date2num(lidar_obj.time_arr[lidar_obj.samp_chan_map[0]*nhori]), 
                   mdates.date2num(lidar_obj.time_arr[ 
                   lidar_obj.samp_chan_map[(lidar_obj.nr)-1]*nhori ]) ]
        xlabel = 'UTC'                   
    else:
        xax_lims = [ 0, lidar_obj.nr ]
        xlabel = 'Record number'
    CPlot, im_width, im_Bbox = make_curtain_plot(samp_chan,nb,vrZ,
                                                 z,canvas_ctrl, color_bar_min,
                                                 color_bar_max, xlabel, ylabel,
                                                 title, yax_type, yax_lims, xax_type, xax_lims,
                                                 mpl_win.get())
    im2 = 'current_curtain.png'
    photos[1] = ImageTk.PhotoImage(Image.open(im2))
    canvas.configure(width=photos[1].width(),height=photos[1].height())
    canvas.create_image(0,0,image=photos[1],anchor=NW)
    canvas_ctrl.pix2prof = lidar_obj.nr/im_width
    
    # Update img_bounds list ---
    
    img_bounds =  [ im_Bbox[0,0],im_Bbox[0,1],im_Bbox[1,0],im_Bbox[1,1] ] 
    for i in range(0,4): 
        img_bounds[i] = int(img_bounds[i])
    canvas_ctrl.img_bounds = img_bounds  
    print('img_bounds: ',img_bounds)      
 
    
def canvas_Lclick(event):
    """What happens when user LEFT clicks on curtain?"""
    
    if canvas_ctrl.var_window_open is True:
        print('Close EM window first please.')
        return
    else:
        if event.x > canvas_ctrl.img_bounds[0] and event.x < canvas_ctrl.img_bounds[2]:
            selected_chan = "".join(chan_sel.get().split())
            prof_num = int( (event.x-canvas_ctrl.img_bounds[0]) * canvas_ctrl.pix2prof ) 
            #print('Prof # ', prof_num, ' of ', lidar_obj.nr)
            if canvas_ctrl.prof_window_open is False:
                # Retreive plot params 
                curtain_type = curt_type.get() # radio button number
                yax_type = yax_drop.get()      # get user-selected axis type
                [opt, lab_tag, wl, ax_lims] = get_plot_parameters_for_current_data_type(curtain_type,
                                               yax_type,get_selected_channels('first'))  
                canvas_ctrl.prof_window_open = True
                fig_title= wl_str[wl]+' nm '+lab_tag+'\n channel '+str(selected_chan).strip()
                canvas_ctrl.create_prof_plot_window(fig2,fig_title)
                pp, = plt.plot( [], [], figure=fig2 )
                plt.title('\n\n') #leave blank lines as space for title
                plt.xlabel(lab_tag)
                ylab = 'Bin number'
                if (yax_type == 'alt'):
                    ylab = 'Altitude (m)'
                plt.ylabel(ylab)
                #plt.tight_layout()
                canvas_ctrl.lines2D_obj = pp
                canvas_ctrl.lines2D_obj.axes.axis(ax_lims)
                # Create protocol for user clicking "X" in prof plot window...
                canvas_ctrl.prof_plot_top.protocol("WM_DELETE_WINDOW",close_prof_plot_window)
            
            datetimeformat = "%Y-%m-%d %H:%M:%S"
            nx0 = lidar_obj.data['meta'].shape
            # The "prof_num*nhori" will account for any averaging of data (any nhori >=1)
            prof_transl = lidar_obj.samp_chan_map[prof_num]*nhori
            prof_plot_title = 'UTC: ' + DT.datetime.strftime(lidar_obj.time_arr[prof_transl],datetimeformat) + \
                              '\nlat: ' + str(lidar_obj.data['meta']['Nav']['GPS_Latitude'][prof_transl])[0:6] + \
                              ', lon: ' + str(lidar_obj.data['meta']['Nav']['GPS_Longitude'][prof_transl])[0:6] + \
                              '\nroll: ' + str(lidar_obj.data['meta']['Nav']['RollAngle'][prof_transl])[0:6] + \
                              ', pitch: ' + str(lidar_obj.data['meta']['Nav']['PitchAngle'][prof_transl])[0:6]                              
            plt.title(prof_plot_title)
            canvas_ctrl.update_any_canvas_ctrl_plot(lidar_obj.samp_chan[prof_num,:],
                                         lidar_obj.z, fig2)
            
def plot_Emon():
    """This function will plot the energy monitor values
       NOTE: It seems figure must be created and closed within
             this function!
    """  
    
    if lidar_obj.ingested is False:
        print('No data have been loaded yet. Please load some data.')
    elif canvas_ctrl.prof_window_open is True:
        print('Close profile plot window first, please.')
    elif canvas_ctrl.var_window_open is True:
        print('Side plot window already open.')        
    else:
        EMs = convert_raw_energy_monitor_values(lidar_obj.data['meta']['Engineering']['LaserEnergyMonitors'],nwl,'CPL',e_flg)
        fig3 = plt.figure(3,figsize=(figW-2,figL-1)) 
        nprofs0 = lidar_obj.counts.shape
        xax_type = xax_drop.get()      # get user-selected axis type
        if xax_type == "time":
            xvals = lidar_obj.time_arr
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            EMp_ax_bounds = [xvals[0],xvals[nprofs0[0]-1]] + EMp_yax_bounds
            EM_xlab = 'UTC'
        else:            
            xvals = np.arange(0,nprofs0[0])
            EMp_ax_bounds = [0,nprofs0[0]] + EMp_yax_bounds
            EM_xlab = 'Record number'
        # 355 as blue, 532 as green, 1064 as red
        EMp1, = plt.plot(xvals, EMs[:,0], 'b', label = '355 nm', figure=fig3)
        EMp2, = plt.plot(xvals, EMs[:,1], 'g', label = '532 nm', figure=fig3)
        EMp3, = plt.plot(xvals, EMs[:,2], 'r', label = '1064 nm',figure=fig3)
        plt.legend( handles = [EMp1, EMp2, EMp3] )
        plt.xlabel(EM_xlab)
        plt.ylabel('Energy (micro-joules per record)')
        plt.title('CPL energy monitors')
        lines2D_obj = EMp1
        lines2D_obj.axes.axis(EMp_ax_bounds)  
        fig3.canvas.draw_idle() 
        canvas_ctrl.create_var_plot_window(fig3,'EM Plot')
        # Create protocol for user clicking "X" in prof plot window...
        canvas_ctrl.var_plot_top.protocol("WM_DELETE_WINDOW",close_var_plot_window)
        plt.savefig(source_dir+'current_EM_plot.png',bbox_inches='tight',pad_inches=CPpad)      
        plt.close(fig3)
        canvas_ctrl.var_window_open = True  
   
def make_gps_plots():
    """ This function will make plots of the GPS data, perhaps
        print out some basic analysis of the data (ie. where's it invalid).
    """
    
    # NOTE:
    # As of 12/18/17, no class has been made for GPS data. The GPS data
    # that are loaded here, are forgotten as soon as this function is
    # exited. Whenever this function is entered, the GPS data are loaded.
    
    if canvas_ctrl.prof_window_open is True:
        print('Close profile plot window first, please.')
    elif canvas_ctrl.var_window_open is True:
        print('Side plot window already open.')        
    else:
        gps_file_list = 'gps_file_list.txt'
        search_str = 'gps*'
        create_a_file_list(gps_file_list,search_str)
        gps_yax_bounds = [-10,10000]
        gps_data = ingest_gps_data(gps_file_list)
        fig3 = plt.figure(3,figsize=(figW-2,figL-1)) 
        ngps = gps_data.shape[0]
        xax_type = xax_drop.get()      # get user-selected axis type
        if xax_type == "time":
            xvals = np.zeros(ngps,dtype=DT.datetime)
            for i in range(0,ngps):
                #pdb.set_trace()
                xvals[i] = weeksecondstoutc(gps_data['GpsWeek'][i]*1.0,
                                    gps_data['Gps_sec'][i],leapsecs )
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            gps_ax_bounds = [xvals[0],xvals[nprofs0[0]-1]] + gps_yax_bounds
            gps_xlab = 'UTC'
        else:            
            xvals = np.arange(0,ngps)
            gps_ax_bounds = [0,ngps] + gps_yax_bounds
            gps_xlab = 'Record number'
        gpsp1, = plt.plot(xvals, gps_data['alt'], 'b', label = 'roll', figure=fig3)
        gpsp2, = plt.plot(xvals, gps_data['lon'], 'g', label = 'pitch', figure=fig3)
        gpsp3, = plt.plot(xvals, gps_data['lat'], 'r', label = 'yaw', figure=fig3)
        plt.legend( handles = [gpsp1] )
        plt.xlabel(gps_xlab)
        plt.ylabel('y')
        plt.title('CPL GPS data')
        lines2D_obj = gpsp1
        lines2D_obj.axes.axis(gps_ax_bounds)  
        fig3.canvas.draw_idle() 
        canvas_ctrl.create_var_plot_window(fig3,'GPS Plot')
        # Create protocol for user clicking "X" in prof plot window...
        canvas_ctrl.var_plot_top.protocol("WM_DELETE_WINDOW",close_var_plot_window)
        plt.savefig(source_dir+'current_GPS_plot.png',bbox_inches='tight',pad_inches=CPpad)      
        plt.close(fig3)
        canvas_ctrl.var_window_open = True  


########################################################################
########################################################################
###################                               ######################
################### PROGRAM EXECUTION BEGINS HERE ######################
###################                               ######################
########################################################################
########################################################################

root = Tk()
root.title(instrument_name+" GUI")
root.bind("<Return>",load_and_plot)

# Set some "global" variables to initial values
# NOTE: It seems you can only hold global data in classes!
# in pixels...[x0,y0,x1,y1]
img_bounds = [0,0,0,0]
lidar_obj = standard_lidar_object(None,None,None,None,None,None,None,None,
       None,None,None,None,None,None,None,None,None,None,None,None,None,
       None,None,None,None,None,)
canvas_ctrl = canvas_control(img_bounds, None, False, None, None, None, 
                             None, False, None, None)    
file_ctrl = file_control([-999])                                                     
fig2 = plt.figure(2,figsize=(pp_figW,pp_figL))
# Save a 1" x 1" blank image for use in computing PPI
fig_ppi = plt.figure(3,figsize=(1,1))
plt.savefig('PPI_ref_img.png')
plt.close(fig_ppi)


########## INITIALIZE THE MAIN WIDGETS HERE ###########################
#                                                                     |
#                                                                     v
# Frame to hold all button on RHS
LHSframe = ttk.Frame(root, borderwidth=2)

# Image frame
imgframe = ttk.Frame(root, width=1350, height=350, borderwidth= 2)
im1 = 'default_GUI_image.png'
photos = [ ImageTk.PhotoImage(Image.open(im1)), None ] # this list stays at main prog. level
canvas_def_W = photos[0].width()
canvas_def_H = photos[0].height()

# canvas to contain curtain plots
canvas = Canvas(imgframe,width=canvas_def_W,height=canvas_def_H)
canvas.create_image(0,0,image=photos[0],anchor=NW)

# Buttons in LHS frame
load_b = Button(LHSframe,text='Load & plot',command=load_and_plot)
savcurt_b = Button(LHSframe,text='Save curtain',command=save_curtain_plot)
Emon_b = Button(LHSframe,text='Plot E mon.',command=plot_Emon)
savwin_b = Button(LHSframe,text='Save window plot',command=save_popup_window_plot)
mpl_win = IntVar()
mpl_chk = Checkbutton(LHSframe, text='Matplotlib', variable=mpl_win)
gps_b = Button(LHSframe, text='GPS plot', command=make_gps_plots)
close_b = Button(LHSframe,text='Exit',command=close_program)

# Text entry boxes in LHS frame
chan_sel_l = ttk.Label(LHSframe, text='Selected channel')
chan_sel = StringVar()
chan_entry = ttk.Entry(LHSframe, width=6, textvariable=chan_sel)
chan_entry.insert(0,'1') # put in a default value

# Other labels and widgets
timeslice_group = LabelFrame(LHSframe,text="Plot files") # Timeslice group
sel_files = StringVar()
files_ent = ttk.Entry(timeslice_group, width=7, textvariable=sel_files)
sel_files.set('1')
ts_switch = IntVar()
ts_all_rad = Radiobutton(timeslice_group,text="All",variable=ts_switch,value=1)
ts_files_rad = Radiobutton(timeslice_group,text="files...",variable=ts_switch,value=0)
ts_switch.set(1)
angle_group = LabelFrame(LHSframe,text="Select angle")    # Angle group
ang_list = ttk.Combobox(angle_group,values=("all"),width=10)
gaps = IntVar()
gaps_cb = Checkbutton(angle_group, text='gaps', variable=gaps)
ylim_group = LabelFrame(LHSframe,text="Y-lims (bins only)")       # Y limit group
ylim_max_l = ttk.Label(ylim_group, text='max')
ylim_max = StringVar()
ylim_max_entry = ttk.Entry(ylim_group, width=7, textvariable=ylim_max)
ylim_max_entry.insert(0,str(maxbin).strip())
ylim_min_l = ttk.Label(ylim_group, text='min')
ylim_min = StringVar()
ylim_min_entry = ttk.Entry(ylim_group, width=7, textvariable=ylim_min)
ylim_min_entry.insert(0,str(minbin).strip())
palt_group = LabelFrame(LHSframe,text="Force platform alt (m)")       # Plane alt group
palt = StringVar()
palt_entry = ttk.Entry(palt_group, width=7, textvariable=palt)
palt_entry.insert(0,'N')
ax_dbox_group = LabelFrame(LHSframe)                     # Axes dropbox group
yax_drop_l = Label(ax_dbox_group,text="Yaxis type")
yax_drop = ttk.Combobox(ax_dbox_group,values=("bins","alt"),width=7)
yax_drop.current(0)
xax_drop_l = Label(ax_dbox_group,text="Xaxis type")
xax_drop = ttk.Combobox(ax_dbox_group,values=("recs","time"),width=7)
xax_drop.current(0)
dtype_radios = LabelFrame(LHSframe,text="Data level")   # Data type radios
curt_type = IntVar()
raw_radio = Radiobutton(dtype_radios, text='raw counts',variable=curt_type,value=1)
bgsub_radio = Radiobutton(dtype_radios, text='BG sub counts',variable=curt_type,value=2)
NRB_radio = Radiobutton(dtype_radios, text='NRB',variable=curt_type,value=3)
curt_type.set(1)
CBar_scale_group = LabelFrame(LHSframe,text="CBar scale") # Color bar scale group
cbar_max_l = ttk.Label(CBar_scale_group, text='max')
cbar_max = StringVar()
cbar_max_entry = ttk.Entry(CBar_scale_group, width=7, textvariable=cbar_max)
cbar_max_entry.insert(0,str(CBar_max).strip())
cbar_min_l = ttk.Label(CBar_scale_group, text='min')
cbar_min = StringVar()
cbar_min_entry = ttk.Entry(CBar_scale_group, width=7, textvariable=cbar_min)
cbar_min_entry.insert(0,str(CBar_min).strip())

############# SET GEOMETERY OF MAIN WIDGETS HERE ######################
#                                                                     |
#                                                                     V

# In grid of root ...
LHSframe.grid(column=0,row=0,stick=(N, W, E, S))
LHSframe.rowconfigure(0,weight=1)
imgframe.grid(column=1, row=0, sticky=(N, W, E, S))
imgframe.columnconfigure(0, weight=1)
imgframe.rowconfigure(0, weight=1)

# In grid of imgframe...
canvas.grid(column=0, row=0, sticky=E)

# In grid of LHSframe...
timeslice_group.grid(column=0, row=0, sticky=W)
angle_group.grid(column=0, row=1, sticky=W)
ylim_group.grid(column=0, row=2, sticky=W)
ax_dbox_group.grid(column=0, row=3, sticky=W)
palt_group.grid(column=0, row=4, sticky=W)
yax_drop_l.grid(column=0, row=5, stick=SW)
yax_drop.grid(column=0, row=6, sticky=NW)
dtype_radios.grid(column=0, row=7, sticky=W)
CBar_scale_group.grid(column=0, row=8, sticky=W)
chan_sel_l.grid(column=0, row=9, sticky=SW)
chan_entry.grid(column=0, row=10, sticky=NW)
load_b.grid(column=0, row=11, sticky=W)
savcurt_b.grid(column=0, row=12, sticky=W)
Emon_b.grid(column=0, row=13, sticky=W)
savwin_b.grid(column=0, row=14, sticky=W)
mpl_chk.grid(column=0, row=15, sticky=W)
gps_b.grid(column=0, row=16, sticky=W)
close_b.grid(column=0, row=17, sticky=W)

# Plane alt group
palt_entry.grid(column=0, row=0, sticky=S)

# Radio buttons grouped in dtype_radios
raw_radio.grid(column=0, row=0, sticky=SW)
bgsub_radio.grid(column=0, row=1, sticky=NW)
NRB_radio.grid(column=0, row=2, sticky=NW)

# All widgets in color bar scale group
cbar_min_l.grid(column=0, row=0, sticky=SW)
cbar_min_entry.grid(column=0, row=1, sticky=NW)
cbar_max_l.grid(column=1, row=0, sticky=W)
cbar_max_entry.grid(column=1, row=1, sticky=W)

# All sorts of widgets grouped in timeslice_group
files_ent.grid(column=0, row=0, columnspan=2, sticky=EW)
ts_all_rad.grid(column=0, row=1, sticky=EW)
ts_files_rad.grid(column=1, row=1,sticky=W)

# All widgets in ax_dbox_group (axes drop box group)
yax_drop_l.grid(column=0, row=0, sticky=EW)
yax_drop.grid(column=0, row=1, sticky=N)
xax_drop_l.grid(column=1, row=0, sticky=EW)
xax_drop.grid(column=1, row=1, sticky=N)

# All widgets in angle group
ang_list.grid(column=0, row=0, sticky=W)
gaps_cb.grid(column=1, row=0, sticky=W)

# All widgets in color bar scale group
ylim_min_l.grid(column=0, row=0, sticky=SW)
ylim_min_entry.grid(column=0, row=1, sticky=NW)
ylim_max_l.grid(column=1, row=0, sticky=W)
ylim_max_entry.grid(column=1, row=1, sticky=W)


############ SET OS EVENT BINDINGS OF MAIN WIDGETS HERE ###############
#                                                                     |
#                                                                     V
canvas.bind("<B1-Motion>",canvas_Lclick)
canvas.bind("<Button-1>",canvas_Lclick)

# The following commands tell the "LHSframe" how to change the shape of its
# columns if the user resizes the window.
for child in LHSframe.winfo_children(): child.grid_configure(padx=4,pady=4)

raw_radio.grid_configure(pady=2)
bgsub_radio.grid_configure(pady=2)
NRB_radio.grid_configure(pady=2)


root.mainloop()

