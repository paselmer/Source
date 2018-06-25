# NOTE: I don't want this library to rely on "initialization.py" [12/19/17]
from subprocess import check_output
import numpy as np
import pdb
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math as math
from scipy import constants as const
import ctypes
import matplotlib.dates as mdates


def average_lidar_data(a,ne,nd):
    """Average lidar data (or other arrays).
       INPUTS: 
               a  -> the input np array. MUST BE 2 DIMENSIONS!
               ne -> the number of elements to average
               nd -> the dimension # to be averaged (start @0)   
               
       LAST MODIFIED: 9/26/17                           
    """
    
    
    # This code will implicitly assume that the larger dimension is supposed
    # to be the x-dimension, or first dimension. In typical lidar data,
    # this corresponds to the number of records. The user should enter an
    # nd value corresponding to the "records" dimension.
    if nd == 1:
        a = a.transpose()
    elif nd != 0:
        print('Cannot subscript array with '+str(nd))
        print('nd must either be == 1 or 0')
        return None
    
    # Compute the number of elements in the output array.
    # Then initialize the output array.
    nelems = a.shape
    nx0 = nelems[0]
    ny0 = nelems[1]
    nx = int(nx0/ne)
    print(str(nx0 % ne), 'leftover')
    ny = ny0
    avgd_a = np.empty([nx,ny])
    
    i=0
    j=0
    while j < ny:
        i = 0
        while i < nx:
            avgd_a[i,j] = np.mean( a[i*ne:(i+1)*ne,j] )
            i += 1
        j += 1            
    
    if nd == 1:
        avgd_a = avgd_a.transpose()
    
    return avgd_a
    

def compute_solar_background(rc,bg_st_bin,bg_ed_bin):
    """This code will compute the solar background of raw counts
       given a 2D counts array and a background region. It will
       return the background as a 2D array, where values are constant
       in the altitude dimension. The altitude dimension should be
       the second dimension (index 1).
    """
    
    # NOTE: This function assumes array of profs x bins [X x Y]
    bg_reg_subset = rc[:,bg_st_bin:bg_ed_bin]
    bg = np.mean(bg_reg_subset, axis=1)
    return bg
    
    
def correct_raw_counts(rc,E,r,nr,nb,bg_st_bin,bg_ed_bin,opt):
    """This function will produce NRB from raw photon counts"""
    
    # INPUTS:
    #         - rc        = raw counts      [nr,nb]
    #         - E         = energy (Joules) [nr]
    #         - r         = range (meters)  [nb]
    #         - nr        = number of records (scalar-int)
    #         - nb        = number of bins    (scalar-int)
    #         - bg_st_bin = first bin of solar background region (scalar-int)
    #         - bg_ed_bin = last bin of solar background region (scalar-int)
    #         - opt       = string = 'bg_sub' or 'NRB' or 'NRB_no_range' 
    #
    # OUTPUS: 
    #         - NRB = NRB, corrected counts [nr,nb]

    # NOTES:
    # [5/8/18] 'NRB_no_range' opt added. Select if counts are already 
    # background-subtracted and range-corrected.
    
    # Compute the solar background
    
    bg = compute_solar_background(rc,bg_st_bin,bg_ed_bin)  
    
    # Compute NRB

    # Broadcast 1D arrays to 2D so computation is vectorized (read: fast)
    if opt is 'NRB': r = np.broadcast_to(r,(nr,nb))
    E = np.broadcast_to(E,(nb,nr)).transpose()
    bg = np.broadcast_to(bg,(nb,nr)).transpose()
    
    if opt == 'NRB':
        NRB =  ((rc - bg) * r**2) / E 
    elif opt == 'bg_sub':
        NRB = rc - bg
    elif opt == 'NRB_no_range':
        NRB = rc / E
    else:
        print('Invalid option entered. stopping execution in lidar.py')
        pdb.set_trace()
    
    return NRB
    


def convert_raw_energy_monitor_values(raw_EMs,nwl,instr='CAMAL',e_flg=None):
    """As of 10/18/17, this routine is only specific to CAMAL"""

    # INPUTS:
    # raw_EMs -> Raw, untouched, energy monitor values. In array form.
    #            Dimension should be something like [nchans x nwl] but
    #            this could depend on instrument. I plan on put multiple
    #            instruments' conversions in this code.
    # nwl -> The number of wavelengths.
    # instr -> String identifying the instrument.
    # e_flg -> A carry-over from CPL-IDL world which tells this code which
    #          equations should be used on CPL conversions.
    #
    # OUTPUT:
    # EMall -> array containing converted energy, in micro-joules, for all
    #          wavelengths.
    
    if instr == 'CAMAL':
    
        sz = raw_EMs.shape
        nr = sz[0]
    
        EMall = np.zeros((nr,nwl),dtype=np.int32) #nwl is in initializations.py
    
        EMall[:,0] = raw_EMs[:,5].astype(np.int32)*16**0 + \
            raw_EMs[:,6].astype(np.int32)*16**2 +     \
            raw_EMs[:,7].astype(np.int32)*16**4 +     \
            raw_EMs[:,8].astype(np.int32)*16**6 +     \
            raw_EMs[:,9].astype(np.int32)*16**8
        EMall[:,1] = raw_EMs[:,0].astype(np.int32)*16**0 + \
            raw_EMs[:,1].astype(np.int32)*16**2 +     \
            raw_EMs[:,2].astype(np.int32)*16**4 +     \
            raw_EMs[:,3].astype(np.int32)*16**6 +     \
            raw_EMs[:,4].astype(np.int32)*16**8 
        EMall[:,2] = raw_EMs[:,10].astype(np.int32)*16**0 + \
            raw_EMs[:,11].astype(np.int32)*16**2 +     \
            raw_EMs[:,12].astype(np.int32)*16**4 +     \
            raw_EMs[:,13].astype(np.int32)*16**6 +     \
            raw_EMs[:,14].astype(np.int32)*16**8           
        
        # After this line EMall is now of dtype float
        EMall = EMall.astype(np.float32) 
        # As of 1/26/18, the following 3 lines of calculations are obsolete
        #EMall[:,0] = 8.63101e-4 + 2.88557e-10*EMall[:,0] - 1.76057e-18*EMall[:,0]**2 + 0.0*EMall[:,0]**3
        #EMall[:,1] = 0.003932537 + 4.11985e-11*EMall[:,1] - 8.82687e-21*EMall[:,1]**2 + 0.0*EMall[:,1]**3
        #EMall[:,2] = 0.012758984 + 2.45103e-10*EMall[:,2] - 1.17876e-19*EMall[:,2]**2 + 0.0*EMall[:,2]**3
        # As of 1/26/18, Andrew Kupchock says convert the energy monitors as follows
        EMall[:,0] = 1.29089e-10*EMall[:,0] + 0.024267116
        EMall[:,1] = 2.71388E-11*EMall[:,1] + 0.005726336
        EMall[:,2] = 1.36981E-10*EMall[:,2] + 0.00206153
        # At this point all EM values are in milli-joules. Convert to micro-joules
        # in return statement.
        EMall = EMall*1e3
    
    elif instr == 'CPL':

        # NOTES:
        # This CPL Energy Monitor conversion section was copied & translated from 
        # Mr. Dennis Hlavka's "engin_corr_v3.pro" code. I even included the "e_flg"
        # and kept the meanings the same.
        #
        # CPL's EMs -> 0=355, 1=532, 2=1064

        sz = raw_EMs.shape
        nr = sz[0]   

        if e_flg == 6:
            # TC4-07 settings . Initially valid 6/19/2006
            emon_t = raw_EMs.astype(np.float64)/500.0
            emon_c = np.copy(emon_t)
            m0_355= 0.0
            m1_355= 0.05
            m2_355= 0.0
            m3_355= 0.0
            emon_c[:,0]= m0_355 + m1_355*emon_t[:,0] + m2_355*emon_t[:,0]**2 + m3_355*emon_t[:,0]**3
            m0_532= -0.25721
            m1_532= 0.038868
            m2_532= -1.533e-6
            m3_532= 0.0
            emon_c[:,1]= m0_532 + m1_532*emon_t[:,1] + m2_532*emon_t[:,1]**2 + m3_532*emon_t[:,1]**3
            m0_1064= -26.806
            m1_1064= 0.19851
            m2_1064= -1.2727e-4
            m3_1064= 0.0
            emon_c[:,2]= m0_1064 + m1_1064*emon_t[:,2] + m2_1064*emon_t[:,2]**2 + m3_1064*emon_t[:,2]**3
        elif e_flg == 7:
            # GloPac10 settings (UAV-CPL) (initially valid 02/16/2010)
            emon_t = raw_EMs.astype(np.float64)/500.0
            emon_c = np.copy(emon_t)
            m0_355= -2.4349
            m1_355= 0.10332
            m2_355= 2.5793e-04
            m3_355= 0.0
            emon_c[:,0]= m0_355 + m1_355*emon_t[:,0] + m2_355*emon_t[:,0]**2 + m3_355*emon_t[:,0]**3
            m0_532= 0.36922
            m1_532= 0.013169
            m2_532= 2.4631e-07
            m3_532= 0.0000
            emon_c[:,1]= m0_532 + m1_532*emon_t[:,1] + m2_532*emon_t[:,1]**2 + m3_532*emon_t[:,1]**3
            m0_1064= 5.2746
            m1_1064= 0.047994
            m2_1064= -7.1374e-06
            m3_1064= 0.000
            emon_c[:,2]= m0_1064 + m1_1064*emon_t[:,2] + m2_1064*emon_t[:,2]**2 + m3_1064*emon_t[:,2]**3
        elif e_flg == 9:
            # new laser box settings (ER2-CPL) 26Feb17	  
            emon_t = raw_EMs.astype(np.float64)/500.0
            emon_c = np.copy(emon_t)
            m0_355= -3.618073187e0
            m1_355= 0.334612636e0
            m2_355= -0.001398674e0
            m3_355= 0.0
            emon_c[:,0]= m0_355 + m1_355*emon_t[:,0] + m2_355*emon_t[:,0]**2 + m3_355*emon_t[:,0]**3
            lowlocs = np.where(emon_c[:,0] < 0.001)
            emon_c[lowlocs[0],0] = 0.001
            highlocs = np.where(emon_c[:,0] > 20.0)
            emon_c[highlocs[0],0] = 20.0
            m0_532= -0.08109323e0
            m1_532= 0.074648283e0
            m2_532= -0.0000079241e0
            m3_532= 0.0
            emon_c[:,1]= m0_532 + m1_532*emon_t[:,1] + m2_532*emon_t[:,1]**2 + m3_532*emon_t[:,1]**3
            m0_1064= 2.172313321e0
            m1_1064= 0.063209038e0
            m2_1064= -0.0000237511e0
            m3_1064= 0.0
            emon_c[:,2]= m0_1064 + m1_1064*emon_t[:,2] + m2_1064*emon_t[:,2]**2 + m3_1064*emon_t[:,2]**3
        elif e_flg == 10:
            # ACT-A, C-130 energy conversion
            # (updated 12/09/2016)  EMON_C in units of micro-Joules
            emon_t = raw_EMs.astype(np.float64)/500.0
            emon_c = np.copy(emon_t)
            m0_355= -1.491632131e0
            m1_355= 1.05116511e-1
            m2_355= 0.0
            m3_355= 0.0
            emon_c[:,0]= m0_355 + m1_355*emon_t[:,0] + m2_355*emon_t[:,0]**2 + m3_355*emon_t[:,0]**3
            #if (emon_c(0) lt .001D) then emon_c[:,0]= 0.001 <- IDL equiv. for ref
            emon_c[np.argwhere(emon_c[:,0] < 0.001),0] = 0.001
            #if (emon_c(0) gt 20.00D) then emon_c(0)= 20.00D <- IDL equiv. for ref
            emon_c[np.argwhere(emon_c[:,0] > 20.0),0] = 20.0
            m0_532= 1.45462874e-1
            m1_532= 7.142194e-3
            m2_532= 3.1427e-8
            m3_532= 0.0
            emon_c[:,1]= m0_532 + m1_532*emon_t[:,1] + m2_532*emon_t[:,1]**2 + m3_532*emon_t[:,1]**3
            m0_1064= 1.090521825e0
            m1_1064= 4.1418361e-2
            m2_1064= -4.2348e-7
            m3_1064= 0.0
            emon_c[:,2]= m0_1064 + m1_1064*emon_t[:,2] + m2_1064*emon_t[:,2]**2 + m3_1064*emon_t[:,2]**3 
        elif e_flg == 12:
            # ER2-CPL settings (initially valid 5/1/2018)
            # Transcribed from email sent by Andrew Kupchock on 4/25/18
            emon_t = raw_EMs.astype(np.float64)/250.0
            emon_c = np.copy(emon_t)
            m0_355= 0.098998716468010
            m1_355= 2.300471133483622
            m2_355= 0.0
            m3_355= 0.0
            emon_c[:,0]= m0_355*emon_t[:,0] + m1_355
            m0_532= 0.091769246955732
            m1_532= -1.671840124168081
            m2_532= 0.0
            m3_532= 0.0
            emon_c[:,1]=  m0_532*emon_t[:,1] + m1_532
            m0_1064= -0.000000000662217
            m1_1064= 0.000001334169843
            m2_1064= -0.000718950067427
            m3_1064= 0.146344091037070
            m4_1064= -3.508899821568728
            emon_c[:,2]= m0_1064*emon_t[:,2]**4 + m1_1064*emon_t[:,2]**3 + m2_1064*emon_t[:,2]**2 + emon_t[:,2]*m3_1064 + m4_1064                      
        else:
             print(str(e_flg)+' is an invalid e_flg value. Stopping in Energy Monitor conversion funtion.')
             pdb.set_trace()
        EMall = emon_c

    else:
        
        print("Instrument name "+instr+" is not valid.")
        print("Stopping in Energy Monitor conversion function.")
        pdb.set_trace()
    
    return EMall   
        


def calculate_off_nadir_angle(phi0,phi0_azi,Y,P,R, **kwargs):
    """ Calculates the off-nadir angle from given inputs.
        Equations used are derived from matrix rotations.
        Fuselage coordinates are rotated into nadir-fixed coordinates.
        Calculation follows what is shown in {G. Cai et al., Unmanned 
        Rotorcraft Systems, Advances in Industrial Control,
        DOI 10.1007/978-0-85729-635-1_2, © Springer-Verlag London 
        Limited 2011} , {Accuracy of Vertical Air Motions from Nadir-
        Viewing Doppler Airborne Radars, Gerald M. Heymsfield, Journal 
        of Atmospheric and Oceanic Technology, December 1989} , and
        {Mapping of Airborn Doppler Radar Data, Wn-Chau Lee et al, 
        Journal of Atmospheric and Oceanic Technology, 1994}.
        !!! ALL ANGLES ARE RADIANS IN CONTEXT OF THIS FUNCTION !!!
        Y, P, and R angles are assumed to be with respect to the
        instrument coordinate system. The instrument and fuselage 
        coordinate systems are assumed to be perfectly aligned.
    """
    
    # ASSUMPTIONS:
    # - All angles (input and output) are in radians.
    # - Y, P, and R angles are with respect to the instrument coordinate
    #   system.
    # - Instrument coordinate origins: +X points perpendicular to right  
    #   of motion, +Y aligns with motion, and +Z points down.
    # - Y is defined as a rotation about the z-axis.
    # - R is defined as a rotation about the y-axis.
    # - P is defined as a rotation about the x-axis.
    # - Rotation angles (Y,P,R) are positive if they move towards a
    #   positive axis origin.
    
    # INPUTS:
    # phi0 -> The off-nadir angle built into the instrument. When the
    #         instrument coordinate system is perfectly level, this is 
    #         off-nadir angle at which the beam points. phi0 should always
    #         be positive. It's not advised that phi0 be greater than or
    #         equal to 90 degrees. If it is, the code will likely run fine;
    #         however, the accuracy of the off-nadir angle calculation
    #         cannot be garanteed. Negative values of phi0 will also produce
    #         unexpected results.
    # phi0_azi -> The counterclockwise angle from the +X instrument axis
    #             (+X defined as perpendicular and to the right of instrument
    #             motion). Counterclockwise from +X obviously defines
    #             a positive phi0_azi.
    # Y -> Yaw, in radians
    # P -> Pitch, in radians
    # R -> Roll, in radians
    
    # OUTPUTS:
    # ONA -> off-nadir angle (in radians)
    
    # Assuming an initial nadir-pointing vector with length equal to 1,
    # compute the x, y, and z components of the intial vector after it's
    # been rotated to the fuselage-ONA (phi0).
    
    x0 = 0.0
    y0 = 0.0
    z0 = 1.0
    
    # Default point angle of beam, broken into components of instrument
    # coordinate system.
    hori_xy_disp = z0 * math.sin(phi0)
    x = hori_xy_disp * math.cos(phi0_azi) # x comp of def. pointing ang.
    y = hori_xy_disp * math.sin(phi0_azi) # y comp of def. pointing ang.
    z = z0 * math.cos(phi0)
    
    # Now compute the nadir-fixed coordinate system components, after
    # rotating through the 3 attitude angles, yaw (Y, about the z-ax), 
    # roll (R, about the y-ax), and pitch (P, about the x-ax), in that 
    # order.
    x1 = x*math.cos(Y)*math.cos(R) + y*(math.sin(Y)*math.cos(P)+math.sin(P)*math.sin(R)*math.cos(Y)) + \
         z*(math.sin(P)*math.sin(Y)-math.sin(R)*math.cos(Y)*math.cos(P))
    y1 = x*(-1*math.sin(Y)*math.cos(R)) + y*(math.cos(Y)*math.cos(P)-math.sin(P)*math.sin(Y)*math.sin(R)) \
         + z*(math.sin(P)*math.cos(Y)+math.cos(P)*math.sin(Y)*math.sin(R))
    z1 = x*math.sin(R) + y*(-1*math.sin(P)*math.cos(R)) + z*math.cos(P)*math.cos(R)
    
    # Now use the dot-product of the initial, nadir-pointing vector, and
    # the rotated vector, to compute the off-nadir angle.
    dotproduct = x1*x0 + y1*y0 +z1*z0
    magvect = math.sqrt(x1**2 + y1**2 +z1**2)
    magnull = math.sqrt(x0**2 + y0**2 +z0**2)
    ONA = math.acos(dotproduct/(magvect*magnull))
    #print('input phi0 is ',phi0*(180./const.pi),' degrees.')
    #print('Off-nadir angle is ',ONA*(180./const.pi),' degrees.')
    
    if ('xz_ang' in kwargs):
        # Compute the angle in the x-z plane (CAMAL purposes in mind)
        dotproduct = x1*x0 + z1*z0
        magvect = math.sqrt(x1**2 +z1**2)
        magnull = math.sqrt(x0**2 +z0**2)
        xz_ang = math.acos(dotproduct/(magvect*magnull))
        return [ONA, xz_ang]
        
    if ('xy_ang' in kwargs):
        # Compute the angle in the x-y plane (CAMAL purposes in mind)

        # The next if block prevents a div-by-zero error.
        # This occurs when R == 0.0
        if x1 == 0.0: 
            if P >= 0.0: xy_ang = const.pi/2.0 # points forward
            if P < 0.0: xy_ang = -1.0*(const.pi/2.0) # points behind
        else:
            xy_ang = np.arctan(y1/x1)
            if x1 < 0.0:
                xy_ang = xy_ang + const.pi
            elif (x1 > 0.0) and (y1 < 0.0):
                xy_ang = xy_ang + 2.0*const.pi
            elif (x1 > 0.0) and (y1 > 0.0):
                pass

        return [ONA, xy_ang]
    
    return ONA

def laser_spot_aircraft(alt0,lat0,lon0,azi,altB,ONA,H):
    """ This code will calculate the lat/lon coordinates at the point
        where the laser beam hits a specified altitude.
    """
    
    # INPUTS:
    #
    # alt0 -> The altitude of the aircraft (m)
    # lat0 -> The nadir latitude in radians
    # lon0 -> The nadir longitude in radians
    # azi -> The azimuth angle measured in the horizontal reference
    #        counterclockwise from X origin (radians)
    # altB -> The altitude at which you'd like the lat/lon coordinates (m)
    # ONA -> The off-nadir angle (radians)
    # H -> The aircraft heading (radians)
    
    
    # OUTPUT:
    #
    # A list of [lat, lon]
    
    # DESCRIPTION & NOTES:
    #
    # This code assumes the altitude (of the aircraft) is close enough to
    # the Earth, that Earth's curvature does not have to be considered.
    # There are two basic conversions used: degrees_lat-to_km and
    # degrees_lon-to_km, the fomer is constant, the latter varies with
    # latitude (distance between lines of lon decreases as you move away
    # from the equator).
    
    # Calculate horizontal distance between point at alt0 and point at altB
    
    # Define equitorial and polar radii of Earth
    a = 6378.137e3
    b = 6356.752e3
    
    D = (alt0 - altB) * np.tan(ONA)
    
    B = azi - H
    
    dx = D * np.sin(B)
    dy = D * np.cos(B)
    
    Clat = 111.0e3              # meters per degree
    Clat = Clat * (180.0/np.pi) # meters per radian
    
    R = math.sqrt(  ( (a**2*np.cos(lat0))**2 + (b**2*np.sin(lat0))**2 ) / 
               ( (a*np.cos(lat0))**2 + (b*np.sin(lat0))**2 )  )
               
    A = R * np.cos(lat0) # distance from Earth's N-S axis
    
    Clon = 2.0*np.pi*A*(1.0/(np.pi*2.0)) # meters per radian
    
    latB = lat0 + dx/Clat
    lonB = lon0 + dy/Clon
    
    #print('ONA: ',ONA*(180.0/np.pi))
    #print(lat0*(180.0/np.pi),lon0*(180.0/np.pi))
    #print(latB*(180.0/np.pi),lonB*(180.0/np.pi))
    
    if abs(lat0*(180.0/np.pi)-latB*(180.0/np.pi)) > 0.5:
        print('\n************************************\n')
        print('The laserspot latitude differs from nadir')
        print('latitude by more than half a degree!')
        print(lat0*(180.0/np.pi),latB*(180.0/np.pi))
        print('\n************************************\n')
        
    
    return [latB, lonB]
               
         
def set_fixed_bin_alt_frame(bot_alt,top_alt,bwid,nb,pdir):
    """ Define  a fixed altitude frame, into which you'd like
        to interpolate data.
        Return array of top edge altitudes of each bin.
    """
    
    if pdir == "Down":
        return np.flip(np.arange(bot_alt,top_alt,bwid,dtype=np.float32),0)
    else:
        return np.arange(bot_alt,top_alt,bwid,dtype=np.float32)


def compute_raw_bin_altitudes(nb,pdir,z0,bsize,ONA):
    """ Compute the altitudes of each raw bin (at single record) """
    
    # INPUTS:
    # nb -> the number of bins
    # pdir -> the pointing direction, must equal "Up" or "Down"
    # z0 -> the altitude of the lidar (meters)
    # bsize -> the size of one bin in meters
    # ONA -> the off-nadir angle in radians
    
    # OUTPUTS:
    # z -> altitudes of each bin. altitudes are the top edge of the
    #      bin.
    
    cosfact = math.cos(ONA)      

    if pdir == "Up":
        z = z0 + np.arange(0, nb)*bsize*cosfact
    elif pdir == "Down":
        z = z0 - np.arange(0, nb)*bsize*cosfact
    else:
        z = z0 - np.arange(0, nb*bsize, bsize)*bsize*cosfact
        
    return z        


def put_bins_onto_fixed_frame(ff,af,orig_counts):
    """ "Interp" counts onto fixed frame. The method of interpolation
        used here is not well-suited for much more than lidar data.
        In fact, it's more appropriately called a bin reassignment, not
        an interpolation. Created based off conversations with Dennis.
    """
    
    # *** IMPORTANT NOTE ***
    # If code crashes within this function, specifically in the innermost
    # loop, it probably means that the actual bin alts are outside the 
    # range of the fixed bin alts (as defined in the initializations.py
    # file by you!)
   
    nc = orig_counts.shape[0] # number of channels, presumably

    af_midpoints = (af[:-1] + af[1:])/2.0
    ff_bin_numbers = np.digitize(af_midpoints,ff)
    u, ui, ncounts = np.unique(ff_bin_numbers,return_index=True,return_counts=True)
    new_counts_frme = np.zeros((nc,ff.shape[0]),dtype=np.float32)

    new_counts_frme[:,ff_bin_numbers-1] = orig_counts[:,:ff_bin_numbers.shape[0]]
    places = np.where(ncounts > 1)[0]
    if places.shape[0] > 0:
        for k in range(0,places.shape[0]):
            orig_indx = np.argwhere(ff_bin_numbers == ff_bin_numbers[places[k]])
            for ch in range(0,nc):
                new_counts_frme[ch,ff_bin_numbers[places[k]]] = np.mean(orig_counts[ch,orig_indx])

    return new_counts_frme
    
def put_bins_onto_fixed_frame_C(np_clib,ff,af,orig_counts,nc):
    """ This function performs the same task as the regular
        put_bins_onto_fixed_frame function; however this version
        uses a C function contained within a library to decrease
        run time.
    """
    
    # NOTES:
    #
    # [3/5/18]
    # ff only gets created once by the calling program. ff stands for 
    # "fixed frame" after all. Therefore I've decided to only use numpy.require
    # on it if its flags aren't set properly.
    #
    #
    
    if not ff.flags['ALIGNED'] or not ff.flags['C_CONTIGUOUS'] or ff.dtype != np.float64: 
        ff = np.require(ff,float,['ALIGNED','C_CONTIGUOUS'])
    if not af.flags['ALIGNED'] or not af.flags['C_CONTIGUOUS'] or af.dtype != np.float64:
        af = np.require(af,float,['ALIGNED','C_CONTIGUOUS'])
    if not orig_counts.flags['ALIGNED'] or not orig_counts.flags['C_CONTIGUOUS'] or orig_counts.dtype != np.float64:
        orig_counts = np.require(orig_counts,float,['ALIGNED','C_CONTIGUOUS'])
    #mult = np.zeros(ff.shape,dtype=np.float32)
    #mult = np.require(mult,float,['ALIGNED','C_CONTIGUOUS'])
    new_counts_frme = np.zeros((nc,ff.shape[0]))
    new_counts_frme = np.require(new_counts_frme,float,['ALIGNED','C_CONTIGUOUS'])

    np_clib.rebin_into_fixed_frame_v2.restype = None
    np_clib.rebin_into_fixed_frame_v2.argtypes = [np.ctypeslib.ndpointer(float,
        ndim=1, flags='aligned'), np.ctypeslib.ndpointer(float, ndim=1,
        flags='aligned'), np.ctypeslib.ndpointer(float, ndim=2, flags='aligned'),
        np.ctypeslib.ndpointer(float, ndim=2, flags='aligned, writeable'),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp)]
    np_clib.rebin_into_fixed_frame_v2(ff, af, orig_counts, new_counts_frme, 
        ff.ctypes.strides,ff.ctypes.shape, af.ctypes.strides, af.ctypes.shape,
        orig_counts.ctypes.strides, orig_counts.ctypes.shape,
        new_counts_frme.ctypes.strides, new_counts_frme.ctypes.shape)
        
    #np_clib.rebin_into_fixed_frame.restype = None
    #np_clib.rebin_into_fixed_frame.argtypes = [np.ctypeslib.ndpointer(float,
        #ndim=1, flags='aligned'), np.ctypeslib.ndpointer(float, ndim=1,
        #flags='aligned'), np.ctypeslib.ndpointer(float, ndim=2, flags='aligned'),
        #np.ctypeslib.ndpointer(float, ndim=2, flags='aligned, writeable'),
        #ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        #ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        #ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        #ctypes.POINTER(np.ctypeslib.c_intp), ctypes.POINTER(np.ctypeslib.c_intp),
        #np.ctypeslib.ndpointer(float, ndim=1, flags='aligned')]        
        
    #np_clib.rebin_into_fixed_frame(ff, af, orig_counts, new_counts_frme, 
        #ff.ctypes.strides,ff.ctypes.shape, af.ctypes.strides, af.ctypes.shape,
        #orig_counts.ctypes.strides, orig_counts.ctypes.shape,
        #new_counts_frme.ctypes.strides, new_counts_frme.ctypes.shape,mult)        

    return new_counts_frme
    
    

    
def get_a_color_map():
    """ Define or pick a color map for use with matplotlib """
    
    rainbow = [
    (     0,         0,         0),
    (     0.0500,     0.02745098,    0.1000),
    (     0.1300,     0.02745098,    0.2000),
    (     0.1800,     0.02745098,    0.3200),
    (     0.2300,     0.02745098,    0.3600),
    (     0.2800,     0.02745098,    0.4000),
    (     0.3300,     0.02745098,    0.4400),
    (     0.3800,     0.02745098,    0.4900),
    (     0.3843,     0.02745098,    0.4941), # rev to purple by here
    (     0.3500,     0.02745098,    0.4941),
    (     0.3400,     0.02745098,    0.5000),
    (     0.3100,     0.02745098,    0.5300),
    (     0.2500,     0.02745098,    0.5700),
    (     0.2000,     0.02745098,    0.6200),
    (     0.1500,     0.02745098,    0.7000),
    (     0.1000,     0.02745098,    0.8500),
    (     0.0500,     0.02745098,    0.9400),
    (     0,         0,    1.0000), ###
    (     0,    0.0500,    1.0000),
    (     0,    0.1000,    1.0000),
    (     0,    0.1500,    1.0000),
    (     0,    0.2000,    1.0000),
    (     0,    0.2500,    1.0000),
    (     0,    0.3000,    1.0000),
    (     0,    0.3500,    1.0000),
    (     0,    0.4000,    1.0000),
    (     0,    0.4500,    1.0000),
    (     0,    0.5000,    1.0000),
    (     0,    0.5500,    1.0000),
    (     0,    0.6000,    1.0000),
    (     0,    0.6500,    1.0000),
    (     0,    0.7000,    1.0000),
    (     0,    0.7500,    1.0000),
    (     0,    0.8000,    1.0000),
    (     0,    0.8500,    1.0000),
    (     0,    0.9000,    1.0000),
    (     0,    0.9500,    1.0000),
    (     0,    1.0000,    1.0000),
    (     0,    1.0000,    1.0000),
    (     0,    1.0000,    0.9500),
    (     0,    1.0000,    0.9000),
    (     0,    1.0000,    0.8500),
    (     0,    1.0000,    0.8000),
    (     0,    1.0000,    0.7500),
    (     0,    1.0000,    0.7000),
    (     0,    1.0000,    0.6500),
    (     0,    1.0000,    0.6000),
    (     0,    1.0000,    0.5500),
    (     0,    1.0000,    0.5000),
    (     0,    1.0000,    0.4500),
    (     0,    1.0000,    0.4000),
    (     0,    1.0000,    0.3500),
    (     0,    1.0000,    0.3000),
    (     0,    1.0000,    0.2500),
    (     0,    1.0000,    0.2000),
    (     0,    1.0000,    0.1500),
    (     0,    1.0000,    0.1000),
    (     0,    1.0000,    0.0500),
    (     0,    1.0000,         0),
    (     0,    1.0000,         0),
    (0.0500,    1.0000,         0),
    (0.1000,    1.0000,         0),
    (0.1500,    1.0000,         0),
    (0.2000,    1.0000,         0),
    (0.2500,    1.0000,         0),
    (0.3000,    1.0000,         0),
    (0.3500,    1.0000,         0),
    (0.4000,    1.0000,         0),
    (0.4500,    1.0000,         0),
    (0.5000,    1.0000,         0),
    (0.5500,    1.0000,         0),
    (0.6000,    1.0000,         0),
    (0.6500,    1.0000,         0),
    (0.7000,    1.0000,         0),
    (0.7500,    1.0000,         0),
    (0.8000,    1.0000,         0),
    (0.8500,    1.0000,         0),
    (0.9000,    1.0000,         0),
    (0.9500,    1.0000,         0),
    (1.0000,    1.0000,         0),
    (1.0000,    1.0000,         0),
    (1.0000,    0.9500,         0),
    (1.0000,    0.9000,         0),
    (1.0000,    0.8500,         0),
    (1.0000,    0.8000,         0),
    (1.0000,    0.7500,         0),
    (1.0000,    0.7000,         0),
    (1.0000,    0.6500,         0),
    (1.0000,    0.6000,         0),
    (1.0000,    0.5500,         0),
    (1.0000,    0.5000,         0),
    (1.0000,    0.4500,         0),
    (1.0000,    0.4000,         0),
    (1.0000,    0.3500,         0),
    (1.0000,    0.3000,         0),
    (1.0000,    0.2500,         0),
    (1.0000,    0.2000,         0),
    (1.0000,    0.1500,         0),
    (1.0000,    0.1000,         0),
    (1.0000,    0.0500,         0),
    (1.0000,         0,         0),
    (1.0000,    1.0000,    1.0000)
    ]
    
    cm = LinearSegmentedColormap.from_list('rainbow',rainbow,N=len(rainbow)*2)
    
    return cm
    
    
def curtain_plot(counts_imgarr, nb, vrZ, z, cb_min, cb_max, hori_cap, pointing_dir,
                      figW, figL, CPpad, xlab, ylab, tit, yax, yax_lims, xax, 
                      xax_lims, scale_alt_OofM, mpl_flg, out_dir):
    """ Function that will create a curtain plot
        This function was basically copied from a function called
        make_curtain_plot in the GUI_function library. The name was changed
        and some things were deleted, but that's about the only 
        differences.
    """
    
    # USE THIS FUNCTION FOR IMAGE
    # The variable counts_imgarr is local to this function, and will not
    # mess up anything outside this scope. Therefore it can be sliced and
    # diced here, without crashing anything outside this scope.
    
    # INPUTS:
    #
    # counts_imgarr -> The 2D array to image. [profs x bins]
    # nb            -> The number of bins. (scalar int/float)
    # vrZ           -> Vertical resolution in meters
    # z             -> The y-axis values, corresponding to "bins" dimension of
    #                  counts_imgarr.
    # cb_min        -> The min value for the color bar.
    # cb_max        -> The max value for the color bar.
    # hori_cap      -> # of profs can't exceed this. Done to prevent code from 
    #                  burdening computer by rendering image of a super massive
    #                  gigantic array.
    # pointing_dir  -> Is the lidar looking "Up" or "Down?"
    # figW          -> The figure width
    # figL          -> The figure length
    # CPpad         -> Parameter that controls padding around the plot (inches).
    # xlab          -> xaxis title (string)
    # ylab          -> yaxis title (string)
    # tit           -> Main title (string)
    # yax           -> String identifying type of yaxis, "alt" or "bins"
    # yax_lims      -> [Bin number of lowest alt, Bin number of greatest alt]
    # xax           -> String identifying type of xaxis, "recs" or "time"
    # xax_lims      -> [record # A, record # B]
    # scale_alt_OofM-> The order of magnitude of z's scale (10 m, 1000 m, etc)
    # mpl_flg       -> If this flag == 1, interactive MPL window appears
    # out_dir       -> Directory where image will be saved; "new_curtain.png"
    
    # OUTPUTS:
    # 

    # expand vertical dimension of image by using np.repeat [retired Oct 2017]
    # subsetting array by user-input bins installed [12/6/17]
    counts_imgarr = counts_imgarr[int(min(yax_lims)):int(max(yax_lims)),:]
    if xax != 'time': counts_imgarr = counts_imgarr[:,int(min(xax_lims)):int(max(xax_lims))]
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
        
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right",size="2%",pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.autoscale
    plt.savefig(out_dir+'new_curtain.png',bbox_inches='tight',pad_inches=CPpad)
    if mpl_flg == 1: plt.show()
    plt.close(fig1)
    
    return None  
    
    
def determine_look_angles(angle_arr, scan_pos_uplim, scan_pos_lowlim, scan_pos_bw ):
    """ 
    This function will determine all the look angles within data.
    Originally wittten with CAMAL in mind because CAMAL has
    programmable, continuously variable scan angles.
    
    This function will produce a list of the float values of the
    angles. 
    """
        
    # NOTE: Any angle offset should already be incorporated into 
    #       angle_arr before entering this function
    
    # The histogram parameters here come from the initialization file
    nhbins = int( (scan_pos_uplim - scan_pos_lowlim)/scan_pos_bw )
    dens, edges = np.histogram(angle_arr,
                  bins=nhbins,range=(scan_pos_lowlim,scan_pos_uplim))
    delta = float(scan_pos_uplim - scan_pos_lowlim) / float(nhbins)
    print("delta is ",delta)
    center = edges[:len(dens)] + delta
    mask = dens != 0
    #print("The bin edges are :",edges)
    print("This many bins have a freq of at least 1: ",dens[mask].shape)
    print("These are the bins: ",center[mask])
    hist_bin_centers = center[mask]
    hist_bin_width = delta
    #plt.plot(center,dens)
    #plt.show()
    
    return [hist_bin_width, hist_bin_centers]
      
      
        
def apply_look_angle_mask(cnt_arr,angle_arr,cs,hist_bin_width,nprofs,nhori,opt):
    """ Code to apply mask to science data based on look angle 
    """   
        
    # This maps the full-length dataset to cnt_arr
    cnt_arr_map = np.arange(0,nprofs).astype('uint32')
    cnt_arr_map = cnt_arr_map # could be overwritten below
        
    # Go no farther if no angle is selected
    if (opt == 'no_mask'): return
    # If you make it to this line, a mask is going to be applied
    print("An angle mask is being applied.")
    
    lim1 = cs - hist_bin_width
    lim2 = cs + hist_bin_width
    print( 'limit1 = '+str(lim1).strip()+' limit2 = '+str(lim2).strip() )
    # Make an array mask to mask "where" angle fall within bin
    # First reduce array based on averaging, ang_mask then will be
    # in "averaged space."
    scan_pos_reduced = angle_arr[::nhori]
    scan_pos_reduced = scan_pos_reduced[0:nprofs]
    ang_mask = ( (scan_pos_reduced >= lim1) & 
         (scan_pos_reduced <= lim2)  )
    #pdb.set_trace()
        
    if opt == 'nogaps':
            
        # Here we make a sample (cnt_arr) with no gaps
        # Remember, ang_mask will be in averaged space.
        cnt_arr = cnt_arr[:,:][cnt_arr_map[ang_mask]]
        print('The shape of cnt_arr is ',cnt_arr.shape)
        # Finalize the cnt_arr_map. The numbers within this array
        # should correspond to full-rez dataset indicies.
        cnt_arr_map = cnt_arr_map[ang_mask]
            
    elif opt == 'gaps':
            
        # Here we make a sample (samp chan) with gaps
        cnt_arr[:,:][cnt_arr_map[(ang_mask==False)]] = 0
        
    # Wrap up this method, overwrite nprofs
    cnt_arr_shape = cnt_arr.shape
    print('The shape of cnt_arr is ',cnt_arr_shape) 
    nprofs = cnt_arr_shape[0]
    print('Length of cnt_arr: ',cnt_arr_shape[0])       
    
    return cnt_arr 
    
    
def make_custom_plot(counts_imgarr, nb, vrZ, z, cb_min, cb_max, hori_cap, pointing_dir,
                      figW, figL, CPpad, xlab, ylab, tit, yax, yax_lims, xax, 
                      xax_lims, scale_alt_OofM, mpl_flg, out_dir, savefname,
                      X1,Y1,X1b,Y1b,X2,Y2):
    """ Function intended to be used to experiment and make custom [1/26/18]
        1-time plots. If plots prove valuable for future reproduction,
        then make a new function with those commands.
    """
    
    # USE THIS FUNCTION FOR IMAGE
    # The variable counts_imgarr is local to this function, and will not
    # mess up anything outside this scope. Therefore it can be sliced and
    # diced here, without crashing anything outside this scope.
    
    # INPUTS:
    #
    
    # OUTPUTS:
    # 

    # expand vertical dimension of image by using np.repeat [retired Oct 2017]
    # subsetting array by user-input bins installed [12/6/17]
    counts_imgarr = counts_imgarr[int(min(yax_lims)):int(max(yax_lims)),:]
    counts_imgarr = counts_imgarr[:,int(min(xax_lims)):int(max(xax_lims))]
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
    plt.subplot(311)
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
        
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right",size="2%",pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.autoscale
    
    # Now that you've plotted the curtain, make other subplots...
    
    plt.subplot(312,sharex=ax)
    plt.plot(X1,Y1,'r',marker='x')
    plt.plot(X1b,Y1b,'g',marker='D')
    plt.ylabel('scan angle (deg)')
    plt.subplot(313,sharex=ax)
    plt.plot(X2,Y2,'b',marker='o')
    plt.ylabel('IWG1 roll (deg)')
   
    plt.savefig(out_dir+savefname,bbox_inches='tight',pad_inches=CPpad)
    
    if mpl_flg == 1: plt.show()
    plt.close(fig1)    
    return None

def stacked_curtains_plot(counts_imgarr, nb, vrZ, z, cb_min, cb_max, hori_cap, pointing_dir,
                      figW, figL, CPpad, xlab, ylab, tit, yax, yax_lims, xax, 
                      xax_lims, scale_alt_OofM, mpl_flg, out_dir, savefname,
                      im1,im2,im3):
    """ Function intended to be used to make stacked curtain plots of
        TC4, 09 Aug 07 data.
    """
    
    # INPUTS:
    #
    interp_param = 'none'#'nearest'
    
    # OUTPUTS:
    # 

    from matplotlib.ticker import StrMethodFormatter

    # expand vertical dimension of image by using np.repeat [retired Oct 2017]
    # subsetting array by user-input bins installed [12/6/17]
    counts_imgarr = counts_imgarr[int(min(yax_lims)):int(max(yax_lims)),:]
    if xax != 'time': counts_imgarr = counts_imgarr[:,int(min(xax_lims)):int(max(xax_lims))]
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
    plt.subplot(411)
    plt.title(tit)
    plt.ylabel(ylab)
    cm = get_a_color_map()
    im = plt.imshow(counts_imgarr, cmap=cm, clim=(cb_min,cb_max),
                   interpolation=interp_param,aspect='auto', 
                   extent=xax_lims+yax_lims)
    ax = plt.gca()                 
                   
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
    ytick_lab = ytick_lab[actual_ticks_mask].astype(np.uint32)
    ytick_ind = ytick_ind[actual_ticks_mask]          	    
    plt.yticks(ytick_ind,ytick_lab)
        
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right",size="1%",pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.autoscale
    ax.tick_params(labelbottom='off')
    
    # Now that you've plotted the curtain, make other subplots...
    
    plt.subplot(412, sharex=ax, sharey=ax)
    # Averaged counts ...
    plt.imshow(im1, cmap=cm, clim=(cb_min,cb_max),
                   interpolation=interp_param,aspect='auto', 
                   extent=xax_lims+yax_lims)
    ax1 = plt.gca()
    plt.ylabel(ylab)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right",size="1%",pad=0.05)
    cax.set_axis_off()    
    ax1.tick_params(labelbottom='off') 
 
    # Sub-sampled counts ...
    plt.subplot(413, sharex=ax, sharey=ax)
    plt.imshow(im2, cmap=cm, clim=(cb_min,cb_max),
                   interpolation=interp_param,aspect='auto', 
                   extent=xax_lims+yax_lims)
    ax2 = plt.gca()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right",size="1%",pad=0.05)
    cax.set_axis_off()
    
    # Avergaged sub-sampled counts ...
    plt.subplot(414, sharex=ax, sharey=ax)
    plt.imshow(im3, cmap=cm, clim=(cb_min,cb_max),
                   interpolation=interp_param,aspect='auto', 
                   extent=xax_lims+yax_lims)
    ax3 = plt.gca()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right",size="1%",pad=0.05)
    cax.set_axis_off()
    plt.tight_layout()
   
    plt.savefig(out_dir+savefname,bbox_inches='tight',pad_inches=CPpad)
    
    if mpl_flg == 1: plt.show()
    plt.close(fig1)   
    return None        

# Put classes below this line <-----------------------------------------
