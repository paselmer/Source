# Experimenting with making 3D line plots for CAMAL profiles

from initializations import *
#
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import colors as mcolors
import numpy as np
import matplotlib.pyplot as plt
import pdb

#mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = Axes3D(fig)
#theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
#z = np.linspace(-2, 2, 100)
#r = z**2 + 1
#x = r * np.sin(theta)
#y = r * np.cos(theta)
#ax.plot(x, y, z, label='parametric curve')
#ax.legend()

in_arr_file = out_dir + 'save_prof_arrays_file26.npz'
npzfile = np.load(in_arr_file)
p = npzfile['p'] # [angles X bins]  , bg-sub count profiles
z = npzfile['z'] # [bins]           , bin-altitude array (fix-frame)
Oang = npzfile['Oang']   # [angles] , mean ONA for each prof
XZang = npzfile['XZang'] # [angles] , mean XY angle for each prof
XZang[0:2] = -1.0*XZang[0:2]

# Original-length arrays ...
ne = z.shape[0]
nang = Oang.shape[0]
z0 = z[0]
zbin1 = z[337]
r_offset = z[0] - zbin1
r = np.arange(0,ne)*30.0 - r_offset # range

## Trimmed-length arrays ...
st = 430
ed = 1050
z = z[st:ed]
ne = z.shape[0]
r = r[st:ed]
p = p[:,st:ed]


def lines_3d_plot():
    
    print('running the 3D line funtion') 
    s = 0 # scan angle index
    fc = ['r','g','b','y','k'] # line colors
    # Loop thru all angles
    for s in range(0,nang):
    
        # Loop thru all count values
        x_val = r * np.tan(XZang[s])
        y_val = np.zeros(ne)
        p_map = np.arange(0,ne)
        for c in range(0,300):
            c_mask = (y_val >= c)
            try:
                y_val[c_mask] = p[s,c_mask]
            except IndexError:
                pass
    
        ax.plot(x_val,y_val,z,color=fc[s])
 
    # Use a hack to put legend on 3D plot
    fake_lines = ['','','','','']
    for s in range(0,nang):
        fake_lines[s] =  mpl.lines.Line2D([0],[0],linestyle="none",c=fc[s],
            marker='o')
    ax.legend(fake_lines,['angle: {0:.1f}'.format(XZang[0]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[1]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[2]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[3]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[4]*180.0/np.pi)],numpoints=1)
                
    ax.set_xlim(-5000,5000)
    ax.set_ylim(-10,30)
    ax.set_zlim(-0.5e3,18e3)
    ax.set_xlabel('distance (m)')
    ax.set_ylabel('BG-sub counts')
    ax.set_zlabel('altitude (m)')
        
def poly_3d_plot():
    """ For this one, use swap Z and Y axes """
    
    low_cut = 0.0
    mask_to = 0.0
    print('running the poly 3D function')
    s = 0 # scan angle index
    fc = ['r','g','b','y','k'] # polygon face colors
    # Loop thru all angles
    for s in range(0,nang):
    
        # Loop thru all count values
        x_val = r * np.tan(XZang[s])
        y_val = np.zeros(ne)
        p_map = np.arange(0,ne)
        for c in range(0,300):
            c_mask = (y_val >= c)
            try:
                y_val[c_mask] = p[s,c_mask]
            except IndexError:
                pass
        # Mask out values below cutoff                
        y_val[ y_val < low_cut ] = mask_to
        v = []
        h = 0.0 # reference height
        for k in range(0, x_val.shape[0]-1):
            X = [ x_val[k], x_val[k+1], x_val[k+1], x_val[k] ]
            Y = [ z[k], z[k+1], z[k+1], z[k] ]
            Z = [ y_val[k], y_val[k+1], h, h ]
            #X = [ x_val[k], x_val[k+1], x_val[k+1], x_val[k] ]
            #Y = [ y_val[k], y_val[k+1], h, h ]
            #Z = [ z[k], z[k+1], z[k+1], z[k] ]
            v.append(list(zip(X,Y,Z)))
        poly3dCollection = Poly3DCollection(v,facecolors=fc[s])
        ax.add_collection3d(poly3dCollection)

    # Use a hack to put legend on 3D plot
    fake_lines = ['','','','','']
    for s in range(0,nang):
        fake_lines[s] =  mpl.lines.Line2D([0],[0],linestyle="none",c=fc[s],
            marker='o')
    ax.legend(fake_lines,['angle: {0:.1f}'.format(XZang[0]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[1]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[2]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[3]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[4]*180.0/np.pi)],numpoints=1)
    
    pdb.set_trace()
    ax.set_xlabel("distance (m)")
    ax.set_ylabel("altitude (m)")
    ax.set_zlabel("BG-sub counts")
    ax.set_xlim3d(-5000,5000)
    ax.set_ylim3d(-500,18000)
    ax.set_zlim3d(-10,20)
    
    
    
       
def contour_3d_plot():
    """ Try a 3D contour plot """
    
    print('running the 3D contour funtion') 
    
    
    s = 0 # scan angle index
    # Loop thru all angles
    for s in range(0,nang):
    
        # Loop thru all count values
        x_val = r * np.tan(XZang[s])
        y_val = np.zeros(ne)
        p_map = np.arange(0,ne)
        for c in range(0,300):
            c_mask = (y_val >= c)
            try:
                y_val[c_mask] = p[s,c_mask]
            except IndexError:
                pass
            
        # Try this...
        x_val_2d = np.zeros((ne,ne))
        x_val
                        
        ##x_val_2d = np.column_stack((x_val,x_val))
        #x_val_2d = np.broadcast_to(x_val,(ne,ne))
        #y0 = np.zeros(ne)
        #y2 = np.column_stack((y0,y_val))
        #y_val_2d = np.repeat(y2,ne,axis=1)
        #z_val_2d = np.broadcast_to(z,(ne,ne))
        ##z_val_2d = np.column_stack((z,z))
        #pdb.set_trace()
    
        ax.contour(x_val_2d,y_val_2d,z_val_2d,extend3d=True)
        print('X',x_val_2d.shape)
        print('Y',y_val_2d.shape)
        print('Z',z_val_2d.shape)
        
def scatter_3d_plot():
    
    print('running the scatter plot function') 
    s = 0 # scan angle index
    fc = ['r','g','b','y','k'] # line colors
    # Loop thru all angles
    for s in range(0,nang):
    
        # Loop thru all count values
        x_val = r * np.tan(XZang[s])
        y_val = np.zeros(ne)
        p_map = np.arange(0,ne)
        for c in range(0,300):
            c_mask = (y_val >= c)
            try:
                y_val[c_mask] = p[s,c_mask]
            except IndexError:
                pass
    
        ax.scatter(x_val,y_val,z,color=fc[s])
        ax.set_ylim(-10,30)
 
    # Use a hack to put legend on 3D plot
    fake_lines = ['','','','','']
    for s in range(0,nang):
        fake_lines[s] =  mpl.lines.Line2D([0],[0],linestyle="none",c=fc[s],
            marker='o')
    ax.legend(fake_lines,['angle: {0:.1f}'.format(XZang[0]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[1]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[2]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[3]*180.0/np.pi),
        'angle: {0:.1f}'.format(XZang[4]*180.0/np.pi)],numpoints=1)
                
    ax.set_zlim(-0.5e3,20e3)
    ax.set_xlabel('distance (m)')
    ax.set_ylabel('BG-sub counts')
    ax.set_zlabel('altitude (m)')        
              
lines_3d_plot()                 
#poly_3d_plot()
#contour_3d_plot()
#scatter_3d_plot()


plt.show()
pdb.set_trace()

