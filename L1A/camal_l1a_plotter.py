from initializations import *
import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime as DT
from read_routines import *
from lidar import *
from time_conversions import *
from mutable_data_structs import define_CLS_structure
import h5py
import matplotlib.dates as mdates

# Set limits to define subset of whole dataset

z0 = -500.0  # meters
z1 = 20000.0
t0 = DT.datetime(2017,12,7,21,20,0) # times
t1 = DT.datetime(2017,12,7,21,50,0)
wl_choice = 1 # 0=355, 1=532, 2=1064
n_avg = 100  # number of profs to avg. 0 = don't average.
cscale = 1e8 #1e9 # color bar scaling

# Load background-substracted counts file (NRB-style format)

h5_file = L1_dir + 'NRB_CAMAL_17_20171207_nav.hdf5'
h5f = h5py.File(h5_file, 'r')
# Keys to data...
#['DEM_laserspot', 'DEM_laserspot_surftype', 'DEM_nadir', 'DEM_nadir_surftype', 
#'EM', 'ONA', 'PGain', 'bin_alt_array', 'laserspot', 'nav', 'nrb', 
#'num_ff_bins', 'num_recs']

z = np.array(h5f['bin_alt_array'])
nr = np.array(h5f['num_recs'])
nr = nr[0]
nb = np.array(h5f['num_ff_bins'])
nb = nb[0]
c = np.array(h5f['nrb'])
nav = np.array(h5f['nav'],dtype=nav_save_struct)

# Initialize some new variables
t = np.zeros(nr,dtype=DT.datetime)

for i in range(0,nr):
    # Decode time into datetime object array
    str_time = str( nav['UTC'][i,:].view('S26') )[3:29]
    try:
        t[i] = DT.datetime.strptime(str_time,"%Y-%m-%dT%H:%M:%S.%f")
    except:
        t[i] = bad_cls_nav_time_value

# Subset the data

indicies = np.arange(0,nr)
testt0 = t >= t0
testt1 = t <= t1
r0 = indicies[testt0][0]
r1 = indicies[testt1][-1]
bin_nums = np.arange(0,nb)
testz0 = z >= z0
testz1 = z >= z1
b1 = bin_nums[testz0].argmax()
b0 = bin_nums[testz1].argmax()
if wl_choice == 0:
    c = c[4,r0:r1+1,b0:b1+1]
elif wl_choice == 1:
    c = c[2,r0:r1+1,b0:b1+1] + c[3,r0:r1+1,b0:b1+1]
else:
    c = c[0,r0:r1+1,b0:b1+1] + c[1,r0:r1+1,b0:b1+1]
z = z[b0:b1+1]
nav = nav[r0:r1+1]
t = t[r0:r1+1]
print('Shape of c is ',c.shape)
if n_avg > 0:
    c = average_lidar_data(c,n_avg,0)


# Close the file once you've extracted the necessary data.
h5f.close()

tit = '532 nm NRB'
print("Subset start time = ",t[0])
print("Subset end time = ",t[-1])
xlimits = [mdates.date2num(t[0]) , mdates.date2num(t[-1])]
ylimits = [z.shape[0]-1,0]
tit='NRB'
xtit = 'UTC'
xax_type = 'time' # time or recs
ytit = 'altitude (km)'
yax_type = 'alt'  # alt or bins
xsize = 25
ysize = 13
im_file_name = 'ACEPOL_stack.png'

curtain_plot(c.transpose(), nb, 30.0/1e3, z/1e3, 0, cscale, hori_cap, pointing_dir,
  xsize, ysize, CPpad, xtit, ytit, tit,  yax_type,[ylimits[0],ylimits[1]], xax_type,
  [xlimits[0],xlimits[1]], 1, 1, out_dir)
pdb.set_trace()




