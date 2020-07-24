""" This code takes the CSV output of make_cls2l1a_map_cpl_l1a_v?.py and
    creates its own CSV mapping the L1A to L2 CPL data.
    
    This code replicates the time-based averaging logic in the 
    'detect_layers_v13_cpl.pro' IDL code used in the CPL L2 process.
    
    Last known date code was last modified: 16 July 2020
"""

import numpy as np
from matplotlib import pyplot as plt
import datetime as DT
import pdb

# --------------> Settings are here <--------------
input_file = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\L1\\CLS2L1A_map_PELIcoe_19_23oct19_cls.csv'
sch = '\\'      # for splitting the above string later by path
output_file_path = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\L1\\'
L1B_offset = 1  # For some reason, Rebecca had a hard-coded offset in L1B code
nhori = 5       # From Level 2 cpl_config
l1b_rate = 1.0  # From Level 2 cpl_config
cutoff_last = 2 # From Level 2 cpl_config
cutoff = 2

with open(input_file,'r') as f_obj:
    raw2L1A = f_obj.readlines()
    
nr_CLS = len(raw2L1A)
CLS_rn = np.zeros(nr_CLS, dtype=np.uint32)
L1A_rn = np.zeros(nr_CLS, dtype=np.uint32)
UTC = np.zeros(nr_CLS, dtype=DT.datetime)
PAlt = np.zeros(nr_CLS, dtype=np.float64)
i = 0
for CLS_rec in raw2L1A:
    split_line = CLS_rec.split(',')
    CLS_rn[i] = int(split_line[0].strip())
    L1A_rn[i] = int(split_line[1].strip())
    try:
        UTC[i] = DT.datetime.strptime(split_line[2].strip(),"%Y-%m-%dT%H:%M:%S.%f")
    except:
        pdb.set_trace()
    PAlt[i] = float(split_line[3].strip())
    i += 1
    
# Apply the L1B_offset
if L1B_offset > 0:
    mask = L1A_rn > L1A_rn[L1B_offset-1]
    CLS_rn = CLS_rn[mask]
    L1A_rn = L1A_rn[mask]
    UTC = UTC[mask]
    PAlt = PAlt[mask]
    
# Convert UTC to Julian Day (JDay)

def djul_day(year, month, day, hour, minute, sec):
    daynum = 0.0
    ldpmon = np.zeros((12,2), dtype=np.uint16)
    ldpmon[:,0] = np.array([0,31,59,90,120,151,181,212,243,273,304,334])
    ldpmon[:,1] = np.array([0,31,60,91,121,152,182,213,244,274,305,335])
    k=0
    if (year % 4) == 0: k=1
    jday = day + ldpmon[(month-1),k]
    daynum = hour/24.0 + minute/1440.0 + sec/86400.0 + jday
    return daynum

# You need the unique L1A record numbers, since they repeat to map to raw
uq_L1A_rn, ui, ncounts = np.unique(L1A_rn, return_index=True, return_counts=True)
djDay = np.zeros(uq_L1A_rn.shape,dtype=np.float64)
for i in range(0,uq_L1A_rn.shape[0]):
    djDay[i] = djul_day(UTC[ui[i]].year, UTC[ui[i]].month, UTC[ui[i]].day, UTC[ui[i]].hour, UTC[ui[i]].minute, UTC[ui[i]].second)    
print('The median # of seconds between records is: {0:6.4f}'.format(np.median(np.diff(djDay))*86400.0))

# ------------------> Begin mimicry of L2 averaging process <------------------

# Compute parameters by which to formulate time-averaging
dt = l1b_rate
expected_dtime = nhori*dt # the expected amount of time in nhori # of profiles
ttol = dt/2.0 # time tolerance. set to the time between 2 profiles as of 5/22/17
daysecs = 3600.0*24.0
print( 'Expected delta time = ', expected_dtime )
#
nx0 = uq_L1A_rn.shape[0]
nhori0 = nhori

i = 0 # averaged profile counter
j = 0 # raw profile counter
already_tried = 0 # Important to avoid potential infinite loop
n_avg = []
L2_rn = np.zeros(uq_L1A_rn.shape[0]+100, dtype=np.int32) - 99
L1A_rn_v2 = np.zeros(uq_L1A_rn.shape[0]+100, dtype=np.int32) - 99
while(j < nx0):
    
    # // This paragraph deals with the leftover profiles at the end of the granule \\
    o = nx0 - j
    if (o < nhori0) or (already_tried == 1):
        if (o < cutoff_last) and (already_tried == 0):
            already_tried = 1 # Is 1 if this end-of-granule block was already tried once
            i = i - 1         # This line and line below will essentially re-do the last...
            j = j - nhori
            nhori = nhori + o # I had nhori0 here, but it should just be nhori
        elif (already_tried == 1) and (o < cutoff_last):
            print(o, ' # of profiles deleted at very end of granule.')
            print('counter at:', j,'| total raw profiles:', nx0)
            j = nx0
            continue
        else:
            nhori = o
    expected_dtime = nhori*dt # expand the expected time window for this last avg'd prof
    #print('nhori set to ', nhori,'::', j, nhori, nx0)
    
    # // Logic to ensure no averaging across time gaps \\    [5/22/17]
    delta_ts = (djDay[j:j+nhori] - djDay[j]) * daysecs
    tol_pl = np.where( delta_ts < (expected_dtime+ttol) )[0]
    ntol_pl = tol_pl.shape[0]
    if (ntol_pl < cutoff):
        print('Less than cutoff profiles within expected time window...',j,ntol_pl,i)
        j = j + ntol_pl
        continue
    n_avg.append(ntol_pl)
    
    # Save the indexes for the mapping
    L2_rn[j:j+nhori] = i
    L1A_rn_v2[j:j+ntol_pl] = np.arange(j,j+ntol_pl,1,dtype=np.int32)
    
    j = j + ntol_pl
    i = i + 1

L2_rn = L2_rn[:j]
L1A_rn_v2 = L1A_rn_v2[:j] + L1A_rn.min() # add offset

# ------------------> End mimicry of L2 averaging process <------------------

# Now write a file that will map CLS to L1A to L2

name_body = input_file.split(sch)[-1][7:]
master_map_file = output_file_path + 'CLS2L2' + name_body
f_obj = open(master_map_file,'w')
f_obj.write('CLS, L1A, L2, L1A_Plane_Alt\n')

for k in range(0,uq_L1A_rn.shape[0]):
    
    L2_space_mask = L1A_rn_v2 == uq_L1A_rn[k]
    
    if L2_space_mask.sum() > 0:
        
        L2_rec_str = str(L2_rn[L2_space_mask][0])      # [0] cuz 1-elem array
        L1A_rec_str = str(L1A_rn_v2[L2_space_mask][0])
        
        for m in range(ui[k],ui[k]+ncounts[k]):
            CLS_rec_str = str(CLS_rn[m])
            PAlt_str = str(PAlt[m])
            outlist = [CLS_rec_str, L1A_rec_str, L2_rec_str, PAlt_str]
            str2write = ','.join(outlist) + '\n'
            f_obj.write(str2write)
            
            
f_obj.close()