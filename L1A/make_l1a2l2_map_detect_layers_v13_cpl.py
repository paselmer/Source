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
input_file = 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\Wallops_12\\L1\\CLS2L1A_map_Wallops_12_09sep12_cls.csv'
L1B_offset = 1 # For some reason, Rebecca had a hard-coded offset in L1B code
nhori = 5      # From Level 2 cpl_config
l1b_rate = 1.0 # From Level 2 cpl_config

with open(input_file,'r') as f_obj:
    raw2L1A = f_obj.readlines()
    
nr_CLS = len(raw2L1A)
CLS_rn = np.zeros(nr_CLS, dtype=np.uint32)
L1A_rn = np.zeros(nr_CLS, dtype=np.uint32)
UTC = np.zeros(nr_CLS, dtype=DT.datetime)
i = 0
for CLS_rec in raw2L1A:
    split_line = CLS_rec.split(',')
    CLS_rn[i] = int(split_line[0].strip())
    L1A_rn[i] = int(split_line[1].strip())
    try:
        UTC[i] = DT.datetime.strptime(split_line[2].strip(),"%Y-%m-%dT%H:%M:%S.%f")
    except:
        pdb.set_trace()
    i += 1
    
pdb.set_trace()    

# ------------------> Begin mimicry of L2 averaging process <------------------

# Compute parameters by which to formulate time-averaging
dt = l1b_rate
expected_dtime = nhori*dt # the expected amount of time in nhori # of profiles
ttol = dt/2.0 # time tolerance. set to the time between 2 profiles as of 5/22/17
print( 'Expected delta time = ', expected_dtime )
#
nx0 = L1A_rn.shape[0]
nhori0 = nhori

i = 0 # averaged profile counter
j = 0 # raw profile counter
already_tried = 0 # Important to avoid potential infinite loop
n_avg = []
while(j < nx0):
    
    # // This paragraph deals with the leftover profiles at the end of the granule \\
    o = nx0 - j
    if (o < nhori0) or (already_tried == 1):
        if (o < cutoff_last) and (already_tried == 0):
            already_tried = 1 # Is 1 if this end-of-granule block was already tried once
            i = i - 1         # This line and line below will essentially re-do the last...
            j = j - nhori
            nhoir = nhori + o # I had nhori0 here, but it should just be nhori
        elif (already_tried == 1) and (o < cutoff_last):
            print(o, ' # of profiles deleted at very end of granule.')
            print('counter at:', j,'| total raw profiles:', nx0)
            j = nx0
            continue
        else:
            nhori = o
    expected_dtime = nhori*dt # expand the expected time window for this last avg'd prof
    print('nhori set to ', nhori,'::', j, nhori, nx0)
    
    # // Logic to ensure no averaging across time gaps \\    [5/22/17]
    delta_ts = (djDay[j:j+nhori] - djDay[j]) * daysecs
    tol_pl = np.where( delta_ts < (expected_dtime+ttol) )[0]
    ntol_pl = tol_pl.sum()
    if (ntol_pl < cutoff):
        print('Less than cutoff profiles within expected time window...',j,ntol_pl,i)
        j = j + ntol_pl
        continue
    n_avg.append(ntol_pl)
    
    j = j + ntol_pl
    i = i + 1
