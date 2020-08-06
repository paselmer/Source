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

def make_L1A_to_L2_map(input_file, output_file_path, sch, L2_configs):
    """ Wanted same function as main-level code make_l1a2l2_map_detect...
        but wanted to be able to automatically call this repeatedly.
    """
    
    L1B_offset, nhori, l1b_rate, cutoff_last, cutoff = L2_configs

    with open(input_file,'r') as f_obj:
        raw2L1A = f_obj.readlines()
    
    nr_CLS = len(raw2L1A)
    CLS_rn = np.zeros(nr_CLS, dtype=np.uint32)
    L1A_rn = np.zeros(nr_CLS, dtype=np.uint32)
    UTC = np.zeros(nr_CLS, dtype=DT.datetime)
    PAlt = np.zeros(nr_CLS, dtype=np.float64)
    ONA = np.zeros(nr_CLS, dtype=np.float64)
    E0 = np.zeros(nr_CLS, dtype=np.float64)
    E1 = np.zeros(nr_CLS, dtype=np.float64)
    E2 = np.zeros(nr_CLS, dtype=np.float64)
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
        ONA[i] = float(split_line[4].strip())
        E0[i] = float(split_line[5].strip())
        E1[i] = float(split_line[6].strip())
        E2[i] = float(split_line[7].strip())
        i += 1
    
    # Apply the L1B_offset
    if L1B_offset > 0:
        mask = L1A_rn > L1A_rn[L1B_offset-1]
        CLS_rn = CLS_rn[mask]
        L1A_rn = L1A_rn[mask]
        UTC = UTC[mask]
        PAlt = PAlt[mask]
        ONA = ONA[mask]
        E0 = E0[mask]
        E1 = E1[mask]
        E2 = E2[mask]
    

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
    f_obj.write('CLS, L1A, L2, L1A_Plane_Alt, L1A_ONA, L1A_E0, L1A_E1, L1A_E2\n')

    for k in range(0,uq_L1A_rn.shape[0]):
    
        L2_space_mask = L1A_rn_v2 == uq_L1A_rn[k]
    
        if L2_space_mask.sum() > 0:
        
            L2_rec_str = str(L2_rn[L2_space_mask][0])      # [0] cuz 1-elem array
            L1A_rec_str = str(L1A_rn_v2[L2_space_mask][0])
        
            for m in range(ui[k],ui[k]+ncounts[k]):
                CLS_rec_str = str(CLS_rn[m])
                PAlt_str = str(PAlt[m])
                ONA_str = str(ONA[m])
                E0_str = str(E0[m])
                E1_str = str(E1[m])
                E2_str = str(E2[m])
                outlist = [CLS_rec_str, L1A_rec_str, L2_rec_str, PAlt_str, ONA_str, E0_str, E1_str, E2_str]
                str2write = ','.join(outlist) + '\n'
                f_obj.write(str2write)
            
            
    f_obj.close()



###############################################################################
###############################################################################
#
#                           MAIN FOR AUTOMATION
#                            
###############################################################################
###############################################################################                            
                              
                            
# --------------> Settings are here <--------------


########################### FIREX ##########################################
###############################################################################

sch = '\\'      # for splitting the above string later by path
L1B_offset = 1  # For some reason, Rebecca had a hard-coded offset in L1B code
nhori = 5       # From Level 2 cpl_config
l1b_rate = 1.0  # From Level 2 cpl_config
cutoff_last = 2 # From Level 2 cpl_config
cutoff = 2
L2_configs = [L1B_offset, nhori, l1b_rate, cutoff_last, cutoff]
                   


#input_files = [ 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-901\\CLS2L1A_map_firex_01aug19_cls.csv',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-902\\CLS2L1A_map_firex_02aug19_cls.csv',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-903\\CLS2L1A_map_firex_06aug19_cls.csv',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-904\\CLS2L1A_map_firex_07aug19_cls.csv',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-905\\CLS2L1A_map_firex_08aug19_cls.csv',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-907\\CLS2L1A_map_firex_13aug19_cls.csv',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-908\\CLS2L1A_map_firex_15aug19_cls.csv',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-909\\CLS2L1A_map_firex_16aug19_cls.csv',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-910\\CLS2L1A_map_firex_19aug19_cls.csv',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-911\\CLS2L1A_map_firex_20aug19_cls.csv',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-912\\CLS2L1A_map_firex_21aug19_cls.csv']
#output_file_paths = [ 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-901\\',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-902\\',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-903\\',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-904\\',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-905\\',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-907\\',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-908\\',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-909\\',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-910\\',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-911\\',
                    #'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-912\\']
# Short flight with primarily clear air ==> 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\firex\\19-906\\',                    
input_files = [ 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-012\\CLS2L1A_map_IMPACTS_20_16dec19_cls.csv',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-013\\CLS2L1A_map_IMPACTS_20_15jan20_iwg1.csv',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-014\\CLS2L1A_map_IMPACTS_20_18jan20_cls.csv',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-016\\CLS2L1A_map_IMPACTS_20_01feb20_cls.csv',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-017\\CLS2L1A_map_IMPACTS_20_05feb20_cls.csv',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-018\\CLS2L1A_map_IMPACTS_20_07feb20_cls.csv',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-019\\CLS2L1A_map_IMPACTS_20_23feb20_cls.csv',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-020\\CLS2L1A_map_IMPACTS_20_25feb20_cls.csv',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-022\\CLS2L1A_map_IMPACTS_20_02mar20_cls.csv']
output_file_paths = [ 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-012\\',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-013\\',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-014\\',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-016\\',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-017\\',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-018\\',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-019\\',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-020\\',
                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-022\\']     
# 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-015\\CLS2L1A_map_IMPACTS_20_25jan20_cls.csv',
# 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\IMPACTS_20\\20-021\\CLS2L1A_map_IMPACTS_20_27feb20_cls.csv',                    




for input_file, output_file_path in zip(input_files, output_file_paths):
    make_L1A_to_L2_map(input_file.strip(), output_file_path.strip(), sch, L2_configs)
    print('Done with ',input_file.split(sch)[-1])                            
    
###############################################################################
###############################################################################



# ############################ PELICOE ##########################################
# ###############################################################################

# sch = '\\'      # for splitting the above string later by path
# L1B_offset = 1  # For some reason, Rebecca had a hard-coded offset in L1B code
# nhori = 5       # From Level 2 cpl_config
# l1b_rate = 1.0  # From Level 2 cpl_config
# cutoff_last = 2 # From Level 2 cpl_config
# cutoff = 2
# L2_configs = [L1B_offset, nhori, l1b_rate, cutoff_last, cutoff]


# input_files = [ 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\20-001\\CLS2L1A_map_PELIcoe_19_21oct19_cls.csv',
#                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\20-002\\CLS2L1A_map_PELIcoe_19_22oct19_cls.csv',
#                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\20-003\\CLS2L1A_map_PELIcoe_19_23oct19_cls.csv',
#                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\20-004\\CLS2L1A_map_PELIcoe_19_29oct19_cls.csv',
#                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\20-005\\CLS2L1A_map_PELIcoe_19_30oct19_cls.csv']
# output_file_paths = [ 'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\20-001\\',
#                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\20-002\\',
#                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\20-003\\',
#                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\20-004\\',
#                    'C:\\Users\\pselmer\\Documents\\CPL_stuff\\PELIcoe_19\\20-005\\']


# for input_file, output_file_path in zip(input_files, output_file_paths):
#     make_L1A_to_L2_map(input_file.strip(), output_file_path.strip(), sch, L2_configs)
#     print('Done with ',input_file.split(sch)[-1])                            
    
# ###############################################################################
# ###############################################################################    
