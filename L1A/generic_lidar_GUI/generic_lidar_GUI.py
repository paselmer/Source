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
from read_routines import * # no longer imports L1A initializations [10/9/19]
from generic_lidar_GUI_initializations import * # !!NEEDS TO BE LAST!!

filepath = 'C:\\Users\\pselmer\\Documents\\Roscoe\\data\\23Mar18-roscoe\\'

upfile = 'dataup_20180323_191345.data'
downfile = 'datadown_20180323_190845.data'
hk_file = 'hk__20180323_193557.data'

updata = read_in_mcs_data_r(filepath+upfile)
downdata = read_in_mcs_data_r(filepath+upfile)

file_list = 'raw_signal_file_list.txt'
create_a_file_list(file_list,'*datadown*',raw_dir)

pdb.set_trace()




