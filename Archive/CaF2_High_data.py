## necessary import statements
from HARP.scan_anal import Scan
from HARP.iris_functions import Iris
from HARP.image_processing import Image
import os
from scipy.optimize import curve_fit
import numpy as np 
import matplotlib.pyplot as plt
from functions import Andor_calib
from functions import Adams_wedge_thickness
import warnings
warnings.filterwarnings('ignore')

# # define your iris calibration that was used on that day (check if the hour is correct)
# power = np.array([4, 160, 450, 570, 630, 670, 720, 720, 31, 35, 63, 110, 150, 230, 280, 330,350])*1e-3 #list is in mw as that is what we wrote down but the multiplication is to change it to w
# iris_pos= np.array([-45, -40, -35, -30, -25, -20, -15, -10,-44,-43, -42, -41, -40, -39, -38,-37,-36]) #make sure you have same number of iris positions and powers
# iris = Iris()
# iris.specify_calib(iris_positions=iris_pos, powers=power)
# iris.specify_params(w0_init=100e-6, f=0.75, wl=1800e-9, m2=1, reprate=1000,pulse_duration=15e-15) #specify your params!
data_folder = '/volumes/qolslc/20120009 - matthews/andorscan-018-20220202-2334 - caf2_power_rot'
h5files = [i for i in os.listdir(data_folder) if i != 'camera.h5']

# caf2 = Scan()
# caf2.set_folder(data_folder)
# # fill out which are available from 
# caf2.set_params(iris=None, wedge=-3250, rotation=None, MCPPos=90000) # here fill out the variables that you are not scanning through. This scan is Rotation and iris so I fill out the other MCPPos and Wedge
# caf2.set_verlim(125, 275) #this is the vertical cropping you want to do
# caf2.set_eVlim((7,23)) #this is the energy limits you want to look at
# f = h5files[3]
# # for f in h5files:
# try:
#     caf2.populate_scan_h5(f, function=Andor_calib)
# except:
#     print('Failed to add {}'.format(f))

# caf2.add_calibration_stage('intensity', iris.get_intensity_TWcm2, 'iris')
# new_folder = os.path.join( '/Volumes/qolslc/20120009 - Matthews/','Processed_CaF2_High' )

# irises= np.unique(caf2.scan_data.iris)
# for i in range(len(irises)):
#     df = caf2.scan_data[caf2.scan_data['iris']==irises[i]]
#     x,y,Z = caf2.return_scan_data('rotation', df = df)
#     np.savez(os.path.join(new_folder,'RotationScan_Irir{}_cycle3.npz'.format(irises[i])), x, y, Z)

# rots = np.unique(caf2.scan_data.rotation)
# for i in range(len(rots)):
#     df = caf2.scan_data[caf2.scan_data['rotation']==rots[i]]
#     x,y,Z = caf2.return_scan_data('iris', df = df)
#     np.savez(os.path.join(new_folder,'IrisScan_Rot{}_cycle3.npz'.format(rots[i])), x, y, Z)

