from HARP.scan_anal import Scan
from HARP.iris_functions import Iris
from HARP.image_processing import Image
import os
import numpy as np 
import matplotlib.pyplot as plt
from functions import Adam_calib
import pandas as pd
import json

# Define where to save the processed scan
save_folder = '/Volumes/KasiaDrive/Processed_Artemis'

#Where to take the scan data from
data_folder = '/Volumes/qolslc/20120009 - Matthews' 
# data_folder =os.path.join(data_folder,'AndorScan-019-20220203-0717 - CaF2_power_rot_lowE')
data_folder = '/Volumes/qolslc/20120009 - Matthews/AndorScan-018-20220202-2334 - CaF2_power_rot'
h5files = [i for i in os.listdir(data_folder) if i != 'Camera.h5']
#
date = '220202'
sample = 'CaF2200',
wl0 = '800long'
set_params = 'HighE'
MCPpos = 80000

save_name = 'Scan_{}_{}um_{}nm_{}_{}MCP.npz'.format(date, sample, wl0, set_params, MCPpos)

power = np.array([4,31, 35, 63, 110, 150,230, 280, 330, 350,  450, 570, 630, 670, 720, 720])*1e-3
iris_pos = np.array([-45, -44, -43, -42, -41, -40,-39, -38, -37, -36, -35, -25, -20, -20, -15, -10 ])

iris = Iris()
iris.specify_calib(iris_positions=iris_pos, powers=power)
iris.specify_params(w0_init=200e-6, f=0.75, wl=800e-9, M2=1, reprate=1000,pulse_duration=35e-15) #specify your params!

#Open the Scan class
scan = Scan()
scan.set_folder(data_folder)
# fill out which are available from 
scan.set_params(iris=None, wedge=0, rotation=None, MCPPos=80000) # here fill out the variables that you are not scanning through. This scan is Rotation and iris so I fill out the other MCPPos and Wedge
scan.set_verlim(100, 30) #this is the vertical cropping you want to do
scan.set_eVlim((1,30)) #this is the energy limits you want to look at


#To add all files or at least multiples uncomment:
# for f in h5files[:3]:
#    try:
#        scan.populate_scan_h5(f, function=Adam_calib)
#    except:
#        print('Failed to add {}'.format(f))


scan.populate_scan_h5(h5files[2], function=Adam_calib)
scan.add_calibration_stage('intensity', iris.get_intensity_TWcm2, 'iris')



scan.save_npz(save_folder,save_name)

