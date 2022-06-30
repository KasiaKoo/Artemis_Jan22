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
data_folder = '/Volumes/KasiaDrive/Data/20220620'
data_folder =os.path.join(data_folder,'AndorScan-011-20220621-0143 - Th-X Scan 84')
h5files = [i for i in os.listdir(data_folder) if i != 'Camera.h5']

date = '220620'
sample = 'MgO200',
wl0 = '800long'
set_params = 'part2_rotfocus'
MCPpos = 70000

save_name = 'Scan_{}_{}um_{}nm_{}_{}MCP.npz'.format(date, sample, wl0, set_params, MCPpos)

# power = np.array([4,31, 35, 63, 110, 150,230, 280, 330, 350,  450, 570, 630, 670, 720, 720])*1e-3
# iris_pos = np.array([-45, -44, -43, -42, -41, -40,-39, -38, -37, -36, -35, -25, -20, -20, -15, -10 ])

# iris = Iris()
# iris.specify_calib(iris_positions=iris_pos, powers=power)
# iris.specify_params(w0_init=200e-6, f=0.75, wl=800e-9, M2=1, reprate=1000,pulse_duration=35e-15) #specify your params!

#Open the Scan class
scan = Scan()
scan.set_folder(data_folder)
# fill out which are available from 
scan.set_params(iris=84, wedge=0, rotation=None, MCPPos=70000) # here fill out the variables that you are not scanning through. This scan is Rotation and iris so I fill out the other MCPPos and Wedge
scan.set_verlim(0, -1) #this is the vertical cropping you want to do
scan.set_eVlim((3,40)) #this is the energy limits you want to look at


#To add all files or at least multiples uncomment:
# for f in h5files[:3]:
#    try:
#        scan.populate_scan_h5(f, function=Adam_calib)
#    except:
#        print('Failed to add {}'.format(f))


scan.populate_scan_h5(h5files[0], function=Adam_calib)
# scan.add_calibration_stage('intensity', iris.get_intensity_TWcm2, 'iris')
scan.scan_data['sampleX'] = scan.scan_data['sample x']
scan.scan_data['rotation'] = scan.scan_data['sample rotation']

test1 = scan.scan_data.copy()
test1['data'] = test1.apply(lambda row: row.Data.data, axis=1)
test1['eV'] = test1.apply(lambda row: row.Data.eV_axis, axis=1)
test1['ver'] = test1.apply(lambda row: row.Data.ver_lim, axis=1)
test1['wl0'] = test1.apply(lambda row: row.Data.wl0, axis=1)
test1 = test1.drop(['Data'], axis=1)
np.savez(os.path.join(save_folder, save_name), 
            iris = test1['iris'].values, 
            rotation=test1['rotation'].values,
            mcppos=test1['MCP Pos'].values, 
            data=test1['data'].values,
            ver=test1['ver'].values,
            wl0=test1['wl0'].values,
            eV = test1['eV'].values, 
            sampleX = test1['sampleX'].values,
        )
