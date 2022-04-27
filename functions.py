from tkinter import W
import numpy as np 

"""
This is the file that specifies functions unique to the Artemis 009 beamline
"""

def pos_to_off(position_value, popt = [1.26493902e-02,-1.25115854e+03]):
    offset_val = popt[0]*position_value +  popt[1]
    return offset_val

def Andor_calib(pix_no, trans_pos=30000):
    coeff = [-0.07732270146188233,1240/(15*1.55)]#[-5.29194209e-02,  56.6780912]
    offset = 1281
    trans_pos_offset = pos_to_off(trans_pos)
    wl = (np.arange(pix_no)+trans_pos_offset-offset) * coeff[0] + coeff[1]
    return wl 

def Adams_wedge_thickness(Step, T0 = 0.1, Lead = 0.8, Rev = 200,theta=2.3):
    return 2 * (T0 + (Step * Lead/Rev) * np.tan(np.radians(theta)))

def gaus(x, a, x0, sig):
    return a*np.exp(-(x-x0)**2/sig**2)