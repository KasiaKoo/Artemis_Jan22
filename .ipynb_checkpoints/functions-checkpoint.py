from tkinter import W
import numpy as np 
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib import cycler
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from mycolorpy import colorlist as mcp
from scipy import ndimage

"""
This is the file that specifies functions unique to the Artemis 009 beamline
"""
c = 300000000.0
h = 6.63e-34

def loading(save_folder, save_name):
    df = pd.DataFrame()
    with np.load(os.path.join(save_folder, save_name), allow_pickle=True) as f:
        for name in ['iris', 'intensity','wedge', 'rotation', 'Lens', 'MCP Pos','data', 'eV', 'ver', 'wl0', 'sampleX', 'wp', 'wp2']:
            try:
                df[name] = f[name] 
            except:
                print(name, 'not added')
                
    return df

def pos_to_off(position_value, popt = [1.26493902e-02,-1.25115854e+03]):
    offset_val = popt[0]*position_value +  popt[1]
    return offset_val

def Andor_calib(pix_no, trans_pos=30000):
    coeff = [-0.07732270146188233,1240/(15*1.55)]#[-5.29194209e-02,  56.6780912]
    offset = 1281
    trans_pos_offset = pos_to_off(trans_pos)
    wl = (np.arange(pix_no)+trans_pos_offset-offset) * coeff[0] + coeff[1]
    return wl 

# def Adam_calib(pix_no, pos):
#     a = -9.32768153e-15*pos+2.52898805e-09
#     b = 1.01644468e-10*pos +2.18827411e-05
#     c = 6.96885173e-07*pos-2.65767535e-01
#     d = -2.78778796e-03*pos+5.28655448e+02
#     return a*np.arange(pix_no)**3 + b*np.arange(pix_no)**2 + c*np.arange(pix_no) +d

def Adam_calib(pix_no, trans_pos):
    p = np.arange(pix_no)
    S = trans_pos
    P0 = -0.607656359704894
    gamma = 0.0126547547726546
    S0 = 30000
    X0 = 80.7242798185754
    Gamma = 0.0231429412719341
    g = 120
    alpha = np.deg2rad(86.9625004329325)
    L = 234.351550449396
    theta = np.deg2rad(-1.06891224023216)
    m = 1 #diffraction order
    P_det = P0 + gamma*(S-S0) +p
    X_det = X0-Gamma*P_det 
    beta_dash = np.arctan(X_det*np.cos(theta)/(L + X_det*np.sin(theta)))
    beta = np.deg2rad(90) - beta_dash
    wl = (np.sin(alpha) - np.sin(beta))/(m*g)
    return (wl*1e-3)*1e9

def d_lambda(pix_no, dp, trans_pos):
    p = pix_no
    S = trans_pos
    P0 = -0.607656359704894
    gamma = 0.0126547547726546
    S0 = 30000
    X0 = 80.7242798185754
    Gamma = 0.0231429412719341
    g = 120
    alpha = np.deg2rad(86.9625004329325)
    L = 234.351550449396
    theta = np.deg2rad(-1.06891224023216)
    m = 1 #diffraction order
    P_det = P0 + gamma*(S-S0) +p
    X_det = X0-Gamma*P_det 
    beta_dash = np.arctan(X_det*np.cos(theta)/(L + X_det*np.sin(theta)))
    beta = np.deg2rad(90) - beta_dash
    dwl =((np.cos(beta))/(m*g))*(Gamma)*((L*np.cos(theta))/((L+X_det*np.sin(theta))**2 +(X_det*np.cos(theta))**2))*dp
    return (dwl*1e-3)*1e9

def Adams_wedge_thickness(Step, T0 = 0.1, Lead = 0.8, Rev = 200,theta=2.3):
    return 2 * (T0 + (Step * Lead/Rev) * np.tan(np.radians(theta)))

def wedgeold_thickness(Step):
    theta = np.deg2rad(11)
    dis = Step+7+0.11
    return dis*np.tan(theta)



def gaus(x, a, x0, sig):
    return a*np.exp(-(x-x0)**2/sig**2)

def add_har_order(scan_df ,stage, energy_lim, vlim, wl, fig, ax, dp=20,minh=1, maxh=30, color='white'):
    scan = scan_df.copy()
    energy_mask = (scan.eV.values[0]>energy_lim[0])*(scan.eV.values[0]<energy_lim[1])
    N = len(scan[stage].unique())
    har_gen = []
    for i in np.arange(minh,maxh)[::2]:
        if i*1240/wl <=energy_lim[1] and i*1240/wl >=energy_lim[0]:
            ax.axvline(i*1240/wl, color=color, alpha=0.3)
            pix_no = np.argmin(abs(scan.eV.values[0][energy_mask]-i*1240/wl))
            d_wl = d_lambda(pix_no, dp, scan['MCPPos'].values[0])
            dE = (i*1240/wl)*(d_wl/(wl/i))
            ax.axvspan(-dE + scan.eV.values[0][energy_mask][pix_no], dE +scan.eV.values[0][energy_mask][pix_no], color=color, alpha=0.1)
            har_gen.append((-dE + scan.eV.values[0][energy_mask][pix_no],dE +scan.eV.values[0][energy_mask][pix_no]))
            ax.text(i*1240/wl, ax.get_ylim()[1], 'H{}'.format(i))
    return har_gen


def plot_linouts(scan_df, stage, energy_lim, vlim, wl, fig, ax, dp=20):
    scan = scan_df.copy()
    energy_mask = (scan.eV.values[0]>energy_lim[0])*(scan.eV.values[0]<energy_lim[1])
    N = len(scan[stage].unique())
    ax.axes.set_prop_cycle(cycler(color = mcp.gen_color(cmap="plasma",n=N+1)))
    l_all = []
    for i in np.arange(N)[scan[stage].unique().argsort()]:
        lineout = scan[scan[stage]==scan[stage].unique()[i]].data.sum()[vlim[0]:vlim[1],energy_mask].mean(axis=0)
        l_all.append(lineout)
        ax.plot(scan.eV.values[0][energy_mask],lineout, label = '%s %.1f' %(stage, scan[stage].unique()[i]), alpha=0.8)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel('Counts')
    ax.set_xlabel('energy [eV]')
    return scan.eV.values[0][energy_mask], l_all


def plot_heatmap(scan_df, stage, energy_lim, vlim, wl, fig, ax, levels = 50, dp=10, norm=LogNorm()):
    scan = scan_df.copy()
    energy_mask = (scan.eV.values[0]>energy_lim[0])*(scan.eV.values[0]<energy_lim[1])
    N = len(scan[stage].unique())
    Z = np.zeros((N, len(scan.eV.values[0][energy_mask])))
    Y = np.ones(N)  
    count=0
    for i in np.arange(N)[scan[stage].unique().argsort()]:
        lineout = scan[scan[stage]==scan[stage].unique()[i]].data.mean()[vlim[0]:vlim[1],energy_mask].mean(axis=0)
        Z[count] = lineout
        Y[count] = scan[stage].unique()[i]
        count += 1
        
    im = ax.pcolormesh(scan.eV.values[0][energy_mask], Y, Z,norm=norm)
    ax.set_ylabel(stage)
    ax.set_xlabel('energy [eV]')
    fig.colorbar(im, ax=ax)
    ax.set_xlim(energy_lim[0], energy_lim[1]) 
    return scan.eV.values[0][energy_mask], Y, Z

def plot_all(scan_df, stage, energy_lim, vlim, wl, fig, ax, dp=10, norm=LogNorm()):
    scan = scan_df.copy()
    energy_mask = (scan.eV.values[0]>energy_lim[0])*(scan.eV.values[0]<energy_lim[1])
    N = len(scan[stage].unique())
    shape = scan.data[0][vlim[0]:vlim[1], energy_mask].shape
    Z = np.zeros((N, shape[0], shape[1]))
    Y = np.zeros(N)
    count = 0
    for i in np.arange(N)[scan[stage].unique().argsort()]:
        Z[count] = scan[scan[stage]==scan[stage].unique()[i]].data.mean()[vlim[0]:vlim[1], energy_mask]
        Y[count] = scan[stage].unique()[i]
        count+=1
    x = scan.eV.values[0][energy_mask]
    y = np.arange(len(Y)*shape[0])
    X_grid, Y_grid = np.meshgrid(x,y)
    im = ax.pcolormesh(X_grid,Y_grid,Z.reshape(N*shape[0], shape[1]),norm=norm)
    ax.set_yticks(y[shape[0]//2::shape[0]])
    ax.set_yticklabels(Y)
    ax.set_ylabel(stage)
    ax.set_xlabel('energy [eV]')
    fig.colorbar(im, ax=ax)
    ax.set_xlim(energy_lim[0], energy_lim[1])
    ax.tick_params(axis='y', which=u'both',length=0)
    return X_grid, Y_grid, Z.reshape(N*shape[0], shape[1])
# defining the colormap cite: Oli and Leszek

cdict =     {'red':       ((0.0,     0.0, 0.0),
           (0.16667, 0.0, 0.0),
           (0.33333, 0.5, 0.5),
           (0.5,     0.0, 0.0),
           (0.66667, 1.0, 1.0),
           (1,       1.0, 1.0)),
          
'green':     ((0.0,     0.0, 0.0),
           (0.16667, 0.0, 0.0),
           (0.33333, 1.0, 1.0), 
           (0.5,     0.5, 0.5),
           (0.66667, 1.0, 1.0),
           (0.83333, 0.0, 0.0),
           (1.0,     1.0, 1.0)),
 
'blue':      ((0.0,     0.0, 0.0),
           (0.33333, 1.0, 1.0),
           (0.5,     0.0, 0.0),
           (0.83333, 0.0, 0.0),
           (1.0,     1.0, 1.0))}

Alps1 = LinearSegmentedColormap('Alps1', cdict, N=10000)

def cmap_map(function, cmap):
    """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
    This routine will break any discontinuous points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(['red','green','blue']):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j,i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector

    return LinearSegmentedColormap('colormap',cdict,1024)


def plot_har_yield(scan_df, stage, energy_lim, vlim, h_lims, wl, fig, ax, dp=20):
    scan = scan_df.copy()
    energy_mask = (scan.eV.values[0]>energy_lim[0])*(scan.eV.values[0]<energy_lim[1])
    N = len(scan[stage].unique())
    ax.axes.set_prop_cycle(cycler(color = mcp.gen_color(cmap="plasma",n=len(h_lims)+2)))
    Z = np.zeros((N, len(h_lims)))
    x = np.zeros(N)
    count=0
    for i in np.arange(N)[scan[stage].unique().argsort()]:
        lineout = scan[scan[stage]==scan[stage].unique()[i]].data.sum()[vlim[0]:vlim[1],energy_mask].mean(axis=0)
        for j in range(len(h_lims)):
            id_min = np.argmin(abs(scan.eV.values[0][energy_mask]-h_lims[j][0]))
            id_max = np.argmin(abs(scan.eV.values[0][energy_mask]-h_lims[j][1]))
            Z[count,j] = lineout[id_min:id_max].sum()
        x[count] = scan[stage].unique()[i]
        count+=1
    
    ax.plot(x, Z/Z[0])
