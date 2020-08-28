# Aidan Patton, 07/2020

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import scipy
from scipy.interpolate import griddata
import os

# original uninterpolated table
QRatio_vals_l = [0.77928, 1.00727,1.0374,1.03768,1.03096,1.0196,1.0051,0.991946,0.981098,0.96872,0.962636,0.957555,0.953856,0.957146,0.958103,0.963894,0.971631,0.978256,0.983962,0.980827,0.968658,0.940343,0.890848,0.703507,0.762225,0.982418,1.03719,1.03791,1.02987,1.01991,1.00747,0.995097,0.982086,0.974063,0.965585,0.960914,0.957779,0.959508,0.962636,0.965191,0.972025,0.97936,0.983152,0.983059,0.972897,0.950605,0.902662,0.694369,0.677243,0.919279,1.03388,1.0371,1.03271,1.02682,1.01285,0.998873,0.987924,0.97936,0.970079,0.966496,0.964149,0.961771,0.968411,0.97126,0.977662,0.984125,0.987568,0.988711,0.981955,0.963199,0.918561,0.725762,9999,0.861375,1.02716,1.03629,1.03425,1.02829,1.01657,1.00604,0.996556,0.986719,0.978611,0.975847,0.973268,0.972666,0.974758,0.978488,0.985607,0.991777,0.995645,0.994139,0.988603,0.971731,0.890123,9999,9999,0.794453,1.00761,1.03484,1.0355,1.03167,1.02371,1.01554,1.00389,0.99573,0.989259,0.98675,0.983553,0.983715,0.986719,0.990572,0.99417,0.999135,0.999583,0.994502,0.987429,0.969909,0.823671,9999,9999,0.7876,0.881011,1.02565,1.03279,1.03394,1.02987,1.02316,1.01511,1.00815,1.00282,0.998602,0.994773,0.995244,0.999668,1.0013,1.00364,1.00525,1.00354,0.99929,0.989468,0.94327,0.762397,9999,9999,9999,0.777064,0.950567,1.02743,1.03208,1.03382,1.02916,1.02441,1.01956,1.01507,1.00965,1.00608,1.00681,1.00832,1.00969,1.01117,1.00788,1.00721,0.999236,0.975422,0.821934,9999,9999,9999,9999,9999,0.810344,0.994471,1.02692,1.0293,1.03196,1.0292,1.02733,1.02452,1.01989,1.01814,1.01621,1.01668,1.01579,1.0151,1.01213,1.00289,0.984379,0.849106,9999,9999,9999,9999,9999,9999,9999,0.813301,0.964041,1.02192,1.02886,1.02819,1.02924,1.02716,1.02445,1.02392,1.02052,1.01817,1.01578,1.0126,1.00276,0.982248,0.845261,9999,9999,9999,9999,9999,9999,9999,9999,9999,0.783002,0.890439,1.00466,1.01838,1.02332,1.02348,1.02247,1.01743,1.01738,1.01306,1.00708,0.99705,0.930236,0.819749,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,0.715659,0.80924,0.880208,0.95303,0.99461,1.00066,0.996008,0.992896,0.974395,0.902639,0.835006,0.758059,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,0.752661,0.793449,0.810552,0.815054,0.801217,0.773636,9999,9999,9999,9999,9999,9999,9999,9999,9999]
QRatio_vals_interp = np.reshape(QRatio_vals_l, (12,24)).T

Rho_vals = [250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 
            750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750, 750,
            1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250,
            1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750,
            2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250, 2250,
            2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750, 2750,
            3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250, 3250,
            3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 3750, 
            4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250, 4250,
            4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750, 4750,
            5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250, 5250,
            5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750, 5750]

Rho_vals_interp = np.array(Rho_vals)

Z_vals = [-5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750,  
          -5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750,
          -5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750,
          -5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750, 
          -5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750,
          -5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750, 
          -5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750, 
          -5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750, 
          -5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750, 
          -5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750, 
          -5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750, 
          -5750, -5250, -4750, -4250, -3750, -3250, -2750, -2250, -1750, -1250, -750, -250, 250, 750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750, 5250, 5750]

Z_vals_interp = np.array(Z_vals)


# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def roundTraditional(val,digits): # round values as expected
   return round(val+10**(-len(str(val))-1), digits)
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def format(l): # to be used to format final list to have 5 numerals after the decimal like in the original table
    for i in range(len(l)):
        if l[i] < 2:
            l[i] = roundTraditional(l[i], 5)
        else:
            l[i] = int(l[i])
    return l
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def interpolate2d_grid(QRatio_vals_l, Rho_vals_interp, Z_vals_interp): 

    ''' 
    Use scipy's griddata method of interpolation to interpolate the charge ratio table 
    Use a mix of linear interpolation values and original points (using nearest type interpolation) where linear interpolation failed.
    
    Parameters:
    QRatio_vals_l : list
        Original table to interpolate
    Rho_vals_interp : numpy array
        The rho values for which the original charge ratio table has values 
    Z_vals_interp : numpy array
        The z values for which the original charge ratio table has values 

    Returns:
    QRatioTableValues : list
        The interpolated charge ratio table values.
        0th -> 24th indices of return goes from z = -5875 to 5875 at rho = 125, next 25 -> 48 is same z's at rho = 375, and so on.

    '''




    np.set_printoptions(threshold=np.inf) # print more digits
    
    QRatio_vals = np.array(QRatio_vals_l) # structured like: go through z from -6000 to 6000 with fixed rho, then repeat for every rho
    
    rho = np.linspace(125, 5875, 24, endpoint=True)
    z = np.linspace(-5875, 5875, 48, endpoint=True)
    rho_pl, z_pl = np.meshgrid(rho, z) # for pic
    z_pr, rho_pr = np.meshgrid(z, rho) # for printing
    rho_pr = rho_pr.flatten()
    z_pr = z_pr.flatten()
    points_new = (rho_pl, z_pl)

    points = np.array([Rho_vals_interp, Z_vals_interp]) # points already have real data on

    QNew_nearest = griddata(points.T, QRatio_vals, points_new, method='nearest')
    QNew_nearest = QNew_nearest.flatten()

    QNew_linear = griddata(points.T, QRatio_vals, points_new, method='linear')
    QNew_linear = QNew_linear.flatten()

    QRatioTableValues = list(QNew_linear)

    for i in range(len(QRatioTableValues)):

        if np.isnan(QRatioTableValues[i]) or QRatioTableValues[i] > 2 : # remove NaN and values in 7000's that aren't useful (since 9999's are more like dummy values)
            QRatioTableValues[i] = QNew_nearest[i]

        if QRatioTableValues[i] < 2: # round to 5 decimal places
            QRatioTableValues[i] = roundTraditional(QRatioTableValues[i], 5)
        else: # just covers values == 9999
            QRatioTableValues[i] = int(QRatioTableValues[i])
            
    # ----------------------------------------------------
    #print(QRatioTableValues) # FINAL LIST
    # ----------------------------------------------------
    return QRatioTableValues
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def QRatioTable(QRatio_mix, outFile):

    ''' 
    Plot the charge ratio table 

    Parameters:
    QRatio_mix : list
        Charge ratio table values
    outFile : str
        What to save the table pdf as
    
    '''

    for i in range(len(QRatio_mix)):
        if QRatio_mix[i] > 2: # will only be 9999 values
            QRatio_mix[i] = np.nan # doing this so appears white/empty in plotted table
    QRatio_mix = np.array([QRatio_mix])
    QRatio_mix = np.reshape(QRatio_mix, (48,24))
    QRatio_mix = np.flipud(QRatio_mix)

    # Plot the heatmap, customize and label the ticks
    fig, ax = plt.subplots()

    rho = np.linspace(125, 5875, 24, endpoint=True)
    z = np.linspace(-5875, 5875, 48, endpoint=True)

    cmap = plt.get_cmap('hsv_r', 80)  # define the colormaps, make identical to original plot of table for consistency
    newcmp = ListedColormap(cmap(np.linspace(0.25, 1, 20)))

    im = ax.imshow(QRatio_mix, cmap=newcmp, vmin=0.6, vmax=1.1, interpolation='none', aspect='auto') # plot cells
    for ri in range(len(rho)): # plot values in each cell
        for zi in range(len(z)):
            if str(QRatio_mix[zi][ri]) != 'nan':
                ax.text(rho[ri]/6000, (z[zi]+6000)/12000, str(QRatio_mix[zi][ri]), va='center', ha='center', fontname='Helvetica', fontsize=2.5, transform=ax.transAxes)

    # formatting stuff for each axis
    ax.tick_params(direction='in')
    plt.style.use('classic')
    ax.yaxis.set_major_locator(ticker.LinearLocator(7))
    ax.yaxis.set_minor_locator(ticker.LinearLocator(49))
    ax.set_yticklabels(['6000', '4000', '2000', '0', '-2000', '-4000', '-6000'], fontname='Helvetica', fontweight='bold')
    ax.set_ylabel('Z [mm]', fontname='Helvetica', fontweight='bold', labelpad=3)
    ax.yaxis.set_tick_params(which='minor', direction = 'in')

    ax.xaxis.set_major_locator(ticker.LinearLocator(7))
    ax.xaxis.set_minor_locator(ticker.LinearLocator(25))
    ax.set_xticklabels(['0', '1000', '2000', '3000', '4000', '5000', '6000'], fontname='Helvetica', fontweight='bold')
    ax.set_xlabel('Rho [mm]', fontname='Helvetica', fontweight='bold', labelpad=8)
    ax.xaxis.set_tick_params(which='minor', direction = 'in')

    # Add a colour bar 
    cbar = fig.colorbar(ax=ax, cmap=newcmp, mappable=im, orientation='vertical')
    cbar.ax.yaxis.set_major_locator(ticker.LinearLocator(11))
    cbar.ax.set_yticklabels(['0.6', '0.65', '0.7', '0.75', '0.8', '0.85', '0.9', '0.95', '1', '1.05', '1.1'], fontname='Helvetica', fontsize=10, fontweight='bold')
    cbar.set_label('Charge-Ratio', fontname='Helvetica', fontweight='bold')

    #plt.show()
    plt.savefig(outFile)
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////