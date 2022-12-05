"""
.. py:module:: TuMAGtools.photon_flux
.. module:: utils
        :platform: Unix
        :synopsis: Function for plotting the mean value as a function of the exposure time.

        Command line execution :

            python/python3 photon_flux.py Path_to_folder -c1 xcoord ycoord -c2 xcoord ycoord -r rad
        
            Params:
                - Path_to_image : Path to a .img file.
                - -c1 xcoord(int) ycoord(int) (Optional) : Coordinates of the center of the box used
                  to compute the mean for camera 1. Default values : 150 150.
                - -c2 xcoord(int) ycoord(int) (Optional) : Coordinates of the center of the box used
                  to compute the mean for camera 2. Default values : 150 150.               
                - -r rad(int) (Optional): Half the size of the box used to compute the mean. 
                  Default value of rad : 50. 

.. moduleauthor:: Pablo Santamarina  (SPG - IAA) <psanta@iaa.es>
"""

# ============================ IMPORTS ====================================== #

# Built-in libs 
import os
import sys
import glob
import numpy as np

import matplotlib.pyplot as plt

# Own libs
from utils_v2 import read_Tumag
 
# ============================= CONFIG ====================================== #

plt.style.use('dark_background')

# =========================================================================== #

def photon_flux(Folder, xc1, yc1, xc2, yc2, rad):
    
    """
    Function that calculates the mean value in counts for a series of images.
    
    Params :
        - Folder : Path to the folder containing the images
        - fp : first image to make the parabola fitting.
        - lp : last point to make the parabola fitting
        
    """
    
    print('Selected Params:')
    print('c1', xc1, yc1)
    print('c2', xc2, yc2)
    print('rad', rad)
    print('Folder : ', Folder)
    
    All_images = sorted(glob.glob(os.path.join(Folder, '*.img')))
    
    Cam1 = []
    Cam2 = []
    
    for img in All_images:
        I, H = read_Tumag(img)
    
        if H['CameraID'] == 0:
            Cam1.append(I[yc1 - rad : yc1 + rad, xc1 - rad : xc1 + rad])
            
        else:
            Cam2.append(I[yc2 - rad : yc2 + rad, xc2 - rad : xc2 + rad])
            
    Cam1 = np.array(Cam1)
    Cam2 = np.array(Cam2)
    
    # Compute mean value
    mean_1 = np.mean(Cam1, axis = (1, 2))
    mean_2 = np.mean(Cam2, axis = (1, 2))

    # Mean per accumulation
    mean_1 /= 16
    mean_2 /= 16

    #ind = np.arange(len(mean_1))
    exp_times = [20, 25, 30, 35, 40, 45, 50]
    
    if len(mean_1) == len(exp_times):
        print('Found', len(mean_1), ' images per camera. Using preasigned exposure times')
        xlabel = 'Exposure Time [ms]'

    else:
        print('Found', len(mean_1), 'images. Expected :', len(exp_times))
        xlabel = 'Nº Images'
        exp_times = np.arange(len(mean_1))

    fig, axs = plt.subplots(1, 2, figsize = (16, 8))
    axs[0].set_title('Camera 1')
    axs[1].set_title('Camera 2')
    
    axs[0].plot(exp_times, mean_1, c = 'gold', ls = '--', lw = 2)
    axs[0].scatter(exp_times, mean_1, marker = 'X', c = 'deeppink', edgecolor = 'w', s = 200, zorder = 10)

    axs[1].plot(exp_times, mean_2, c = 'gold', ls = '--', lw = 2)
    axs[1].scatter(exp_times, mean_2, marker = 'X', c = 'deeppink', edgecolor = 'w', s = 200, zorder = 10)

    axs[0].grid(True, c = 'w', alpha = 0.2)
    axs[1].grid(True, c = 'w', alpha = 0.2)
    axs[0].set_title('Camera 1')
    axs[1].set_title('Camera 2')
    axs[0].set_xlabel(xlabel)
    axs[1].set_xlabel(xlabel)
    axs[0].set_ylabel('Mean value / 16')
    fig.suptitle('Mean value vs Exposure time')

    plt.tight_layout()
    
    plt.show()

def parser(args, params):
    
    xc1 = params[0]
    yc1 = params[1]
    xc2 = params[2]
    yc2 = params[3]
    rad = params[4]

    for ind, ar in enumerate(args):

        if ar == '-c1':
            xc1 = int(args[ind + 1])
            yc1 = int(args[ind + 2])
            
        if ar == '-c2':
            xc2 = int(args[ind + 1])
            yc2 = int(args[ind + 2])
            
        if ar == '-r':
            rad = int(args[ind + 1])
                            
    params = [xc1, yc1, xc2, yc2, rad]
         
    return params
    
if __name__ == "__main__":
    
    # Default Values
    xc1 = 150
    yc1 = 150
    xc2 = 150
    yc2 = 150
    rad = 50
    
    params = [xc1, yc1, xc2, yc2, rad]
    
    args = sys.argv
    
    [xc1, yc1, xc2, yc2, rad] = parser(args, params)

    Folder = args[1]
    
    photon_flux(Folder, xc1, yc1, xc2, yc2, rad)
    
    
    
    
    







