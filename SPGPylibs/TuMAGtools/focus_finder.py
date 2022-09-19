"""
.. py:module:: TuMAGtools.focus_finder
.. module:: utils
        :platform: Unix
        :synopsis: Function for finding the focus of a series of images.

        Command line execution :

            python/python3 photon_flux.py Path_to_folder -c1 xcoord ycoord -c2 xcoord ycoord -r rad -fp firstpoint -lp lastpoint
        
            Params:
                - Path_to_image : Path to a .img file.
                - -c1 xcoord(int) ycoord(int) (Optional) : Coordinates of the center of the box used
                  to compute the mean for camera 1. Default values : 150 150.
                - -c2 xcoord(int) ycoord(int) (Optional) : Coordinates of the center of the box used
                  to compute the mean for camera 2. Default values : 150 150.               
                - -r rad(int) (Optional): Half the size of the box used to compute the mean. 
                  Default value of rad :50. 
                - -fp firstpoint(int). First point to use for the fitting (cardinal 0-n). Default : 0
                - -lp lastpoint(int). Last point to use for the fitting (cardinal 0-n). Default : -1
.. moduleauthor:: Pablo Santamarina  (SPG - IAA) <psanta@iaa.es>
"""

# ============================ IMPORTS ====================================== #

# Built-in libs 
import os
import sys
import time
import glob
import numpy as np

from shutil import copyfile

# Only needed if Fits_Flag is set to true
from astropy.io import fits

from scipy.optimize import minimize

import matplotlib.pyplot as plt

# Own libs
from utils_v2 import read_Tumag
 
# ============================= CONFIG ====================================== #

plt.style.use('dark_background')

# =========================================================================== #

def focus_finder(Folder, fp, lp, xc1, yc1, xc2, yc2, rad):
    
    """
    Function that calculates the focus of a series of images, prints the
    result in terminal and plots the analysis.
    
    Params :
        - Folder : Path to the folder containing the images
        - fp : first image to make the parabola fitting.
        - lp : last point to make the parabola fitting
        
    """
    
    print('Selected Params:')
    print('fp', fp)
    print('lp', lp)
    print('c1', xc1, yc1)
    print('c2', xc2, yc2)
    print('rad', rad)
    print('Folder : ', Folder)
    
    All_images = sorted(glob.glob(os.path.join(Folder, '*')))
    
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
    
    ctr_1 = np.std(Cam1, axis = (1, 2)) / np.mean(Cam1, axis = (1, 2)) * 100
    ctr_2 = np.std(Cam2, axis = (1, 2)) / np.mean(Cam2, axis = (1, 2)) * 100
    
    ind = np.arange(len(ctr_1))
    
    params_cam1 = np.polyfit(ind[fp :lp], ctr_1[fp :lp], 2)
    params_cam2 = np.polyfit(ind[fp:lp], ctr_2[fp :lp], 2)
    
    parabola1 = np.poly1d(params_cam1)
    parabola2 = np.poly1d(params_cam2)
    
    xx = np.linspace(ind[fp], ind[lp], 100)
    
    def p1_min(x):
        return - parabola1(x)
    
    def p2_min(x):
        return - parabola2(x)
    
    max_c1 = minimize(p1_min, x0 = ind[int(len(ind) / 2)], method = 'Powell').x
    max_c2 = minimize(p2_min, x0 = ind[int(len(ind) / 2)], method = 'Powell').x
    
    print('---------------------------------------------------------------')
    print('Maximum for Camara 1:', max_c1[0])
    print('Maximum for Camara 2:', max_c2[0])
    print('---------------------------------------------------------------')
    
    fig, axs = plt.subplots(2, 2, figsize = (10, 7))
    axs[0, 0].set_title('Camera 1')
    axs[0, 1].set_title('Camera 2')
    axs[0, 0].plot(ind, ctr_1, c = 'crimson', linewidth = 2, ls = '--', zorder = 9)
    axs[0, 1].plot(ind, ctr_2, c = 'crimson', linewidth = 2, ls = '--', zorder = 9)
    axs[0, 0].scatter(ind, ctr_1, c = 'crimson', marker = 'D', s = 100,  label = 'Measures', zorder = 10)
    axs[0, 1].scatter(ind, ctr_2, c = 'crimson', marker = 'D', s = 100, zorder = 10)
    axs[0, 0].plot(xx, parabola1(xx), color = 'dodgerblue', linewidth = 2, label = 'Parabolic fit')
    axs[0, 1].plot(xx, parabola2(xx), color = 'dodgerblue', linewidth = 2)
    
    axs[0, 0].scatter(max_c1, parabola1(max_c1), s = 150, marker = 'x', color = 'w', zorder = 10,label =  'Maximum')
    axs[0, 1].scatter(max_c2, parabola2(max_c2), s = 150, marker = 'x', color = 'w', zorder = 10)
    
    axs[0, 0].legend()
    
    axs[0, 1].yaxis.tick_right()
    axs[0, 0].set_ylabel('Contrast [%]')
    axs[0, 0].grid(True, color = 'w', alpha = 0.2)
    axs[0, 0].set_xticks(ind)
    axs[0, 0].set_xlabel('Nº images')
    axs[0, 1].set_xlabel('Nº images')
    
    axs[0, 1].grid(True, color = 'w', alpha = 0.2)
    axs[0, 1].set_xticks(ind)
    
    im1 = Cam1[round(max_c1[0])]
    im2 = Cam2[round(max_c2[0])]
    
    axs[1, 0].set_title('Best focus Cam1 - Image : ' + str(round(max_c1[0])))
    axs[1, 1].set_title('Best focus Cam2 - Image : ' + str(round(max_c2[0])))
    im = axs[1, 0].imshow(im1, cmap = 'inferno')
    plt.colorbar(im, ax = axs[1, 0])
    im = axs[1, 1].imshow(im2, cmap = 'inferno')
    plt.colorbar(im, ax = axs[1, 1])
    axs[1, 0].set_xticks([])
    axs[1, 0].set_yticks([])
    axs[1, 1].set_xticks([])
    axs[1, 1].set_yticks([])
    
    plt.tight_layout()
    
    plt.show()

def parser(args, params):
    
    xc1 = params[0]
    yc1 = params[1]
    xc2 = params[2]
    yc2 = params[3]
    rad = params[4]
    fp  = params[5]
    lp  = params[6]
    
    for ind, ar in enumerate(args):

        if ar == '-c1':
            xc1 = int(args[ind + 1])
            yc1 = int(args[ind + 2])
            
        if ar == '-c2':
            xc2 = int(args[ind + 1])
            yc2 = int(args[ind + 2])
            
        if ar == '-r':
            rad = int(args[ind + 1])
            
        if ar == '-fp':
            fp = int(args[ind + 1])
            
        if ar == '-lp':
            lp = int(args[ind + 1])
                    
    params = [xc1, yc1, xc2, yc2, rad, fp, lp]
         
    return params
    
if __name__ == "__main__":
    
    # Default Values
    xc1 = 150
    yc1 = 150
    xc2 = 150
    yc2 = 150
    rad = 50
    fp = 0
    lp = -1
    
    params = [xc1, yc1, xc2, yc2, rad, fp, lp]
    
    args = sys.argv
    
    [xc1, yc1, xc2, yc2, rad, fp, lp] = parser(args, params)

    Folder = args[1]
    
    focus_finder(Folder, fp, lp, xc1, yc1, xc2, yc2, rad)
    
    
    
    
    







