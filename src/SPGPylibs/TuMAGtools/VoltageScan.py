"""
.. py:module:: TuMAGtools.VoltageScan
.. module:: utils
        :platform: Unix
        :synopsis: Function to plot the mean value of a series of images as a
        function of the Read HVPS voltage. 

        Command line execution :

            python/python3 VoltageScan.py Path_to_folder -c1 xcoord ycoord -r rad
        
            Params:
                - Path_to_image : Path to a .img file.
                - -c1 xcoord(int) ycoord(int) (Optional) : Coordinates of the center of the box used
                  to compute the mean. Default values : 150 150
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

def voltage_scan(Folder, xc1, yc1, rad):
    
    """   
    Params :
        - Folder : Path to the folder containing the images
        - fp : first image to make the parabola fitting.
        - lp : last point to make the parabola fitting
    """
    
    print('Selected Params:')
    print('c1', xc1, yc1)
    print('rad', rad)
    print('Folder : ', Folder)
    
    All_images = sorted(glob.glob(os.path.join(Folder, '*')))
    
    Intensity = []
    Read_volts = []
    Commanded_volts = []

    for img in All_images:
        
        I, H = read_Tumag(img)

        if H['CameraID'] == 0:
            Intensity.append(I[yc1 - rad : yc1 + rad, xc1 - rad : xc1 + rad])
            Read_volts.append(5000 * ( ( ( 2 * H['EtalonVoltsReading'] ) / 4095 ) - 1 ))
            Commanded_volts.append((-1) ** (int(H['EtalonSign']) + 1) \
                        * (int(H['EtalonDN']) * 4999 / (2 ** 12 - 1)))
            
    print('Images size : ', np.shape(I))
    Intensity  = np.array(Intensity)
    Read_volts = np.array(Read_volts)
    
    # Compute mean value
    Intensity = np.mean(Intensity, axis = (1, 2))

    plt.figure(figsize = (15, 10))
    plt.title('Voltage Scan - ' + str(len(Intensity)) + ' Images' )

    plt.plot(Read_volts, Intensity, color = 'gold', lw =2)
    plt.scatter(Read_volts, Intensity, marker = 'x', c = 'w', s = 100)

    plt.xticks(Read_volts)

    plt.ylabel('Intensity')
    plt.xlabel('Voltage [V]')

    plt.grid(True, color = 'w', alpha = 0.2)
    plt.xticks(rotation=60)
    plt.tight_layout()
    plt.show()


def parser(args, params):
    
    xc1 = params[0]
    yc1 = params[1]
    rad = params[2]

    for ind, ar in enumerate(args):

        if ar == '-c1':
            xc1 = int(args[ind + 1])
            yc1 = int(args[ind + 2])
                        
        if ar == '-r':
            rad = int(args[ind + 1])
                            
    params = [xc1, yc1, rad]
         
    return params
    
if __name__ == "__main__":
    
    # Default Values
    xc1 = 150
    yc1 = 150
    rad = 50
    
    params = [xc1, yc1, rad]
    
    args = sys.argv
    
    [xc1, yc1, rad] = parser(args, params)

    Folder = args[1]
    
    voltage_scan(Folder, xc1, yc1, rad)
    
    
    
    
    







