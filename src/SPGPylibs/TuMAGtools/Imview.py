"""
.. py:module:: TuMAGtools.imview
.. module:: utils
        :platform: Unix
        :synopsis: Function to visualize TuMag (and SCIP) images in .img format.

        Command line execution :

            python/python3 Imview.py "Path_to_image" -c1 xcoord ycoord -r rad
        
            Params:
                - Path_to_image : Path to a .img file.
                - -c1 xcoord(int) ycoord(int) (Optional) : Coordinates of the center of the box used
                  to compute statistics. Default values : Center of image.
                - -r rad(int) (Optional): Half the size of the box used to compute statistics. 
                  Default value of rad : 1/5 the size of the image. 

.. moduleauthor:: Pablo Santamarina  (SPG - IAA) <psanta@iaa.es>
"""
# ============================ IMPORTS ====================================== #

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
# Own libs
from utils_v2 import read_Tumag
 
# ============================= CONFIG ====================================== #

plt.style.use('dark_background')

# =========================================================================== #

def parser(args, params):
    """
    Function used to fix the parameters given in the command line.
    """

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
        
    args = sys.argv
    
    Image = args[1]

    # Read the image
    I, H = read_Tumag(Image)

    # Statistics positions - DEFAULT
    x_stats   = int(H['Roi_x_size'] / 2)
    y_stats   = int(H['Roi_y_size'] / 2)
    rad = int(H['Roi_x_size'] / 5)

    # Read command line arguments
    params = [x_stats, y_stats,  rad]       
    [x_stats, y_stats, rad] = parser(args, params)

    print('HEADER : ')
    print(H)

    # Statistics
    mean = np.mean(I[x_stats - rad : x_stats + rad, y_stats - rad : y_stats + rad])
    std =  np.std(I[ x_stats - rad : x_stats + rad, y_stats - rad : y_stats + rad])

    # Representation of image
    fig, ax = plt.subplots(figsize = (12, 9))
    fig.suptitle('Mean : ' + str(round(mean)) + ' - STD : ' + str(round(std)) + ' \n Min Value: ' + str(round(np.min(I))) + ' - Max Value: ' + str(round(np.max(I))))
    im = ax.imshow(I, cmap = 'gray')
    ax.xaxis.tick_top()
    cbar = plt.colorbar(im, ax = ax)
    cbar_ticks = np.linspace(np.min(I), np.max(I), num=8, endpoint=True)
    cbar_min = np.min(I)
    cbar_max = np.max(I)
    cbar.set_ticks(cbar_ticks)

    # Statistics Square
    ax.plot([x_stats - rad, x_stats - rad], [y_stats - rad, y_stats + rad], ls = '--', c = 'crimson', alpha = 0.4)
    ax.plot([x_stats + rad, x_stats + rad], [y_stats - rad, y_stats + rad], ls = '--', c = 'crimson', alpha = 0.4)
    ax.plot([x_stats - rad, x_stats + rad], [y_stats - rad, y_stats - rad], ls = '--', c = 'crimson', alpha = 0.4)
    ax.plot([x_stats - rad, x_stats + rad], [y_stats + rad, y_stats + rad],  ls = '--', c = 'crimson', alpha = 0.4)

    # Sliders Config
    min_slider = plt.axes([0.05, 0.1, 0.02, 0.8], facecolor='k')
    max_slider = plt.axes([0.1, 0.1, 0.02, 0.8], facecolor='k')
    slider_min = Slider(min_slider, 'MIN', 0, np.max(I), valinit=np.min(I), orientation = 'vertical', color = 'crimson')
    slider_max = Slider(max_slider, 'MAX', 0, np.max(I), valinit=np.max(I), orientation = 'vertical', color = 'crimson')

    def update_min(val):

        global cbar
        cbar.remove()       
        cbar_min = val
        im = ax.imshow(I, cmap = 'gray', vmin = cbar_min, vmax =cbar_max )
        cbar = plt.colorbar(im, ax = ax)
        cbar_ticks = np.linspace(cbar_min, cbar_max, num=8, endpoint=True)
        cbar.set_ticks(cbar_ticks) 
        cbar.draw_all() 
        plt.draw()
        plt.show()

    def update_max(val):

        global cbar
        cbar.remove()
        cbar_max = val
        im = ax.imshow(I, cmap = 'gray', vmin = cbar_min, vmax =cbar_max )
        cbar = plt.colorbar(im, ax = ax)
        cbar_ticks = np.linspace(cbar_min, cbar_max, num=8, endpoint=True)
        cbar.set_ticks(cbar_ticks) 
        cbar.draw_all() 
        plt.draw()
        plt.show()

    slider_min.on_changed(update_min)
    slider_max.on_changed(update_max)

    plt.show()