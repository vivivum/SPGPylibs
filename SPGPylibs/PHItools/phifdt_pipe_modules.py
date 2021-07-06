### phifdt_pipe_modules

import numpy as np 
from .tools import printc,bcolors
from .phi_fits import fits_get

import SPGPylibs.GENtools.plot_lib as plib

def phi_correct_dark(dark_f,data_f,header,data,verbose = None,get_dark = False):
    #-----------------
    # READ AND CORRECT DARK FIELD
    #-----------------

    printc('-->>>>>>> Reading Darks                   ',color=bcolors.OKGREEN)
    printc('          Input should be [y-dim,x-dim].',color=bcolors.OKGREEN)
    printc('          DARK IS DIVIDED by 256.   ',color=bcolors.OKGREEN)

    try:
        dark,h = fits_get(dark_f)
        dark = dark / 256.
    except Exception:
        printc("ERROR, Unable to open darks file: {}",dark_f,color=bcolors.FAIL)

    dark_scale = fits_get(dark_f,scaling = True)
    data_scale = fits_get(data_f,scaling = True)
    if dark_scale["Present"][0] == data_scale["Present"][0]:
        scaling = dark_scale["scaling"][0] / data_scale["scaling"][0]
    else:
        scaling = dark_scale["scaling"][1] / data_scale["scaling"][1] * dark_scale["scaling"][0]

    if scaling != 1:
        printc('          checking scalling and correcting for it in the dark.',dark_scale,data_scale,scaling,color=bcolors.WARNING)
        dark = dark * scaling

    if get_dark:  #for kll
        printc('-->>>>>>> Dark is output in phi_correct_dark()',color=bcolors.OKGREEN)
        return dark, scaling

    printc('-->>>>>>> Correcting dark current.',color=bcolors.OKGREEN)
    PXBEG1  = int(header['PXBEG1']) - 1           
    PXEND1  = int(header['PXEND1']) - 1          
    PXBEG2  = int(header['PXBEG2']) - 1           
    PXEND2  = int(header['PXEND2']) - 1   
    data = data - dark[np.newaxis,np.newaxis,PXBEG2:PXEND2+1,PXBEG1:PXEND1+1]
    data = np.abs(data)

    if verbose:
        plib.show_one(dark,vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Dark',cbarlabel='DN',save=None,cmap='gray')
        plib.show_one(data[0,0,:,:],vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Data after dark',cbarlabel='DN',save=None,cmap='gray')

    return data,dark_scale,data_scale