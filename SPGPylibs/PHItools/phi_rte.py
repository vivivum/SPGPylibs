#=============================================================================
# Project: SoPHI
# File:    phi_rte.py
# Author:  David Orozco Suárez (orozco@iaa.es)
# Contributors: 
#-----------------------------------------------------------------------------
# Description: programs for accesing data and fits files
#-----------------------------------------------------------------------------
from .tools import *
from .cmilos import pymilos
import subprocess
import numpy as np

@timeit
def phi_rte(data,wave_axis,rte,cmilos = None):
    ''' For the moment this is just isolated from the main pipeline'''

    # Inv Iterations                                   = 15
    # Initial Model Continuum Absoprtion               = 12.0000
    # Initial Model Vector Magnetic Field Strength     = 1200.00000
    # Initial Model Line-Of-Sight Velocity             = 0.05000
    # Initial Model Doppler Width Of Line              = 0.05000
    # Initial Model Damping Parameter                  = 0.05000
    # Initial Model Vector Magnetic Field Inclination  = 170.00000
    # Initial Model Vector Magnetic Field Azimuth      = 25.00000
    # Initial Model Source Function Ordinate At Origin = 0.30000
    # Initial Model Source Function Slope              = 0.80000
    # Wavelength Setting Initial Sampling Wavelength   = 6173.20117
    # Wavelength Setting Spectral Line Step            = 0.07000
    # Wavelength Setting Wavelength Of Line Continuum  = 6173.04117 (blue) 6173.64117 (red)

    l,p,x,y = data.shape
    print(l,p,x,y)

    filename = 'dummy_in.txt'
    with open(filename,"w") as f:
        for i in range(x):
            for j in range(y):
                for k in range(l):
                    f.write('%e %e %e %e %e \n' % (wave_axis[k],data[k,0,j,i],data[k,1,j,i],data[k,2,j,i],data[k,3,j,i]))

    printc('  ---- >>>>> Inverting data.... ',color=bcolors.OKGREEN)

    if rte == 'RTE':
        rte_on = subprocess.call(cmilos+" 6 15 0 0 dummy_in.txt  >  dummy_out.txt",shell=True)
    if rte == 'CE':
        rte_on = subprocess.call(cmilos+" 6 15 2 0 dummy_in.txt  >  dummy_out.txt",shell=True)
    if rte == 'CE+RTE':
        rte_on = subprocess.call(cmilos+" 6 15 1 0 dummy_in.txt  >  dummy_out.txt",shell=True)

    print(rte_on)
    printc('  ---- >>>>> Finishing.... ',color=bcolors.OKGREEN)

    printc('  ---- >>>>> Reading results.... ',color=bcolors.OKGREEN)
    res = np.loadtxt('dummy_out.txt')

    del_dummy = subprocess.call("rm dummy_in.txt",shell=True)
    print(del_dummy)
    del_dummy = subprocess.call("rm dummy_out.txt",shell=True)
    print(del_dummy)

    npixels = res.shape[0]/12.
    print(npixels)
    print(npixels/x)
    result = np.zeros((12,y*x)).astype(float)
    for i in range(y*x):
        result[:,i] = res[i*12:(i+1)*12]
    result = result.reshape(12,y,x)
    result = np.einsum('ijk->ikj', result)

    return result

