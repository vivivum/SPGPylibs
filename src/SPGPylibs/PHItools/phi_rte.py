"""
=============================================================================
# Project: SoPHI
# File:    phi_rte.py
# Author:  David Orozco Suárez (orozco@iaa.es)
# Contributors: Alex Feller (feller@mps.mpg.de)
#-----------------------------------------------------------------------------
# Description: programs for accesing data and fits files
-----------------------------------------------------------------------------
"""

import concurrent.futures
import os
import inspect
from .tools import printc, bcolors, timeit

try:
    from .cmilos import pmilos
except ImportError:
    print("unable to import pmilos in phi_rte.py (this is o.k.)")

import subprocess
import numpy as np

# although not necessary these are the dtype to be passed to C
DTYPE_INT = np.intc
DTYPE_DOUBLE = np.float_


# @timeit
def phi_rte(
    data: np.ndarray, wave_axis: np.ndarray, rte_mode: str, output_dir: str = None,
    cmilos=None, options: list = None,
    parallel=False, num_workers=10
):
    ''' For the moment this is just isolated from the main pipeline
    input should be:
            l,p,x,y = data.shape  -cmilos  (DEFAULT)
    '''

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

    # the next input to pymilos is the options
    # for the moment no PSF is included
    # and classical estimates are deactivated.
    # these will be added in following versions
    try:
        if options is None:
            options = np.zeros((4))
            options[0] = len(wave_axis)  # NLAMBDA wave axis dimension
            options[1] = 30  # MAX_ITER max number of iterations
            options[2] = 1  # CLASSICAL_ESTIMATES [0,1] classical estimates ON or OFF
            options[3] = 0  # RFS [0,1,2] 0.-> Inversion, 1-> synthesis 2-> RFS
            print('No options')
            # options[4] = 100 #FWHM = atof(argv[5]);
            # options[5] = 40  #DELTA = atof(argv[6]);
            # options[6] = 25  #NMUESTRAS_G = atoi(argv[7]);
        else:
            print(len(options) == 4)
    except TypeError:
        print('ups')

    if rte_mode == 'RTE':
        options[2] = 0
    elif rte_mode == 'CE':
        options[2] = 2
    elif rte_mode == 'CE+RTE':
        options[2] = 1
    elif rte_mode == 'CE+RTE+PSF':
        options = np.zeros((7))
        options[0] = len(wave_axis)  # NLAMBDA wave axis dimension
        options[1] = 30  # MAX_ITER max number of iterations
        options[2] = 1  # CLASSICAL_ESTIMATES 0: disabled[0,1] classical estimates ON or OFF
        options[3] = 0  # RFS [0,1,2] 0.-> Inversion, 1-> synthesis 2-> RFS
        options[4] = 105  # FWHM = atof(argv[5]);
        options[5] = 70  # DELTA = atof(argv[6]);
        options[6] = 25  # NMUESTRAS_G = atoi(argv[7]);

    else:
        printc('RET option not recognized: ', rte_mode, color=bcolors.FAIL)
        return

    if cmilos is None:  # use python wrapper
        print('shape input', data.shape)

        # Adjust and check data dimensions
        data = np.einsum('ijkl->klji', data)
        print('Data shape', data.shape, "should be (y, x, pol, wave) for C")

        # the pipeline has (for the moment being) the data in
        # (4, 6, 298, 1176) (pol,wave, y,x)
        # This has to be changed to (y,x,pol,wave) for C

        global phi_rte_map_func

        def phi_rte_map_func(args):
            stripe, data = args
            return stripe, pmilos(options, data, wave_axis)

        ny, nx, npol, nwave = data.shape
        data = np.reshape(data, (num_workers, ny // num_workers, nx, npol, nwave))  # split data spatially into row-wise stripes
        args_list = [(stripe, data[stripe, :]) for stripe in range(num_workers)]

        if parallel:
            with concurrent.futures.ProcessPoolExecutor(num_workers) as executor:
                results = executor.map(phi_rte_map_func, args_list)
        else:
            results = map(phi_rte_map_func, args_list)

        stripes = []
        data = []
        for result in results:
            stripes.append(result[0])
            data.append(result[1])
        data = np.array(data)
        data = np.reshape(data, (ny, nx, 12))

        return np.einsum('ijk->kij', data)

        # outputdata[cnt_model] = counter;
        # outputdata[cnt_model+1] = iter;
        # outputdata[cnt_model+2] = initModel.B;
        # outputdata[cnt_model+3] = initModel.gm;
        # outputdata[cnt_model+4] = initModel.az;
        # outputdata[cnt_model+5] = initModel.eta0;
        # outputdata[cnt_model+6] = initModel.dopp;
        # outputdata[cnt_model+7] = initModel.aa;
        # outputdata[cnt_model+8] = initModel.vlos;
        # outputdata[cnt_model+9] = initModel.S0;
        # outputdata[cnt_model+10] = initModel.S1;
        # outputdata[cnt_model+11] = chisqrf;

    else:  # use CMILOS executable via subprocess call

        wave, p, y, x = data.shape
        printc('   saving data into dummy_in.txt for RTE input. dimensions (l,p,y,x):', wave, p, y, x)

        if output_dir is None:
            # get file path from calling script
            output_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

        file_dummy_in = os.path.join(output_dir, 'dummy_in.txt')
        file_dummy_out = os.path.join(output_dir, 'dummy_out.txt')

        filename = file_dummy_in
        with open(filename, "w") as f:
            for i in range(x):
                for j in range(y):
                    for k in range(wave):
                        f.write('%e %e %e %e %e \n' % (wave_axis[k], data[k, 0, j, i], data[k, 1, j, i], data[k, 2, j, i], data[k, 3, j, i]))

        printc('  ---- >>>>> Inverting data.... ', color=bcolors.OKGREEN)

        if rte_mode == 'CE+RTE+PSF':
            trozo = " " \
                + str(options[0].astype(int)) + " " \
                + str(options[1].astype(int)) + " " \
                + str(options[2].astype(int)) + " " \
                + str(options[3].astype(int)) + " " \
                + str(options[4].astype(int)) + " " \
                + str(options[5].astype(int)) + " " \
                + str(options[6].astype(int))
        else:
            trozo = " " \
                + str(options[0].astype(int)) + " " \
                + str(options[1].astype(int)) + " " \
                + str(options[2].astype(int)) + " " \
                + str(options[3].astype(int))

        cmd = cmilos + trozo + " " + file_dummy_in + " > " + file_dummy_out
        printc(cmd, color=bcolors.OKGREEN)
        rte_on = subprocess.call(cmd, shell=True)
        printc(rte_on, color=bcolors.OKGREEN)

        print(rte_on)
        printc('  ---- >>>>> Finishing.... ', color=bcolors.OKGREEN)

        printc('  ---- >>>>> Reading results.... ', color=bcolors.OKGREEN)
        result = np.loadtxt(file_dummy_out)
        result = np.reshape(result, (x, y, 12))
        result = np.einsum('jik->kij', result)

        return result
