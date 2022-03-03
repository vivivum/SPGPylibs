#=============================================================================
# Project: SoPHI
# File:    phi_rte.py
# Author:  David Orozco Suárez (orozco@iaa.es)
# Contributors: 
#-----------------------------------------------------------------------------
# Description: programs for accesing data and fits files
#-----------------------------------------------------------------------------
from .tools import *
#from .cmilos import pymilos
import subprocess
import numpy as np
#although not necessary these are the dtype to be passed to C
DTYPE_INT = np.intc
DTYPE_DOUBLE = np.float_

@timeit
def phi_rte(data,wave_axis,rte_mode,cmilos = None,options = None,loopthis=0):
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

    if cmilos == None: #meaning you will use python wrapper
        print('shape input',data.shape) 
        # CHECK DATA DIMENSIONS
        data = np.einsum('ijkl->klji',data)
        y,x,p,l = data.shape
        print('Data shape',data.shape, "should be (y,x,pol,wave) for C")                                                                                                                                                                                                                      

        # the pipeline has (for the moment being) the data in
        # (4, 6, 298, 1176) (pol,wave, y,x)
        # This has to be changed to (y,x,pol,wave) for C

        #the next input to pymilos is the options
        #for the moment no PSF is included 
        # and classical estimates are deactivated.
        # these will be added in following versions
        if options == None:
            options = np.zeros((4))#,dtype=DTYPE_INT)
            options[0] = len(wave_axis) #NLAMBDA wave axis dimension
            options[1] = 15 #MAX_ITER max number of iterations
            options[2] = 0 #CLASSICAL_ESTIMATES [0,1] classical estimates ON or OFF
            options[3] = 0 #RFS [0,1,2] 0.-> Inversion, 1-> synthesis 0-> RFS

        result =  pymilos.pmilos(options,data,wave_axis)

        return np.einsum('ijk->kij',result)                                                                                                                                                                                     

			# outputdata[cnt_model] = contador;
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

    else:

        l,p,x,y = data.shape
        print(l,p,x,y)
        if loopthis == 5:
            print(data[:,0,0,0])
            datad = np.copy(data)
            datad[0:5,:,:,:] = data[1:,:,:,:]
            datad[5,:,:,:] = data[0,:,:,:]
            data = np.copy(datad)
            del datad
            print(data[:,0,0,0])
            printc('RUNING TEST LOOP 5 - BAD INPUT TO RTE',color=bcolors.FAIL)
        else:
            printc('pass',color=bcolors.FAIL)

        filename = 'dummy_in.txt'
        with open(filename,"w") as f:
            for i in range(x):
                for j in range(y):
                    for k in range(l):
                        f.write('%e %e %e %e %e \n' % (wave_axis[k],data[k,0,j,i],data[k,1,j,i],data[k,2,j,i],data[k,3,j,i]))

        printc('  ---- >>>>> Inverting data.... ',color=bcolors.OKGREEN)

        if rte_mode == 'RTE':
            rte_on = subprocess.call(cmilos+" 6 15 0 0 dummy_in.txt  >  dummy_out.txt",shell=True)
        if rte_mode == 'CE':
            rte_on = subprocess.call(cmilos+" 6 15 2 0 dummy_in.txt  >  dummy_out.txt",shell=True)
        if rte_mode == 'CE+RTE':
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

