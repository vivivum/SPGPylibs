#=============================================================================
# Project: SoPHI
# File:    phi_rte.py
# Author:  David Orozco Suárez (orozco@iaa.es)
# Contributors: 
#-----------------------------------------------------------------------------
# Description: programs for accesing data and fits files
#-----------------------------------------------------------------------------
from .tools import *
try:
    from .cmilos import pmilos
except:
    print("unable to import pmilos in phi_rte.py (this is o.k.)")

import subprocess
import numpy as np
#although not necessary these are the dtype to be passed to C
DTYPE_INT = np.intc
DTYPE_DOUBLE = np.float_

@timeit
def phi_rte(data: np.ndarray,wave_axis: np.ndarray,rte_mode:str,output_dir:str,cmilos = None,options:list = None):
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

    #the next input to pymilos is the options
    #for the moment no PSF is included 
    # and classical estimates are deactivated.
    # these will be added in following versions
    try:
        if options is None:
            options = np.zeros((4))#,dtype=DTYPE_INT)
            options[0] = len(wave_axis) #NLAMBDA wave axis dimension
            options[1] = 30 #MAX_ITER max number of iterations
            options[2] = 1 #CLASSICAL_ESTIMATES [0,1] classical estimates ON or OFF
            options[3] = 0 #RFS [0,1,2] 0.-> Inversion, 1-> synthesis 2-> RFS
            print('No options')
            # options[4] = 100 #FWHM = atof(argv[5]);
            # options[5] = 40  #DELTA = atof(argv[6]);
            # options[6] = 25  #NMUESTRAS_G = atoi(argv[7]);
        else:
            print(len(options) == 4)
    except:
        print('ups')
    # options = np.zeros((7))#,dtype=DTYPE_INT)
    # options[0] = len(wave_axis) #NLAMBDA wave axis dimension
    # options[1] = 15 #MAX_ITER max number of iterations
    # options[2] = 0 #CLASSICAL_ESTIMATES [0,1] classical estimates ON or OFF
    # options[3] = 0 #RFS [0,1,2] 0.-> Inversion, 1-> synthesis 0-> RFS

    if rte_mode == 'RTE':
        options[2] = 0
    elif rte_mode == 'CE':
        options[2] = 2
    elif rte_mode == 'CE+RTE':
        options[2] = 1
    elif rte_mode == 'CE+RTE+PSF':
        options = np.zeros((7))#,dtype=DTYPE_INT)
        options[0] = len(wave_axis) #NLAMBDA wave axis dimension
        options[1] = 30 #MAX_ITER max number of iterations
        options[2] = 1 #CLASSICAL_ESTIMATES [0,1] classical estimates ON or OFF
        options[3] = 0 #RFS [0,1,2] 0.-> Inversion, 1-> synthesis 2-> RFS
        options[4] = 105 #FWHM = atof(argv[5]);
        options[5] = 70  #DELTA = atof(argv[6]);
        options[6] = 25  #NMUESTRAS_G = atoi(argv[7]);

    else:
        printc('RET option not recognized: ',rte_mode,color=bcolors.FAIL)
        return

    if cmilos == None: #meaning you will use python wrapper
        print('shape input',data.shape) 
        # CHECK DATA DIMENSIONS
        data = np.einsum('ijkl->klji',data)
        y,x,p,l = data.shape
        print('Data shape',data.shape, "should be (y,x,pol,wave) for C")                                                                                                                                                                                                                      

        # the pipeline has (for the moment being) the data in
        # (4, 6, 298, 1176) (pol,wave, y,x)
        # This has to be changed to (y,x,pol,wave) for C

        result =  pmilos(options,data,wave_axis)

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
        printc('   saving data into dummy_in.txt for RTE input. dimensions:',l,p,x,y)
        # GV output_dir was include in the chain (generate_level_2, pft_pipe_modules_ phi_rte) to have it availabel here
        file_dummy_in=output_dir+'dummy_in.txt'
        file_dummy_out=output_dir+'dummy_out.txt'
        filename=file_dummy_in
        #filename = 'dummy_in.txt'
        with open(filename,"w") as f:
            for i in range(x):
                for j in range(y):
                    for k in range(l):
                        f.write('%e %e %e %e %e \n' % (wave_axis[k],data[k,0,j,i],data[k,1,j,i],data[k,2,j,i],data[k,3,j,i]))

        printc('  ---- >>>>> Inverting data.... ',color=bcolors.OKGREEN)

        if rte_mode == 'CE+RTE+PSF':
            trozo = " "+str(options[0].astype(int))+" "+str(options[1].astype(int))+" "+str(options[2].astype(int))+" "+str(options[3].astype(int))+" "+str(options[4].astype(int))+" "+str(options[5].astype(int))+" "+str(options[6].astype(int))
        else:
            trozo = " "+str(options[0].astype(int))+" "+str(options[1].astype(int))+" "+str(options[2].astype(int))+" "+str(options[3].astype(int))
        #GV adding dir to filein/out
        #printc(cmilos+trozo+" dummy_in.txt  >  dummy_out.txt",color=bcolors.OKGREEN)
        #rte_on = subprocess.call(cmilos+trozo+" dummy_in.txt  >  dummy_out.txt",shell=True)
        printc(cmilos+trozo+ " " +file_dummy_in+" > "+file_dummy_out ,color=bcolors.OKGREEN)
        rte_on = subprocess.call(cmilos+trozo+" "+ file_dummy_in+" > "+file_dummy_out ,shell=True)
        printc(rte_on,color=bcolors.OKGREEN)

        print(rte_on)
        printc('  ---- >>>>> Finishing.... ',color=bcolors.OKGREEN)

        printc('  ---- >>>>> Reading results.... ',color=bcolors.OKGREEN)
        #GV res = np.loadtxt('dummy_out.txt')
        res = np.loadtxt(file_dummy_out)

        #GV del_dummy = subprocess.call("rm dummy_in.txt",shell=True)
        del_dummy = subprocess.call("rm "+file_dummy_in,shell=True)
        print(del_dummy)
        #GV del_dummy = subprocess.call("rm dummy_out.txt",shell=True)
        del_dummy = subprocess.call("rm "+file_dummy_out,shell=True)
        print(del_dummy)

        npixels = res.shape[0]/12.
        result = np.zeros((12,y*x)).astype(float)

        for i in range(y*x):
            result[:,i] = res[i*12:(i+1)*12]
        result = result.reshape(12,y,x)
        result = np.einsum('ijk->ikj', result)

        return result

