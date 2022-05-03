### phifdt_pipe_modules

from sys import exit
import numpy as np 
import random, statistics
from scipy.ndimage import gaussian_filter
from matplotlib import pyplot as plt
from .tools import printc,bcolors,fix_path
from .phi_fits import fits_get
from .phi_rte import phi_rte
from .phi_utils import find_string,azimutal_average,newton,limb_darkening,genera_2d
from .phi_gen import bin_annulus,shift,find_center,apod,rebin,generate_circular_mask
from .phi_reg import shift_subp,moments

import SPGPylibs.GENtools.plot_lib as plib

# from platform import node
import os

def phi_correct_dark(dark_f,data,header,data_scale,verbose = False,get_dark = False):

    #-----------------
    # READ AND CORRECT DARK FIELD
    #-----------------

    printc('-->>>>>>> Reading Darks                   ',color=bcolors.OKGREEN)
    printc('          Input should be [y-dim,x-dim].',color=bcolors.OKGREEN)
    printc('          DARK IS DIVIDED by 256.   ',color=bcolors.OKGREEN)

    try:
        dark,dark_header = fits_get(dark_f)
        dark = dark / 256.
    except Exception:
        printc("ERROR, Unable to open darks file: {}",dark_f,color=bcolors.FAIL)
        raise 

    # locations = find_string(dark_f,'_')
    # try:
    #     DID = dark_f[locations[-1]+1:locations[-1]+10]
    #     print('DID: ',np.float(DID))
    # except Exception:
    #     printc("Unable to get DID: {}",dark_f,color=bcolors.FAIL)         
    #     printc('DID: ',DID,' -->> NOT A NUMBER',color=bcolors.FAIL)
    #     raise 

    DID = dark_header['PHIDATID']
    printc('Dark DID: ',DID,color=bcolors.OKBLUE)
    dark_scale = fits_get(dark_f,scaling = True)

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
    #CHECK NACC
    acc = int(header['ACCACCUM']) * int(header['ACCCOLIT'])
    acc_dark = int(dark_header['ACCACCUM']) * int(dark_header['ACCCOLIT'])
    if acc != acc_dark:
        printc('WARNING - NACC NOT IDENTICAL DURING DARK CORRECTION',color=bcolors.FAIL)
        printc('DARK NACC ',acc_dark,' DATA NACC ',acc,color=bcolors.FAIL)
        
    if verbose:
        dummy = data[0,0,:,:]
    data = data - dark[np.newaxis,np.newaxis,PXBEG2:PXEND2+1,PXBEG1:PXEND1+1]
    data = np.abs(data)

    if 'CAL_DARK' in header:  # Check for existence
        header['CAL_DARK'] = DID

    if verbose:
        md = np.mean(dark)
        if verbose != True:
            plib.show_three(dark,dummy,data[0,0,:,:],vmin=[-md,-md,-md],vmax=[md*2,md*2,md*2],block=True,pause=0.1,title=['Dark','Data','Data after dark correction'],
                xlabel='Pixel',ylabel='Pixel',cmap='gray',save=verbose)
        else:
            plib.show_three(dark,dummy,data[0,0,:,:],vmin=[-md,-md,-md],vmax=[md*2,md*2,md*2],block=True,pause=0.1,title=['Dark','Data','Data after dark correction'],
                xlabel='Pixel',ylabel='Pixel',cmap='gray')

    return data,header

def interpolateImages(image1, image2, dist1I, distI2):
    ''' interpolate 2D images - 
    '''
    imageInterp = (image1 * distI2 + image2 * dist1I) / (dist1I + distI2)
    return imageInterp

def phi_correct_prefilter(prefilter_fits,header,data,voltagesData,verbose = False):

    printc('-->>>>>>> Read prefilter and correct for it ')
    printc('          ',prefilter_fits,'   ')
    PXBEG1  = int(header['PXBEG1']) - 1           
    PXEND1  = int(header['PXEND1']) - 1          
    PXBEG2  = int(header['PXBEG2']) - 1           
    PXEND2  = int(header['PXEND2']) - 1   
    xd  = int(header['NAXIS1']) 
    yd  = int(header['NAXIS2'])    
    zd  = int(header['NAXIS3'])    

    prefdata,h = fits_get(prefilter_fits)
    prefdata = prefdata.astype(float)
    prefdata = prefdata[:,PXBEG2:PXEND2+1,PXBEG1:PXEND1+1]
    #PREFILTER INFO
    #——————
    #I call the prefilter correction from the python pipeline like this:
    #data_corr = np.reshape(data_corr,(6,4,2048, 2048)) (6,4,y,x)
    prefscale = 8388608./2.
    # the scale of the prefilter is irrelevant because it is a relative variation
    # data = data / prefilter. Scale = 2 makes the prefilter to be aroung 1
    # later, everything is normalized wrt the continuum so this is not important.
    scale = 1.
    prefVoltages = [-1300, -1234, -1169, -1103, -1038, -972, -907, -841, -776,\
        -710, -645, -579, -514, -448, -383, -317, -252, -186,-121, -56,9,74,\
        140,205,271,336,402,467,533,598,664,729,795,860,926,991,1056,1122,1187,\
        1253,1318,1384,1449,1515,1580,1646,1711,1777,1842]
    if verbose:
        datap = np.copy(data)
        data = applyPrefilter_dos(data, voltagesData, prefdata, prefscale, prefVoltages, -1,scaledown = scale, verbose = verbose)
    if verbose:
        plt.plot(data[:,0,yd//2,xd//2],'o-',label='corrected')
        plt.plot(datap[:,0,yd//2,xd//2],'--',label='original')
        plt.legend()
        plt.show()
        del datap
    slash,nothing = find_string(prefilter_fits,'/')
    if nothing == -1:
        slash = [0]
    if 'CAL_PRE' in header:  # Check for existence
        header['CAL_PRE'] = prefilter_fits[slash[-1]+1:-4]
    else:
        header.set('CAL_PRE', prefilter_fits[slash[-1]+1:-4], 'prefilter file',after='CAL_DARK')
    return data,header

def applyPrefilter(data, wvltsData, prefilter, prefScale, wvltsPref, direction, scaledown=8,verbose = False):
    '''PHI prefilter. Version from K. Albert.
    '''
    prefToApply = np.zeros((6,prefilter.shape[1],prefilter.shape[2]))
    #prefilter = prefilter/prefScale
    for i in range(0,6):
        wvlCurr = wvltsData[i]
        valueClosest = min(wvltsPref, key=lambda x:abs(x-wvlCurr))
        if verbose:
            print("iter", i, "wvlCurr", wvlCurr)
            print("iter", i, "valueClosest", valueClosest)
        indexClosest = wvltsPref.index(valueClosest)
        if verbose:
            print("iter", i, "indexClosest", indexClosest)
        if (valueClosest < wvlCurr):
            indexBefore = indexClosest
            indexAfter = indexClosest + 1
        else:
            indexAfter = indexClosest
            indexBefore = indexClosest - 1

        dist1I = abs(wvltsPref[indexBefore] - wvltsData[i])
        distI2 = abs(wvltsPref[indexAfter] - wvltsData[i])
        prefToApply[i,:,:] = interpolateImages(prefilter[indexBefore], prefilter[indexAfter], dist1I, distI2)

        if verbose:
            print("mean prefValue Before:", np.mean(prefilter[indexBefore])*256)
            print("mean prefValue After:", np.mean(prefilter[indexAfter])*256)
            print("distance1:", dist1I)
            print("distance2:", distI2)
            print("percentage:", distI2 / (dist1I + distI2))

        #Remove scale factor from prefilter
        if verbose:
            print("mean prefilter:", np.mean(prefToApply[i,:,:])*256)
        prefToApply[i,:,:] = prefToApply[i,:,:] / prefScale
        if verbose:
            print("mean prefilter:", np.mean(prefToApply[i,:,:]))

    if verbose:
        print("Reshaping prefilter:")
        print(prefToApply.shape)
        print(data.shape)
    if(data.shape[2] != prefToApply.shape[1]):
        FOV_Start_y = int(prefToApply.shape[1]/2 - data.shape[2]/2)
        FOV_End_y = int(prefToApply.shape[1]/2 + data.shape[2]/2)
        prefToApply = prefToApply[:,FOV_Start_y:FOV_End_y,:]
    if verbose:
         print(prefToApply.shape)
    if(data.shape[3] != prefToApply.shape[2]):
        FOV_Start_x = int(prefToApply.shape[2]/2 - data.shape[3]/2)
        FOV_End_x = int(prefToApply.shape[2]/2 + data.shape[3]/2)
        prefToApply = prefToApply[:,:,FOV_Start_x:FOV_End_x]
    if verbose:
        print(prefToApply.shape)

    dataPrefApplied = np.zeros(data.shape)
    for i in range(0,4):
        if(direction == 1):
            dataPrefApplied[:,i,:,:] = data[:,i,:,:] * prefToApply
        elif(direction == -1):
            dataPrefApplied[:,i,:,:] = data[:,i,:,:] / prefToApply / scaledown
        else:
            print("Ivnalid direction! Must be 1 (mult) or -1 (div).")
    return dataPrefApplied

    #/**
    # * This is the maximum range scaled from the division of the prefilter.
    # * The smallest number in the prefilter is 0.13977 -> 1/0.13977 = 7.3 ~= 8
    # */
    #define RNG_RES_APPL_PREF 8  -> reason for division by 8.

def applyPrefilter_dos(data, wvltsData, prefilter, prefScale, wvltsPref, direction, scaledown=8,verbose = False):
    '''PHI prefilter. Modified version from K. Albert.
    '''
    prefToApply = np.zeros((6,prefilter.shape[1],prefilter.shape[2]))
    prefilter = prefilter / prefScale  #dos
    
    for i in range(0,6):
        wvlCurr = wvltsData[i]
        valueClosest = min(wvltsPref, key=lambda x:abs(x-wvlCurr))
        if verbose:
            print("iter", i, "wvlCurr", wvlCurr)
            print("iter", i, "valueClosest", valueClosest)
        indexClosest = wvltsPref.index(valueClosest)
        if verbose:
            print("iter", i, "indexClosest", indexClosest)
        if (valueClosest < wvlCurr):
            indexBefore = indexClosest
            indexAfter = indexClosest + 1
        else:
            indexAfter = indexClosest
            indexBefore = indexClosest - 1

        dist1I = abs(wvltsPref[indexBefore] - wvltsData[i])
        distI2 = abs(wvltsPref[indexAfter] - wvltsData[i])
        prefToApply[i,:,:] = interpolateImages(prefilter[indexBefore], prefilter[indexAfter], dist1I, distI2)

        if verbose:
            print("mean prefValue Before:", np.mean(prefilter[indexBefore])*256)
            print("mean prefValue After:", np.mean(prefilter[indexAfter])*256)
            print("distance1:", dist1I)
            print("distance2:", distI2)
            print("percentage:", distI2 / (dist1I + distI2))

        #Remove scale factor from prefilter
        if verbose:
            print("mean prefilter:", np.mean(prefToApply[i,:,:])*256)
        #prefToApply[i,:,:] = prefToApply[i,:,:] / prefScale   #dos
        if verbose:
            print("mean prefilter:", np.mean(prefToApply[i,:,:]))


    if verbose:
        print("Reshaping prefilter:")
        print(prefToApply.shape)
        print(data.shape)
    if(data.shape[2] != prefToApply.shape[1]):
        FOV_Start_y = int(prefToApply.shape[1]/2 - data.shape[2]/2)
        FOV_End_y = int(prefToApply.shape[1]/2 + data.shape[2]/2)
        prefToApply = prefToApply[:,FOV_Start_y:FOV_End_y,:]
    if verbose:
         print(prefToApply.shape)
    if(data.shape[3] != prefToApply.shape[2]):
        FOV_Start_x = int(prefToApply.shape[2]/2 - data.shape[3]/2)
        FOV_End_x = int(prefToApply.shape[2]/2 + data.shape[3]/2)
        prefToApply = prefToApply[:,:,FOV_Start_x:FOV_End_x]
    if verbose:
        print(prefToApply.shape)

    dataPrefApplied = np.zeros(data.shape)
    for i in range(0,4):
        if(direction == 1):
            dataPrefApplied[:,i,:,:] = data[:,i,:,:] * prefToApply
        elif(direction == -1):
            dataPrefApplied[:,i,:,:] = data[:,i,:,:] / prefToApply # / scaledown  #dos
        else:
            print("Ivnalid direction! Must be 1 (mult) or -1 (div).")
    return dataPrefApplied

def check_pmp_temp(hdr_arr):
    """
    check science scans have same PMP temperature set point
    From HRT software (Jonas et al)
    """
    first_pmp_temp = hdr_arr['HPMPTSP1']
    # result = all(hdr['HPMPTSP1'] == first_pmp_temp for hdr in hdr_arr)
    if first_pmp_temp:
        print(f"The scan have a PMP Temperature Set Point: {first_pmp_temp}")
        pmp_temp = str(first_pmp_temp)
        return pmp_temp
    else:
        print("The scans have different PMP Temperatures! Please fix \n Ending Process")

        exit(1)

def phi_apply_demodulation(data,instrument,header = False,demod=False,verbose = False,modulate = False):
    '''
    Use demodulation matrices to demodulate data size ls,ps,ys,xs  (n_wave*S_POL,N,M)
    ATTENTION: FDT40 is fixed to the one Johann is using!!!!
    '''

    if instrument == 'FDT40': #MODEL FIT  INTA April 2022
        mod_matrix = np.array( [[ 0.99913 , -0.69504 , -0.38074 , -0.60761 ],\
                                [ 1.0051  ,  0.41991 , -0.73905 ,  0.54086 ],\
                                [ 0.99495 ,  0.44499 ,  0.36828 , -0.8086  ],\
                                [ 1.0008  , -0.38781 ,  0.91443 ,  0.13808 ]] )
        demodM = np.linalg.inv(mod_matrix)
    elif instrument == 'FDT40_old': # Johanns (it is the average in the central area of the one onboard) 
        mod_matrix = np.array([[1.0006,-0.7132, 0.4002,-0.5693],
                                  [1.0048, 0.4287,-0.7143, 0.5625],
                                  [0.9963, 0.4269,-0.3652,-0.8229],
                                  [0.9983,-0.4022, 0.9001, 0.1495]])
        demodM = np.linalg.inv(mod_matrix)
    elif instrument == 'FDT40_dos': #NUMERICAL FIT  (0,360,10)-90   *** David Orozco March 2022
        mod_matrix = np.array( [[ 0.99493621 ,-0.69068949 ,-0.36619221 ,-0.60735698 ],\
                                   [ 1.00582274 , 0.41415974 ,-0.73663873 , 0.54241813 ],\
                                   [ 0.99684194 , 0.4458513  , 0.36831856 , -0.8134466 ],\
                                   [ 1.00239911 ,-0.38416442 , 0.9165126  , 0.13481876]] )
        demodM = np.linalg.inv(mod_matrix)
    elif instrument == 'FDT40_mod': #ANALITICA FIT   (0,360,10)-90    *** David Orozco March 2022 - based on PHI / FDT model
        mod_matrix = np.array([[ 1.    ,     -0.69634929 ,-0.37330967 ,-0.61297435],\
                                    [ 1.  ,        0.41324362 ,-0.73330607,  0.53989992],\
                                    [ 1.  ,        0.44598252 , 0.36860311, -0.81561715],\
                                    [ 1.   ,      -0.3846079  , 0.91317683,  0.13485121]])
        demodM = np.linalg.inv(mod_matrix)
    elif instrument == 'FDT45':  #MODEL FIT  INTA April 2022
        mod_matrix = np.array([[1.0023,-0.64814, -0.56202,-0.51859],
                                  [1.0041, 0.54693, -0.55299, 0.633],
                                  [0.99523, 0.46132,0.54165,-0.69603],
                                  [0.99838,-0.61944,0.66189, 0.42519]])
        demodM = np.linalg.inv(mod_matrix)
    elif instrument == 'FDT45_HREW':  #MODEL FIT  INTA April 2022
        mod_matrix = np.array([[1.0022,-0.6543, -0.57471,-0.49363],
                                  [1.0038, 0.54924, -0.53817, 0.64209],
                                  [0.99438, 0.45867,0.515,-0.71314],
                                  [0.99964,-0.61277,0.67894, 0.40843]])
        demodM = np.linalg.inv(mod_matrix)
    elif instrument == 'FDT45_E2E': #E2E doc
        mod_matrix = np.array([[1.0035,-0.6598, 0.5817,-0.4773],
                                  [1.0032, 0.5647, 0.5275, 0.6403],
                                  [0.9966, 0.4390,-0.5384,-0.7150],
                                  [0.9968,-0.6169,-0.6443, 0.4425]])
        demodM = np.linalg.inv(mod_matrix)
        
    elif instrument == 'HRT40_old': #E2E doc
        demodM = np.array([[ 0.26450154,  0.2839626,   0.12642948,  0.3216773 ],
                            [ 0.59873885,  0.11278069, -0.74991184,  0.03091451],
                            [ 0.10833212, -0.5317737,  -0.1677862,   0.5923593 ],
                            [-0.46916953,  0.47738808, -0.43824592,  0.42579797]])
    elif instrument == 'HRT40': #MODEL FIT  INTA April 2022
        mod_matrix = np.array([[ 0.99816  ,0.61485 , 0.010613 ,-0.77563 ], 
                               [ 0.99192 , 0.08382 , 0.86254 , 0.46818],
                               [ 1.0042 , -0.84437 , 0.12872 ,-0.53972],
                               [ 1.0057 , -0.30576 ,-0.87969 , 0.40134]])
        demodM = np.linalg.inv(mod_matrix)
    elif instrument == 'HRT40_dos': #NUMERICAL FIT  (0,360,10)-90   *** David Orozco March 2022
        mod_matrix = np.array([[ 0.98889418  ,0.63311402 , 0.01490015 ,-0.7816713 ],  
                            [ 0.99603854 , 0.07763719 , 0.89156538 , 0.46243721],
                            [ 1.00427941 ,-0.84080523 , 0.12554703 ,-0.55159669],
                            [ 1.01078787 ,-0.29564517 ,-0.88365522 , 0.41235959]])
        demodM = np.linalg.inv(mod_matrix)
    elif instrument == 'HRT45': #E2E doc
        mod_matrix = np.array([[1.00159,-0.50032, 0.7093,-0.4931],
                                      [1.0040, 0.6615, 0.3925, 0.6494],
                                      [0.9954, 0.3356,-0.6126,-0.7143],
                                      [0.9989,-0.7474,-0.5179, 0.4126]]) #MIA
        demodM = np.linalg.inv(mod_matrix)
    elif instrument == 'HRT50_E2E': #E2E doc
        demodM = np.array([[ 0.28037298,  0.18741922,  0.25307596,  0.28119895],
                     [ 0.40408596,  0.10412157, -0.7225681,   0.20825675],
                     [-0.19126636, -0.5348939,   0.08181918,  0.64422774],
                     [-0.56897295,  0.58620095, -0.2579202,   0.2414017 ]])
    elif instrument == 'HRT50': #MODEL FIT  INTA April 2022
        mod_matrix = np.array([[ 1.0014  ,0.56715 , 0.3234 ,-0.74743 ], 
                               [ 1.0007 , 0.0037942 , 0.69968 , 0.71423],
                               [ 1.0002 , -0.98937 , 0.04716 ,-0.20392],
                               [ 0.99769 , 0.27904 ,-0.86715 , 0.39908]])
        demodM = np.linalg.inv(mod_matrix)
    else:
        printc('No demod available in demod_phi.py',color = bcolors.FAIL)
        raise SystemError()

    printc('Demodulation matrix for ', instrument,color = bcolors.WARNING)
    printc(demodM,color = bcolors.WARNING)

    if demod:
        return demodM

    ls,ps,ys,xs = data.shape

    if modulate:
        for i in range(ls):
            data[i,:,:,:] = np.reshape(np.matmul(mod_matrix,np.reshape(data[i,:,:,:],(ps,xs*ys))),(ps,ys,xs))
        return data


    for i in range(ls):
        data[i,:,:,:] = np.reshape(np.matmul(demodM,np.reshape(data[i,:,:,:],(ps,xs*ys))),(ps,ys,xs))

    if header != False:
        if 'CAL_IPOL' in header:  # Check for existence
            header['CAL_IPOL'] = instrument
        else:
            header.set('CAL_IPOL', instrument, 'Onboard calibrated for instrumental polarization',after='CAL_DARK')
        return data, header
    else:
        return data

def crosstalk_ItoQUV(data_demod,verbose=False,npoints=2000):
    
    limit=0.2
    PLT_RNG = 3
    my = []
    sy = []

    x = data_demod[:,0,:,:].flatten()
    ids = x > limit
    x = x[ids].flatten()

    N = x.size
    idx = random.sample(range(N),npoints)
    mx = x[idx].mean() 
    sx = x[idx].std() 
    xp = np.linspace(x.min(), x.max(), 100)

    A = np.vstack([x, np.ones(len(x))]).T

    # I to Q
    yQ = data_demod[:,1,:,:].flatten()
    yQ = yQ[ids].flatten()
    my.append(yQ[idx].mean())
    sy.append(yQ[idx].std())

    # cQ = np.polyfit( x , yQ , 1)
    # pQ = np.poly1d(cQ)
    m, c = np.linalg.lstsq(A, yQ, rcond=None)[0]
    cQ = [m,c]
    pQ = np.poly1d(cQ)

    # I to U
    yU = data_demod[:,2,:,:].flatten()
    yU = yU[ids].flatten()
    my.append(yU[idx].mean())
    sy.append(yU[idx].std())

    # cU = np.polyfit( x , yU , 1)
    # pU = np.poly1d(cU)
    m, c = np.linalg.lstsq(A, yU, rcond=None)[0]
    cU = [m,c]
    pU = np.poly1d(cU)

    # I to V
    yV = data_demod[:,3,:,:].flatten()
    yV = yV[ids].flatten()
    my.append(yV[idx].mean())
    sy.append(yV[idx].std())

    # cV = np.polyfit( x , yV , 1)
    # pV = np.poly1d(cV)
    m, c = np.linalg.lstsq(A, yV, rcond=None)[0]
    cV = [m,c]
    pV = np.poly1d(cV)

    if verbose:
        
        plt.figure(figsize=(8, 8))
        plt.scatter(x[idx],yQ[idx],color='red',alpha=0.6,s=10)
        plt.plot(xp, pQ(xp), color='red', linestyle='dashed',linewidth=3.0)

        plt.scatter(x[idx],yU[idx],color='blue',alpha=0.6,s=10)
        plt.plot(xp, pU(xp), color='blue', linestyle='dashed',linewidth=3.0)

        plt.scatter(x[idx],yV[idx],color='green',alpha=0.6,s=10)
        plt.plot(xp, pV(xp), color='green', linestyle='dashed',linewidth=3.0)

        plt.xlim([mx - PLT_RNG * sx,mx + PLT_RNG * sx])
        plt.ylim([min(my) - 1.8*PLT_RNG * statistics.mean(sy),max(my) + PLT_RNG * statistics.mean(sy)])
        plt.xlabel('Stokes I')
        plt.ylabel('Stokes Q/U/V')
        plt.text(mx - 0.9*PLT_RNG * sx, min(my) - 1.4*PLT_RNG * statistics.mean(sy), 'Cross-talk from I to Q: slope = {: {width}.{prec}f} ; off-set = {: {width}.{prec}f} '.format(cQ[0],cQ[1],width=8,prec=4), style='italic',bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 1}, fontsize=15)
        plt.text(mx - 0.9*PLT_RNG * sx, min(my) - 1.55*PLT_RNG * statistics.mean(sy), 'Cross-talk from I to U: slope = {: {width}.{prec}f} ; off-set = {: {width}.{prec}f} '.format(cU[0],cU[1],width=8,prec=4), style='italic',bbox={'facecolor': 'blue', 'alpha': 0.1, 'pad': 1}, fontsize=15)
        plt.text(mx - 0.9*PLT_RNG * sx, min(my) - 1.7*PLT_RNG * statistics.mean(sy), 'Cross-talk from I to V: slope = {: {width}.{prec}f} ; off-set = {: {width}.{prec}f} '.format(cV[0],cV[1],width=8,prec=4), style='italic',bbox={'facecolor': 'green', 'alpha': 0.1, 'pad': 1}, fontsize=15)
        plt.show()

    print('Cross-talk from I to Q: slope = {: {width}.{prec}f} ; off-set = {: {width}.{prec}f} '.format(cQ[0],cQ[1],width=8,prec=4))
    print('Cross-talk from I to U: slope = {: {width}.{prec}f} ; off-set = {: {width}.{prec}f} '.format(cU[0],cU[1],width=8,prec=4))
    print('Cross-talk from I to V: slope = {: {width}.{prec}f} ; off-set = {: {width}.{prec}f} '.format(cV[0],cV[1],width=8,prec=4))

    return cQ,cU,cV

def cross_talk_QUV(data,nran = 2000,nlevel=0.3,verbose=False,block=True):
    
    limit=0.1
    PLT_RNG = 3
    my = []
    sy = []

    #para el limite en I
    lx = data[:,0,:,:].flatten()
    lv = np.abs(data[:,3,:,:]).flatten()
    if verbose:
        print(lv)
        print(nlevel/100.)
    ids = ((lx > limit) & (lv > nlevel/100.))
    ids = (lv > nlevel/100.) #and (lv > nlevel/100.)

    x = data[:,3,:,:].flatten()
    x = x[ids].flatten()
    if verbose:
        print(x.size)
    N = x.size
    if N > nran:
        nran = N
    idx = random.sample(range(N),nran)
    mx = x[idx].mean() 
    sx = x[idx].std() 
    xp = np.linspace(x.min(), x.max(), 100)

    A = np.vstack([x, np.ones(len(x))]).T

    # V to Q
    yQ = data[:,1,:,:].flatten()
    yQ = yQ[ids].flatten()
    my.append(yQ[idx].mean())
    sy.append(yQ[idx].std())

    m, c = np.linalg.lstsq(A, yQ, rcond=None)[0]
    cVQ = [m,c]
    pQ = np.poly1d(cVQ)

    # V to U
    yU = data[:,2,:,:].flatten()
    yU = yU[ids].flatten()
    my.append(yU[idx].mean())
    sy.append(yU[idx].std())

    m, c = np.linalg.lstsq(A, yU, rcond=None)[0]
    cVU = [m,c]
    pU = np.poly1d(cVU)

    if verbose:
        plt.figure(figsize=(8, 8))
        plt.scatter(x[idx],yQ[idx],color='red',alpha=0.6,s=10)#,c=dS[idx],cmap=cm.Paired)
        plt.plot(xp, pQ(xp), color='red', linestyle='dashed',linewidth=3.0)

        plt.scatter(x[idx],yU[idx],color='blue',alpha=0.6,s=10)
        plt.plot(xp, pU(xp), color='blue', linestyle='dashed',linewidth=3.0)

        plt.xlim([mx - PLT_RNG * sx,mx + PLT_RNG * sx])
        plt.ylim([min(my) - 1.8*PLT_RNG * statistics.mean(sy),max(my) + PLT_RNG * statistics.mean(sy)])
        plt.xlabel('Stokes V')
        plt.ylabel('Stokes Q/U')

        plt.text(mx - 0.9*PLT_RNG * sx, min(my) - 1.4*PLT_RNG * statistics.mean(sy), 'Cross-talk from V to Q: slope = {: {width}.{prec}f} ; off-set = {: {width}.{prec}f} '.format(cVQ[0],cVQ[1],width=8,prec=4), style='italic',bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 1}, fontsize=15)
        plt.text(mx - 0.9*PLT_RNG * sx, min(my) - 1.55*PLT_RNG * statistics.mean(sy), 'Cross-talk from V to U: slope = {: {width}.{prec}f} ; off-set = {: {width}.{prec}f} '.format(cVU[0],cVU[1],width=8,prec=4), style='italic',bbox={'facecolor': 'blue', 'alpha': 0.1, 'pad': 1}, fontsize=15)
        plt.show(block=block)

    print('Cross-talk from V to Q: slope = {: {width}.{prec}f} ; off-set = {: {width}.{prec}f} '.format(cVQ[0],cVQ[1],width=8,prec=4))
    print('Cross-talk from V to U: slope = {: {width}.{prec}f} ; off-set = {: {width}.{prec}f} '.format(cVU[0],cVU[1],width=8,prec=4))

    return cVQ,cVU

def crosstalk_ItoQUV2d(data_demod,size=4):
    
    from sklearn.feature_extraction import image
    PLT_RNG = 3
    limit=0.1

    #datos son lambda, stokes
    #hago primero I to Q 
    #miro cuandos cuadraditos deben ser
    nwave,npol,yy,xx = data_demod.shape
    wave = 0
    pol = 0
    max_patches = xx//size * yy//size - 1
    print('check1',nwave,npol,yy,xx,max_patches)
    patches = image.extract_patches_2d(data_demod[wave,pol,:,:], (size, size))#,max_patches=max_patches)
    npatches,yy,xx = patches.shape
    print('check2',npatches,yy,xx)
    patchy = np.zeros((npatches,nwave,npol,size,size),dtype='float32')
    patchy[:,wave,pol,:,:] = patches
    print('check3',patchy.shape)
    for i in range(nwave):
        print('doing ',i)
        for j in range(npol):
            print('doing ',i, j)
            patchy[:,i,j,:,:] = image.extract_patches_2d(data_demod[i,j,:,:], (size, size))#,max_patches=max_patches)
    
    print('check4',patchy.shape,npatches)
    #loop over npatches
    #cV = np.zeros((yy,xx),dtype='float32')
    cV0 = np.zeros((npatches),dtype='float32')
    cV1 = np.zeros((npatches),dtype='float32')

    print('check5',' before loop')
    for i in range(npatches):
        if np.remainder(i, npatches//100) == 0:
            print(i, ' of ', npatches,end='', flush=True)
        x = patchy[i,:,0,:,:].flatten()
        ids = x > limit
        x = x[ids].flatten()
        A = np.vstack([x, np.ones(len(x))]).T

        y = patchy[i,:,3,:,:].flatten()
        y = y[ids].flatten()

        m, c = np.linalg.lstsq(A, y, rcond=None)[0]
        cV1[i] = c
        cV0[i] = m

    print(cV0.shape,xx,yy)
    #cV0 = image.reconstruct_from_patches_2d(cV0, (yy,xx))
    #cV1 = image.reconstruct_from_patches_2d(cV1, (yy,xx))
    return cV0,cV1

def phi_correct_ghost_single(data,header,rad,verbose=False): 
    '''
    Startup version on April 2022. Based on phi_correct_ghost but for just one image
    '''
    version = 'phi_correct_ghost_single V1.0 April 2022'

    center = np.array([header['CRPIX1'],header['CRPIX2']]).astype(int)
    printc('        Read center from header (updated): x=',center[0],' y=',center[1],color=bcolors.OKBLUE)
    xd  = int(header['NAXIS1'])    
    yd  = int(header['NAXIS2'])    
    zd  = int(header['NAXIS3'])    
    PXBEG1  = int(header['PXBEG1']) - 1           
    PXEND1  = int(header['PXEND1']) - 1          
    PXBEG2  = int(header['PXBEG2']) - 1           
    PXEND2  = int(header['PXEND2']) - 1   

    printc('-->>>>>>> Correcting ghost image ',color=bcolors.OKGREEN)

    #common part
    coef = [-1.98787669,1945.28944245] #empirically (first version)
    coef = [-1.9999,1942.7] #empirically (updated by trial and error)

    #correct center:
    center_c = np.copy(center)
    center_c[0] += PXBEG1 
    center_c[1] += PXBEG2 
    poly1d_fn = np.poly1d(coef)
    sh = poly1d_fn(center_c).astype(int) 
    sh_float = poly1d_fn(center_c)
    printc('          image center: x: ',center[0],' y: ',center[1],color=bcolors.OKGREEN)
    printc('          image center [for 2048]: x: ',center_c[0],' y: ',center_c[1],color=bcolors.OKGREEN)
    printc('          ghost displacements: x: ',sh_float[0],' y: ',sh_float[1],color=bcolors.OKGREEN)

    #generate a ring mask out of the solar disk to see how much is the ghost image
    #we will be using complex histrograms....
    mask_anulus = bin_annulus([yd,xd],rad + 20, 10, full = False)
    mask_anulus = shift(mask_anulus, shift=(center[0]-xd//2,center[1]-yd//2), fill_value=0)

    idx = np.where(mask_anulus == 1)
    mask_anulus_big = bin_annulus([yd,xd],rad - 150, 100, full = False)
    mask_anulus_big = shift(mask_anulus_big, shift=(center[0]-xd//2,center[1]-yd//2), fill_value=0)
    idx_big = np.where(data*mask_anulus_big == 1)

    printc('          computing azimuthal averages  ',color=bcolors.OKGREEN)

    centers = np.zeros((2))
    radius = 0
    ints = np.zeros((int(np.sqrt(xd**2+yd**2))))
    ints_rad = np.zeros((int(np.sqrt(xd**2+yd**2))))
    ints_fit = np.zeros((int(np.sqrt(xd**2+yd**2))))
    ints_syn = np.zeros((int(np.sqrt(xd**2+yd**2))))
    ints_fit_pars = np.zeros((5))
    factor = 0
    mean_intensity = 0

    # STEP --->>> Find center 
    centers[1],centers[0],radius = find_center(data)  #cy,cx....

    # STEP --->>> Generate CLV from averaged data
    intensity, rad = azimutal_average(data,[centers[0],centers[1]])
    ints[0:len(intensity)] = intensity
    ints_rad[0:len(intensity)] = rad

    # STEP --->>> FIT LIMB DATA
    rrange = int(radius + 2) #2
    clv = ints[0:rrange]
    clv_r = ints_rad[0:rrange]
    mu = np.sqrt( (1 - clv_r**2/clv_r[-1]**2) )

    if verbose:
        plt.plot(clv_r,clv)
        plt.xlabel('Solar radious [pixel]')
        plt.ylabel('Intensity [DN]')
        plt.show()

    u = 0.5
    I0 = 100
    ande = np.where(mu > 0.1)
    pars = newton(clv[ande],mu[ande],[I0,u,0.2,0.2,0.2],limb_darkening)
    fit, _ = limb_darkening(mu,pars)
    ints_fit[0:len(fit)] = fit
    ints_fit_pars[:] = pars

    ints_syn = np.copy(ints)
    ints_syn[0:len(fit)] = fit

    # STEP --->>> NORMALIZE
    ints_syn = ints_syn / ints_fit_pars[0]
    ints_fit = ints_fit / ints_fit_pars[0]
    ints = ints / ints_fit_pars[0]

    if verbose:
        plt.plot(ints_fit,label='fitted clv')
        plt.plot(ints,'.',label='real clv')
        plt.plot(ints_syn,'--',label='synt clv')
        plt.xlabel('Heliocentric angle ['+r'$\theta$]')
        plt.ylabel('Intensity [DN]')
        plt.legend()
        plt.show()

    # STEP --->>> GENERATE GHOST

    
    nc = (PXEND2-PXBEG2+1)//2 #  center of frame
    limb_2d = np.zeros((PXEND2-PXBEG2+1,PXEND1-PXBEG1+1)) #generate image for ghost
    #fill the image with the revoluted fit
    s_of_gh = int(radius*1.1)
    limb_2d[ nc - s_of_gh:nc + s_of_gh + 1, nc - s_of_gh:nc + s_of_gh + 1] = genera_2d(ints_syn[0:s_of_gh])
    xl,yl = limb_2d.shape
    if verbose:
        plt.imshow(limb_2d)
        plt.show()

    # STEP --->>> Smooth and SHIFT GHOST

    limb_2d = gaussian_filter(limb_2d, sigma=(8, 8)) 

    # #shift ghost to center of image
    # limb_2d = shift(limb_2d, shift=(int(centers[1,i])-yd//2,int(centers[0,i])-xd//2), fill_value=0)
    limb_2d = shift_subp(limb_2d, shift=[centers[1]-yd//2,centers[0]-xd//2])
    #OJO, shift_subp coge como parametros los de shift pero al reves!!!!!!
    if verbose:
        plib.show_one(limb_2d,vmax=1,vmin=0,xlabel='pixel',ylabel='pixel',title='limb 2D',cbarlabel=' ',cmap='gray')
    # #shift to the position of the ghost
    reflection = shift(limb_2d, shift=(sh[0],sh[1]), fill_value=0) 

    # STEP --->>> Correct

    mean_intensity = np.mean(data[idx_big])

    # FORMA 1 - con anillo
    values = data[idx].flatten() #Take the ring
    #show the histogram
    meanv = np.mean(values)
    idx_l = np.where(values <= meanv)
    m_l = np.mean(values[idx_l])
    idx_r = np.where(values >= meanv)
    m_r = np.mean(values[idx_r])
    factor = (m_r - m_l) * 100. / ints_fit_pars[0] 
    print("factor",factor)

    if verbose:
        plt.hist(values, bins=40)
        plt.title('signal')
        plt.axvline(meanv, lw=2, color='yellow', alpha=0.4)
        plt.axvline(m_l, lw=2, color='red', alpha=0.4)
        plt.axvline(m_r, lw=2, color='blue', alpha=0.4)
        plt.axvline(factor*ints_fit_pars[0] / 100., lw=2, color='green', alpha=0.4)
        plt.show()

    data = data - reflection * factor / 100. * ints_fit_pars[0]  

    if 'CAL_GHST' in header:  # Check for existence
        header['CAL_GHST'] = version
    else:
        header.set('CAL_GHST', version, 'ghost correction version py module (phifdt_pipe_modules.py)',after='CAL_DARK')
    
    return data, header

def phi_correct_ghost(data,header,rad,verbose=False): 
    '''
    Startup version on Jun 2021
    '''
    version = 'phi_correct_ghost V1.0 Jun 2021'

    only_one_vorbose = 1
    center = np.array([header['CRPIX1'],header['CRPIX2']]).astype(int)
    printc('        Read center from header (updated): x=',center[0],' y=',center[1],color=bcolors.OKBLUE)
    xd  = int(header['NAXIS1'])    
    yd  = int(header['NAXIS2'])    
    zd  = int(header['NAXIS3'])    
    PXBEG1  = int(header['PXBEG1']) - 1           
    PXEND1  = int(header['PXEND1']) - 1          
    PXBEG2  = int(header['PXBEG2']) - 1           
    PXEND2  = int(header['PXEND2']) - 1   

    # Center at: x= 1005  y= 953  radius= 517
    # printc('          Center at: x=',c[1],' y=',c[0],' radius=',radius,color=bcolors.OKBLUE)

    if verbose and only_one_vorbose:
        datap = np.copy(data)

    printc('-->>>>>>> Correcting ghost image ',color=bcolors.OKGREEN)

    #common part
    coef = [-1.98787669,1945.28944245] #empirically (first version)
    coef = [-1.9999,1942.7] #empirically (updated by trial and error)
    #coef = [-1.999,1944] #empirically (updated by trial and error)
    #TBD FINETUNING!!!!!

    #correct center:
    center_c = np.copy(center)
    center_c[0] += PXBEG1 
    center_c[1] += PXBEG2 
    poly1d_fn = np.poly1d(coef)
    sh = poly1d_fn(center_c).astype(int) 
    sh_float = poly1d_fn(center_c)
    printc('          image center: x: ',center[0],' y: ',center[1],color=bcolors.OKGREEN)
    printc('          image center [for 2048]: x: ',center_c[0],' y: ',center_c[1],color=bcolors.OKGREEN)
    printc('          ghost displacements: x: ',sh_float[0],' y: ',sh_float[1],color=bcolors.OKGREEN)

    #generate a ring mask out of the solar disk to see how much is the ghost image
    #we will be using complex histrograms....
    mask_anulus = bin_annulus([yd,xd],rad + 20, 10, full = False)
    mask_anulus = shift(mask_anulus, shift=(center[0]-xd//2,center[1]-yd//2), fill_value=0)

    idx = np.where(mask_anulus == 1)
    mask_anulus_big = bin_annulus([yd,xd],rad - 150, 100, full = False)
    mask_anulus_big = shift(mask_anulus_big, shift=(center[0]-xd//2,center[1]-yd//2), fill_value=0)
    idx_big = np.where(data[0,0,:,:]*mask_anulus_big == 1)

    printc('          computing azimuthal averages  ',color=bcolors.OKGREEN)

    centers = np.zeros((2,6))
    radius = np.zeros((6))
    ints = np.zeros((6,int(np.sqrt(xd**2+yd**2))))
    ints_rad = np.zeros((6,int(np.sqrt(xd**2+yd**2))))
    ints_fit = np.zeros((6,int(np.sqrt(xd**2+yd**2))))
    ints_syn = np.zeros((6,int(np.sqrt(xd**2+yd**2))))
    ints_fit_pars = np.zeros((6,5))
    factor = np.zeros((6,4))
    mean_intensity = np.zeros((6,4))

    # LOOP in wavelengths!
    for i in range(zd//4):

        # STEP --->>> average data in polarization for fitting the limb
        dummy_data = np.mean(data[i,:,:,:],axis=0)

        # STEP --->>> Find center of average data
        centers[1,i],centers[0,i],radius[i] = find_center(dummy_data)  #cy,cx....

        # STEP --->>> Generate CLV from averaged data
        intensity, rad = azimutal_average(dummy_data,[centers[0,i],centers[1,i]])
        ints[i,0:len(intensity)] = intensity
        ints_rad[i,0:len(intensity)] = rad

        # STEP --->>> FIT LIMB DATA
        rrange = int(radius[i] + 2) #2
        clv = ints[i,0:rrange]
        clv_r = ints_rad[i,0:rrange]
        mu = np.sqrt( (1 - clv_r**2/clv_r[-1]**2) )

        if verbose and only_one_vorbose:
            plt.plot(clv_r,clv)
            plt.xlabel('Solar radious [pixel]')
            plt.ylabel('Intensity [DN]')
            plt.show()

        u = 0.5
        I0 = 100
        ande = np.where(mu > 0.1)
        pars = newton(clv[ande],mu[ande],[I0,u,0.2,0.2,0.2],limb_darkening)
        fit, _ = limb_darkening(mu,pars)
        ints_fit[i,0:len(fit)] = fit
        ints_fit_pars[i,:] = pars

        # plt.plot(clv_r,clv,label='real clv')
        # plt.plot(clv_r,fit,label='fit order = 4')
        # plt.xlabel('Heliocentric angle ['+r'$\theta$]')
        # plt.ylabel('Intensity [DN]')
        # plt.legend()
        # plt.show()

        # STEP --->>> REPLICATE LIMB FIT WITH AVERAGE INTENSITY
        #Now there are two options. Either we use the fitted CLV or the azimutally averged CLV.
        #it is the second here but can be easily changed

        ints_syn[i,:] = ints[i,:]
        ints_syn[i,0:len(fit)] = fit

        # STEP --->>> NORMALIZE
        ints_syn[i,:] = ints_syn[i,:] / ints_fit_pars[i][0]
        ints_fit[i,:] = ints_fit[i,:] / ints_fit_pars[i][0]
        ints[i,:] = ints[i,:] / ints_fit_pars[i][0]

        # plt.plot(ints_fit[i,:],label='fitted clv')
        # plt.plot(ints[i,:],'.',label='real clv')
        # plt.plot(ints_syn[i,:],'--',label='synt clv')
        # plt.xlabel('Heliocentric angle ['+r'$\theta$]')
        # plt.ylabel('Intensity [DN]')
        # plt.legend()
        # plt.show()

        # STEP --->>> GENERATE GHOST

    
        nc = (PXEND2-PXBEG2+1)//2 #  center of frame
        limb_2d = np.zeros((PXEND2-PXBEG2+1,PXEND1-PXBEG1+1)) #generate image for ghost
        #fill the image with the revoluted fit
        s_of_gh = int(radius[i]*1.1)
        limb_2d[ nc - s_of_gh:nc + s_of_gh + 1, nc - s_of_gh:nc + s_of_gh + 1] = genera_2d(ints_syn[i,0:s_of_gh])
        xl,yl = limb_2d.shape
        # plt.imshow(limb_2d)
        # plt.show()

        # old down here
        # limb_real = genera_2d(intensity)
        # xl,yl = limb_real.shape
        # limb_2d = np.zeros((PXEND2-PXBEG2+1,PXEND1-PXBEG1+1))
        # limb_2d[ : , : ] = limb_real[xl//2 - yd//2 : xl//2 + yd//2 , yl//2 - xd//2:yl//2 + xd//2] #VERY SAME DATA

        # nuevorango = np.linspace(0,rad[-1]+2,len(rad))
        # f = interp1d(nuevorango, intensity)
        # intensity = f(rad)

        # STEP --->>> Smooth and SHIFT GHOST

        limb_2d = gaussian_filter(limb_2d, sigma=(8, 8)) 

        # #shift ghost to center of image
        # limb_2d = shift(limb_2d, shift=(int(centers[1,i])-yd//2,int(centers[0,i])-xd//2), fill_value=0)
        limb_2d = shift_subp(limb_2d, shift=[centers[1,i]-yd//2,centers[0,i]-xd//2])
        #OJO, shift_subp coge como parametros los de shift pero al reves!!!!!!
        if verbose and only_one_vorbose:
            plib.show_one(limb_2d,vmax=1,vmin=0,xlabel='pixel',ylabel='pixel',title='limb 2D',cbarlabel=' ',cmap='gray')
        # #shift to the position of the ghost
        reflection = shift(limb_2d, shift=(sh[0],sh[1]), fill_value=0) 

        # reflection = shift_subp(limb_2d, shift=[sh_float[1],sh_float[0]])
        if verbose and only_one_vorbose:
            plib.show_one(reflection,vmax=1,vmin=0,xlabel='pixel',ylabel='pixel',title='reflection',cbarlabel=' ',cmap='gray')
        #s_x,s_y,_ = PHI_shifts_FFT(data[i,:,:,:],prec=500,verbose=True,norma=False)

        # STEP --->>> Correct each modulation

        for j in range(4):

            dummy = data[i,j,:,:] 
            mean_intensity[i,j] = np.mean(dummy[idx_big])

            # FORMA 1 - con anillo
            values = dummy[idx].flatten() #Take the ring
            #show the histogram
            meanv = np.mean(values)
            idx_l = np.where(values <= meanv)
            m_l = np.mean(values[idx_l])
            idx_r = np.where(values >= meanv)
            m_r = np.mean(values[idx_r])
            factor[i,j] = (m_r - m_l) * 100. / ints_fit_pars[i][0] 
            print("factor",factor[i,j])

            # FORMA 2 - ajustando histogramas

            # y,x = np.histogram(values, bins=80)
            # x = x[:-1]
            # # weighted arithmetic mean (corrected - check the section below)
            # # m_l = sum(x * y) / sum(y)
            # # sigma = np.sqrt(sum(y * (x - m_l)**2) / sum(y))

            # sigma_l = m_l * 0.1
            # sigma_r = m_r * 0.1
            # popt,pcov = curve_fit(Gauss2, x, y, p0=[max(y),m_l,sigma_l,max(y),m_r,sigma_r])
            # _, m_l, _, _ ,m_r, _ = popt #unpack
            # factor[i,j] = (m_r - m_l) * 100. / ints_fit_pars[i][0] 
            # print("factor-n",factor[i,j])

            # plt.plot(x, y, 'b+:', label='data')
            # plt.plot(x, Gauss2(x, *popt), 'r-', label='fit')
            # plt.legend()
            # plt.title('Histo Fit wl: '+str(i)+' pol: '+str(j)+' factor: '+str(factor[i,j]))
            # plt.xlabel('DN')
            # plt.ylabel('Ad')
            # plt.show()

            # FORMA 3 - con cajitas
            # bsize = 20
            # box_right = data[i,j,int(centers[0,i]-bsize/2):int(centers[0,i]+bsize/2),int(centers[1,i]-radius[i]-bsize/2-20 ):int(centers[1,i]-radius[i]+bsize/2-20)]
            # box_left = data[i,j,int(centers[0,i]-bsize/2):int(centers[0,i]+bsize/2),int(centers[1,i]+radius[i]-bsize/2+20):int(centers[1,i]+radius[i]+bsize/2+20)]
            # # plib.show_one(box_left)
            # # plib.show_one(box_right)
            # m_l = np.mean(box_left)
            # m_r = np.mean(box_right)
            # factor[i,j] = (m_r - m_l) * 100. / ints_fit_pars[i][0] 
            # print("factor",factor[i,j])

            if verbose and only_one_vorbose:
                plt.hist(values, bins=40)
                plt.title('signal')
                plt.axvline(meanv, lw=2, color='yellow', alpha=0.4)
                plt.axvline(m_l, lw=2, color='red', alpha=0.4)
                plt.axvline(m_r, lw=2, color='blue', alpha=0.4)
                plt.axvline(factor[i,j]*ints_fit_pars[i][0] / 100., lw=2, color='green', alpha=0.4)
                plt.show()

            #sub-pixel shift to the position of the ghost
            #reflection = shift_subp(reflection, shift=[s_x[j],s_y[j]])
            #MINIMIZA EL ANILLO!!!
            #The problem is that V works but QU there is a residual at the disk. I use 0.9 though 

            # rms = np.zeros((100))
            # for k in range(100):
            #     dummy = data[i,j,:,:] - reflection * factor[i,j] / 100. * ints_fit_pars[i][0]  * (k+50)/100.#0.9 # * mean_intensity[i,j] / mean_intensity[i,0]
            #     rms[k] = np.std(dummy[idx])
            # plt.plot((np.arange(100)+50)/100.,rms)
            # plt.show()
            # cf = np.where(rms == np.min(rms))
            # data[i,j,:,:] = data[i,j,:,:] - reflection * factor[i,j] / 100. * ints_fit_pars[i][0]  * float((cf[0] + 50)/100.)#0.9 # * mean_intensity[i,j] / mean_intensity[i,0]
            # print('new',factor[i,j]*float((cf[0] + 50)/100.))

            data[i,j,:,:] = data[i,j,:,:] - reflection * factor[i,j] / 100. * ints_fit_pars[i][0]  # * mean_intensity[i,j] / mean_intensity[i,0]
            #data[i,j,:,:] = data[i,j,:,:] - reflection * factor[i,j] / 100. * ints_fit_pars[i][0]  * 0.9 # * mean_intensity[i,j] / mean_intensity[i,0]
            if verbose and only_one_vorbose:
                plib.show_two(datap[i,j,:,:],data[i,j,:,:],vmin=[0,0],vmax=[1,1],block=True,pause=0.1,title=['Before','After'],xlabel='Pixel',ylabel='Pixel')
                plt.plot(datap[0,0,0:200,200])
                plt.plot(data[0,0,0:200,200])
                plt.ylim([0, 5])
                plt.show()
                plt.plot(datap[0,0,200,0:200])
                plt.plot(data[0,0,200,0:200])
                plt.ylim([0, 5])
                plt.show()

            only_one_vorbose = 0
        
    if 'CAL_GHST' in header:  # Check for existence
        header['CAL_GHST'] = version
    else:
        header.set('CAL_GHST', version, 'ghost correction version py module (phifdt_pipe_modules.py)',after='CAL_DARK')
    
    return data, header

def phi_correct_ghost_dm(data,header,rad,verbose=False):
    '''
    Startup version on Jun 2021
    '''
    version = 'phi_correct_ghost_dm V1.0 Sep 2021 - appied to demodulated images'

    only_one_vorbose = 1
    center = np.array([header['CRPIX1'],header['CRPIX2']]).astype(int)
    printc('        Read center from header (updated): x=',center[0],' y=',center[1],color=bcolors.OKBLUE)
    xd  = int(header['NAXIS1'])    
    yd  = int(header['NAXIS2'])    
    zd  = int(header['NAXIS3'])    
    PXBEG1  = int(header['PXBEG1']) - 1           
    PXEND1  = int(header['PXEND1']) - 1          
    PXBEG2  = int(header['PXBEG2']) - 1           
    PXEND2  = int(header['PXEND2']) - 1   

    if verbose and only_one_vorbose:
        datap = np.copy(data)

    printc('-->>>>>>> Correcting ghost image ',color=bcolors.OKGREEN)

    #common part
    coef = [-1.98787669,1945.28944245] #empirically (first version)

    #find center of ghost (empirically determined - TBD Newton minimization):
    center_c = np.copy(center)
    center_c[0] += PXBEG1 
    center_c[1] += PXBEG2 
    poly1d_fn = np.poly1d(coef)
    sh = poly1d_fn(center_c).astype(int) 
    sh_float = poly1d_fn(center_c)
    printc('          image center: x: ',center[0],' y: ',center[1],color=bcolors.OKGREEN)
    printc('          image center [for 2048]: x: ',center_c[0],' y: ',center_c[1],color=bcolors.OKGREEN)
    printc('          ghost displacements: x: ',sh_float[0],' y: ',sh_float[1],color=bcolors.OKGREEN)

    #generate a ring mask out of the solar disk to see how much is the ghost image
    #we will be using complex histrograms....
    mask_anulus = bin_annulus([yd,xd],rad - 20, 10, full = False)
    mask_anulus = shift(mask_anulus, shift=(center[0]-xd//2,center[1]-yd//2), fill_value=0)

    idx = np.where(mask_anulus == 1)
    # mask_anulus_big = bin_annulus([yd,xd],rad - 150, 100, full = False)
    # mask_anulus_big = shift(mask_anulus_big, shift=(center[0]-xd//2,center[1]-yd//2), fill_value=0)
    # idx_big = np.where(data[0,0,:,:]*mask_anulus_big == 1)

    printc('          computing azimuthal averages  ',color=bcolors.OKGREEN)

    centers = np.zeros((2,6))
    radius = np.zeros((6))
    ints = np.zeros((6,int(np.sqrt(xd**2+yd**2))))
    ints_rad = np.zeros((6,int(np.sqrt(xd**2+yd**2))))
    ints_fit = np.zeros((6,int(np.sqrt(xd**2+yd**2))))
    ints_syn = np.zeros((6,int(np.sqrt(xd**2+yd**2))))
    ints_fit_pars = np.zeros((6,5))
    factor = np.zeros((6,4))
    mean_intensity = np.zeros((6,4))

    #do for QUV in continuum (first wavelength)
    dummy = data[0,1,:,:] 
    mean_intensity[0,1] = np.mean(dummy[idx])

    # FORMA 1 - con anillo
    values = dummy[idx].flatten() #Take the ring
    #show the histogram
    meanv = np.mean(values)
    idx_l = np.where(values <= meanv)
    m_l = np.mean(values[idx_l])
    idx_r = np.where(values >= meanv)
    m_r = np.mean(values[idx_r])
    factor[0,1] = (m_r - m_l)# * 100. / ints_fit_pars[i][0] 
    print("factor",factor[0,1])

    plt.hist(values, bins=40)
    plt.title('signal')
    plt.axvline(meanv, lw=2, color='yellow', alpha=0.4)
    plt.axvline(m_l, lw=2, color='red', alpha=0.4)
    plt.axvline(m_r, lw=2, color='blue', alpha=0.4)
    # plt.axvline(factor[0,1]*ints_fit_pars[0][0] / 100., lw=2, color='green', alpha=0.4)
    plt.axvline(factor[0,1], lw=2, color='green', alpha=0.4)
    plt.show()

    stop

    # LOOP in wavelengths!
    for i in range(zd//4):

        # STEP --->>> average data in polarization for fitting the limb
        dummy_data = np.mean(data[i,:,:,:],axis=0)

        # STEP --->>> Find center of average data
        centers[1,i],centers[0,i],radius[i] = find_center(dummy_data)  #cy,cx....

        # STEP --->>> Generate CLV from averaged data
        intensity, rad = azimutal_average(dummy_data,[centers[0,i],centers[1,i]])
        ints[i,0:len(intensity)] = intensity
        ints_rad[i,0:len(intensity)] = rad

        # STEP --->>> FIT LIMB DATA
        rrange = int(radius[i] + 2) #2
        clv = ints[i,0:rrange]
        clv_r = ints_rad[i,0:rrange]
        mu = np.sqrt( (1 - clv_r**2/clv_r[-1]**2) )

        if verbose and only_one_vorbose:
            plt.plot(clv_r,clv)
            plt.xlabel('Solar radious [pixel]')
            plt.ylabel('Intensity [DN]')
            plt.show()

        u = 0.5
        I0 = 100
        ande = np.where(mu > 0.1)
        pars = newton(clv[ande],mu[ande],[I0,u,0.2,0.2,0.2],limb_darkening)
        fit, _ = limb_darkening(mu,pars)
        ints_fit[i,0:len(fit)] = fit
        ints_fit_pars[i,:] = pars

        # STEP --->>> REPLICATE LIMB FIT WITH AVERAGE INTENSITY
        #Now there are two options. Either we use the fitted CLV or the azimutally averged CLV.
        #it is the second here but can be easily changed

        ints_syn[i,:] = ints[i,:]
        ints_syn[i,0:len(fit)] = fit

        # STEP --->>> NORMALIZE
        ints_syn[i,:] = ints_syn[i,:] / ints_fit_pars[i][0]
        ints_fit[i,:] = ints_fit[i,:] / ints_fit_pars[i][0]
        ints[i,:] = ints[i,:] / ints_fit_pars[i][0]

        # STEP --->>> GENERATE GHOST

    
        nc = (PXEND2-PXBEG2+1)//2 #  center of frame
        limb_2d = np.zeros((PXEND2-PXBEG2+1,PXEND1-PXBEG1+1)) #generate image for ghost
        #fill the image with the revoluted fit
        s_of_gh = int(radius[i]*1.1)
        limb_2d[ nc - s_of_gh:nc + s_of_gh + 1, nc - s_of_gh:nc + s_of_gh + 1] = genera_2d(ints_syn[i,0:s_of_gh])
        xl,yl = limb_2d.shape
        # plt.imshow(limb_2d)
        # plt.show()

        # old down here
        # limb_real = genera_2d(intensity)
        # xl,yl = limb_real.shape
        # limb_2d = np.zeros((PXEND2-PXBEG2+1,PXEND1-PXBEG1+1))
        # limb_2d[ : , : ] = limb_real[xl//2 - yd//2 : xl//2 + yd//2 , yl//2 - xd//2:yl//2 + xd//2] #VERY SAME DATA

        # nuevorango = np.linspace(0,rad[-1]+2,len(rad))
        # f = interp1d(nuevorango, intensity)
        # intensity = f(rad)

        # STEP --->>> Smooth and SHIFT GHOST

        #limb_2d = gaussian_filter(limb_2d, sigma=(8, 8)) 

        # #shift ghost to center of image
        # limb_2d = shift(limb_2d, shift=(int(centers[1,i])-yd//2,int(centers[0,i])-xd//2), fill_value=0)
        limb_2d = shift_subp(limb_2d, shift=[centers[1,i]-yd//2,centers[0,i]-xd//2])
        #OJO, shift_subp coge como parametros los de shift pero al reves!!!!!!
        if verbose and only_one_vorbose:
            plib.show_one(limb_2d,vmax=1,vmin=0,xlabel='pixel',ylabel='pixel',title='limb 2D',cbarlabel=' ',cmap='gray')
        # #shift to the position of the ghost
        reflection = shift(limb_2d, shift=(sh[0],sh[1]), fill_value=0) 

        # reflection = shift_subp(limb_2d, shift=[sh_float[1],sh_float[0]])
        if verbose and only_one_vorbose:
            plib.show_one(reflection,vmax=1,vmin=0,xlabel='pixel',ylabel='pixel',title='reflection',cbarlabel=' ',cmap='gray')
        #s_x,s_y,_ = PHI_shifts_FFT(data[i,:,:,:],prec=500,verbose=True,norma=False)

        # STEP --->>> Correct each modulation

        for j in range(4):

            dummy = data[i,j,:,:] 
            mean_intensity[i,j] = np.mean(dummy[idx_big])

            # FORMA 1 - con anillo
            values = dummy[idx].flatten() #Take the ring
            #show the histogram
            meanv = np.mean(values)
            idx_l = np.where(values <= meanv)
            m_l = np.mean(values[idx_l])
            idx_r = np.where(values >= meanv)
            m_r = np.mean(values[idx_r])
            factor[i,j] = (m_r - m_l) * 100. / ints_fit_pars[i][0] 
            print("factor",factor[i,j])


            if verbose and only_one_vorbose:
                plt.hist(values, bins=40)
                plt.title('signal')
                plt.axvline(meanv, lw=2, color='yellow', alpha=0.4)
                plt.axvline(m_l, lw=2, color='red', alpha=0.4)
                plt.axvline(m_r, lw=2, color='blue', alpha=0.4)
                plt.axvline(factor[i,j]*ints_fit_pars[i][0] / 100., lw=2, color='green', alpha=0.4)
                plt.show()

            #sub-pixel shift to the position of the ghost
            #reflection = shift_subp(reflection, shift=[s_x[j],s_y[j]])
            #MINIMIZA EL ANILLO!!!
            #The problem is that V works but QU there is a residual at the disk. I use 0.9 though 

            # rms = np.zeros((100))
            # for k in range(100):
            #     dummy = data[i,j,:,:] - reflection * factor[i,j] / 100. * ints_fit_pars[i][0]  * (k+50)/100.#0.9 # * mean_intensity[i,j] / mean_intensity[i,0]
            #     rms[k] = np.std(dummy[idx])
            # plt.plot((np.arange(100)+50)/100.,rms)
            # plt.show()
            # cf = np.where(rms == np.min(rms))
            # data[i,j,:,:] = data[i,j,:,:] - reflection * factor[i,j] / 100. * ints_fit_pars[i][0]  * float((cf[0] + 50)/100.)#0.9 # * mean_intensity[i,j] / mean_intensity[i,0]
            # print('new',factor[i,j]*float((cf[0] + 50)/100.))

            data[i,j,:,:] = data[i,j,:,:] - reflection * factor[i,j] / 100. * ints_fit_pars[i][0]  # * mean_intensity[i,j] / mean_intensity[i,0]
            #data[i,j,:,:] = data[i,j,:,:] - reflection * factor[i,j] / 100. * ints_fit_pars[i][0]  * 0.9 # * mean_intensity[i,j] / mean_intensity[i,0]
            if verbose and only_one_vorbose:
                plib.show_two(datap[i,j,:,:],data[i,j,:,:],vmin=[0,0],vmax=[1,1],block=True,pause=0.1,title=['Before','After'],xlabel='Pixel',ylabel='Pixel')
                plt.plot(datap[0,0,0:200,200])
                plt.plot(data[0,0,0:200,200])
                plt.ylim([0, 5])
                plt.show()
                plt.plot(datap[0,0,200,0:200])
                plt.plot(data[0,0,200,0:200])
                plt.ylim([0, 5])
                plt.show()

            only_one_vorbose = 1
        
        stop
    if 'CAL_GHST' in header:  # Check for existence
        header['CAL_GHST'] = version
    else:
        header.set('CAL_GHST', version, 'ghost correction version py module (phifdt_pipe_modules.py)',after='CAL_DARK')
    
    return data, header

def phi_correct_fringes(data,header,option,verbose=False):
    '''
    Startup version on Jun 2021
    24 Feb 2022: change int in freq by round since we we were floor rounding the pixel positions 
    '''
    version = 'phi_correct_fringes V1.0 Jun 2021'
    
    xd  = int(header['NAXIS1'])    
    yd  = int(header['NAXIS2'])    
    zd  = int(header['NAXIS3'])    

    if option == 'auto':
        printc('-->>>>>>> Looking for fringes and removing them --',color=bcolors.OKGREEN)

        freq_x = np.zeros((zd//4,3,50))  #Frecuencias Maximo de 10 ventanas
        freq_y = np.zeros((zd//4,3,50))
        freq_x2 = np.zeros((zd//4,3,50))  #Frecuencias Maximo de 10 ventanas
        freq_y2 = np.zeros((zd//4,3,50))
        rad_min = 10
        rad_max = 30
        wsize = 50
        wbin = 1

        win_halfw = 2 # for 0

        win = apod(win_halfw*2+1,0.6)
        x,y = np.ogrid[0:win_halfw*2 + 1, 0:win_halfw*2 + 1]
        level_theshold = [1.5,1.5,2]
        plt.ion()

        for i in range(zd//4):
            for j in np.arange(1,4):
                print('Wavelengh ',i,' pol state: ',j)
                data_fringes = rebin(data[i,j,:,:], [yd//wbin,xd//wbin])
                F=np.fft.fft2(data_fringes)
                F=np.fft.fftshift(F)
                h  = F.shape[0]
                w  = F.shape[1]
                #First compute FFT of single image
                power2d = np.log10( np.abs( (F*np.conj(F)).astype(np.float) ) ) 
                power2d = gaussian_filter(power2d, sigma=(1, 1)) 

                im = power2d[w//2-wsize:w//2+wsize+1,h//2-wsize:h//2+wsize+1]
                imc = im[2:-2,2:-2]
                # mean = np.mean(imc[wsize//2+5:,wsize//2+5:])
                # rms = np.std(imc[wsize//2+5:,wsize//2+5:])
                minimum = np.min(imc[wsize-rad_max:wsize+rad_max+1,wsize-rad_max:wsize+rad_max+1])
                mean = np.mean(imc[wsize-rad_max:wsize+rad_max+1,wsize-rad_max:wsize+rad_max+1] - minimum)
                rms = np.std(imc[wsize-rad_max:wsize+rad_max+1,wsize-rad_max:wsize+rad_max+1])

                stack = ( 
                    (im[2:-2,2:-2] > shift(im,[-2,-2])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[-2,-1])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[-2,0 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[-2,1 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[-2,2 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[-1,-2])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[-1,-1])[2:-2,2:-2])
                    * (im[2:-2,2:-2] > shift(im,[-1,0 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[-1,1 ])[2:-2,2:-2])
                    * (im[2:-2,2:-2] > shift(im,[-1,2 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[0 ,-2])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[0 ,-1])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[0 ,1 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[0 ,2 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[1 ,-2])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[1 ,-1])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[1 ,0 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[1 ,1 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[1 ,2 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[2 ,-2])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[2 ,-1])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[2 ,0 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[2 ,1 ])[2:-2,2:-2]) 
                    * (im[2:-2,2:-2] > shift(im,[2 ,2 ])[2:-2,2:-2]) 
                    )

                # mask = ( 
                #     (im[2:-2,2:-2] > shift(im,[-1,-1])[2:-2,2:-2])
                #     * (im[2:-2,2:-2] > shift(im,[-1,0 ])[2:-2,2:-2]) 
                #     * (im[2:-2,2:-2] > shift(im,[-1,1 ])[2:-2,2:-2])
                #     * (im[2:-2,2:-2] > shift(im,[0 ,-1])[2:-2,2:-2]) 
                #     * (im[2:-2,2:-2] > shift(im,[0 ,1 ])[2:-2,2:-2]) 
                #     * (im[2:-2,2:-2] > shift(im,[1 ,-1])[2:-2,2:-2]) 
                #     * (im[2:-2,2:-2] > shift(im,[1 ,0 ])[2:-2,2:-2]) 
                #     * (im[2:-2,2:-2] > shift(im,[1 ,1 ])[2:-2,2:-2]) 
                #     )

                idx = np.where(stack == 1)
                sm = imc.shape
                plt.imshow(imc)
                if len(idx[0]) > 0:
                    loop = 0
                    for idx_i in range(len(idx[0])):
                        if (imc[idx[0][idx_i],idx[1][idx_i]] - minimum ) > level_theshold[j-1]*mean:
                            if (np.abs(np.sqrt((idx[0][idx_i]-sm[0]//2)**2+(idx[1][idx_i]-sm[1]//2)**2)) > rad_min) and np.abs(np.sqrt((idx[0][idx_i]-sm[0]//2)**2+(idx[1][idx_i]-sm[1]//2)**2)) < rad_max:
                                plt.plot(idx[1][idx_i],idx[0][idx_i],"og",markersize=3)
                                subm = imc[idx[0][idx_i]-win_halfw:idx[0][idx_i]+win_halfw+1,idx[1][idx_i]-win_halfw:idx[1][idx_i]+win_halfw + 1]
                                if np.max(subm < 0):
                                    subm = 1 - subm
                                height, xcoor, ycoor, width_x, width_y = moments(subm)
                                freq_x2[i,j-1,loop] = (idx[0][idx_i] - win_halfw + xcoor - wsize + 2  )/h
                                freq_y2[i,j-1,loop] = (idx[1][idx_i] - win_halfw + ycoor - wsize + 2  )/w
                                freq_x[i,j-1,loop] = (idx[0][idx_i] - wsize + 2)/h
                                freq_y[i,j-1,loop] = (idx[1][idx_i] - wsize + 2)/w
                                # f_gauss = height * np.exp(-((x-xcoor)**2/(2*width_x**2) + (y-ycoor)**2/(2*width_y**2)))
                                f_gauss = 1 - np.exp(-((x-xcoor)**2/(2*(width_x*3)**2) + (y-ycoor)**2/(2*(width_y*3)**2)))
                                #plib.show_four_row(subm,f_gauss,subm - f_gauss,imc)
                                F[ idx[0][idx_i] + (h//2 - wsize + 2) - win_halfw : idx[0][idx_i]  + (h//2 - wsize + 2) + win_halfw + 1 , idx[1][idx_i]  + (w//2 - wsize + 2) - win_halfw : idx[1][idx_i]  + (w//2 - wsize + 2) + win_halfw + 1] *= f_gauss #win 1e-6#
                                power2d[ idx[0][idx_i] + (h//2 - wsize + 2) - win_halfw : idx[0][idx_i]  + (h//2 - wsize + 2) + win_halfw + 1, idx[1][idx_i]  + (w//2 - wsize + 2) - win_halfw : idx[1][idx_i]  + (w//2 - wsize + 2) + win_halfw + 1 ] *=  f_gauss #win 1e-6#
                                print(freq_x[i,j-1,loop],freq_y[i,j-1,loop])
                                print(i,j,level_theshold[j-1]*mean,3.*level_theshold[j-1]*mean, rms, 3*rms, imc[idx[0][idx_i],idx[1][idx_i]] - minimum,freq_x[i,j-1,loop],freq_y[i,j-1,loop])
                                loop += 1
                    plt.colorbar()
                    plt.show(block=True)
                    plt.pause(1)
                    plt.clf()
                    dum = np.copy(data_fringes)
                    data_fringes = np.fft.ifft2(np.fft.fftshift(F)).astype(np.float)
                    #plib.show_four_row(data_fringes,dum,dum-data_fringes,power2d,svmin=[-0.002,-0.002,-0.0002,-3],svmax=[0.002,0.002,0.0002,3])
                    data[i,j,:,:] = np.fft.ifft2(np.fft.fftshift(F)).astype(np.float)
        plt.ioff()

    #  fx = 0.012500000 ;16   PX/Size
    #  fy = 0.0070312503 ;9
    #  px = round(fx * sx)
    #  py = round(fy * sy)
        for i in range(zd//4):
            for j in np.arange(1,3):
                print(i,j,freq_y[i,j,:6],freq_x[i,j,:6])

        if 'CAL_FRIN' in header:  # Check for existence
            header['CAL_FRIN'] = version
        else:
            header.set('CAL_FRIN', version, 'Fringe correction ( name+version of py module if True )',after='CAL_DARK')

    #-----------------
    # FINGING - Need to include Ajusta senos contras!!!!
    #-----------------
    elif option == 'manual':
        printc('-->>>>>>> Removing fringes with fixed freq. --',color=bcolors.OKGREEN)
        printc('          ',version,'--',color=bcolors.OKGREEN)

        # printc('Freq. provided 9-July-2021 (H. Strecker and D. Orozco Suarez',color=bcolors.WARNING)
        # freq_x_Q = np.array([0.01328125,0.01328125]) 
        # freq_y_Q = np.array([0.00234375,0.00703125])
        # freq_x_U = np.array([0.01328125,0.01328125]) 
        # freq_y_U = np.array([0.00234375,0.00703125])
        # freq_x_V = np.array([0.01328125,0.01328125,0.0078125,0.01015625]) 
        # freq_y_V = np.array([0.00234375,0.00703125,0.0109375,0.00859375])
        printc('Freq. updated on 11-August-2021 (H. Strecker and D. Orozco Suarez',color=bcolors.WARNING)
        freq_x_Q = np.array([0.01318359375 ,0.01318359375]) 
        freq_y_Q = np.array([0.00195312500 ,0.00732421875])
        freq_x_U = np.array([0.01318359375 ,0.01318359375]) 
        freq_y_U = np.array([0.00195312500 ,0.00732421875])
        freq_x_V = np.array([0.01318359375 ,0.01318359375, 0.00976562500, 0.0078125000]) 
        freq_y_V = np.array([0.00195312500 ,0.00732421875, 0.00830078125, 0.0107421875])

        #freq to pixel yd,xd
        px_x_Q = freq_x_Q*xd
        px_y_Q = freq_y_Q*yd
        px_x_U = freq_x_U*xd
        px_y_U = freq_y_U*yd
        px_x_V = freq_x_V*xd 
        px_y_V = freq_y_V*yd
        #reflection
        # printc(px_x_Q,xd - px_x_Q,color=bcolors.OKBLUE)
        # printc( round(px_x_Q,(xd - px_x_Q) ),color=bcolors.OKBLUE)

        px_x_Q = np.append(px_x_Q,xd - px_x_Q - 1)
        px_y_Q = np.append(px_y_Q,yd - px_y_Q - 1)
        px_x_U = np.append(px_x_U,xd - px_x_U - 1)
        px_y_U = np.append(px_y_U,yd - px_y_U - 1)
        px_x_V = np.append(px_x_V,xd - px_x_V - 1)
        px_y_V = np.append(px_y_V,yd - px_y_V - 1)

        px_x_Q = np.round(px_x_Q).astype(int)
        px_y_Q = np.round(px_y_Q).astype(int)
        px_x_U = np.round(px_x_U).astype(int)
        px_y_U = np.round(px_y_U).astype(int)
        px_x_V = np.round(px_x_V).astype(int)
        px_y_V = np.round(px_y_V).astype(int)

        wsize = 50
        win_halfw = 2 

        printc('freq_x_Q [f,px] ',freq_x_Q,px_x_Q,color=bcolors.OKBLUE)
        printc('freq_y_Q [f,px] ',freq_y_Q,px_y_Q,color=bcolors.OKBLUE)
        printc('freq_x_U [f,px] ',freq_x_U,px_x_U,color=bcolors.OKBLUE)
        printc('freq_y_U [f,px] ',freq_y_U,px_y_U,color=bcolors.OKBLUE)
        printc('freq_x_V [f,px] ',freq_x_V,px_x_V,color=bcolors.OKBLUE)
        printc('freq_y_V [f,px] ',freq_y_V,px_y_V,color=bcolors.OKBLUE)
        printc('win_halfw ',win_halfw,color=bcolors.OKBLUE)

        mask_QUV = np.ones((3,yd,xd))
        maski,coords = generate_circular_mask([2*win_halfw,2*win_halfw],win_halfw,win_halfw)
        print(maski)
        print(KeyboardInterrupt)
        for k in range(len(px_x_Q)):
            # mask_QUV[0,px_y_Q[k]-win_halfw:px_y_Q[k]+win_halfw+1,px_x_Q[k]-win_halfw:px_x_Q[k]+win_halfw+1] = 1e-9
            print(k,px_y_Q[k]-win_halfw,px_y_Q[k]+win_halfw+1,px_x_Q[k]-win_halfw,px_x_Q[k]+win_halfw+1)
            mask_QUV[0,px_y_Q[k]-win_halfw:px_y_Q[k]+win_halfw+1,px_x_Q[k]-win_halfw:px_x_Q[k]+win_halfw+1] *= (1 - maski)
        for k in range(len(px_x_U)):
            mask_QUV[1,px_y_U[k]-win_halfw:px_y_U[k]+win_halfw+1,px_x_U[k]-win_halfw:px_x_U[k]+win_halfw+1] *= (1 - maski)
        for k in range(len(px_x_V)):
            mask_QUV[2,px_y_V[k]-win_halfw:px_y_V[k]+win_halfw+1,px_x_V[k]-win_halfw:px_x_V[k]+win_halfw+1] *= (1 - maski)
        # if verbose:
        # plib.show_one(mask_QUV[0,0:30,0:30])
        # plib.show_one(mask_QUV[1,0:30,0:30])
        # plib.show_one(mask_QUV[2,0:30,0:30])
        for i in range(zd//4):
            for j in np.arange(1,4):
                F = np.fft.fft2(data[i,j,:,:])
                # power2d = np.abs( (F*np.conj(F)).astype(np.float) ) 
                # power2d = gaussian_filter(power2d, sigma=(1, 1)) 
                # im = np.log10( power2d[0:30,0:30] )
                # plt.imshow(im,cmap='Greys')
                # plt.show()
                F *= mask_QUV[j - 1 , :, :]
                data[i,j,:,:]= np.fft.ifft2(F)

        if 'CAL_FRIN' in header:  # Check for existence
            header['CAL_FRIN'] = version
        else:
            header.set('CAL_FRIN', version, 'Fringe correction ( name+version of py module if True )',after='CAL_DARK')

    else:
        print('No fringe correction')
        return data, header
    #     PLT_RNG = 2
    #     plib.show_four_row(data_d[1,0,:,:],data_d[1,1,:,:],data_d[1,2,:,:],data_d[1,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003],block=False)
    #     plib.show_four_row(data[1,0,:,:],data[1,1,:,:],data[1,2,:,:],data[1,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003])
    #     plib.show_four_row(data_d[3,0,:,:],data_d[3,1,:,:],data_d[3,2,:,:],data_d[3,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003],block=False)
    #     plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003])
    #     plib.show_four_row(data_d[5,0,:,:],data_d[5,1,:,:],data_d[5,2,:,:],data_d[5,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003],block=False)
    #     plib.show_four_row(data[5,0,:,:],data[5,1,:,:],data[5,2,:,:],data[5,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003])

    # np.save('data_f',data)
    # return

    #     PLT_RNG = 2
    #     plib.show_four_row(data_d[1,0,:,:],data_d[1,1,:,:],data_d[1,2,:,:],data_d[1,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003],block=False)
    #     plib.show_four_row(data[1,0,:,:],data[1,1,:,:],data[1,2,:,:],data[1,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003])
    #     plib.show_four_row(data_d[3,0,:,:],data_d[3,1,:,:],data_d[3,2,:,:],data_d[3,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003],block=False)
    #     plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003])
    #     plib.show_four_row(data_d[5,0,:,:],data_d[5,1,:,:],data_d[5,2,:,:],data_d[5,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003],block=False)
    #     plib.show_four_row(data[5,0,:,:],data[5,1,:,:],data[5,2,:,:],data[5,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003])

    # np.save('data_f',data)
    # return

    return data, header

#GV need to pass outiput_dir too
def generate_level2(data,wave_axis,rte_mode,output_dir,ref_wavelenth = 6173.3354, milos_executable = None,options = None,center_line = True):

    '''
    generate_level2(data,wave_axis,rte_mode,milos_executable = MILOS_EXECUTABLE,options = None)
    return phi_rte(data,wave_axis,rte_mode,cmilos = cmd,options = options)

    '''
    if milos_executable != 'pmilos':
        try:
            print(milos_executable)
            CMILOS_LOC = os.path.realpath(__file__) 
            CMILOS_LOC = CMILOS_LOC[:CMILOS_LOC.rfind('/')-len(CMILOS_LOC)+1] + 'cmilos/'
            if os.path.isfile(CMILOS_LOC+milos_executable):
                printc("Cmilos executable located at:", CMILOS_LOC,color=bcolors.WARNING)
            else:
                raise ValueError('Cannot find cmilos:', CMILOS_LOC)
        except ValueError as err:
            printc(err.args[0],color=bcolors.FAIL)
            printc(err.args[1],color=bcolors.FAIL)
            return        

        cmd = CMILOS_LOC+"./"+milos_executable
        cmd = fix_path(cmd)
    elif milos_executable == 'pmilos':
        cmd = None
    else:
        printc("Using python milos wrapper",color=bcolors.WARNING)
        return 0
    ref_wavelenth = 6173.3354
    #OJO, REMOVE. NEED TO CHECK THE REF WAVE FROM S/C-PHI H/K
    if center_line !=False:
        if center_line == True:
            shift_w =  wave_axis[3] - ref_wavelenth
            wave_axis = wave_axis - shift_w
        elif type(center_line) is int:
            shift_w =  wave_axis[center_line] - ref_wavelenth
            wave_axis = wave_axis - shift_w

    #wave_axis = np.array([-300,-160,-80,0,80,160])/1000.+wavelength
    #wave_axis = np.array([-300,-140,-70,0,70,140])
    printc('   It is assumed the wavelength is given by the header info ')
    printc('         wave axis: ', wave_axis,color = bcolors.WARNING)
    printc('         wave axis (step):  ',(wave_axis - ref_wavelenth)*1000.,color = bcolors.WARNING)
        
    #GV need to pass outiput_dir too
    return phi_rte(data,wave_axis,rte_mode,output_dir,cmilos = cmd,options = options)
    

def create_output_filenames(filename, DID, version = '01'):
    """
    creating the L2 output filenames from the input, assuming L1
    - comes from HRT pipeline! - 
    - it does not use write_output_inversion for writing the opuput 
        since it heavely modified the header and FDT does it in every step
    """
    try:
        file_start = filename.split('solo_')[1] 
        file_start = 'solo_' + file_start
        L2_str = file_start.replace('L1', 'L2')
        versioned = L2_str.split('V')[0] + 'V' + version + '_' + DID + '.fits'
        stokes_file = versioned.replace('ilam', 'stokes')
        icnt_file = versioned.replace('ilam', 'icnt')
        bmag_file = versioned.replace('ilam', 'bmag')
        bazi_file = versioned.replace('ilam', 'bazi')
        binc_file = versioned.replace('ilam', 'binc')
        blos_file = versioned.replace('ilam', 'blos')
        vlos_file = versioned.replace('ilam', 'vlos')

        return stokes_file, icnt_file, bmag_file, bazi_file, binc_file, blos_file, vlos_file

    except Exception:
        print("The input file: {file_path} does not contain 'L1'")
        raise KeyError
