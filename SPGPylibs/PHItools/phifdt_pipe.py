# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% 
import numpy as np 
import os.path
#from astropy.io import fits as pyfits
from scipy.ndimage import gaussian_filter, rotate
import random
import statistics

from .tools import *
from .phi_fits import *
from .phi_gen import *
from .phi_reg import *
import SPGPylibs.GENtools.plot_lib as plib

# import scipy.integrate
# from scipy import rotate
# import scipy.ndimage as ndimage

def interpolateImages(image1, image2, dist1I, distI2):
    ''' interpolate 2D images - 
    '''
    imageInterp = (image1 * distI2 + image2 * dist1I) / (dist1I + distI2)
    return imageInterp

def applyPrefilter(data, wvltsData, prefilter, prefScale, wvltsPref, direction, scaledown=8,verbose = 0):
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

def applyPrefilter_dos(data, wvltsData, prefilter, prefScale, wvltsPref, direction, scaledown=8,verbose = 0):
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

def demod_phi(data,instrument,demod=False,verbose = 0):
    '''
    Use demodulation matrices to demodulate data size (n_wave*S_POL,N,M)
    ATTENTION: FDT40 is fixed to the one Johann is using!!!!
    '''

    if instrument == 'FDT40':
        mod_matrix_40 = np.array([[1.0006,-0.7132, 0.4002,-0.5693],
                    [1.0048, 0.4287,-0.7143, 0.5625],
                    [0.9963, 0.4269,-0.3652,-0.8229],
                    [0.9983,-0.4022, 0.9001, 0.1495]])
        demodM = np.linalg.inv(mod_matrix_40)
    # Johanns (it is the average in the central area of the one onboard)
        demodM  = np.array([[0.168258,      0.357277,     0.202212,     0.273266],\
            [-0.660351,     0.314981,     0.650029,    -0.299685],\
            [ 0.421242,     0.336994,    -0.183068,    -0.576202],\
            [-0.351933,     0.459820,    -0.582167,     0.455458]])
    elif instrument == 'FDT45':
        mod_matrix_45 = np.array([[1.0035,-0.6598, 0.5817,-0.4773],
                        [1.0032, 0.5647, 0.5275, 0.6403],
                        [0.9966, 0.4390,-0.5384,-0.7150],
                        [0.9968,-0.6169,-0.6443, 0.4425]])
        demodM = np.linalg.inv(mod_matrix_45)
    elif instrument == 'HRT40':
        mod_matrix_40 = np.array([[1.0040,-0.6647, 0.5928,-0.4527],
                        [1.0018, 0.5647, 0.5093, 0.6483],
                        [0.9964, 0.4348,-0.5135,-0.7325],
                        [0.9978,-0.6128,-0.6567, 0.4283]]) #HREW
        demodM = np.linalg.inv(mod_matrix_40)
    elif instrument == 'HRT45':
        mod_matrix_45_dos = np.array([[1.00159,-0.50032, 0.7093,-0.4931],
                        [1.0040, 0.6615, 0.3925, 0.6494],
                        [0.9954, 0.3356,-0.6126,-0.7143],
                        [0.9989,-0.7474,-0.5179, 0.4126]]) #MIA
        demodM = np.linalg.inv(mod_matrix_45_dos)
    else:
        printc('No demod available in demod_phi.py',color = bcolors.FAIL)
        raise SystemError()
    printc('Demodulation matrix for ', instrument,color = bcolors.WARNING)
    printc(demodM,color = bcolors.WARNING)
    if demod:
        return demodM
    ls,ps,ys,xs = data.shape
    for i in range(ls):
        data[i,:,:,:] = np.reshape(np.matmul(demodM,np.reshape(data[i,:,:,:],(ps,xs*ys))),(ps,ys,xs))
    return data

def crosstalk_ItoQUV(data_demod,verbose=0,npoints=2000):
    
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

def cross_talk_QUV(data,nran = 2000,nlevel=0.3,verbose=0,block=True):
    
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
    p = np.poly1d(cVU)

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
        print(i, ' of ', len(npatches),end='', flush=True)
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

def phifdt_pipe(data_f,dark_f,flat_f,instrument = 'FDT40',flat_c = True,dark_c = True,
    inner_radius = 250, outer_radius = 800, steps = 100, normalize = 0., flat_n = 1.,
    index = None, prefilter = 1, prefilter_fits = '0000990710_noMeta.fits',
    realign = False, verbose = True, outfile=None, mask_margin = 7, 
    individualwavelengths = False,correct_ghost = False,putmediantozero=False,
    vqu = False, do2d = 0, rte = False, debug = False,nlevel = 0.3,center_method=None):

    '''
    PHI-FDT naive data reduction pipeline (HRT pipeline is rather similar)
    The steps in the "naive" procesing are:
    1- Read data
    2- Check dimensions
    3- Read flats
    4- Read and correct dark field (taking into account the scaling)
    5- Find center of the Sun in the data
    6- get wavelength samppling from header
    7- move the continuum to the blue (if needed) in flat and data
    8- Correct flatfield
    9- Correct prefilter - needs prefilter data!
    10- realign data (not activated)
    11- Demodulate data using appropriate dem matrix
    12- correct cross-talk from I to QUV

    Parameters
    ----------
        Input:
    data_f : string
        Fits file of the raw FDT data  
    dark_f : string
        Fits file of a Valid dark file (processed dark).   
    flat_f : string
        Fits file of a Valid FDT flatfield.  
    ** Options:
    outfile = 'data_red.fits' : string
        File to store final processed data
    instrument = 'FDT40' : select the instrument and PMP temperature (for demod)
        -> implemented cases: -- 'FDT40','FDT45' --
    flat_c = True : default is to apply flat field correction to the data
    dark_c = True : default is to apply dark field correction to the data
    inner_radius = 250, outer_radius = 600, steps = 100 : initial values for finding sun center
    normalize = 0 : To normalize flats internally to the mean value of 5% of the disk (central) intensity  
    flat_n = 1.0 : flat scaling (flat = flat / flat_n) 
    index = None : in case you want a particular flat to be applied at another wave, e.g.,
        index = [5,1,2,3,4,0] exchange the first and last wave flats
        This is for testing stuff, mainly. 
    prefilter = 1 : To correct for the prefilter 
    prefilter_fits = '../RSW1/0000990710_noMeta.fits' : User should provide prefilter data fits file location
    realign = False : bool
        Realign all images before demodulating using FFT 
    individualwavelengths = False : bool
        Correct crosstalk from I to QUV for individual wavelengths

    Returns
    -------
    None 

    Raises
    ------

    References
    ----------
    
    Examples
    --------
    >>> import SPGPylibs as spg

    Notes
    -----
    This program is not optimized for speed. It assumes that input data 
        is 6 wavelength 
    '''

    PLT_RNG = 5

    printc('--------------------------------------------------------------',bcolors.OKGREEN)
    printc('PHI FDT data reduction software (for develping purposes only) ',bcolors.OKGREEN)
    printc('--------------------------------------------------------------',bcolors.OKGREEN)

    #-----------------
    # READ DATA
    #-----------------
    printc('-->>>>>>> Reading Data              ',color=bcolors.OKGREEN)
    printc('          DATA IS DIVIDED by 256.   ',color=bcolors.OKGREEN)
    #
    # PXBEG1  =                  385 ; First read-out pixel in dimension 1            
    # PXEND1  =                 1664 ; Last read-out pixel in dimension 1             
    # PXBEG2  =                  385 ; First read-out pixel in dimension 2            
    # PXEND2  =                 1664 ; Last read-out pixel in dimension 2             

    #  
    try:
        data, header = fits_get(data_f)
        printc('-->>>>>>> Reshaping data to [wave,Stokes,y-dim,x-dim] ',color=bcolors.OKGREEN)
        zd,yd,xd = data.shape
        data = np.reshape(data,(zd//4,4,yd, xd))
        data = data / 256.

    except Exception:
        printc("ERROR, Unable to open fits file: {}",data_f,color=bcolors.FAIL)

    if verbose:
        plib.show_one(data[0,0,:,:],vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Data first wave',cbarlabel='DN',save=None,cmap='gray')

    #-----------------
    # TAKE DATA DIMENSIONS 
    #-----------------
    PXBEG1  = int(header['PXBEG1']) - 1           
    PXEND1  = int(header['PXEND1']) - 1          
    PXBEG2  = int(header['PXBEG2']) - 1           
    PXEND2  = int(header['PXEND2']) - 1   
    if xd != (PXEND1 - PXBEG1 + 1) or yd != (PXEND2 - PXBEG2 + 1):
        printc('ERROR, Keyword dimensions and data array dimensions dont match ',color=bcolors.FAIL)
        raise SystemExit
    if xd < 2047:    
        printc('         data cropped to: [',PXBEG1,',',PXEND1,'],[',PXBEG2,',',PXEND2,']',color=bcolors.WARNING)
    
    #-----------------
    # READ FLAT FIELDS
    #-----------------

    if flat_c:
        printc('-->>>>>>> Reading Flats                    ',color=bcolors.OKGREEN)
        printc('          Assumes they are already normalized to ONE ',color=bcolors.OKGREEN)
        printc('          input should be [wave X Stokes,y-dim,x-dim].',color=bcolors.OKGREEN)
        try:
            flat,h = fits_get(flat_f)
            fz,fy,fx = flat.shape
            flat = np.reshape(flat,(fz//4,4,fy,fx))
            printc('-->>>>>>> Reshaping Flat to [wave,Stokes,y-dim,x-dim] ',color=bcolors.OKGREEN)
        except Exception:
            printc("ERROR, Unable to open flats file: {}",flat_f,color=bcolors.FAIL)
        if verbose:
            plib.show_one(flat[0,0,:,:],vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Flat first wave',cbarlabel='DN',save=None,cmap='gray')

    else:
        printc('-->>>>>>> No flats mode                    ',color=bcolors.WARNING)


    #-----------------
    # READ AND CORRECT DARK FIELD
    #-----------------
    if dark_c:

        printc('-->>>>>>> Reading Darks                   ',color=bcolors.OKGREEN)
        printc('          Input should be [y-dim,x-dim].',color=bcolors.OKGREEN)
        printc('          DARK IS DIVIDED by 256.   ',color=bcolors.OKGREEN)
        try:
            dark,h = fits_get(dark_f)
            dark = dark / 256.
        except Exception:
            printc("ERROR, Unable to open darks file: {}",dark_f,color=bcolors.FAIL)

        dy,dx = dark.shape
        dark_scale = fits_get(dark_f,scaling = True)
        data_scale = fits_get(data_f,scaling = True)
        if dark_scale["Present"][0] == data_scale["Present"][0]:
            scaling = dark_scale["scaling"][0] / data_scale["scaling"][0]
        else:
            scaling = dark_scale["scaling"][1] / data_scale["scaling"][1] * dark_scale["scaling"][0]

        if scaling != 1:
            printc('          checking scalling and correcting for it in the dark.',dark_scale,data_scale,scaling,color=bcolors.WARNING)
            dark = dark * scaling
        printc('-->>>>>>> Correcting dark current.',color=bcolors.OKGREEN)
        data = data - dark[np.newaxis,np.newaxis,PXBEG2:PXEND2+1,PXBEG1:PXEND1+1]
        data = np.abs(data)

        if verbose:
            plib.show_one(dark,vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Dark',cbarlabel='DN',save=None,cmap='gray')
            plib.show_one(data[0,0,:,:],vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Data after dark',cbarlabel='DN',save=None,cmap='gray')
    else:
        printc('-->>>>>>> No darks mode                    ',color=bcolors.WARNING)


    #-----------------
    # FIND DATA CENTER 
    #-----------------

    printc('-->>>>>>> finding the center of the solar disk (needed for masking) ',color=bcolors.OKGREEN)
    if center_method == 'Hough':
        c, r,threshold = find_circle_hough(data[0,0,:,:],inner_radius,outer_radius,steps,threshold = 0.01,normalize=False,verbose=False)
        c = np.roll(c,1)
    else:
        cx,cy,r=find_center(data[0,0,:,:])
        c = np.array([int(cx),int(cy)])
        r = int(r)

    printc('          Center at: x=',c[1],' y=',c[0],' radius=',r,color=bcolors.OKBLUE)
    #OJO.
    # find_circle_hough devuelve c = c[0] = x and c[1] = y !!!!!!!!!!!!!!
    # Esto viene porque en el KLL esta definido así (al reves) en la rutina votes()
 
    #-----------------
    # TAKE ONLY DISK WITH MARGIN
    #-----------------
    printc('-->>>>>>> Creating a mask for RTE with 7 px margin')
    size_of_mask = r - mask_margin
    rx = [int(c[1]-size_of_mask),int(c[1]+size_of_mask)]
    ry = [int(c[0]-size_of_mask),int(c[0]+size_of_mask)]
    mask,coords = generate_circular_mask([yd-1,xd-1],size_of_mask,size_of_mask)
    mask = shift(mask, shift=(c[1]-yd//2,c[0]-xd//2), fill_value=0)

    #-----------------
    # GET INFO ABOUT VOLTAGES/WAVELENGTHS, determine continuum and new flat
    #-----------------
    printc('-->>>>>>> Obtaining voltages from data ',color=bcolors.OKGREEN)
    wave_axis,voltagesData,tunning_constant,cpos = fits_get_sampling(data_f) 
    printc('          Data FG voltages: ',voltagesData,color=bcolors.OKBLUE)
    printc('          Continuum position at wave: ', cpos,color=bcolors.OKBLUE)
    printc('          Data wave axis [mA]: ',wave_axis,color=bcolors.OKBLUE)
    rubish = (voltagesData-np.roll(voltagesData,-1))/tunning_constant/1e4
    sampling = np.mean(np.sort(np.abs(rubish[0:-2])))
    printc('          Data sampling [mA]: ',sampling,color=bcolors.OKBLUE)

    #-----------------
    # ROLL DATA IF CONTINUUM IS IN DIFFERENT POSITION 
    # Probably before demodulation and flat is better!!!!
    #-----------------
    if cpos != 0:
        datar = np.copy(data)
        voltagesDatar = np.copy(voltagesData)
        wave_axisr = np.copy(wave_axis)
        if voltagesData[cpos] < voltagesData[0]:
        #continuum is on the right but it is in the blue!!!!
            printc('          Rolling data to move continuum from right (red) to left (blue)',color=bcolors.WARNING)
            for i in range(zd//4):
                #print((i+1)%(zd//4),i%(zd//4),i,(zd//4))
                datar[(i+1)%(zd//4),:,:,:] = data[i%(zd//4),:,:,:] # np.roll(data, 4, axis=0)
                voltagesDatar[(i+1)%(zd//4)] = voltagesData[i%(zd//4)] # np.roll(data, 4, axis=0)
                wave_axisr[(i+1)%(zd//4)] = wave_axis[i%(zd//4)] # np.roll(data, 4, axis=0)
        data = np.copy(datar)
        voltagesData = np.copy(voltagesDatar)
        wave_axis = np.copy(wave_axisr)
        del datar
        del voltagesDatar
        del wave_axisr
        cpos = 0
        printc('          New FG voltages: ',voltagesData,color=bcolors.OKBLUE)
        printc('          NEW continuum position at wave: ', cpos,color=bcolors.OKBLUE)
        printc('          NEW data wave axis [mA]: ',wave_axis,color=bcolors.OKBLUE)
    
    #READ FLAT ORIGINAL FILE TODO: mantener la cabecera de los flats!!!!!!
    #info = phi.fits_read('../RSW1/Add-data/solo_L0_phi-fdt-ilam_20200618T035946_V202007101227C_0066180100.fits',head=3) 
  
    if flat_c:
        printc('-->>>>>>> Obtaining voltages from flats ',color=bcolors.OKGREEN)
        #ff =  '../Nov-2020-STP122/solo_L0_phi-fdt-flat_0645767986_V202012091123I_0066181100.fits'

        wave_axis_f,voltagesFlat,tunning_constant_f,cpos_f = fits_get_sampling(flat_f) 
        printc('          FLAT FG voltages: ',voltagesFlat,color=bcolors.OKBLUE)
        printc('          FLAT Continuum position at wave: ', cpos_f,color=bcolors.OKBLUE)
        printc('          FLAT wave axis [mA]: ',wave_axis_f,color=bcolors.OKBLUE)
        rubish = (voltagesFlat-np.roll(voltagesFlat,-1))/tunning_constant/1e4
        sampling_f = np.mean(np.sort(np.abs(rubish[0:-2])))
        printc('          FLAT sampling [mA]: ',sampling_f,color=bcolors.OKBLUE)

        printc('-->>>>>>> Reshaping flat to [wave,Stokes,y-dim,x-dim]',color=bcolors.OKGREEN)
        flat = np.reshape(flat,(fz//4,4,fy, fx))

    #-----------------
    # ROLL FLAT IF CONTINUUM IS IN DIFFERENT POSITION 
    # Probably before demodulation and flat is better!!!!
    #-----------------
    if flat_c:
        if cpos_f != 0:
            flatr = np.copy(flat)
            voltagesFlatr = np.copy(voltagesFlat)
            wave_axis_fr = np.copy(wave_axis_f)
            if voltagesFlat[cpos_f] < voltagesFlat[0]:
            #continuum is on the right but it is in the blue!!!!
                printc('          Rolling flat to move continuum from right (red) to left (blue)',color=bcolors.WARNING)
                for i in range(fz//4):
                    #print((i+1)%(fz//4),i%(fz//4),i)
                    flatr[(i+1)%(fz//4),:,:,:] = flat[i%(fz//4),:,:,:] # np.roll(data, 4, axis=0)
                    voltagesFlatr[(i+1)%6] = voltagesFlat[i%(fz//4)] # np.roll(data, 4, axis=0)
                    wave_axis_fr[(i+1)%(zd//4)] = wave_axis_f[i%(zd//4)] # np.roll(data, 4, axis=0)
            flat = np.copy(flatr)
            voltagesFlat = np.copy(voltagesFlatr)
            wave_axis_f = np.copy(wave_axis_fr)
            del flatr
            del voltagesFlatr
            del wave_axis_fr
            cpos_f = 0
            printc('          New Flat FG voltages: ',voltagesFlat,color=bcolors.OKBLUE)
            printc('          NEW Flat continuum position at wave: ', cpos_f,color=bcolors.OKBLUE)
            printc('          NEW Flat data wave axis [mA]: ',wave_axis_f,color=bcolors.OKBLUE)

    # TODO: INTERPOLATE THE FLAT TO MATCH THE WAVELENG (CAVITY)

    # from scipy.interpolate import RegularGridInterpolator
    # x = np.linspace(0,2047,2048).astype(int)
    # y = np.linspace(0,2047,2048).astype(int)
    # z = np.array([-300.,-140.,-70.,0.,70.,140.,300.]) #ojo con el -300
    # zn = np.array([-175.,-140.,-105.,-70.,-35.,0.,35.,70.,105.,140.,175.,300.])

    # flat_rsw1 = np.concatenate(((flat_rsw1[5,:,:,:])[np.newaxis,:,:,:],flat_rsw1))
    # fn = RegularGridInterpolator((z,y,x), flat_rsw1[:,0,:,:])
    # pts = np.array([-40,10,10])
    # print(fn(pts))

    #     pts = np.meshgrid(-40.,y,x)
    # pts = np.array([m.flatten() for m in pts])
    # flat_n = fn(pts.T)
    # result = flat_n.reshape((2048,2048))
    # plt.imshow(result,vmin=0.9,vmax=1.1)

    # flat_n = np.zeros((12,4,2048,2048))

    # for i in range(4):
    #   fn = RegularGridInterpolator((z,y,x), flat_rsw1[:,i,:,:],bounds_error=False)
    #   for j in range(12):
    #     print(i,zn[j])
    #     pts_list = np.meshgrid(zn[j],y,x)
    #     pts = np.array([m.flatten() for m in pts_list])
    #     flat_n[j,i,:,:] = fn(pts.T).reshape((2048,2048))

    #-----------------
    # APPLY FLAT CORRECTION 
    # TODO: TAKE THE REAL FLAT 
    #-----------------
    factor = 0.05
    rrx = [int(c[1]-r*factor),int(c[1]+r*factor)]
    rry = [int(c[0]-r*factor),int(c[0]+r*factor)]

    if flat_c:
        printc('-->>>>>>> Correcting Flatfield',color=bcolors.OKGREEN)
        try:
            if (len(index)) == 6:
                print('          Changing flat index to ',index)
        except:
            index = [0,1,2,3,4,5]
        for p in range(4):
            for l in range(int(zd//4)):
                print('          ... pol: ',p,' wave: ',l,' index: ',index[l])
                dummy_flat = (flat[index[l],p,PXBEG2:PXEND2+1,PXBEG1:PXEND1+1]/flat_n)
                if normalize ==1:
                    print('          normalizing flats using region x = [',rrx[0],':',rrx[1],'] y = ]',rry[0],':',rry[1],']')
                    mm = np.mean(dummy_flat[rry[0]:rry[1],rrx[0]:rrx[1]])
                    dummy_flat = dummy_flat / mm
                data[l,p,:,:] = data[l,p,:,:]/dummy_flat
        del dummy_flat

        if verbose:
            plib.show_one(data[cpos,0,:,:],vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Data / flat at continuum',cbarlabel='DN',save=None,cmap='gray')

    #-----------------
    # CORRECT PREFILTER 
    #-----------------
    if prefilter == 1:
        printc('-->>>>>>> Read prefilter and correct for it ')
        printc('          ',prefilter_fits,'   ')

        prefdata,h = fits_get(prefilter_fits)
        prefdata = prefdata.astype(fmt)
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
        if verbose == 1:
            datap = np.copy(data)
        data = applyPrefilter_dos(data, voltagesData, prefdata, prefscale, prefVoltages, -1,scaledown = scale, verbose = verbose)
        if verbose == 1:
            plt.plot(data[:,0,yd//2,xd//2],'o-',label='corrected')
            plt.plot(datap[:,0,yd//2,xd//2],'--',label='original')
            plt.legend()
            plt.show()
            del datap
        printc('\n')

    #-----------------
    # REALIGN DATA BEFORE DEMODULATION
    #-----------------

    if realign:
        printc('-->>>>>>> Realigning data...           ',color=bcolors.OKGREEN)
        for i in range(zd//4):
            s_x,s_y,_ = PHI_shifts_FFT(data[i,:,:,:],prec=500,verbose=True,norma=False)
            for j in range(4):
                data[i,j,:,:] = shift_subp(data[i,j,:,:], shift=[s_x[j],s_y[j]])

    #-----------------
    # APPLY DEMODULATION 
    #-----------------
    printc('-->>>>>>> Demodulating data...         ',color=bcolors.OKGREEN)
    if debug:
        datan = np.copy(data)
        ds = np.copy(data)
        demodM  = np.array([[0.168258,      0.357277,     0.202212,     0.273266],\
            [-0.660351,     0.314981,     0.650029,    -0.299685],\
            [ 0.421242,     0.336994,    -0.183068,    -0.576202],\
            [-0.351933,     0.459820,    -0.582167,     0.455458]])
        for i in range(zd//4):
            for l in range(xd):
                for m in range(yd):
                    datan[i,:,m,l] = np.matmul(demodM, ds[i,:,m,l] )
        plib.show_four_row(datan[3,0,:,:],datan[3,1,:,:],datan[3,2,:,:],datan[3,3,:,:],svmin=[0,-0.2,-0.2,-1.],svmax=[100,0.1,0.1,0.1])

    data = demod_phi(data,instrument)

    if verbose == 1:
        plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'])

    #-----------------
    # APPLY NORMALIZATION 
    #-----------------
    printc('-->>>>>>> Applying normalization --',color=bcolors.OKGREEN)
    nrm = np.mean(data[cpos,0,rry[0]:rry[1],rrx[0]:rrx[1]])
    print('          Norma is: ',nrm,' evaluated in x = [',rrx[0],':',rrx[1],'] y = [',rry[0],':',rry[1],']')
    data = data/nrm
    if debug:
        datan = datan/nrm

    if debug:
        plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:])
        plib.show_four_row(datan[3,0,:,:],datan[3,1,:,:],datan[3,2,:,:],datan[3,3,:,:])

    if verbose == 1:
        plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'])

    # %%

    #-----------------
    # CROSS-TALK CALCULATION 
    #-----------------
    printc('-->>>>>>> Cross-talk correction from Stokes I to Stokes Q,U,V --',color=bcolors.OKGREEN)
    factor = 0.80 # 80% of the disk
    printc('          Using ',factor*100,'% of the disk                     ',color=bcolors.OKGREEN)

    rrx = [int(c[1]-r*factor),int(c[1]+r*factor)]
    rry = [int(c[0]-r*factor),int(c[0]+r*factor)]

    if individualwavelengths:
        for i in range(zd//4):
            printc('          Crosstalk evaluated in x = [',rrx[0],':',rrx[1],'] y = [',rry[0],':',rry[1],']',' using ',factor*100,"% of the disk",color=bcolors.OKBLUE)
            printc('          Individual wavelengths....',color=bcolors.OKBLUE)
            broadcastd = data[i,:,rry[0]:rry[1],rrx[0]:rrx[1]]
            data_dummy = data[:,:,rry[0]:rry[1],rrx[0]:rrx[1]]*0. + broadcastd[np.newaxis,:,:,:]
            cQ,cU,cV = crosstalk_ItoQUV(data_dummy[:,:,rry[0]:rry[1],rrx[0]:rrx[1]],npoints=10000,verbose=verbose)
            #-----------------
            # CROSS-TALK CORRECTION 
            #-----------------
            data[i,1,:,:] = data[i,1,:,:] - cQ[0]*data[i,0,:,:] - cQ[1]
            data[i,2,:,:] = data[i,2,:,:] - cU[0]*data[i,0,:,:] - cU[1]
            data[i,3,:,:] = data[i,3,:,:] - cV[0]*data[i,0,:,:] - cV[1]
        if verbose:
            plt.hist(data[4,1,900:1100,900:1100].flatten(), bins='auto')
            plt.hist(data[4,2,900:1100,900:1100].flatten(), bins='auto')
            plt.hist(data[4,3,900:1100,900:1100].flatten(), bins='auto')
            plt.show()
    else:
        printc('          Crosstalk evaluated in x = [',rrx[0],':',rrx[1],'] y = [',rry[0],':',rry[1],']',' using ',factor*100,"% of the disk",color=bcolors.OKBLUE)
        cQ,cU,cV = crosstalk_ItoQUV(data[:,:,rry[0]:rry[1],rrx[0]:rrx[1]],verbose=True,npoints=10000)
    #-----------------
    # CROSS-TALK CORRECTION 
    #-----------------
    printc('          Applying cross-talk correction...',color=bcolors.OKGREEN)
    data[:,1,:,:] = data[:,1,:,:] - cQ[0]*data[:,0,:,:] - cQ[1]
    data[:,2,:,:] = data[:,2,:,:] - cU[0]*data[:,0,:,:] - cU[1]
    data[:,3,:,:] = data[:,3,:,:] - cV[0]*data[:,0,:,:] - cV[1]
    if verbose:
        plt.hist(data[4,1,900:1100,900:1100].flatten(), bins='auto')
        plt.hist(data[4,2,900:1100,900:1100].flatten(), bins='auto')
        plt.hist(data[4,3,900:1100,900:1100].flatten(), bins='auto')
        plt.show()
        plib.show_four_row(data[2,0,:,:],data[2,1,:,:],data[2,2,:,:],data[2,3,:,:],title=['I','Q','U','V'])

    if correct_ghost:
        printc('-->>>>>>> Correcting ghost image ',color=bcolors.OKGREEN)

        coef = [-1.98787669,1945.28944245] #empirically
        poly1d_fn = np.poly1d(coef)
        sh = poly1d_fn(c).astype(int) 

        #fit limb data!!!!!
        # generate image with fit
        # shift fit
        # remove fit

        # # reflection = image[4] - phi_f.shift(image[4], shift=sh) * 0.004
        # reflection,coords = generate_circular_mask([yd-1,xd-1],r-15,r-15) #mascara del mismo radio
        # #la centro en la imagen
        # reflection = phi.shift(reflection, shift=(c[1]-yd//2,c[0]-xd//2), fill_value=0)
        # #la desplazo para que coja la imagen fantasma
        # reflection = phi.shift(reflection, shift=(sh[1]-1024,sh[0]-1024), fill_value=0)
        #     
        # corner1 = np.mean(data[0,3,int(c[1]+r-5):int(c[1]+r),int(c[0]+r-5):int(c[0]+r)])
        # corner2 = np.mean(data[0,3,int(c[1]+r-5):int(c[1]+r),int(c[0]-r-5):int(c[0]-r)])
        # corner3 = np.mean(data[0,3,int(c[1]-r-5):int(c[1]-r),int(c[0]+r-5):int(c[0]+r)])
        # corner4 = np.mean(data[0,3,int(c[1]-r-5):int(c[1]-r),int(c[0]-r-5):int(c[0]-r)])
        # cor = np.array([corner1,corner2,corner3,corner4])
        # max_index = np.where(cor == np.amax(cor))
        # min_index = np.where(cor == np.amin(cor))
        # factor = cor[max_index] - cor[min_index]
        # factor = 0.0025
        # dd = data[0,3,:,:] - reflection * factor
        # plt.imshow(dd,vmin=-0.001,vmax=0.003)
        # plt.colorbar()
        # plt.imshow(data[0,3,:,:] - reflection*0.0025,vmin=-0.05,vmax=0.05,cmap="gray")

        blurred_Q = gaussian_filter(data[0,1,:,:], sigma=(5, 5))
        blurred_U = gaussian_filter(data[0,2,:,:], sigma=(5, 5))
        blurred_V = gaussian_filter(data[0,3,:,:], sigma=(5, 5))
        for i in range(zd//4):
            data[i,1,:,:] = data[i,1,:,:] - blurred_Q
            data[i,2,:,:] = data[i,2,:,:] - blurred_U
            data[i,3,:,:] = data[i,3,:,:] - blurred_V

        if verbose:
            plib.show_one(data[2,3,:,:],vmin=-0.006,vmax=0.006)

    if putmediantozero:
        printc('-->>>>>>> Putting median to zero ',color=bcolors.OKGREEN)
        PQ = np.median(data[:,1,rry[0]:rry[1],rrx[0]:rrx[1]],axis=(1,2))
        PU = np.median(data[:,2,rry[0]:rry[1],rrx[0]:rrx[1]],axis=(1,2))
        PV = np.median(data[:,3,rry[0]:rry[1],rrx[0]:rrx[1]],axis=(1,2))
        data[:,1,:,:] = data[:,1,:,:] - PQ[:,np.newaxis,np.newaxis]
        data[:,2,:,:] = data[:,2,:,:] - PU[:,np.newaxis,np.newaxis]
        data[:,3,:,:] = data[:,3,:,:] - PV[:,np.newaxis,np.newaxis]
        printc(PQ,PU,PV)
    if verbose == 1:
        plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'])

    #-----------------
    # CROSS-TALK CALCULATION FROM V TO QU (Interactive)
    #-----------------
    
    if vqu:
        printc('-->>>>>>> Cross-talk correction from Stokes V to Stokes Q,U ',color=bcolors.OKGREEN)

        factor = 0.3 # 30% of the disk
        rrx = [int(c[1]-r*factor),int(c[1]+r*factor)]
        rry = [int(c[0]-r*factor),int(c[0]+r*factor)]
        print(' Cross-talk evaluated in x = [',rrx[0],':',rrx[1],'] y = [',rry[0],':',rry[1],']',' using ',factor*100,"% of the disk")
        cVQ,cVU = cross_talk_QUV(data[:,:,rry[0]:rry[1],rrx[0]:rrx[1]],nran = 2000,nlevel=nlevel, show = 1,block=False)

        option = input('Do you want to apply the correction (y/n) [n]: ')
        if option == 'y':
            datao = np.copy(data)
            print('Applying V to QU cross-talk correction...')
            datao[:,2,:,:] = data[:,2,:,:] - cVQ[0]*data[:,3,:,:] - cVQ[1]
            datao[:,3,:,:] = data[:,3,:,:] - cVU[0]*data[:,3,:,:] - cVU[1]
            plib.show_two(data[3,1,ry[0]:ry[1],rx[0]:rx[1]],datao[3,1,ry[0]:ry[1],rx[0]:rx[1]],block=False,title=['Stokes Q','Stokes Q corrected'])
            plib.show_two(data[3,2,ry[0]:ry[1],rx[0]:rx[1]],datao[3,2,ry[0]:ry[1],rx[0]:rx[1]],block=False,title=['Stokes U','Stokes U corrected'])
            option2 = input('Do you wnat to continue (y/n) [n]: ')
            if option2 == 'y':
                data = np.copy(datao)
                del datao
        plt.close()

    if do2d >= 2:
        printc('---------------------------------------------------------',color=bcolors.OKGREEN)
        printc('-- IN 2-Dimensions                                     --')
        printc('-- Cross-talk correction from Stokes I to Stokes Q,U,V --')
        printc('---------------------------------------------------------',color=bcolors.OKGREEN)
        size = do2d
        cV0,cV1 = crosstalk_ItoQUV2d(data[:,:,ry[0]:ry[1],rx[0]:rx[1]],size=size)
        nsize = size-1
        dim = ry[1]-ry[0]-nsize
        cV0 = cV0.reshape(dim,dim)
        cV1 = cV1.reshape(dim,dim)
        plib.show_one(cV0,vmin=-0.005,vmax=0.005)
        data[:,3,ry[0]+nsize//2:ry[1]-nsize//2-1,rx[0]+nsize//2:rx[1]-nsize//2-1] = \
            data[:,3,ry[0]+nsize//2:ry[1]-nsize//2-1,rx[0]+nsize//2:rx[1]-nsize//2-1] -\
            cV0*data[:,0,ry[0]+nsize//2:ry[1]-nsize//2-1,rx[0]+nsize//2:rx[1]-nsize//2-1] #- 0.95*cV1

    #-----------------
    #CHECK FOR INFs
    #-----------------

    data[np.isinf(data)] = 0
    data[np.isnan(data)] = 0

    #-----------------
    # SAVE DATA TODO: CMILOS FORMAT AND FITS
    #-----------------
    printc('---------------------------------------------------------',color=bcolors.OKGREEN)
    if outfile == None:
        outfile = data_f[:-4]+'_red.fits'
    printc(' Saving data to: ',outfile)

    hdu = pyfits.PrimaryHDU(data)
    hdul = pyfits.HDUList([hdu])
    hdul.writeto(outfile, overwrite=True)

    #-----------------
    # INVERSION OF DATA WITH CMILOS
    #-----------------

    if rte == True:
        printc('---------------------------------------------------------',color=bcolors.OKGREEN)

        try:
            CMILOS_LOC = os.path.realpath(__file__) 
            CMILOS_LOC = CMILOS_LOC[:-14] + 'cmilos/'
            if os.path.isfile(CMILOS_LOC+'milos'):
                printc("Cmilos executable located at:", CMILOS_LOC,color=bcolors.WARNING)
            else:
                raise ValueError('Cannot find cmilos:', CMILOS_LOC)
        except ValueError as err:
            printc(err.args[0],color=bcolors.FAIL)
            printc(err.args[1],color=bcolors.FAIL)
            return        

        wavelength = 6173.3356
        # wave_axis = np.array([-300,-160,-80,0,80,160])/1000.+wavelength
        # wave_axis = np.array([-300,-140,-70,0,70,140])
        printc('   It is assumed the wavelength is given by the header info ')
        printc(wave_axis,color = bcolors.WARNING)
        printc((wave_axis - wavelength)*1000.,color = bcolors.WARNING)
        printc('   saving data into dummy_in.txt for RTE input')

        sdata = data[:,:,ry[0]:ry[1],rx[0]:rx[1]]
        l,p,x,y = sdata.shape
        filename = 'dummy_in.txt'
        with open(filename,"w") as f:
            for i in range(x):
                for j in range(y):
                    for k in range(l):
                        f.write('%e %e %e %e %e \n' % (wave_axis[k],sdata[k,0,j,i],sdata[k,1,j,i],sdata[k,2,j,i],sdata[k,3,j,i]))
        del sdata

        printc('  ---- >>>>> Inverting data.... ',color=bcolors.OKGREEN)
        umbral = 3.

        import subprocess
        cmd = CMILOS_LOC+"./milos"
        cmd = fix_path(cmd)
        rte_on = subprocess.call(cmd+" 6 15 0 0 dummy_in.txt  >  dummy_out.txt",shell=True)
        print(rte_on)
        printc('  ---- >>>>> Finishing.... ',color=bcolors.OKGREEN)
        printc('  ---- >>>>> Reading results.... ',color=bcolors.OKGREEN)
        #del_dummy = subprocess.call("rm dummy_in.txt",shell=True)
        #print(del_dummy)

        res = np.loadtxt('dummy_out.txt')
        npixels = res.shape[0]/12.
        print(npixels)
        print(npixels/x)
        result = np.zeros((12,y*x)).astype(float)
        rte_invs = np.zeros((12,yd,xd)).astype(float)
        for i in range(y*x):
            result[:,i] = res[i*12:(i+1)*12]
        result = result.reshape(12,y,x)
        result = np.einsum('ijk->ikj', result)
        rte_invs[:,ry[0]:ry[1],rx[0]:rx[1]] = result
        del result
        rte_invs_noth = np.copy(rte_invs)

        noise_in_V =  np.mean(data[0,3,rry[0]:rry[1],rrx[0]:rrx[1]])
        low_values_flags = np.max(np.abs(data[:,3,:,:]),axis=0) < noise_in_V*umbral  # Where values are low
        rte_invs[2,low_values_flags] = 0
        rte_invs[3,low_values_flags] = 0
        rte_invs[4,low_values_flags] = 0
        for i in range(12):
            rte_invs[i,:,:] = rte_invs[i,:,:] * mask
        #save plots!!!!
        if verbose:
            plib.show_four_row(rte_invs_noth[2,:,:],rte_invs_noth[3,:,:],rte_invs_noth[4,:,:],rte_invs_noth[8,:,:],svmin=[0,0,0,-3.],svmax=[1200,180,180,+3.],title=['Field strengh [Gauss]','Field inclination [degree]','Field azimuth [degree]','LoS velocity [km/s]'],xlabel='Pixel',ylabel='Pixel',save=outfile+'_VLoS.png')
            plib.show_four_row(rte_invs[2,:,:],rte_invs[3,:,:],rte_invs[4,:,:],rte_invs[8,:,:],svmin=[0,0,0,-3.],svmax=[1200,180,180,+3.],title=['Field strengh [Gauss]','Field inclination [degree]','Field azimuth [degree]','LoS velocity [km/s]'],xlabel='Pixel',ylabel='Pixel',save=outfile+'BLoS.png')
        np.savez_compressed(outfile+'_RTE', rte_invs=rte_invs, rte_invs_noth=rte_invs_noth)
        
        #del_dummy = subprocess.call("rm dummy_out.txt",shell=True)
        #print(del_dummy)
        printc('--------------------- END  ----------------------------',color=bcolors.FAIL)

        b_los_cropped = rte_invs_noth[2,:,:]*np.cos(rte_invs_noth[3,:,:]*np.pi/180.)*mask
        b_los = np.zeros((2048,2048))
        b_los[PXBEG2:PXEND2+1,PXBEG1:PXEND1+1] = b_los_cropped

        v_los_cropped = rte_invs_noth[8,:,:] * mask
        v_los = np.zeros((2048,2048))
        v_los[PXBEG2:PXEND2+1,PXBEG1:PXEND1+1] = v_los_cropped
        if verbose:
            plib.show_one(v_los,vmin=-2.5,vmax=2.5,title='LoS velocity')
            plib.show_one(b_los,vmin=-30,vmax=30,title='LoS magnetic field')

        hdu = pyfits.PrimaryHDU(b_los)
        hdul = pyfits.HDUList([hdu])
        hdul.writeto(outfile+'_blos.fits', overwrite=True)

        hdu = pyfits.PrimaryHDU(v_los)
        hdul = pyfits.HDUList([hdu])
        hdul.writeto(outfile+'_vlos.fits', overwrite=True)

        plib.show_one(rte_invs_noth[3,:,:]*mask,vmin=10,vmax=170)

    return


#%%
  # TODO Level for correct inversions is (*noise)

    #-----------------
    # APPLY COG
    #-----------------
    dr = np.copy(data)
    # co = dr[0,3,:,:]
    # co = ndimage.gaussian_filter(co, sigma=[15,15], order=0)
    # dr[:,3,:,:] = dr[:,3,:,:] - co[np.newaxis,:,:]
    print('Calculating COG.....')
    # vl,bl = COG(dr[:,:,ry[0]:ry[1],rx[0]:rx[1]],wave_axis,glande=3,cpos = 0)
    vl,bl = COG(dr,wave_axis,glande=3,cpos = 0)

    #-----------------
    # CORRECT CAVITY
    #-----------------

    try:
        if cavity == 0:
            cavity=np.zeros_like(vl)
            print("zerooooooooooo")
    except:
        # cavity = cavity[ry[0]:ry[1],rx[0]:rx[1]]
        pass

    factor = 0.5
    rrx = [int(c[1]-r*factor),int(c[1]+r*factor)]
    rry = [int(c[0]-r*factor),int(c[0]+r*factor)]
    print(rrx,rry,' check these for los vel calib')
    off = np.mean(vl[rry[0]:rry[1],rrx[0]:rrx[1]])
    vl = vl - off #- cavity
    print('velocity offset ',off)

    #-----------------
    # MASK DATA AND SAVE
    #-----------------

    vl = vl * mask#[ry[0]:ry[1],rx[0]:rx[1]]
    bl = bl * mask#[ry[0]:ry[1],rx[0]:rx[1]]
    # plib.show_one(vl[rry[0]:rry[1],rrx[0]:rrx[1]],vmin=-1.5,vmax=1.5)
    if verbose == 1:
        plib.show_one(vl,vmin=-2.5,vmax=2.5)
        plib.show_one(bl,vmin=-150,vmax=150)
    # plib.show_four_row(vl*mask[ry[0]:ry[1],rx[0]:rx[1]],vl*masks[ry[0]:ry[1],rx[0]:rx[1]],vl*maskp[ry[0]:ry[1],rx[0]:rx[1]],vl*maskq[ry[0]:ry[1],rx[0]:rx[1]],svmin=[-2,-2,-2,-2],svmax=[2,2,2,2],cagonento = 1)
    # quit()
    hdu = pyfits.PrimaryHDU(vl)
    hdul = pyfits.HDUList([hdu])
    hdul.writeto('data-red-tests/vlos_'+outfile, overwrite=True)
    hdu = pyfits.PrimaryHDU(bl)
    hdul = pyfits.HDUList([hdu])
    hdul.writeto('data-red-tests/blos_'+outfile, overwrite=True)
    
    #-----------------
    # PLOTS VLOS
    #-----------------
    Zm = np.ma.masked_where(mask == 1, mask)

    plt.figure(figsize=(10, 10))
    ax = plt.gca()
    plt.title('PHI-FDT LoS velocity',size=20)

    # Hide grid lines
    ax.grid(False)
    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])

    im = ax.imshow(np.fliplr(ndimage.rotate(vl, 50, reshape=False)), cmap='bwr',vmin=-3.,vmax=3.)

    divider = plib.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plib.plt.colorbar(im, cax=cax)
    cbar.set_label('[km/s]')
    cbar.ax.tick_params(labelsize=16)

    # ax.imshow(Zm[rx[0]:rx[1],ry[0]:ry[1]], cmap='gray')
    plt.savefig('imagenes/velocity-map-'+outfile+'.png',dpi=300)
    plt.close()

    #-----------------
    # PLOTS
    #-----------------

    esto = (dr[1,3,ry[0]:ry[1],rx[0]:rx[1]]-dr[3,3,ry[0]:ry[1],rx[0]:rx[1]])/2
    # im = ax.imshow(np.fliplr(ndimage.rotate(lacosa, 45, reshape=False)*100), cmap='gray',vmin=-0.2,vmax=0.2) 
    # plib.show_one(esto)
    # plib.show_one(dr[0,3,ry[0]:ry[1],rx[0]:rx[1]])
    # plib.show_one(bl - dr[0,3,ry[0]:ry[1],rx[0]:rx[1]])
    
    plt.figure(figsize=(10, 10))
    ax = plt.gca()
    plt.title('PHI-FDT Magnetogram',size=20)

    # Hide grid lines
    ax.grid(False)
    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
                
    lacosa = bl
    im = ax.imshow(np.fliplr(ndimage.rotate(lacosa, 45, reshape=False)), cmap='gray',vmin=-100,vmax=100) 

    divider = plib.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plib.plt.colorbar(im, cax=cax)
    # cbar.set_label('Stokes V amplitude [%]')
    cbar.set_label('LoS magnetic field [Mx/cm$^2$]')
    cbar.ax.tick_params(labelsize=16)

    #ax.imshow(Zm[rx[0]:rx[1],ry[0]:ry[1]], cmap='gray')
    plt.savefig('imagenes/magnetogram-map-'+outfile+'.png',dpi=300)
    plt.close()

    plt.figure(figsize=(10, 10))
    ax = plt.gca()
    plt.title('PHI-FDT Magnetogram',size=20)

    # Hide grid lines
    ax.grid(False)
    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])

    im = ax.imshow(np.fliplr(ndimage.rotate(lacosa, 45, reshape=False)), cmap='BrBG',vmin=-150,vmax=150) 

    divider = plib.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plib.plt.colorbar(im, cax=cax)
    cbar.set_label('LoS magnetic field [Mx/cm$^2$]')
    cbar.ax.tick_params(labelsize=16)

    #ax.imshow(Zm[rx[0]:rx[1],ry[0]:ry[1]], cmap='gray')
    plt.savefig('imagenes/magnetogram-map2-'+outfile+'.png',dpi=300)
    plt.close()

    # from datetime import datetime, timedelta
    # when = [datetime(2020, 5, 21,15,0,0),datetime(2020, 6, 18,14,00,0)]

    return vl