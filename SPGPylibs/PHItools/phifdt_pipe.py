import numpy as np 
import os.path, datetime, subprocess
#from astropy.io import fits as pyfits
#from time import sleep
#from scipy.ndimage import gaussian_filter#, rotate
#from scipy.interpolate import interp1d
#from scipy.optimize import curve_fit
from .tools import *
from .phi_fits import *
from .phi_gen import *
from .phi_reg import *
#from .phi_utils import newton,azimutal_average,limb_darkening,genera_2d,find_string
from .phifdt_pipe_modules import phi_correct_dark,phi_correct_prefilter,phi_apply_demodulation,\
    crosstalk_ItoQUV,cross_talk_QUV,crosstalk_ItoQUV2d,phi_correct_ghost,phi_correct_fringes

import SPGPylibs.GENtools.plot_lib as plib
import SPGPylibs.GENtools.cog as cog

def phifdt_pipe(data_f,dark_f,flat_f,instrument = 'FDT40',flat_c = True,dark_c = True,
    inner_radius = 250, outer_radius = 800, steps = 100, normalize_flat = 0., flat_n = 1.,
    index = None, prefilter = True, prefilter_fits = '0000990710_noMeta.fits',
    realign = False, verbose = True, outfile=None, mask_margin = 2, correct_fringes = False,
    individualwavelengths = False,correct_ghost = False,putmediantozero=True,directory = './',
    rte = False, debug = False,nlevel = 0.3,center_method='circlefit',loopthis=0,
    cross_talk_IQUV = False, cross_talk_VQU = False, do2d = 0):

    '''
	FDT pipeline steps:

    1- Read data
    2- Check dimensions (open to bugs since only few data were tested)
    3- Read flats
    4- Read and correct dark field (taking into account the scaling)
    5- Find center of the Sun in the data for masking and ghost correction
    6- get wavelength sampling from header
    7- move the continuum to the blue (if needed) in flat and data
    8- interpolate flats (Not implemented yet - in case of deviations in voltage)
    9- Correct flat-field
    9- Correct prefilter - needs prefilter data file!
    9- Correct ghost (complex step) [NEEDED in FM - high priority]
    10- realign data before demodulation [NEEDED in FM - low priority]
    11- Demodulate data using appropriate dem matrix
    12- Normalize to average continuum [NEEDED in FM - low priority - determine the Icont automatically]
    13- correct cross-talk from I to QUV [NEEDED in FM - evaluation of these automatically - mid priority]
    14- correct cross-talk from V to QU (interactive) [NEEDED in FM - evaluation of these automatically - mid priority]
    15- correct cross-talk from I to QUV in 2D (do not use this, your PC will have a hangover)
    16- Fringes correction [NEEDED in FM -TBD]
    17- median to zero [NEEDED in FM]
    18- save
    19- RTE (RTE or CE or CE+RTE)
    20- plots


    Parameters
    ----------
        Input:
    data_f : string
        Fits file of the raw FDT data  (if differe t path use directory keyword)
    dark_f : string
        Fits file of a Valid dark file (processed dark) (including path, if necessary)
    flat_f : string
        Fits file of a Valid FDT flatfield (including path, if necessary)
    
    IMPORTANT: dark_f, flat_f, and prefilter file must be provided with the FULL PATH. 
               the data has to be provided as a list of files (fits) and a directory: "directory="
        The output directories (Depending on RTE on or off) are
        A)  directory + Level 2: reduced raw data L2+ilam plus RTE output (so far, BLOS, VLOS and SO1: continuum) 
        B)  directory + Level 2 + png: png figures of the RTE output (so far, BLOS and VLOS) 
        B)  directory + Level 2 + npz: NPZ (python) reduced raw data L2
        
    ** OPTIONAL ARGUMENTS **

    instrument = 'FDT40' : select the instrument and PMP temperature (for demod)
        -> implemented cases: -- 'FDT40','FDT45' --
    flat_c = True : default is to apply flat field correction to the data
    dark_c = True : default is to apply dark field correction to the data
    inner_radius = 250, outer_radius = 600, steps = 100 : initial values for finding sun center
    normalize_flat = 0 : To normalize flats internally to the mean value of 5% of the disk (central) intensity  
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
    vervose: True prints a lot of stuff (and plots)
    mask_margin = 2: 'Number of pixels to contract the sun mask for output of RTE'
    correct_fringes = False: Fringe correction
        'manual': first FM version. Freq, mask applied to all images with fixed frequencies
        'auto' : calculate fringes freq. automatically (in development).
    correct_ghost = False; Correcto ghost images
    putmediantozero=True; puts median value to zero before RTE
    rte = False: Run RTE     if rte == 'RTE' or rte == 'CE' or rte == 'CE+RTE':
        'RTE': RTE only
        'CE+RTE': RTE with classical estiamtes
        'CE': Only classical estimates
        'cog': Only Center of gravity (in testing)
    cross_talk_IQUV= False: apply crostalk correction from Stokes I to Stokes Q, U, and V.
    cross_talk_VQU= False: apply crostalk correction from Stokes V to Stokes Q and U.
    nlevel = 0.3: Noise level above which to evaluate cross_talk_VQU (In testing)
    center_method='circlefit' or 'hough'

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
    This program is not optimized for speed. It assumes that input data is 6 wavelength.
    C-MILOS must be compiled in each specific machine (C)
    The software update some of the information in the keyword:

    Keywords in the header (modified or added) within this program:
 
    CAL_DARK =             26181001 / Onboard calibrated for dark field ! Dark correction ( DID/file of dark if True) 
    CAL_FLAT =             26181101 / Onboard calibrated for gain table ! Dark correction ( DID/file of flat if True)             
    CAL_PRE  =             Prefilter / Prefilter correction ( DID/file of flat if True)                    
    CAL_GHST=              Prefilter / Ghost correction ( name+version of py module if True )           
    CAL_REAL=              Prefilter / Prealigment of images before demodulation ( name+version of py module if True )             
    CAL_IPOL=               990510 / Onboard calibrated for instrumental polarizatio ! demodulation ( DID of demod matrix if True ) - demod matrix may be 4x4 or 2048x2048x4x4
    CAL_CRT0=               float / cross-talk from I to Q (slope value, wrt normalized data in python) 
    CAL_CRT1=               float / cross-talk from I to Q (off-set value, wrt normalized data in python) 
    CAL_CRT2=               float / cross-talk from I to U (slope value, wrt normalized data in python) 
    CAL_CRT3=               float / cross-talk from I to U (off-set value, wrt normalized data in python) 
    CAL_CRT4=               float / cross-talk from I to V (slope value, wrt normalized data in python) 
    CAL_CRT5=               float / cross-talk from I to V (off-set value, wrt normalized data in python) 
    CAL_CRT6=               float / cross-talk from V to Q (slope value, wrt normalized data in python) 
    CAL_CRT7=               float / cross-talk from V to Q (off-set value, wrt normalized data in python) 
    CAL_CRT8=               float / cross-talk from V to U (slope value, wrt normalized data in python) 
    CAL_CRT9=               float / cross-talk from V to U (off-set value, wrt normalized data in python) 
    CAL_NORM=               990510 / Normalization constant PROC_Ic)
    CAL_FRIN=               990510 / Fringe correction ( name+version of py module if True )  TBD (posibly we need the freqs.)  
   * CAL_PSF=                990510 / Onboard calibrated for instrumental PSF  ! TBD
    CAL_RTE=                990510 / ok
   * CAL_SCIP= 'None'               / Onboard scientific data analysis
   * RTE_ITER=           4294967295 / Number RTE inversion iterations

    (*) are not touched in this software.

    Keywords CRPIX1 and CRPIX2 are updated following the new center calculation within the pipeline. Old values are stored in the history.
    Keywords CRVAL1 and CRVAL2 are NOT updated but should be SET to zero!!!!

    '''

    version = 'V1.0 July 2021'
    version_cmilos = 'CMILOS v0.91 (July - 2021)'

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
    data_filename = directory + data_f
    print(data_filename)
    if os.path.isfile(data_filename):
        print("File exist")
    else:
        print("File not exist")

    try:
        data, header = fits_get(data_filename)
        DID = header['PHIDATID']
        printc('-->>>>>>> data DID '+DID,color=bcolors.OKGREEN)
        printc('-->>>>>>> Reshaping data to [wave,Stokes,y-dim,x-dim] ',color=bcolors.OKGREEN)
        zd,yd,xd = data.shape
        data = np.reshape(data,(zd//4,4,yd, xd))
        data = data / 256. #from fix to 32
        data = np.ascontiguousarray(data)

    except Exception:
        printc("ERROR, Unable to open fits file: {}",data_filename,color=bcolors.FAIL)
        return

    header['history'] = ' Data processed with phifdt_pipe.py '+ version
    header['history'] = '      and time '+ str(datetime.datetime.now())
    header['history'] = ' Parameters normalize_flat: '+ str(normalize_flat)
    header['history'] = ' Parameters flat_scaling: '+ str(flat_n)
    header['history'] = ' Parameters mask_margin: '+ str(mask_margin)
    header['history'] = ' Parameters center_method: '+ str(center_method)

    if verbose:
        plib.show_one(data[0,0,:,:],vmin=0,xlabel='pixel',ylabel='pixel',title='Data first wave',cbarlabel='DN',save=None,cmap='gray')
    
#    * CAL_RTE=                990510 / ok
#    * CAL_SCIP= 'None'               / Onboard scientific data analysis
#    * RTE_ITER=           4294967295 / Number RTE inversion iterations
#    * PHIDATID= '142010402'          / PHI dataset Id

    #-----------------
    # TAKE DATA DIMENSIONS 
    #-----------------
    PXBEG1  = int(header['PXBEG1']) - 1           
    PXEND1  = int(header['PXEND1']) - 1          
    PXBEG2  = int(header['PXBEG2']) - 1           
    PXEND2  = int(header['PXEND2']) - 1   
    printc('Dimensions: ',PXBEG1, PXEND1, PXBEG2, PXEND2,color=bcolors.OKGREEN)

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
            dummy,flat_header = fits_get(flat_f)
            flat = np.zeros([24,2048,2048]).astype(np.float32)
            PXBEG1_f  = int(flat_header['PXBEG1']) - 1           
            PXEND1_f  = int(flat_header['PXEND1']) - 1          
            PXBEG2_f  = int(flat_header['PXBEG2']) - 1           
            PXEND2_f  = int(flat_header['PXEND2']) - 1 
            flat[:,PXBEG1_f:PXEND1_f+1,PXBEG2_f:PXEND2_f+1] = dummy
            del dummy
            fz,fy,fx = flat.shape
            flat = np.reshape(flat,(fz//4,4,fy,fx))
            printc('-->>>>>>> Reshaping Flat to [wave,Stokes,y-dim,x-dim] ',color=bcolors.OKGREEN)
        except Exception:
            printc("ERROR, Unable to open flats file: {}",flat_f,color=bcolors.FAIL)
        if verbose:
            plib.show_one(flat[0,0,:,:],xlabel='pixel',ylabel='pixel',title='Flat first wave',cbarlabel='DN',save=None,cmap='gray')

    else:
        printc('-->>>>>>> No flats mode                    ',color=bcolors.WARNING)

    #-----------------
    # READ AND CORRECT DARK FIELD
    #-----------------
    if dark_c:
        data,header  = phi_correct_dark(dark_f,data_filename,header,data,verbose = verbose)

        # printc('-->>>>>>> Reading Darks                   ',color=bcolors.OKGREEN)
        # printc('          Input should be [y-dim,x-dim].',color=bcolors.OKGREEN)
        # printc('          DARK IS DIVIDED by 256.   ',color=bcolors.OKGREEN)
        # try:
        #     dark,h = fits_get(dark_f)
        #     dark = dark / 256.
        # except Exception:
        #     printc("ERROR, Unable to open darks file: {}",dark_f,color=bcolors.FAIL)

        # dy,dx = dark.shape
        # dark_scale = fits_get(dark_f,scaling = True)
        # data_scale = fits_get(data_f,scaling = True)
        # if dark_scale["Present"][0] == data_scale["Present"][0]:
        #     scaling = dark_scale["scaling"][0] / data_scale["scaling"][0]
        # else:
        #     scaling = dark_scale["scaling"][1] / data_scale["scaling"][1] * dark_scale["scaling"][0]

        # if scaling != 1:
        #     printc('          checking scalling and correcting for it in the dark.',dark_scale,data_scale,scaling,color=bcolors.WARNING)
        #     dark = dark * scaling
        # printc('-->>>>>>> Correcting dark current.',color=bcolors.OKGREEN)
        # data = data - dark[np.newaxis,np.newaxis,PXBEG2:PXEND2+1,PXBEG1:PXEND1+1]
        # data = np.abs(data)

        # if verbose:
        #     plib.show_one(dark,vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Dark',cbarlabel='DN',save=None,cmap='gray')
        #     plib.show_one(data[0,0,:,:],vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Data after dark',cbarlabel='DN',save=None,cmap='gray')
    else:
        printc('-->>>>>>> No darks mode                    ',color=bcolors.WARNING)

    #-----------------
    # FIND DATA CENTER 
    #-----------------

    printc('-->>>>>>> finding the center of the solar disk (needed for masking) ',color=bcolors.OKGREEN)
    try:
        if center_method == 'Hough':
            c, radius,threshold = find_circle_hough(data[0,0,:,:],inner_radius,outer_radius,steps,threshold = 0.01,normalize=False,verbose=False)
            #c = np.roll(c,1)
            cx = c[0]
            cy = c[1]
            #TBE PUT IN CORRECT UNITS
        elif center_method == 'circlefit':
            cy,cx,radius=find_center(data[0,0,:,:])  #OJO Cy... Cx
            c = np.array([int(cx),int(cy)])   #El vector es [0,1,2,...] == [x,y,z,...] == [cx,cy,cz,...] Pero esto ultimo esta al reves
            radius = int(radius)
        else:  
            raise ValueError("ERROR in center determination method") 
    except ValueError as err:
        print(err.args)

    #Uptade header with new centers

    printc('          Uptade header with new center:',color=bcolors.OKBLUE)
    printc('          OLD center:',color=bcolors.OKBLUE)
    printc('                  at: CRPIX1[x]=',header['CRPIX1'],' CRPIX2[y]=',header['CRPIX2'],' radius=',radius,color=bcolors.OKBLUE)
    header['history'] = ' CRPIX 1 and CRPIX2 uptated from ' + str(header['CRPIX1'])+ ' and ' + str(header['CRPIX2'])
    header['CRPIX1'] = (round(cx, 2))
    header['CRPIX2'] = (round(cy, 2))
    printc('          NEW center:',color=bcolors.OKBLUE)
    printc('                  at: CRPIX1[x]=',header['CRPIX1'],' CRPIX2[y]=',header['CRPIX2'],' radius=',radius,color=bcolors.OKBLUE)
    printc('ATTENTION: Keywords CRVAL1 and CRVAL2 are NOT updated but should be SET to zero!!!!',color=bcolors.FAIL)

    #OJO.
    # find_circle_hough devuelve c = c[0] = x and c[1] = y !!!!!!!!!!!!!!
    # Esto viene porque en el KLL esta definido así (al reves) en la rutina votes()
 
    #-----------------
    # TAKE ONLY DISK WITH MARGIN
    #-----------------
    printc('-->>>>>>> Creating a mask for RTE with ',mask_margin,' px margin')
    size_of_mask = radius - mask_margin
    rx = [int(c[0]-size_of_mask),int(c[0]+size_of_mask)]
    ry = [int(c[1]-size_of_mask),int(c[1]+size_of_mask)]
    mask,coords = generate_circular_mask([xd-1,yd-1],size_of_mask,size_of_mask)
    mask = shift(mask, shift=(c[0]-xd//2,c[1]-yd//2), fill_value=0)
    printc('   RX = ', rx, 'RY = ', ry, color=bcolors.WARNING)

    #-----------------
    # GET INFO ABOUT VOLTAGES/WAVELENGTHS, determine continuum and new flat
    #-----------------
    printc('-->>>>>>> Obtaining voltages from data ',color=bcolors.OKGREEN)
    wave_axis,voltagesData,tunning_constant,cpos,ref_wavelength = fits_get_sampling(data_filename)  
    printc('          Data FG voltages: ',voltagesData,color=bcolors.OKBLUE)
    printc('          Continuum position at wave: ', cpos,color=bcolors.OKBLUE)
    printc('          Data ref_wavelength [mA]: ',ref_wavelength,color=bcolors.OKBLUE)
    printc('          Data wave axis [mA]: ',wave_axis,color=bcolors.OKBLUE)
    printc('          Data wave axis - axis[0] [mA]: ',wave_axis - wave_axis[0],color=bcolors.OKBLUE)
    dummy_1 = (voltagesData-np.roll(voltagesData,-1))*(tunning_constant*1000)
    dummy_2 = np.sort(np.abs(dummy_1))
    sampling = np.mean(dummy_2[0:-2])
    printc('          Data average sampling [mA]: ',sampling,' using tunning constant: ',(tunning_constant*1000),color=bcolors.OKBLUE)
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
            printc('            * the continuum was TAKEN in the BLUE, stored in the RED, but we want it in the BLUE *',color=bcolors.WARNING)
            for i in range(zd//4):
                #print((i+1)%(zd//4),i%(zd//4),i,(zd//4))
                datar[(i+1)%(zd//4),:,:,:] = data[i%(zd//4),:,:,:] # np.roll(data, 4, axis=0)
                voltagesDatar[(i+1)%(zd//4)] = voltagesData[i%(zd//4)] # np.roll(data, 4, axis=0)
                wave_axisr[(i+1)%(zd//4)] = wave_axis[i%(zd//4)] # np.roll(data, 4, axis=0)
        if voltagesData[cpos] > voltagesData[0]:
        #continuum is on the right but it is in the blue!!!!
            printc('          Rolling data to move continuum from right (red) to left (blue)',color=bcolors.WARNING)
            printc('            * the continuum was TAKEN in the RED but we want it in the BLUE',color=bcolors.WARNING)
            printc('            * Notice that this is necessary for the FLAT correction',color=bcolors.FAIL)
            printc('            *  but the wavelength axis has changed the continuum point *',color=bcolors.FAIL)
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

        wave_axis_f,voltagesFlat,tunning_constant_f,cpos_f,ref_wavelength_f = fits_get_sampling(flat_f) 
        printc('          FLAT FG voltages: ',voltagesFlat,color=bcolors.OKBLUE)
        printc('          FLAT Continuum position at wave: ', cpos_f,color=bcolors.OKBLUE)
        printc('          FLAT wave axis [mA]: ',wave_axis_f,color=bcolors.OKBLUE)
        printc('          FLAT ref_wavelength [mA]: ',ref_wavelength_f,color=bcolors.OKBLUE)
        printc('          FLAT wave axis [mA]: ',wave_axis_f,color=bcolors.OKBLUE)
        printc('          FLAT wave axis - ref_wavelength [mA]: ',wave_axis_f - ref_wavelength_f,color=bcolors.OKBLUE)
        dummy_1 = (voltagesFlat-np.roll(voltagesFlat,-1))*(tunning_constant*1000)
        dummy_2 = np.sort(np.abs(dummy_1))
        sampling_f = np.mean(dummy_2[0:-2])
        printc('          FLAT average sampling [mA]: ',sampling_f,color=bcolors.OKBLUE)

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
                printc('            * the continuum was TAKEN in the BLUE, stored in the RED, but we want it in the BLUE *',color=bcolors.WARNING)
                for i in range(fz//4):
                    #print((i+1)%(fz//4),i%(fz//4),i)
                    flatr[(i+1)%(fz//4),:,:,:] = flat[i%(fz//4),:,:,:] # np.roll(data, 4, axis=0)
                    voltagesFlatr[(i+1)%6] = voltagesFlat[i%(fz//4)] # np.roll(data, 4, axis=0)
                    wave_axis_fr[(i+1)%(zd//4)] = wave_axis_f[i%(zd//4)] # np.roll(data, 4, axis=0)
            if voltagesFlat[cpos_f] > voltagesFlat[0]:
            #continuum is on the right but it is in the blue!!!!
                printc('          Rolling flat to move continuum from right (red) to left (blue)',color=bcolors.WARNING)
                printc('            * the continuum was TAKEN in the RED but we want it in the BLUE',color=bcolors.WARNING)
                printc('            * Notice that this is necessary for the FLAT correction',color=bcolors.FAIL)
                printc('            *  but the wavelength axis has changed the continuum point *',color=bcolors.FAIL)
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
    rrx = [int(c[0]-radius*factor),int(c[0]+radius*factor)]
    rry = [int(c[1]-radius*factor),int(c[1]+radius*factor)]

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
                if normalize_flat ==1:
                    print('          normalizing flats using region x = [',rrx[0],':',rrx[1],'] y = ]',rry[0],':',rry[1],']')
                    mm = np.mean(dummy_flat[rry[0]:rry[1],rrx[0]:rrx[1]])
                    dummy_flat = dummy_flat / mm
                data[l,p,:,:] = data[l,p,:,:]/dummy_flat
        del dummy_flat

        # locations = find_string(flat_f,'_')
        # try:
        #     DID_flat = flat_f[locations[-1]+1:locations[-1]+10]
        #     print('DID: ',np.float(DID_flat))
        # except:
        #     DID_flat = flat_f[:-4]
        #     printc("Unable to get DID from: {}",flat_f,color=bcolors.WARNING)         
        #     locations = find_string(dark_f,'/')
        #     DID_flat = flat_f[locations[-1]+1:]
        #     printc('DID: ',DID_flat,' -->> WILL NOT BE A NUMBER',color=bcolors.WARNING)

        DID_flat = flat_header['PHIDATID']

        if 'CAL_FLAT' in header:  # Check for existence
            header['CAL_FLAT'] = DID_flat
        else:
            header.set('CAL_FLAT', DID_flat, 'Onboard calibrated for gain table',after='CAL_DARK')

        if verbose:
            plib.show_one(data[cpos,0,:,:],vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Data / flat at continuum',cbarlabel='DN',save=None,cmap='gray')

    #-----------------
    # CORRECT PREFILTER 
    #-----------------
    if prefilter == 1:
        data,header  = phi_correct_prefilter(prefilter_fits,header,data,voltagesData,verbose = verbose)


    #-----------------
    # GHOST CORRECTION  
    #-----------------

    if correct_ghost:
        data,header = phi_correct_ghost(data,header,radius,verbose = verbose)

    #-----------------
    # REALIGN DATA BEFORE DEMODULATION
    #-----------------
    if realign:
        printc('-->>>>>>> Realigning data...           ',color=bcolors.OKGREEN)
        for i in range(zd//4):
            s_x,s_y,_ = PHI_shifts_FFT(data[i,:,:,:],prec=500,verbose=verbose,norma=False)
            for j in range(4):
                data[i,j,:,:] = shift_subp(data[i,j,:,:], shift=[s_x[j],s_y[j]]) #estra y,z asi que esta al reves FFT
        if 'CAL_REAL' in header:  # Check for existence
            header['CAL_REAL'] = 'FFT'
        else:
            header.set('CAL_REAL', 'FFT', 'Realigment of data (phifdt_pipe_modules.py)',after='CAL_DARK')

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

    data = phi_apply_demodulation(data,header,instrument)

    if verbose == 1:
        plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'])

    #-----------------
    # APPLY NORMALIZATION 
    #-----------------
    printc('-->>>>>>> Applying normalization --',color=bcolors.OKGREEN)
    nrm = np.mean(data[cpos,0,rry[0]:rry[1],rrx[0]:rrx[1]])
    print('          Norma is: ',nrm,' evaluated in x = [',rrx[0],':',rrx[1],'] y = [',rry[0],':',rry[1],']')
    data = data/nrm

    if 'CAL_NORM' in header:  # Check for existence
        header['CAL_NORM'] = np.round(nrm,6)
    else:
        header.set('CAL_NORM', np.round(nrm,6), 'Normalization constant PROC_Ic',after='CAL_DARK')

    if debug:
        datan = datan/nrm

    if debug:
        plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:])
        plib.show_four_row(datan[3,0,:,:],datan[3,1,:,:],datan[3,2,:,:],datan[3,3,:,:])

    if verbose == 1:
        plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'])

    #-----------------
    # CROSS-TALK CALCULATION 
    #-----------------
    if cross_talk_IQUV:
        printc('-->>>>>>> Cross-talk correction from Stokes I to Stokes Q,U,V --',color=bcolors.OKGREEN)
        factor = 0.80 # 80% of the disk
        printc('          Using ',factor*100,'% of the disk                     ',color=bcolors.OKGREEN)

        rrx = [int(c[0]-radius*factor),int(c[0]+radius*factor)]
        rry = [int(c[1]-radius*factor),int(c[1]+radius*factor)]
        printc('          Crosstalk evaluated in x = [',rrx[0],':',rrx[1],'] y = [',rry[0],':',rry[1],']',' using ',factor*100,"% of the disk",color=bcolors.OKBLUE)

        if individualwavelengths:
            for i in range(zd//4):
                printc('          Individual wavelengths....',color=bcolors.OKBLUE)
                broadcastd = data[i,:,rry[0]:rry[1],rrx[0]:rrx[1]]
                data_dummy = data[:,:,rry[0]:rry[1],rrx[0]:rrx[1]]*0. + broadcastd[np.newaxis,:,:,:]
                cQ,cU,cV = crosstalk_ItoQUV(data_dummy[:,:,rry[0]:rry[1],rrx[0]:rrx[1]],npoints=10000,verbose=verbose)
                #-----------------
                # CROSS-TALK CORRECTION 
                #-----------------
                printc('          Applying cross-talk correction...',color=bcolors.OKGREEN)
                data[i,1,:,:] = data[i,1,:,:] - cQ[0]*data[i,0,:,:] - cQ[1]
                data[i,2,:,:] = data[i,2,:,:] - cU[0]*data[i,0,:,:] - cU[1]
                data[i,3,:,:] = data[i,3,:,:] - cV[0]*data[i,0,:,:] - cV[1]
            if verbose:
                plt.hist(data[4,1,900:1100,900:1100].flatten(), bins='auto')
                plt.hist(data[4,2,900:1100,900:1100].flatten(), bins='auto')
                plt.hist(data[4,3,900:1100,900:1100].flatten(), bins='auto')
                plt.show()
        else:
            cQ,cU,cV = crosstalk_ItoQUV(data[:,:,rry[0]:rry[1],rrx[0]:rrx[1]],verbose=verbose,npoints=10000)
        # rrx = [int(c[1]),int(c[1]+r*factor)]
        # rry = [int(c[0]-r*factor),int(c[0])]
        # cQ,cU,cV = crosstalk_ItoQUV(data[:,:,rry[0]:rry[1],rrx[0]:rrx[1]],verbose=verbose,npoints=10000)
        # rrx = [int(c[1]-r*factor),int(c[1])]
        # rry = [int(c[0]),int(c[0]+r*factor)]
        # cQ,cU,cV = crosstalk_ItoQUV(data[:,:,rry[0]:rry[1],rrx[0]:rrx[1]],verbose=verbose,npoints=10000)
        # return
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
        # PLT_RNG = 2
        # plib.show_four_row(data[1,0,:,:],data[1,1,:,:],data[1,2,:,:],data[1,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003])#,save='t1_'+str(loopthis)+'.png')
        # plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003])#,save='t3_'+str(loopthis)+'.png')
        # plib.show_four_row(data[5,0,:,:],data[5,1,:,:],data[5,2,:,:],data[5,3,:,:],title=['I','Q','U','V'],svmin=[0.1,-0.002,-0.002,-0.003],svmax=[1.1,0.002,0.002,0.003])#,save='t5_'+str(loopthis)+'.png')

        # np.save('data_dummy',data)
        if 'CAL_CRT0' in header:  # Check for existence
            header['CAL_CRT0'] = np.round(cQ[0]*100,3)
        else:
            header.set('CAL_CRT0', np.round(cQ[0]*100,3), 'cross-talk I to Q (slope in %), wrt CAL_NROM ',after='CAL_DARK')

        if 'CAL_CRT1' in header:  # Check for existence
            header['CAL_CRT1'] = np.round(cQ[1]*100,3)
        else:
            header.set('CAL_CRT1', np.round(cQ[1]*100,3), 'cross-talk I to Q (off-set in %),  wrt CAL_NROM ',after='CAL_CRT0')

        if 'CAL_CRT2' in header:  # Check for existence
            header['CAL_CRT2'] = np.round(cU[0]*100,3)
        else:
            header.set('CAL_CRT2', np.round(cU[0]*100,3), 'cross-talk I to U (slope in %) alue,  wrt CAL_NROM ',after='CAL_CRT1')

        if 'CAL_CRT3' in header:  # Check for existence
            header['CAL_CRT3'] = np.round(cU[1]*100,3)
        else:
            header.set('CAL_CRT3', np.round(cU[1]*100,3), 'cross-talk I to U (off-set in %),  wrt CAL_NROM ',after='CAL_CRT2')

        if 'CAL_CRT4' in header:  # Check for existence
            header['CAL_CRT4'] = np.round(cV[0]*100,3)
        else:
            header.set('CAL_CRT4', np.round(cV[0]*100,3), 'cross-talk I to V (slope in %),  wrt CAL_NROM',after='CAL_CRT3')

        if 'CAL_CRT5' in header:  # Check for existence
            header['CAL_CRT5'] = np.round(cV[1]*100,3)
        else:
            header.set('CAL_CRT5', np.round(cV[1]*100,3), 'cross-talk I to V (off-set in %),  wrt CAL_NROM ',after='CAL_CRT4')

    #-----------------
    # CROSS-TALK CALCULATION FROM V TO QU (Interactive)
    #-----------------
    
    if cross_talk_VQU:
        printc('-->>>>>>> Cross-talk correction from Stokes V to Stokes Q,U ',color=bcolors.OKGREEN)

        factor = 0.3 # 30% of the disk
        rrx = [int(c[0]-radius*factor),int(c[0]+radius*factor)]
        rry = [int(c[1]-radius*factor),int(c[1]+radius*factor)]
        print(' Cross-talk evaluated in x = [',rrx[0],':',rrx[1],'] y = [',rry[0],':',rry[1],']',' using ',factor*100,"% of the disk")
        cVQ,cVU = cross_talk_QUV(data[:,:,rry[0]:rry[1],rrx[0]:rrx[1]],nran = 2000,nlevel=nlevel,block=False)

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

        if 'CAL_CRT6' in header:  # Check for existence
            header['CAL_CRT6'] = np.round(cVQ[0]*100,3)
        else:
            header.set('CAL_CRT6', np.round(cVQ[0]*100,3), 'cross-talk V to Q  (slope in %),  wrt CAL_NROM ',after='CAL_CRT5')

        if 'CAL_CRT7' in header:  # Check for existence
            header['CAL_CRT7'] = np.round(cVQ[1]*100,3)
        else:
            header.set('CAL_CRT7', np.round(cVQ[1]*100,3), 'cross-talk V to Q (off-set in %),  wrt CAL_NROM ',after='CAL_CRT6')

        if 'CAL_CRT8' in header:  # Check for existence
            header['CAL_CRT8'] = np.round(cVU[0]*100,3)
        else:
            header.set('CAL_CRT8', np.round(cVU[0]*100,3), 'cross-talk V to U (slope in %),  wrt CAL_NROM ',after='CAL_CRT7')

        if 'CAL_CRT9' in header:  # Check for existence
            header['CAL_CRT9'] = np.round(cVU[1]*100,3)
        else:
            header.set('CAL_CRT9', np.round(cVU[1]*100,3), 'cross-talk V to U (off-set in %),  wrt CAL_NROM ',after='CAL_CRT8')

    #-----------------
    # CROSS-TALK CALCULATION FROM I TO QUV (2D)
    #-----------------

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
    # FRINGING - 
    #-----------------

    data, header = phi_correct_fringes(data,header,option=correct_fringes,verbose=verbose)
    if verbose:
        plib.show_four_row(data[2,0,:,:],data[2,1,:,:],data[2,2,:,:],data[2,3,:,:],title=['I','Q','U','V'])

    #-----------------
    # MEDIAN TO CERO
    #-----------------

    if putmediantozero:
        factor = 0.8
        rrx = [int(c[0]-radius*factor),int(c[0]+radius*factor)]
        rry = [int(c[1]-radius*factor),int(c[1]+radius*factor)]
        printc('-->>>>>>> Putting median to zero ',color=bcolors.OKGREEN)
        printc('          Median evaluated in x = [',rrx[0],':',rrx[1],'] y = [',rry[0],':',rry[1],']',' using ',factor*100,"% of the disk",color=bcolors.OKBLUE)
        maski,coords = generate_circular_mask([xd-1,yd-1],radius*factor,radius*factor)
        maski = shift(maski, shift=(c[0]-xd//2,c[1]-yd//2), fill_value=0).astype(int)
        for i in range(zd//4):
            PQ = np.median( data[i,1, maski > 0])#,axis=(1,2))
            PU = np.median( data[i,2, maski > 0])#,axis=(1,2))
            PV = np.median( data[i,3, maski > 0])#,axis=(1,2))
        # PQ = np.median(data[:,1,rry[0]:rry[1],rrx[0]:rrx[1]],axis=(1,2))
        # PU = np.median(data[:,2,rry[0]:rry[1],rrx[0]:rrx[1]],axis=(1,2))
        # PV = np.median(data[:,3,rry[0]:rry[1],rrx[0]:rrx[1]],axis=(1,2))
            data[i,1,:,:] = data[i,1,:,:] - PQ#[:,np.newaxis,np.newaxis]
            data[i,2,:,:] = data[i,2,:,:] - PU#[:,np.newaxis,np.newaxis]
            data[i,3,:,:] = data[i,3,:,:] - PV#[:,np.newaxis,np.newaxis]
            printc(PQ,PU,PV)
    if verbose == 1:
        plt.hist(data[4,1,900:1100,900:1100].flatten(), bins='auto')
        plt.hist(data[4,2,900:1100,900:1100].flatten(), bins='auto')
        plt.hist(data[4,3,900:1100,900:1100].flatten(), bins='auto')
        plt.show()
        plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'])

        header['history'] = ' Parameters putmediantozero [%]: '+ str(np.round(PQ*100,6))+ ' '+ str(np.round(PU*100,6))+ ' '+ str(np.round(PV*100,6))

    #-----------------
    #CHECK FOR INFs
    #-----------------

    data[np.isinf(data)] = 0
    data[np.isnan(data)] = 0

    #-----------------
    # SAVE DATA TODO: CMILOS FORMAT AND FITS
    #-----------------

    #check if npz,pngs and level2 exist
    dirs = ['npz','pngs','level2']
        
    for checkit in dirs:
        check_dir = os.path.isdir(directory+checkit)
        if not check_dir:
            os.makedirs(directory+checkit)
            print("created folder : ", directory+checkit)
        else:
            print(directory+checkit, "folder already exists.")

    printc('---------------------------------------------------------',color=bcolors.OKGREEN)
    if outfile == None:
        #basically replace L1 by L1.5
        try:
            outfile_L2 = set_level(data_f,'L1','L2')
        except:
            outfile_L2 = set_level(data_f,'L0','L2')

    printc(' Saving data to: ',directory+'level2/'+outfile_L2)

    # hdu = pyfits.PrimaryHDU(data)
    # hdul = pyfits.HDUList([hdu])
    # hdul.writeto(outfile, overwrite=True)

    with pyfits.open(data_filename) as hdu_list:
        hdu_list[0].data = data
#        header = hdu_list[0].header
        hdu_list[0].header = header

        # Add a new key to the header
        # header.insert(20, ('NEWKEY', 'OMIT', 'test'))
        #header.set('NEWKEY','50.5')
        hdu_list.writeto(directory+'level2/'+outfile_L2, clobber=True)
#        hdu_list.writeto(directory+outfile+'_L1.fits', clobber=True)

    # with pyfits.open(data_f) as hdu_list:
    #     hdu_list[0].data = mask
    #     hdu_list.writeto(directory+outfile+'_red-mask.fits', clobber=True)

    #-----------------
    # INVERSION OF DATA WITH CMILOS
    #-----------------

    if rte == 'RTE' or rte == 'CE' or rte == 'CE+RTE':
        printc('---------------------RUNNING CMILOS --------------------------',color=bcolors.OKGREEN)

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

        wavelength = 6173.3354
        #OJO, REMOVE. NEED TO CHECK THE REF WAVE FROM S/C-PHI H/K
        shift_w =  wave_axis[3] - wavelength
        wave_axis = wave_axis - shift_w
        #wave_axis = np.array([-300,-160,-80,0,80,160])/1000.+wavelength
        #wave_axis = np.array([-300,-140,-70,0,70,140])
        printc('   It is assumed the wavelength is given by the header info ')
        printc('         wave axis: ', wave_axis,color = bcolors.WARNING)
        printc('         wave axis (step):  ',(wave_axis - wavelength)*1000.,color = bcolors.WARNING)
        printc('   saving data into dummy_in.txt for RTE input')

        sdata = data[:,:,ry[0]:ry[1],rx[0]:rx[1]]
        l,p,x,y = sdata.shape
        print(l,p,x,y)

        filename = 'dummy_in.txt'
        with open(filename,"w") as f:
            for i in range(x):
                for j in range(y):
                    for k in range(l):
                        f.write('%e %e %e %e %e \n' % (wave_axis[k],sdata[k,0,j,i],sdata[k,1,j,i],sdata[k,2,j,i],sdata[k,3,j,i]))
        del sdata

        printc('  ---- >>>>> Inverting data.... ',color=bcolors.OKGREEN)
        umbral = 3.

        cmd = CMILOS_LOC+"./milos"
        cmd = fix_path(cmd)
        if rte == 'RTE':
            rte_on = subprocess.call(cmd+" 6 15 0 0 dummy_in.txt  >  dummy_out.txt",shell=True)
        if rte == 'CE':
            rte_on = subprocess.call(cmd+" 6 15 2 0 dummy_in.txt  >  dummy_out.txt",shell=True)
        if rte == 'CE+RTE':
            rte_on = subprocess.call(cmd+" 6 15 1 0 dummy_in.txt  >  dummy_out.txt",shell=True)

        print(rte_on)
        printc('  ---- >>>>> Finishing.... ',color=bcolors.OKGREEN)
        printc('  ---- >>>>> Reading results.... ',color=bcolors.OKGREEN)
        del_dummy = subprocess.call("rm dummy_in.txt",shell=True)
        print(del_dummy)

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
            plib.show_four_row(rte_invs_noth[2,:,:],rte_invs_noth[3,:,:],rte_invs_noth[4,:,:],rte_invs_noth[8,:,:],svmin=[0,0,0,-6.],svmax=[1200,180,180,+6.],title=['Field strengh [Gauss]','Field inclination [degree]','Field azimuth [degree]','LoS velocity [km/s]'],xlabel='Pixel',ylabel='Pixel')#,save=outfile+'_VLoS.png')
            plib.show_four_row(rte_invs[2,:,:],rte_invs[3,:,:],rte_invs[4,:,:],rte_invs[8,:,:],svmin=[0,0,0,-6.],svmax=[1200,180,180,+6.],title=['Field strengh [Gauss]','Field inclination [degree]','Field azimuth [degree]','LoS velocity [km/s]'],xlabel='Pixel',ylabel='Pixel')#,save=outfile+'BLoS.png')
        rte_invs_noth[8,:,:] = rte_invs_noth[8,:,:] - np.mean(rte_invs_noth[8,rry[0]:rry[1],rrx[0]:rrx[1]])
        rte_invs[8,:,:] = rte_invs[8,:,:] - np.mean(rte_invs[8,rry[0]:rry[1],rrx[0]:rrx[1]])

        np.savez_compressed(directory+'npz/'+outfile_L2+'_RTE', rte_invs=rte_invs, rte_invs_noth=rte_invs_noth,mask=mask)

        del_dummy = subprocess.call("rm dummy_out.txt",shell=True)
        print(del_dummy)

        b_los = rte_invs_noth[2,:,:]*np.cos(rte_invs_noth[3,:,:]*np.pi/180.)*mask
        # b_los = np.zeros((2048,2048))
        # b_los[PXBEG2:PXEND2+1,PXBEG1:PXEND1+1] = b_los_cropped

        v_los = rte_invs_noth[8,:,:] * mask
        # v_los = np.zeros((2048,2048))
        # v_los[PXBEG2:PXEND2+1,PXBEG1:PXEND1+1] = v_los_cropped
        if verbose:
            plib.show_one(v_los,vmin=-2.5,vmax=2.5,title='LoS velocity')
            plib.show_one(b_los,vmin=-30,vmax=30,title='LoS magnetic field')

        header['history'] = ' RTE CMILOS INVERTER: '+ rte
        header['history'] = ' CMILOS VER: '+ version_cmilos
        
        if 'RTE_ITER' in header:  # Check for existence
            header['RTE_ITER'] = 'FFT'
        else:
            header.set('RTE_ITER', str(15), 'Number RTE inversion iterations',after='CAL_SCIP')

        with pyfits.open(directory+data_f) as hdu_list:
            hdu_list[0].data = b_los
#            header = hdu_list[0].header
            hdu_list[0].header = header
            writeto = set_level(outfile_L2,'ilam','blos')
            hdu_list.writeto(directory+'level2/'+writeto, clobber=True)

        with pyfits.open(directory+data_f) as hdu_list:
            hdu_list[0].data = v_los
#            header = hdu_list[0].header
            hdu_list[0].header = header
            writeto = set_level(outfile_L2,'ilam','vlos')
            hdu_list.writeto(directory+'level2/'+writeto, clobber=True)

        with pyfits.open(directory+data_f) as hdu_list:
            hdu_list[0].data = rte_invs[9,:,:]+rte_invs[10,:,:]
#            header = hdu_list[0].header
            hdu_list[0].header = header
            writeto = set_level(outfile_L2,'ilam','icont')
            hdu_list.writeto(directory+'level2/'+writeto, clobber=True)

        printc('  ---- >>>>> Saving plots.... ',color=bcolors.OKGREEN)

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

    #        im = ax.imshow(np.fliplr(rotate(v_los[PXBEG2:PXEND2+1,PXBEG1:PXEND1+1], 52, reshape=False)), cmap='bwr',vmin=-3.,vmax=3.)
        im = ax.imshow(v_los, cmap='bwr',vmin=-3.,vmax=3.)

        divider = plib.make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plib.plt.colorbar(im, cax=cax)
        cbar.set_label('[km/s]')
        cbar.ax.tick_params(labelsize=16)

        #ax.imshow(Zm, cmap='gray')
        writeto = set_level(outfile_L2,'ilam','vlos')
        writeto = set_level(writeto,'.fits','.png')
        plt.savefig(directory+'pngs/'+writeto,dpi=300)
        plt.close()

        #-----------------
        # PLOTS BLOS
        #-----------------

        plt.figure(figsize=(10, 10))
        ax = plt.gca()
        plt.title('PHI-FDT Magnetogram',size=20)

        # Hide grid lines
        ax.grid(False)
        # Hide axes ticks
        ax.set_xticks([])
        ax.set_yticks([])
                    
        #im = ax.imshow(np.fliplr(rotate(b_los[PXBEG2:PXEND2+1,PXBEG1:PXEND1+1], 52, reshape=False)), cmap='gray',vmin=-100,vmax=100) 
        im = ax.imshow(b_los, cmap='gray',vmin=-100,vmax=100)

        divider = plib.make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plib.plt.colorbar(im, cax=cax)
        # cbar.set_label('Stokes V amplitude [%]')
        cbar.set_label('LoS magnetic field [Mx/cm$^2$]')
        cbar.ax.tick_params(labelsize=16)
        #ax.imshow(Zm, cmap='gray')

        writeto = set_level(outfile_L2,'ilam','blos')
        writeto = set_level(writeto,'.fits','.png')
        plt.savefig(directory+'pngs/'+writeto,dpi=300)
        plt.close()

        printc('--------------------- END  ----------------------------',color=bcolors.FAIL)

    if rte == 'cog':
        printc('---------------------RUNNING COG --------------------------',color=bcolors.OKGREEN)
        wavelength = 6173.3356
        v_los,b_los = cog(data,wavelength,wave_axis,lande_factor=3,cpos = cpos)

        #-----------------
        # MASK DATA AND SAVE
        #-----------------

        v_los = v_los * mask
        b_los = b_los * mask
        plib.show_one(v_los,vmin=-1.5,vmax=1.5)
        plib.show_one(b_los,vmin=-150,vmax=150)
        if verbose == 1:
            plib.show_one(v_los,vmin=-2.5,vmax=2.5)
            plib.show_one(b_los,vmin=-150,vmax=150)

        with pyfits.open(data_f) as hdu_list:
            hdu_list[0].data = v_los
            hdu_list[0].header = header
            writeto = set_level(outfile_L2,'ilam','vlos-cog')
            writeto = set_level(writeto,'.fits','.png')
            plt.savefig(directory+writeto,dpi=300)

        with pyfits.open(data_f) as hdu_list:
            hdu_list[0].data = b_los
            hdu_list[0].header = header
            writeto = set_level(outfile_L2,'ilam','blos-cog')
            writeto = set_level(writeto,'.fits','.png')
            plt.savefig(directory+'pngs/'+writeto,dpi=300)

    return
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
    rrx = [int(c[0]-r*factor),int(c[0]+r*factor)]
    rry = [int(c[1]-r*factor),int(c[1]+r*factor)]
    print(rrx,rry,' check these for los vel calib')
    off = np.mean(vl[rry[0]:rry[1],rrx[0]:rrx[1]])
    vl = vl - off #- cavity
    print('velocity offset ',off)

    # cavity,h = phi.fits_read('HRT_cavity_map_IP5.fits')
    # cavity = cavity * 0.3513e-3/6173.*300000. #A/V 
    # cavity = cavity - np.median(cavity[800:1200,800:1200])
    
    return vl
