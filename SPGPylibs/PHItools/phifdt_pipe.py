import json
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
from .phi_rte import *
#from .phi_utils import newton,azimutal_average,limb_darkening,genera_2d,find_string
from .phifdt_pipe_modules import phi_correct_dark,phi_correct_prefilter,phi_apply_demodulation,\
    crosstalk_ItoQUV,cross_talk_QUV,crosstalk_ItoQUV2d,phi_correct_ghost,phi_correct_fringes,\
    generate_level2,check_pmp_temp

import SPGPylibs.GENtools.plot_lib as plib
import SPGPylibs.GENtools.cog as cog

from platform import node
FIGUREOUT = '.png'
SIX_FLATS = 1

#global variables 
PLT_RNG = 5

def phifdt_pipe(json_input = None, 
    data_f: str = None,  dark_f: str = None,  flat_f: str = None,
    input_data_dir: str = './',   output_dir:str = './',
    instrument: str = 'FDT40',
    flat_c:bool = True, dark_c:bool = True, ItoQUV:bool = False, VtoQU:bool = False, ind_wave:bool = False,        #correction options
    hough_params:list = [250, 800, 100], 
    norm_f:bool = False, flat_scaling:float = 1., flat_index:list = None,                                          #flatfield options
    prefilter:bool = True, prefilter_fits:str = '0000990710_noMeta.fits',
    realign:bool = False, verbose:bool = True, shrink_mask:int = 2, correct_fringes:str = False,
    correct_ghost:bool = False, putmediantozero:bool = True,
    debug:bool = False, nlevel:float = 0.3, center_method:str = 'circlefit',
    vers = '01',
    rte: str = False, RTE_code: str = 'cmilos',
    do2d = 0, outfile = None #not in use
    ) -> int: 
    
    '''
    Parameters
    ----------
    :param str data_f: 

    Input parameters
    ----------    
    json_input = json input (for convenience). All parameters are then described there).
    data_f = data_f : string
        Fits file of the raw FDT data  (for path use input_data_dir keyword)
    dark_f = dark_f : string
        Fits file of a Valid dark file (processed dark) (including path, if necessary)
    flat_f = flat_f : string
        Fits file of a Valid FDT flatfield (including path, if necessary)

    input_data_dir: directory where input data is located. Default is local directory
    output_dir: output directory. If default, takes local './'  

    IMPORTANT: dark_f, flat_f, and prefilter file must be provided with the FULL PATH. 
               the data has to be provided as a list of files (fits) and the directory via "input_data_dir = "
        The output directories (Depending on RTE on or off) are
        A)  directory + Level 2: reduced raw data L2+ilam plus RTE output (so far, BLOS, VLOS and SO1: continuum) 
        B)  directory + Level 2 + png: png figures of the RTE output (so far, BLOS and VLOS) 
        B)  directory + Level 2 + npz: NPZ (python) reduced raw data L2
        
    ** OPTIONAL ARGUMENTS **

    instrument = 'FDT40' : select the instrument and PMP temperature (for demod)
        -> implemented cases: -- 'FDT40','FDT45' --
    flat_c = True : default is to apply flat field correction to the data
    dark_c = True : default is to apply dark field correction to the data
    norm_f = False : To normalize flats internally to the mean value of 5% of the disk (central) intensity  
    flat_scaling = 1.0 : flat scaling (flat = flat / flat_scaling) 
    flat_index = None : in case you want a particular flat to be applied at another wave, e.g.,
        flat_index = [5,1,2,3,4,0] exchange the first and last wave flats
        This is for testing stuff, mainly. 
    prefilter = 1 : To correct for the prefilter 
    prefilter_fits = '../RSW1/0000990710_noMeta.fits' : User should provide prefilter data fits file location
    realign = False : bool
        Realign all images before demodulating using FFT 
    ind_wave = False : bool
        Correct crosstalk from I to QUV for individual wavelengths
    vervose: True prints a lot of stuff (and plots)
    shrink_mask = 2: 'Number of pixels to contract the sun mask for output of RTE'
    correct_fringes = False: Fringe correction
        'manual': first FM version. Freq, mask applied to all images with fixed frequencies
        'auto' : calculate fringes freq. automatically (in development).
    correct_ghost = False; Correcto ghost images
    putmediantozero=True; puts median value to zero before RTE
    rte = False: Run RTE     if rte == 'RTE' or rte == 'CE' or rte == 'CE+RTE':
        'RTE': RTE only
        'CE+RTE': RTE with classical estiamtes
        'CE': Only classical estimates
        'cog': Only Center of gravity (To be implemented)
    ItoQUV= False: apply crostalk correction from Stokes I to Stokes Q, U, and V.
    VtoQU= False: apply crostalk correction from Stokes V to Stokes Q and U.
    nlevel = 0.3: Noise level above which to evaluate cross_talk_VQU (To be implemented)

    center_method = ['circlefit','hough']
        Default is 'circlefit'. If set to 'hough' uses the given find_center parameters 
        If find_center is set to None then uses header information, in any case
    hough_params = [250, 800, 100]; inner_radius = 250, outer_radius = 600, steps = 100 : initial values for finding sun center

    verbose: increase the verbosity (many plots here) - default False

    Returns
    -------
    0 if fail, 1 any other case 

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
    The software update some of the information in the fits keyword:

    TODO:
    # data_f -> input data (single file for FDT - ADD MULTIPLE FILES!!!! )
    keyword to provide fixed cross-talk coefficients 
    keyword to provide fixed data normalization (based on a first obs) 
    # pending to add class stile (for module development)

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
    version = 'V1.0 13th September 2021'
    version = 'V1.0 3th November 2021'
    #added json configuration and modify all keyword names to be consistent with HRT pipe
    version = 'V1.1 13th April 2022'

    version_cmilos = 'CMILOS v0.91 (July - 2021)'

    printc('--------------------------------------------------------------',bcolors.OKGREEN)
    printc('PHI FDT data reduction software (for develping purposes only) ',bcolors.OKGREEN)
    printc('  version: '+ version,bcolors.OKGREEN)
    printc('  version_cmilos: '+ version_cmilos,bcolors.OKGREEN)
    printc('--------------------------------------------------------------',bcolors.OKGREEN)

    if json_input:

        # =========================================================================== #
        # READING CONFIG FILE AND PRINTING
        printc('--------------------------------------------------------------',bcolors.OKGREEN)
        printc(' Reading config json file '+json_input,bcolors.OKGREEN)
        with open(json_input) as j:
            CONFIG = json.load(j)

        verbose = CONFIG['verbose']
        input_data_dir = CONFIG['input_data_dir']
        data_f = CONFIG['data_f']
        shrink_mask = CONFIG['shrink_mask']
        center_method = CONFIG['center_method'] 
        hough_params = CONFIG['hough_params']
        instrument = CONFIG['instrument']
        flat_f = CONFIG['flat_f']
        dark_f = CONFIG['dark_f']
        dark_c = CONFIG['dark_c']
        flat_c = CONFIG['flat_c']
        flat_index = CONFIG['flat_index']
        norm_f = CONFIG['norm_f']
        flat_scaling = CONFIG['flat_scaling']
        prefilter_fits = CONFIG['prefilter_fits']
        prefilter = CONFIG['prefilter']
        output_dir = CONFIG['output_dir']
        rte = CONFIG['rte']
        correct_fringes = CONFIG['correct_fringes']
        correct_ghost = CONFIG['correct_ghost']
        putmediantozero = CONFIG['putmediantozero']
        debug = CONFIG['debug']
        vers = CONFIG['vers']
        RTE_code = CONFIG['RTE_code']  #specify which code will be used in the inversion. Default is cmilos (uses .txt ASCII input files)
        ItoQUV = CONFIG['ItoQUV']
        VtoQU = CONFIG['VtoQU']
        realign = CONFIG['realign']
        ind_wave = CONFIG['ind_wave']
        nlevel = CONFIG['nlevel']

        import pprint
        # Prints the nicely formatted dictionary
        pprint.pprint(CONFIG)#, sort_dicts=False)
    else:
        printc(' Using sequencial mode ',bcolors.OKGREEN)
        printc('    (hopefully with the right inputs since ERROR handling is not yet fully in place) ',bcolors.OKGREEN)

    #CHECK IF input is FITS OR FITS.GZ
    if  data_f.endswith('.fits'):
        filetype = '.fits'
    elif data_f.endswith('.fits.gz'):
        filetype = '.fits.gz'
    else:
        raise ValueError("input data type nor .fits neither .fits.gz")

    #TODO: check if inversions are requested: RTE and then find cmilos in such a case.
    if rte == 'RTE' or rte == 'CE' or rte == 'CE+RTE':
        printc(' RTE on. Looking for milos...',bcolors.OKGREEN)

        if RTE_code == 'cmilos': #TODO add .x after first orbit is done
            #MILOS_EXECUTABLE = 'milos.'+node().split('.')[0]+'.x'
            MILOS_EXECUTABLE = 'milos.'+node().split('.')[0]
            
        elif RTE_code == 'pmilos':
            MILOS_EXECUTABLE = 'pmilos'

        else: #TODO This is nonsense stuff
            #check if RTE is on
            if (rte != 'RTE') and (rte != 'RTE') and (rte != 'RTE'):
                print('Error setting RTE. RTE_code: ',RTE_code,' RTE option: ',rte)
                return

    #-----------------
    # READ DATA
    #-----------------
    
    data_filename = input_data_dir + data_f

    if os.path.isfile(data_filename):
        print("File exist")
    else:
        print("File not exist")

    try:
        data, header = fits_get(data_filename)
        printc('-->>>>>>> Reading Data file: '+data_filename,color=bcolors.OKGREEN)
        #
        # PXBEG1  =                  385 ; First read-out pixel in dimension 1            
        # PXEND1  =                 1664 ; Last read-out pixel in dimension 1             
        # PXBEG2  =                  385 ; First read-out pixel in dimension 2            
        # PXEND2  =                 1664 ; Last read-out pixel in dimension 2             

        DID = header['PHIDATID']
        ACC = header['ACCACCUM']
        printc('-->>>>>>> data DID '+DID,color=bcolors.OKGREEN)
        printc('          DATA IS DIVIDED by 256.   ',color=bcolors.OKGREEN)
        printc('-->>>>>>> Reshaping data to [wave,Stokes,y-dim,x-dim] ',color=bcolors.OKGREEN)
        zd,yd,xd = data.shape
        data = np.reshape(data,(zd//4,4,yd, xd))
        data = data / 256. #from fix to 32
        data = np.ascontiguousarray(data)

        #/ PHI_FITS_FPA_settings 
        # FPIMGCMD= 8 / FPA image command 
        # FPA_SROW= 0 / FPA start row setting FPA_EROW= 1022 / FPA end row setting 
        # FPA_NIMG= 20 / FPA number of images set FPEXPSTC= 1592452786 / [s] FPA exposure start time coarse 
        # FPEXPSTF= 699245 / [us] FPA exposure start time fine 
        # INTTIME = 0.01 / [s] Exposure time of single readout 
        # TELAPSE = 58.1974400877953 / [s] 
        # Elapsed time between start and end of obser
        # NSUMEXP = 480 / Number of detector readouts 
        # XPOSURE = 4.8 / [s] Total effective exposure time 
        # ACCLENGT= 4194304 / ACCU number of pixel set 
        # ACCNROWS= 6 / ACCU number of rows set 
        # ACCROWIT= 1 / ACCU number of row iterations set 
        # ACCNCOLS= 4 / ACCU number of columns set 
        # ACCCOLIT= 1 / ACCU number of column iterations set 
        # ACCACCUM= 20 / ACCU number of accumulations set 
        # ACCADDR = 0 / ACCU readout address (start) 

    except Exception:
        printc("ERROR, Unable to open fits file: {}",data_filename,color=bcolors.FAIL)
        return 0

    header['history'] = ' Data processed with phifdt_pipe.py '+ version
    header['history'] = '      and time '+ str(datetime.datetime.now())
    header['history'] = ' Parameters normalize_flat: '+ str(norm_f)
    header['history'] = ' Parameters flat_scaling: '+ str(flat_scaling)
    header['history'] = ' Parameters shrink_mask: '+ str(shrink_mask)
    header['history'] = ' Parameters center_method: '+ str(center_method)
    header['history'] = ' Parameters Hough: '+ str(hough_params)

    if verbose:
        plib.show_one(data[0,0,:,:],vmin=0,xlabel='pixel',ylabel='pixel',title='Data first image raw (1 of 24)',cbarlabel='DN',save=None,cmap='gray')
    
    #    * CAL_RTE=                990510 / ok
    #    * CAL_SCIP= 'None'               / Onboard scientific data analysis
    #    * RTE_ITER=           4294967295 / Number RTE inversion iterations
    #    * PHIDATID= '142010402'          / PHI dataset Id

    #-----------------
    # TAKE DATA DIMENSIONS AND SCALING
    #-----------------
    PXBEG1  = int(header['PXBEG1']) - 1           
    PXEND1  = int(header['PXEND1']) - 1          
    PXBEG2  = int(header['PXBEG2']) - 1           
    PXEND2  = int(header['PXEND2']) - 1   
    printc('Dimensions: ',PXBEG1, PXEND1, PXBEG2, PXEND2,color=bcolors.OKGREEN)

    if xd != (PXEND1 - PXBEG1 + 1) or yd != (PXEND2 - PXBEG2 + 1):
        printc('ERROR, Keyword dimensions and data array dimensions dont match ',color=bcolors.FAIL)
        return 0
    if xd < 2047:    
        printc('         data cropped to: [',PXBEG1,',',PXEND1,'],[',PXBEG2,',',PXEND2,']',color=bcolors.WARNING)
    
    data_scale = fits_get(data_filename,scaling = True)

    #-----------------
    # CHECK PMP Temperature
    #-----------------
    if instrument == 'auto':
        instrument = check_pmp_temp(header)
        instrument = 'FDT'+instrument
        printc('Instrument: ',instrument,color=bcolors.OKGREEN)

    else:
        pass

    #-----------------
    # READ FLAT FIELDS
    #-----------------

    if flat_c:
        printc('-->>>>>>> Reading flat file'+flat_f,color=bcolors.OKGREEN)
        printc('          Assumes they are already normalized to ONE ',color=bcolors.OKGREEN)
        printc('          input should be [wave X Stokes,y-dim,x-dim].',color=bcolors.OKGREEN)

        try:
            dummy,flat_header = fits_get(flat_f)
            fz_d,fy_d,fx_d = dummy.shape
        except Exception:
            printc("ERROR, something happened while reading the file: {}",flat_f,color=bcolors.FAIL)
            return 0

        flat = np.zeros([24,2048,2048]).astype(np.float32)
        PXBEG1_f  = int(flat_header['PXBEG1']) - 1           
        PXEND1_f  = int(flat_header['PXEND1']) - 1          
        PXBEG2_f  = int(flat_header['PXBEG2']) - 1           
        PXEND2_f  = int(flat_header['PXEND2']) - 1 
        if fx_d < 2047:    
            printc('         input flat was cropped to: [',PXBEG1_f,',',PXEND1_f,'],[',PXBEG2_f,',',PXEND2_f,']',color=bcolors.WARNING)

        flat[:,PXBEG1_f:PXEND1_f+1,PXBEG2_f:PXEND2_f+1] = dummy
        del dummy

        printc('-->>>>>>> Reshaping Flat to [wave,Stokes,y-dim,x-dim] ',color=bcolors.OKGREEN)
        fz,fy,fx = flat.shape
        flat = np.reshape(flat,(fz//4,4,fy,fx))

        # if SIX_FLATS:
        #     for i in range(6):
        #         mm = np.mean(flat[i,:,:,:],axis = 0)
        #         flat[i,:,:,:] = mm[np.newaxis,:,:]

        if verbose:
            plib.show_one(flat[0,0,:,:],xlabel='pixel',ylabel='pixel',title='Flat first image raw (1 of 24)',cbarlabel='Any (as input)',save=None,cmap='gray')

    else:
        printc('-->>>>>>> No flats mode                    ',color=bcolors.WARNING)

    #-----------------
    # READ AND CORRECT DARK FIELD
    #-----------------
    if dark_c:
        data,header  = phi_correct_dark(dark_f,data,header,data_scale,verbose = verbose)
    else:
        printc('-->>>>>>> No darks mode                    ',color=bcolors.WARNING)

    #-----------------
    # FIND DATA CENTER 
    #-----------------

    printc('-->>>>>>> finding the center of the solar disk (needed for masking) ',color=bcolors.OKGREEN)
    try:
        if center_method == 'Hough':
            inner_radius,outer_radius,steps = hough_params
            c, radius,threshold = find_circle_hough(data[0,0,:,:],inner_radius,outer_radius,steps,threshold = 0.01,normalize=False,verbose=False)
            #c = np.roll(c,1)
            cx = c[0]
            cy = c[1]
            #TBE PUT IN CORRECT UNITS
        elif center_method  == 'circlefit':
            cy,cx,radius=find_center(data[0,0,:,:])  #OJO Cy... Cx
            c = np.array([int(cx),int(cy)])   #El vector es [0,1,2,...] == [x,y,z,...] == [cx,cy,cz,...] Pero esto ultimo esta al reves
            radius = int(radius)
        elif center_method == None:
            #get from header
            cx = header['CRPIX1']
            cy = header['CRPIX2']
            c = np.array([int(cx),int(cy)])   #El vector es [0,1,2,...] == [x,y,z,...] == [cx,cy,cz,...] Pero esto ultimo esta al reves
            radius = header['RSUN_ARC']/header['CDELT1']
        else:  
            raise ValueError("ERROR in center determination method - check input 'circlefit','Hough',null/None") 
    except ValueError as err:
        print(err.args)
        return 0

    #Uptade header with new centers
    if center_method == 'Hough' or center_method  == 'circlefit':
        printc('          Uptade header with new center:',color=bcolors.OKBLUE)
        printc('          OLD center:',color=bcolors.OKBLUE)
        printc('                  at: CRPIX1[x]=',header['CRPIX1'],' CRPIX2[y]=',header['CRPIX2'],' radius=',radius,color=bcolors.OKBLUE)
        header['history'] = ' CRPIX 1 and CRPIX2 uptated from ' + str(header['CRPIX1'])+ ' and ' + str(header['CRPIX2'])
        header['CRPIX1'] = (round(cx, 2))
        header['CRPIX2'] = (round(cy, 2))
        printc('          NEW center:',color=bcolors.OKBLUE)
        printc('                  at: CRPIX1[x]=',header['CRPIX1'],' CRPIX2[y]=',header['CRPIX2'],' radius=',radius,color=bcolors.OKBLUE)
        printc('ATTENTION: Keywords CRVAL1 and CRVAL2 are NOT updated but should be SET to zero!!!!',color=bcolors.FAIL)
    else:
        printc('          Using header image center:',color=bcolors.OKBLUE)
        printc('                  at: CRPIX1[x]=',header['CRPIX1'],' CRPIX2[y]=',header['CRPIX2'],' radius=',radius,color=bcolors.OKBLUE)

    #OJO.
    # find_circle_hough devuelve c = c[0] = x and c[1] = y !!!!!!!!!!!!!!
    # Esto viene porque en el KLL esta definido así (al reves) en la rutina votes()
 
    #-----------------
    # TAKE ONLY DISK WITH MARGIN
    #-----------------
    printc('-->>>>>>> Creating a mask for RTE with ',shrink_mask,' px margin')
    size_of_mask = radius - shrink_mask
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

        # printc('-->>>>>>> Reshaping flat to [wave,Stokes,y-dim,x-dim]',color=bcolors.OKGREEN)
        # flat = np.reshape(flat,(fz//4,4,fy, fx))

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
            if (len(flat_index)) == 6:
                print('          Changing flat index to ',flat_index)
        except:
            flat_index = [0,1,2,3,4,5]
        for p in range(4):
            for l in range(int(zd//4)):
                print('          ... pol: ',p,' wave: ',l,' index: ',flat_index[l])
                dummy_flat = (flat[flat_index[l],p,PXBEG2:PXEND2+1,PXBEG1:PXEND1+1]/float(flat_scaling))
                if norm_f:
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
            plib.show_one(data[cpos,0,:,:],vmax=None,vmin=0,xlabel='pixel',ylabel='pixel',title='Data / flat at continuum',cbarlabel='DN',save=None,cmap='gray')

    #-----------------
    # CORRECT PREFILTER 
    #-----------------

    if prefilter:
        data,header  = phi_correct_prefilter(prefilter_fits,header,data,voltagesData,verbose = verbose)


    #-----------------
    # GHOST CORRECTION  
    #-----------------

    if correct_ghost:
        data,header = phi_correct_ghost(data,header,radius,verbose = True)

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

    data, header = phi_apply_demodulation(data,instrument,header=header)

    if verbose == 1:
        plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'],zoom = 3,svmin=[0,-0.004,-0.004,-0.004],svmax=[1.2,0.004,0.004,0.004])

    # with pyfits.open(data_filename) as hdu_list:
    #     hdu_list[0].data = data
    #     hdu_list[0].header = header
    #     hdu_list.writeto('dummy.fits', clobber=True)

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
        plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'],zoom = 3,svmin=[0,-0.004,-0.004,-0.004],svmax=[1.2,0.004,0.004,0.004])

    #-----------------
    # GHOST CORRECTION  AFTER DEMODULATION
    #-----------------

    # if correct_ghost:
    #     data,header = phi_correct_ghost_dm(data,header,radius,verbose = verbose)

    #-----------------
    # CROSS-TALK CALCULATION 
    #-----------------
    if ItoQUV or putmediantozero:
            factor_media = 0.8  # 80% of the disk
            rrx_m = [int(c[0]-radius*factor_media),int(c[0]+radius*factor_media)]
            rry_m = [int(c[1]-radius*factor_media),int(c[1]+radius*factor_media)]
            maski,coords = generate_circular_mask([xd-1,yd-1],radius*factor_media,radius*factor_media)
            maski = shift(maski, shift=(c[0]-xd//2,c[1]-yd//2), fill_value=0).astype(int)

    if ItoQUV:
        printc('-->>>>>>> Cross-talk correction from Stokes I to Stokes Q,U,V --',color=bcolors.OKGREEN)
        printc('          Using ',factor_media*100,'% of the disk                     ',color=bcolors.OKGREEN)
        printc('          Crosstalk evaluated in x = [',rrx_m[0],':',rrx_m[1],'] y = [',rry_m[0],':',rry_m[1],']',' using ',factor_media*100,"% of the disk",color=bcolors.OKBLUE)

        if ind_wave:
            for i in range(zd//4):
                printc('          Individual wavelengths....',color=bcolors.OKBLUE)
                broadcastd = data[i,:,rry_m[0]:rry_m[1],rrx_m[0]:rrx_m[1]]
                data_dummy = data[:,:,rry_m[0]:rry_m[1],rrx_m[0]:rrx_m[1]]*0. + broadcastd[np.newaxis,:,:,:]
                cQ,cU,cV = crosstalk_ItoQUV(data_dummy[:,:,rry_m[0]:rry_m[1],rrx_m[0]:rrx_m[1]],npoints=10000,verbose=verbose)
                #-----------------
                # CROSS-TALK CORRECTION 
                #-----------------
                printc('          Applying cross-talk correction...',color=bcolors.OKGREEN)
                data[i,1,:,:] = data[i,1,:,:] - cQ[0]*data[i,0,:,:] - cQ[1]
                data[i,2,:,:] = data[i,2,:,:] - cU[0]*data[i,0,:,:] - cU[1]
                data[i,3,:,:] = data[i,3,:,:] - cV[0]*data[i,0,:,:] - cV[1]
            if verbose:
                plt.hist(data[0,1,maski > 0].flatten(), bins='auto')
                plt.title('Stokes Q')
                plt.hist(data[0,2,maski > 0].flatten(), bins='auto')
                plt.title('Stokes U')
                plt.hist(data[0,3,maski > 0].flatten(), bins='auto')
                plt.title('Stokes V')

                plt.show()
        else:
            cQ,cU,cV = crosstalk_ItoQUV(data[:,:,rry_m[0]:rry_m[1],rrx_m[0]:rrx_m[1]],verbose=verbose,npoints=10000)

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
            plib.show_hist(data[0,1, maski > 0].flatten(), bins='auto',title=' ',leave = 'open',color='green')
            plib.show_hist(data[0,2, maski > 0].flatten(), bins='auto',title=' ',leave = 'open',color='red')
            plib.show_hist(data[0,3, maski > 0].flatten(), bins='auto',title='Stokes Q/U/V - no zero',color='blue')
            plib.show_four_row(data[2,0,:,:],data[2,1,:,:],data[2,2,:,:],data[2,3,:,:],title=['I - corr','Q - corr','U - corr','V - corr'],zoom=3,svmin=[0,-0.004,-0.004,-0.004],svmax=[1.2,0.004,0.004,0.004])
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
    
    if VtoQU:
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
            plib.show_two(data[3,1,ry[0]:ry[1],rx[0]:rx[1]],datao[3,1,ry[0]:ry[1],rx[0]:rx[1]],block=False,title=['Stokes Q','Stokes Q corrected'],zoom=3)
            plib.show_two(data[3,2,ry[0]:ry[1],rx[0]:rx[1]],datao[3,2,ry[0]:ry[1],rx[0]:rx[1]],block=False,title=['Stokes U','Stokes U corrected'],zoom=3)
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

    if correct_fringes == 'auto' or correct_fringes == 'manual':
        if verbose:
            plib.show_four_row(data[2,0,:,:],data[2,1,:,:],data[2,2,:,:],data[2,3,:,:],title=['I - before fringe','Q','U','V'],zoom=3,svmin=[0,-0.004,-0.004,-0.004],svmax=[1.2,0.004,0.004,0.004])
        data, header = phi_correct_fringes(data,header,option=correct_fringes,verbose=verbose)
        if verbose:
            plib.show_four_row(data[2,0,:,:],data[2,1,:,:],data[2,2,:,:],data[2,3,:,:],title=['I - after fringe','Q','U','V'],zoom=3,svmin=[0,-0.004,-0.004,-0.004],svmax=[1.2,0.004,0.004,0.004])
    elif correct_fringes == False:
        pass
    else:
        printc('Error in option finge correction. Options are "manual", "auto" or false. Given: ',color=bcolors.WARNING)
        print(correct_fringes)
    #-----------------
    # MEDIAN TO CERO
    #-----------------

    if putmediantozero:
        printc('-->>>>>>> Putting median to zero ',color=bcolors.OKGREEN)
        printc('          Median evaluated in x = [',rrx_m[0],':',rrx_m[1],'] y = [',rry_m[0],':',rry_m[1],']',' using ',factor*100,"% of the disk",color=bcolors.OKBLUE)
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
            
    # plib.show_four_row(data[0,0,:,:],data[0,1,:,:],data[0,2,:,:],data[0,3,:,:],title=['I','Q','U','V'],zoom=3,svmin=[0,-0.004,-0.004,-0.004],svmax=[1.2,0.004,0.004,0.004])
    # ctnd = data[0,1:3,:,:]
    # for i in range(zd//4):
    #     data[i,1:3,:,:] = data[i,1:3,:,:] - ctnd
    # plib.show_four_row(data[0,0,:,:],data[0,1,:,:],data[0,2,:,:],data[0,3,:,:],title=['I','Q','U','V'],zoom=3,svmin=[0,-0.004,-0.004,-0.004],svmax=[1.2,0.004,0.004,0.004])

    if verbose == 1:
        plib.show_hist(data[0,1, maski > 0].flatten(), bins='auto',title=' ',leave='open',color='green')
        plib.show_hist(data[0,2, maski > 0].flatten(), bins='auto',title=' ',leave='open',color='red')
        plib.show_hist(data[0,3, maski > 0].flatten(), bins='auto',title='Stokes Q/U/V - zero',color='blue')
        plib.show_four_row(data[3,0,:,:],data[3,1,:,:],data[3,2,:,:],data[3,3,:,:],title=['I','Q','U','V'],zoom=3,svmin=[0,-0.004,-0.004,-0.004],svmax=[1.2,0.004,0.004,0.004])

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
    dirs = ['pngs','level2']
        
    for checkit in dirs:
        check_dir = os.path.isdir(output_dir+checkit)
        if not check_dir:
            os.makedirs(output_dir+checkit)
            print("created folder : ", output_dir+checkit)
        else:
            print(output_dir+checkit, "folder already exists.")

    printc('---------------------------------------------------------',color=bcolors.OKGREEN)
    if outfile == None:
        #basically replace L1 by L1.5
        try:
            outfile_L2 = set_level(data_f,'L1','L2')
            outfile_L2 = set_level(outfile_L2,'ilam','stokes')
            outfile_L2 = append_id(outfile_L2,filetype,vers,DID) 
            # outfile_L2 = outfile_L2.split('V')[0] + 'V' + vers+ '_' + DID  + filetype
        except:
            outfile_L2 = set_level(data_f,'L0','L2')
            outfile_L2 = set_level(outfile_L2,'ilam','stokes')
            outfile_L2 = append_id(outfile_L2,filetype,vers,DID) 
            # outfile_L2 = outfile_L2.split('V')[0] + 'V' + vers+ '_' + DID  + filetype

    else:
        outfile_L2 = outfile

    printc(' Saving data to: ',output_dir+'level2/'+outfile_L2)

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
        hdu_list.writeto(output_dir+'level2/'+outfile_L2, clobber=True)
    #        hdu_list.writeto(directory+outfile+'_L1.fits', clobber=True)

    # with pyfits.open(data_f) as hdu_list:
    #     hdu_list[0].data = mask
    #     hdu_list.writeto(directory+outfile+'_red-mask.fits', clobber=True)

    #-----------------
    # INVERSION OF DATA WITH CMILOS
    #-----------------
    if rte == 'RTE' or rte == 'CE' or rte == 'CE+RTE':

        printc('---------------------RUNNING CMILOS --------------------------',color=bcolors.OKGREEN)
        
        rte_invs = np.zeros((12,yd,xd)).astype(float)
        rte_invs[:,ry[0]:ry[1],rx[0]:rx[1]] = generate_level2(data[:,:,ry[0]:ry[1],rx[0]:rx[1]],wave_axis,rte,output_dir,milos_executable = MILOS_EXECUTABLE)

        rte_invs_noth = np.copy(rte_invs)
        umbral = 3.

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

        # writeto = set_level(outfile_L2,'stokes','_RTE')
        # np.savez_compressed(output_dir+outfile_L2, rte_invs=rte_invs, rte_invs_noth=rte_invs_noth,mask=mask)

        b_los = rte_invs_noth[2,:,:]*np.cos(rte_invs_noth[3,:,:]*np.pi/180.)*mask
        # b_los = rte_invs_noth[2,:,:]*np.cos(rte_invs_noth[3,:,:]*np.pi/180.)*mask
        # b_los = rte_invs_noth[2,:,:]*np.cos(rte_invs_noth[3,:,:]*np.pi/180.)*mask
        # b_los = rte_invs_noth[2,:,:]*np.cos(rte_invs_noth[3,:,:]*np.pi/180.)*mask
        # b_los = np.zeros((2048,2048))
        # b_los[PXBEG2:PXEND2+1,PXBEG1:PXEND1+1] = b_los_cropped

        v_los = rte_invs_noth[8,:,:] * mask
        # v_los = np.zeros((2048,2048))
        # v_los[PXBEG2:PXEND2+1,PXBEG1:PXEND1+1] = v_los_cropped
        if verbose:
            plib.show_one(v_los,vmin=-2.5,vmax=2.5,title='LoS velocity')
            plib.show_one(b_los,vmin=-30,vmax=30,title='LoS magnetic field')

        printc('  ---- >>>>> Updating L2 header.... ',color=bcolors.OKGREEN)

        header['history'] = ' RTE CMILOS INVERTER: '+ rte
        header['history'] = ' CMILOS VER: '+ version_cmilos

        if 'RTE_ITER' in header:  # Check for existence
            header['RTE_ITER'] = str(15)
        else:
            header.set('RTE_ITER', str(15), 'Number RTE inversion iterations',after='CAL_SCIP')

        printc('  ---- >>>>> Saving L2 data.... ',color=bcolors.OKGREEN)

        #HRT version

        with pyfits.open(data_filename) as hdu_list:
            hdu_list[0].data = rte_invs_noth[2,:,:] * mask
    #            header = hdu_list[0].header
            hdu_list[0].header = header
            writeto = set_level(outfile_L2,'stokes','bmag')
            hdu_list.writeto(output_dir+'level2/'+writeto, clobber=True)

        # with pyfits.open(data_filename) as hdu_list:
            hdu_list[0].data = rte_invs_noth[3,:,:] * mask
            # hdu_list[0].header = header
            writeto = set_level(outfile_L2,'stokes','binc')
            hdu_list.writeto(output_dir+'level2/'+writeto, clobber=True)

        # with pyfits.open(data_filename) as hdu_list:
            hdu_list[0].data = rte_invs[4,:,:] * mask
            # hdu_list[0].header = header
            writeto = set_level(outfile_L2,'stokes','bazi')
            hdu_list.writeto(output_dir+'level2/'+writeto, clobber=True)

        # with pyfits.open(data_filename) as hdu_list:
            hdu_list[0].data = b_los
            # hdu_list[0].header = header
            writeto = set_level(outfile_L2,'stokes','blos')
            hdu_list.writeto(output_dir+'level2/'+writeto, clobber=True)

        # with pyfits.open(data_filename) as hdu_list:
            hdu_list[0].data = v_los
            # hdu_list[0].header = header
            writeto = set_level(outfile_L2,'stokes','vlos')
            hdu_list.writeto(output_dir+'level2/'+writeto, clobber=True)

        # with pyfits.open(data_filename) as hdu_list:
            hdu_list[0].data = rte_invs[9,:,:]+rte_invs[10,:,:]
            # hdu_list[0].header = header
            writeto = set_level(outfile_L2,'stokes','icnt')
            hdu_list.writeto(output_dir+'level2/'+writeto, clobber=True)

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
        writeto = set_level(outfile_L2,'stokes','vlos')
        writeto = set_level(writeto,filetype,FIGUREOUT) 

        plt.savefig(output_dir+'pngs/'+writeto,dpi=300)
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

        writeto = set_level(outfile_L2,'stokes','blos')
        writeto = set_level(writeto,filetype,FIGUREOUT) 
        plt.savefig(output_dir+'pngs/'+writeto,dpi=300)
        plt.close()

        #-----------------
        # PLOTS AZIMUTH
        #-----------------
        Zm = np.ma.masked_where(mask == 1, mask)

        plt.figure(figsize=(10, 10))
        ax = plt.gca()
        plt.title('PHI-FDT field azimuth',size=20)

        # Hide grid lines
        ax.grid(False)
        # Hide axes ticks
        ax.set_xticks([])
        ax.set_yticks([])

        im = ax.imshow(rte_invs[4,:,:] * mask, cmap='bwr',vmin=0,vmax=180)

        divider = plib.make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plib.plt.colorbar(im, cax=cax)
        cbar.set_label('[km/s]')
        cbar.ax.tick_params(labelsize=16)

        writeto = set_level(outfile_L2,'stokes','bazi')
        writeto = set_level(writeto,filetype,FIGUREOUT) 

        plt.savefig(output_dir+'pngs/'+writeto,dpi=300)
        plt.close()

        printc('--------------------- END  ----------------------------',color=bcolors.FAIL)

    return wave_axis
    # if rte == 'cog':
    #     printc('---------------------RUNNING COG --------------------------',color=bcolors.OKGREEN)
    #     wavelength = 6173.3356
    #     v_los,b_los = cog(data,wavelength,wave_axis,lande_factor=3,cpos = cpos)

    #     #-----------------
    #     # MASK DATA AND SAVE
    #     #-----------------

    #     v_los = v_los * mask
    #     b_los = b_los * mask
    #     plib.show_one(v_los,vmin=-1.5,vmax=1.5)
    #     plib.show_one(b_los,vmin=-150,vmax=150)
    #     if verbose == 1:
    #         plib.show_one(v_los,vmin=-2.5,vmax=2.5)
    #         plib.show_one(b_los,vmin=-150,vmax=150)

    #     with pyfits.open(data_f) as hdu_list:
    #         hdu_list[0].data = v_los
    #         hdu_list[0].header = header
    #         writeto = set_level(outfile_L2,'ilam','vlos-cog')
    #         writeto = set_level(writeto,'.fits','.png')
    #         plt.savefig(directory+writeto,dpi=300)

    #     with pyfits.open(data_f) as hdu_list:
    #         hdu_list[0].data = b_los
    #         hdu_list[0].header = header
    #         writeto = set_level(outfile_L2,'ilam','blos-cog')
    #         writeto = set_level(writeto,'.fits','.png')
    #         plt.savefig(directory+'pngs/'+writeto,dpi=300)

    # return
    #-----------------
    # CORRECT CAVITY
    #-----------------

    # try:
    #     if cavity == 0:
    #         cavity=np.zeros_like(vl)
    #         print("zerooooooooooo")
    # except:
    #     # cavity = cavity[ry[0]:ry[1],rx[0]:rx[1]]
    #     pass

    # factor = 0.5
    # rrx = [int(c[0]-r*factor),int(c[0]+r*factor)]
    # rry = [int(c[1]-r*factor),int(c[1]+r*factor)]
    # print(rrx,rry,' check these for los vel calib')
    # off = np.mean(vl[rry[0]:rry[1],rrx[0]:rrx[1]])
    # vl = vl - off #- cavity
    # print('velocity offset ',off)

    # # cavity,h = phi.fits_read('HRT_cavity_map_IP5.fits')
    # # cavity = cavity * 0.3513e-3/6173.*300000. #A/V 
    # # cavity = cavity - np.median(cavity[800:1200,800:1200])
    
    # return vl


