#=============================================================================
# Project: SoPHI
# File:    phifdt_flat_fixpoint.py
# Author:  David Orozco Suárez (orozco@iaa.es)
# Contributors: Tobias Lange and Nestor Albelo (albelo@mps.mpg.de)
#-----------------------------------------------------------------------------
# Description: Pipeline implementation for calculating the FDT flat
#              Includes three algorithms for flat calculation and 
#              the circular Hough Transform in the frequency domain.
# This version was modified from phifdt_flat.py to simulate with fix point precision
#-----------------------------------------------------------------------------

from math import nan
import numpy as np
import matplotlib.pyplot as plt
from .tools import printc,bcolors,timeit
from .phi_gen import * 
from .phi_utils import *
from .phi_fits import *
from .phi_reg import *
from .phifdt_pipe_modules import phi_correct_dark
from SPGPylibs.GENtools import *

BIT_DEPTH: str = '>i4'   #"images are 32-bit big-endian integer"
def fixpoint(imagen):
    if imagen.dtype != BIT_DEPTH:
        imagen = imagen.astype(BIT_DEPTH)       
    print('image bit depth ',imagen.dtype)
    return imagen

class phidata():

    bit_depth: str = '>i4'   #"images are 32-bit big-endian integer"
    dark_offset: float = 1   #"images are 32-bit big-endian integer"

    import reprlib
    r = reprlib.Repr()
    r.maxlist = 4        # max elements displayed for lists
    r.maxstring = 100    # max characters displayed for strings

    def __init__(self, file=''):
        self.file = file
        self.darkc: bool = False
        self.flatc: bool = False
        self.image = None
        self.header = None
        self.imageSummary = None
        self.imageSummary_head = None
        self.scaling = {"Present": [False,True], "scaling": [0,0], 'bit-depth': None}
        self.DID = None

    @staticmethod    
    def mprint(what,*args, **kw):
        print(phidata.r.repr(what),*args, **kw)
    
    def dscale(self):
        if self.image.dtype != phidata.bit_depth:
            self.image = self.image.astype(phidata.bit_depth)       
        printc('image bit depth ',self.image.dtype,color=bcolors.OKGREEN)

    def set_file(self,file):
        self.file = file

    def info(self):
        options = vars(self)
        print('Image info: ')
        for item in options.items():
            extension = len(item[0]) 
            spaces = 4 - extension
            string_val = " " * spaces
            self.mprint("%s: %s %s" % (item[0],string_val,item[1]))
    
    def load(self,info=False):
        if info == True:
            try:
                return pyfits.info(self.file)
            except Exception:
                print("Unable to open fits file: {}",self.file)        
                raise
        else:
            try:
                with pyfits.open(self.file) as hdu_list:
                    head = 0
                    self.image = hdu_list[head].data.astype(phidata.bit_depth)
                    self.header = hdu_list[head].header
                    self.DID = self.header['PHIDATID']

                    #get scaling
                    index = 1
                    bear_moved = False
                    while bear_moved == False:
                        try:
                            dummy_head = hdu_list[index].header
                            if dummy_head['EXTNAME'] == 'PHI_FITS_imageSummary':
                                self.imageSummary_head = dummy_head
                                self.imageSummary = hdu_list[index].data
                                #case 1 if that there is only ONE scaling (untouched data)
                                if len(self.imageSummary) == 1:
                                    self.scaling["Present"][0] = False
                                    self.scaling["Present"][1] = True
                                    self.scaling["scaling"][0] = 0.
                                    self.scaling["scaling"][1] = float(self.imageSummary[0][12])
                                #case 2 if that there is more than TWO scaling data
                                if len(self.imageSummary) > 2:
                                    #check the first one from below and if it is IMGFMT_16_0_S store it and continue
                                    if self.imageSummary[-1][3] == 'IMGFMT_16_0_S':
                                        self.scaling['bit-depth'] = '16'
                                        self.scaling["Present"][0] = True
                                        self.scaling["Present"][1] = True
                                        self.scaling["scaling"][0] = float(self.imageSummary[-1][12])
                                        self.scaling["scaling"][1] = float(self.imageSummary[-3][12])
                                    if self.imageSummary[-1][3] == 'IMGFMT_24_8':
                                        self.scaling['bit-depth'] = '24.8'
                                        self.scaling["Present"][0] = False
                                        self.scaling["Present"][1] = True
                                        self.scaling["scaling"][0] = 0.
                                        self.scaling["scaling"][1] = float(self.imageSummary[-1][12])
                        except:
                            bear_moved = True
                        index += 1
                    #set image to 24.8 ALWAYS
                    if self.scaling['bit-depth'] == '24.8':
                        printc('image scale '+self.scaling['bit-depth'],color=bcolors.OKGREEN)
                    if self.scaling['bit-depth'] == '16':
                        self.image = self.image * self.scaling["scaling"][1] / (self.scaling["scaling"][0] + 1)   #81920/128 para 24.8 
                        printc('image scale '+self.scaling['bit-depth']+' --->>> 24.8',color=bcolors.OKGREEN)
                        self.scaling['bit-depth'] = '24.8'
                    self.dscale()

            except Exception:
                print("Unable to open fits file: {}",self.file)
                raise

    def apply_dark(self,dark,verbose: bool =False):

        #def phi_correct_dark(dark_f,data,header,data_scale,verbose = False,get_dark = False):
        printc('Dark correction                   ',color=bcolors.OKGREEN)
        #check if dark was loaded
        if not(dark.DID):
            printc('     reading out dark                   ',color=bcolors.OKGREEN)
            dark.load()
        #check scaling
        printc('Data scaling ' + self.scaling['bit-depth'],color=bcolors.OKGREEN)
        printc('Dark scaling ' + dark.scaling['bit-depth'],color=bcolors.OKGREEN)

        printc('-->>>>>>> Correcting dark current.',color=bcolors.OKGREEN)
        PXBEG1  = int(self.header['PXBEG1']) - 1           
        PXEND1  = int(self.header['PXEND1']) - 1          
        PXBEG2  = int(self.header['PXBEG2']) - 1           
        PXEND2  = int(self.header['PXEND2']) - 1   
        PXBEG1d  = int(dark.header['PXBEG1']) - 1           
        PXEND1d  = int(dark.header['PXEND1']) - 1          
        PXBEG2d  = int(dark.header['PXBEG2']) - 1           
        PXEND2d  = int(dark.header['PXEND2']) - 1   

        if self.darkc == False:
            # self.image = self.image.astype(phidata.bit_depth) - \
                # dark.image[np.newaxis,PXBEG2:PXEND2+1,PXBEG1:PXEND1+1].astype(phidata.bit_depth)
            if verbose:
                dummy = self.image[0,:,:]
            self.image = self.image - (dark.image[np.newaxis,PXBEG2:PXEND2+1,PXBEG1:PXEND1+1]*phidata.dark_offset).astype(BIT_DEPTH)
            #update data header
            self.header['CAL_DARK'] = dark.DID
            self.darkc = True

        else:
            printc('-->>>>>>> data already corrected.',color=bcolors.OKGREEN)

        print(PXBEG1d,PXEND1d,PXBEG2d,PXEND2d)

        if verbose:
            md = np.mean(dark.image)
            plib.show_three(dark.image,dummy,self.image[0,:,:],vmin=[0,0,0],vmax=[md*1.2,md*1.2,md*1.2],block=True,pause=0.1,title=['Dark','Data','Data after dark correction'],
                xlabel='Pixel',ylabel='Pixel',cmap='gray')
            print(np.mean(dark.image[0:100,0:100]),np.mean(self.image[0,0:100,0:100]))

def centers_flat(n_images,inner_radius,outer_radius,steps,r_width,binmask,imsize,verbose=None):
    ############################
    #FIND CENTERS - PART OD DO_HOUGH
    ############################
    centers = []
    radius = []
    radii = np.linspace(inner_radius, outer_radius, steps + 1)
    printc('Analizing ', n_images, ' images',color = bcolors.OKGREEN)

    for i in range(n_images): 
    #acc_conv = find_Circles(
    #    binmask[i], radii_coarse, r_width_coarse, verbose=verbose, full=True)
        acc_conv = find_Circles_ida(binmask[i], radii, r_width)
        center,rad,c,d = votes(acc_conv, radii)

        centers.append(center)
        radius.append(rad)
        printc('Found center: ', centers[i], ' and radius: ', radius[i],color = bcolors.WARNING)
        if verbose == True:
            fig = plt.figure(frameon=False)
            im1 = plt.imshow(binmask[i], cmap=plt.cm.gray, alpha=.5)
            circle_fit = bin_annulus(
                imsize, radius[i], 1, full=False).astype(float)
            dd = np.array(centers[i])
            dx = dd[0] - imsize[0]//2
            dy = dd[1] - imsize[1]//2
            circle_fit = shift(circle_fit, shift=[dx,dy])
            im2 = plt.imshow(circle_fit, cmap=plt.cm.gray, alpha=.5)
            plt.show()

    return centers,radius

@timeit
def do_hough(image,inner_radius, outer_radius, steps, org_centers=None,method='prewitt',save=False,
            dhtr=10,normalize = False,verbose=False,Otsu = None,threshold = 0.15):
    '''
    Calculates the position and radious of the solar disk in a set of input images using the Hough transform.

    Parameters
    ----------
    image : (K, N, M) ndarray
        List or numpy array of K Grayscale images of NxM size.
    inner_radius : int
        Minimum search radious 
    outer_radius : int
        Maximum search radious 
    steps: int
        Number of steps to look for solar radius. 
        step is used to generate:
            (1): coarse find jumps: np.linspace(inner_radius, outer_radius, steps)
            (2): width of the ring for crosscorrelating the disk: (outer_radius - inner_radius)//steps * 2
            (3):if step is a negative number then uses FM find model
                -#-
                4 iterations
                1] inner_radius = 152;  outer_radius = 1048; steps = 64; 15 iterations 
                152_____________600_____________1048
                --|---|---|---|---|---|---|---|---|---|---|---|---|---|---|--
                2] inner_radius = Prev.Radius-32;  outer_radius = Prev.Radius+32; steps = 16; 5 iterations 
                ---------|---------------|---------------|---------------|---------------|--------
                3] inner_radius = Prev.Radius-8;  outer_radius = Prev.Radius+8; steps = 4; 5 iterations 
                -----------|---------------|---------------|---------------|---------------|-----------
                4] inner_radius = Prev.Radius-2;  outer_radius = Prev.Radius+2; steps = 1; 5 iterations 
                -----------|---------------|---------------|---------------|---------------|-----------
                -#-
    org_centers = org_centers: numpy array [K,2] centers for comparison (they are not used)
    method = method: method for finding the limb boundary. default = 'prewitt'
        more info look FindEdges()
    save = False: save the centers as 'hough_centers.txt' -> ASCII (centers_fine,radii_fine)
    dhtr = 10:
    normalize = False:
    verbose = False:
    Otsu = None:
    threshold = 0.15:

    Returns
    -------
    centers : numpy int array of [K,2] elements where [i,0] = x-centers and [i,1] = y-centers
    radious : numpy int array of [K] elements containing the radious of the K-images in pixels 

    Raises
    ------

    References
    ----------
    [1] C. Hollitt, Machine Vision and Applications (2013) 24:683–694 DOI 10.1007/s00138-012-0420-x

    Examples
    --------
    >>> import SPGPylibs as spg

    Notes
    -----
    '''
    imsize = image[0].shape
    n_images = len(image)
    if org_centers is None:
        org_centers = np.tile(np.array([0., 0.], dtype=np.int16), (n_images, 1))
  
    ############################
    #Normalize images (using a box 100x100 in the central image)
    ############################

    if normalize == True:
        norma = np.mean(image[0][imsize[0]//2-100:imsize[0]//2 +
                            100, imsize[0]//2-100:imsize[0]//2+100])
        if verbose == True:
            print('Normalization constant: ', norma, '[calculated with first image assumed to be central one]')

        for i in range(n_images):
            image[i] = image[i]/norma

    ############################
    #CALCULATE THE MASK GRADIENT FOR EACH IMAGE
    ############################

    binmask = []
    image_dummy, threshold = FindEdges(
        image[0], threshold, method=method, dthr=dhtr, verbose=verbose,Otsu=Otsu)
    binmask.append(image_dummy)

    for i in range(1, n_images):
        image_dummy = FindEdges(
            image[i], threshold, method=method, verbose=verbose,Otsu=Otsu)
        binmask.append(image_dummy)

    ############################
    #FIND CENTERS - COARSE SEARCH
    ############################
    #Coarse and fine compressed in one call

    # centers = []
    # radius = []
    # r_width_coarse = (outer_radius - inner_radius)//steps * 2
    # radii_coarse = np.linspace(inner_radius, outer_radius, steps)
    # print('Analizing ', n_images, ' images (coarse search)')

    # for i in range(n_images): 
    # #acc_conv = find_Circles(
    # #    binmask[i], radii_coarse, r_width_coarse, verbose=verbose, full=True)
    #     acc_conv = find_Circles_ida(binmask[i], radii_coarse, r_width_coarse)
    #     center,rad,c,d = votes(acc_conv, radii_coarse)

    #     centers.append(center)
    #     radius.append(rad)
    #     print('Found center: ', centers[i], ' and radius: ', radius[i])
    #     if verbose == True:
    #         fig = plt.figure(frameon=False)
    #         im1 = plt.imshow(binmask[i], cmap=plt.cm.gray, alpha=.5)
    #         circle_fit = bin_annulus(
    #             imsize, radius[i], 1, full=False).astype(float)
    #         dd = np.array(centers[i])
    #         dx = dd[0] - imsize[0]//2
    #         dy = dd[1] - imsize[1]//2
    #         circle_fit = shift(circle_fit, shift=[dx,dy])
    #         im2 = plt.imshow(circle_fit, cmap=plt.cm.gray, alpha=.5)
    #         plt.show()

    # print('Image |   Original  |  Inferred   |   Radius')
    # for i in range(n_images):
    #     print("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
    #         (i, org_centers[i, 0], org_centers[i, 1],
    #         centers[i][1], centers[i][0], radius[i]))

    ############################
    #FIND CENTERS - FINE SEARCH
    ############################

    # centers_fine = []
    # radius_fine = []

    # mean_r = np.mean(radius)
    # print('pp',mean_r)
    # inner_radius = mean_r-20
    # outer_radius = mean_r+20
    # steps = 20
    # r_width_fine = 5
    # radii_fine = np.linspace(inner_radius, outer_radius, steps)
    # print('Analizing ', n_images, ' images (fine case)')

    # for i in range(n_images):
    #     acc_conv = find_Circles_ida(binmask[i], radii_fine, r_width_fine,verbose=False)
    #     center,rad,c,d = votes(acc_conv, radii_fine)
    #     centers_fine.append(center)
    #     radius_fine.append(rad)
    #     print('Found center: ', centers_fine[i],
    #         ' and radius: ', radius_fine[i])
    #     if verbose == True:
    #         fig = plt.figure(frameon=False)
    #         im1 = plt.imshow(binmask[i], cmap=plt.cm.gray, alpha=.5)
    #         circle_fit = bin_annulus(
    #             imsize, radius_fine[i], 1, full=False).astype(float)
    #         dd = np.array(center)
    #         dx = dd[0] - imsize[0]//2
    #         dy = dd[1] - imsize[1]//2
    #         circle_fit = shift(circle_fit, shift=[dx,dy])
    #         im2 = plt.imshow(circle_fit, cmap=plt.cm.gray, alpha=.5)
    #         plt.show()

    # print('Method  |  Image |   Original  |  Inferred   |   Radius')
    # for i in range(n_images):
    #     print(" Coarse  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
    #         (i, org_centers[i, 0], org_centers[i, 1],
    #         centers[i][1], centers[i][0], radius[i]))
    #     print(" Fine    %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
    #         (i, org_centers[i, 0], org_centers[i, 1],
    #         centers_fine[i][1], centers_fine[i][0], radius_fine[i]))
    
    if steps > 0:
        ############################# 
        #FIND CENTERS - COARSE SEARCH
        #############################
        r_width = (outer_radius - inner_radius)//steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))

        printc('from: ',inner_radius,' to: ',outer_radius,' steps: ', steps,' width: ',r_width,color = bcolors.OKGREEN)
        centers, radius = centers_flat(n_images,inner_radius,outer_radius,steps,r_width,binmask,imsize,verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            printc("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                (i, org_centers[i, 0], org_centers[i, 1],
                centers[i][1], centers[i][0], radius[i]),color = bcolors.FAIL)
        ###########################
        #FIND CENTERS - FINE SEARCH
        ###########################
        mean_r = np.int(np.mean(radius))
        inner_radius = mean_r - 32
        outer_radius = mean_r + 32
        steps = 16
        r_width = (outer_radius - inner_radius)//steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))

        printc('from: ',inner_radius,' to: ',outer_radius,' steps: ', steps,' width: ',r_width,color = bcolors.OKGREEN)
        centers, radius = centers_flat(n_images,inner_radius,outer_radius,steps,r_width,binmask,imsize,verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            printc("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                (i, org_centers[i, 0], org_centers[i, 1],
                centers[i][1], centers[i][0], radius[i]),color = bcolors.FAIL)
        ################################
        #FIND CENTERS - VERY FINE SEARCH
        ################################
        mean_r = np.int(np.mean(radius))
        inner_radius = mean_r - 4
        outer_radius = mean_r + 4
        steps = 8
        r_width = (outer_radius - inner_radius)//steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))

        printc('from: ',inner_radius,' to: ',outer_radius,' steps: ', steps,' width: ',r_width,color = bcolors.OKGREEN)
        centers, radius = centers_flat(n_images,inner_radius,outer_radius,steps,r_width,binmask,imsize,verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            printc("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                (i, org_centers[i, 0], org_centers[i, 1],
                centers[i][1], centers[i][0], radius[i]),color = bcolors.FAIL)
    elif steps < 0:
        ##################################
        #FIND CENTERS - FM SEARCH STRATEGY
        ##################################
        r_width = 2
        inner_radius = 128
        outer_radius = 1024
        steps = 32
        r_width = (outer_radius - inner_radius)//steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))
        printc('from: ',inner_radius,' to: ',outer_radius,' steps: ', steps,' width: ',r_width,color = bcolors.OKGREEN)
        centers, radius = centers_flat(n_images,inner_radius,outer_radius,steps,r_width,binmask,imsize,verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            print("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                (i, org_centers[i, 0], org_centers[i, 1],
                centers[i][1], centers[i][0], radius[i]))
        mean_r = np.int(np.mean(radius))
        inner_radius = mean_r - 32
        outer_radius = mean_r + 32
        steps = 16 
        r_width = (outer_radius - inner_radius)//steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))
        printc('from: ',inner_radius,' to: ',outer_radius,' steps: ', steps,' width: ',r_width,color = bcolors.OKGREEN)
        centers, radius = centers_flat(n_images,inner_radius,outer_radius,steps,r_width,binmask,imsize,verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            print("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                (i, org_centers[i, 0], org_centers[i, 1],
                centers[i][1], centers[i][0], radius[i]))
        mean_r = np.int(np.mean(radius))
        inner_radius = mean_r - 8
        outer_radius = mean_r + 8
        steps = 8 
        r_width = (outer_radius - inner_radius)//steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))
        printc('from: ',inner_radius,' to: ',outer_radius,' steps: ', steps,' width: ',r_width,color = bcolors.OKGREEN)
        centers, radius = centers_flat(n_images,inner_radius,outer_radius,steps,r_width,binmask,imsize,verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            print("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                (i, org_centers[i, 0], org_centers[i, 1],
                centers[i][1], centers[i][0], radius[i]))
        mean_r = np.int(np.mean(radius))
        inner_radius = mean_r - 2
        outer_radius = mean_r + 2
        steps = 4
        r_width = (outer_radius - inner_radius)//steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))
        printc('from: ',inner_radius,' to: ',outer_radius,' steps: ', steps,' width: ',r_width,color = bcolors.OKGREEN)
        centers, radius = centers_flat(n_images,inner_radius,outer_radius,steps,r_width,binmask,imsize,verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            print("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                (i, org_centers[i, 0], org_centers[i, 1],
                centers[i][1], centers[i][0], radius[i]))
    else:
        print('NO HOUGH **** WRONG')

    if save == True:
        status = write_shifts('hough_centers.txt', (centers,radius))
        if status != 1:
            print('Error in write_shifts')

    return centers, radius

@timeit  
def fdt_flat_gen_fixpoint(image, rel_centers, radious = 0, thrd=0.05, iter=15, \
    bit_trun = 0,verbose = 0, expand=0,imasize=[2048,2048]):
    '''
    Khun-Lin-Lorantz algorithm ()
    Input: 
    image -> [n_images][y,x]
    rel_centers -> [n_images,2] where [:,0]=dx and [:,1]=dy 
        Displacements are given with respect to image origin (0,0)
    radious = 0 : int
        radious of circular mask. Default = 0. In this case, the code uses the thrd to create the mask
    thrd = 0.05 : float
        threshold above which pixels are valid. Default = 0.05 (assuming image is normalized to one)
    iter = 15 : int
        maximum number of iterations in the kll algorithm.  Default = 15
    expand = 0 : int
        how much the circular mask is expanded (positive = schrinks the mask) in %.
    verbose = 0 : int
        0 = do nothing, 1 = ...., 2 = ...., 
    bit_trun = 0: int
        Do not touch
    imasize = [2048,2048] : Image size     
    '''

    imsize = image[0].shape
    n_images = len(image)

    ############################
    # set displacements of observed images (A) wth respect Object image (centered) 
    ############################
    xyshifts = np.empty([n_images,2],dtype=int)
    xyshifts[:,0], xyshifts[:,1] = rel_centers[:,0] - imsize[0]//2 , rel_centers[:,1] - imsize[1]//2

    ############################
    # calculate masks
    ############################
    mask = np.zeros([n_images, imsize[0], imsize[1]], dtype=np.int8)

    if radious != 0:   # In case radius of solar disk is provided....
        maskn,dummy = generate_circular_mask([imsize[0]-1, imsize[1]-1],radious - expand,radious - expand)  
        print('Using circular mask')
        for i in range(n_images): 
            mask[i] = shift(maskn, shift=[xyshifts[i,0],xyshifts[i,1]], fill_value = 0)
    else:
        # find pixel coordinates with solar information (> thrd given by user, default = 0.05)
        # This step assumes input data has a mean value of one.
        for i in range(n_images): 
            x, y = np.where(image[i] > thrd)
            mask[i][x, y] = 1

    ############################
    # DO LOG
    ############################

    D = np.zeros([n_images, imsize[0], imsize[1]], dtype=BIT_DEPTH)

    for i in range(n_images):
         D[i,:,:] = (np.log10(image[i])*2**26).astype(BIT_DEPTH)

    # replace NaNs and Infs by 0
    D[np.isneginf(D)] = 0
    D[np.isnan(D)] = 0

    ############################
    # CALCULATE CONSTANT
    ############################

    n = np.zeros([imsize[0], imsize[1]], dtype=BIT_DEPTH)
    sum_image = np.zeros([imsize[0], imsize[1]], dtype=BIT_DEPTH)

    print('Rel centers: ',rel_centers)
    #  for [iq, ir] in itertools.combinations(range(n_images), 2): # overall 36 combinations
    for iq in range(1, n_images):  #loop in iq
        for ir in range(iq):  #loop in ir
        # shift of iq with respect ir
            dx = rel_centers[iq, 0] - rel_centers[ir, 0] 
            dy = rel_centers[iq, 1] - rel_centers[ir, 1] 
            if verbose == 2:
                print('dx,dy',dx,dy,iq,ir)

            t_mask_1 = mask[ir] & shift(mask[iq], [-dx, -dy])
            t_mask_2 = mask[iq] & shift(mask[ir], [dx, dy])
            t_mask = t_mask_1 & t_mask_2 # compound mask only used for mean

            t_image_1 = shift(D[iq], [-dx, -dy])
            t_image_2 = shift(D[ir], [dx, dy])
            aa = (D[iq] - t_image_2) * t_mask_2  #add _2
            bb = (D[ir] - t_image_1) * t_mask_1  #add _1 
            image_pair = aa + bb

            sum_image += image_pair 
            
            n += t_mask_1 # accumulate valid pixels first sumatorio
            n += t_mask_2 # accumulate valid pixels second sumatorio

    K = (sum_image / n).astype(BIT_DEPTH) 
    # replace NaNs and Infs by 0
    K[np.isneginf(K)] = 0
    K[np.isnan(K)] = 0
    if verbose == 1:
        plt.imshow(K,cmap='gray')
        plt.clim(vmax=0.02,vmin=-0.02)
        plt.colorbar()
        plt.show()

    G = np.copy(K)

    for itera in range(iter):
        r_res = np.zeros(imasize, dtype=BIT_DEPTH)
        for iq in range(1,n_images):
            for ir in range(iq):
                # shift of iq with respect ir
                dx = rel_centers[iq, 0] - rel_centers[ir, 0]
                dy = rel_centers[iq, 1] - rel_centers[ir, 1]
                if verbose == 2:
                    print('dx,dy',dx,dy,iq,ir)
                t_mask_1 = mask[ir] & shift(mask[iq], [-dx, -dy])
                t_mask_2 = mask[iq] & shift(mask[ir], [dx, dy])

                t_image_1 = shift(G, [-dx, -dy]) * t_mask_1
                t_image_2 = shift(G, [ dx,  dy]) * t_mask_2
                correction = (t_image_1 + t_image_2) 
                r_res +=  correction 

        G = (K + r_res / n).astype(BIT_DEPTH) 
        # replace NaNs and Infs by 0
        G[np.isneginf(G)] = 0
        G[np.isnan(G)] = 0

        idx = np.where(n > 0)
        s = G[idx]
        #  calculate average of gain table for normalization

        #sm = np.mean(s)
        #sm2 = np.mean(s**2)
        #five_sigma = 5*np.sqrt(sm2-sm*sm)
        #idx2 = np.where(np.abs(s-sm) < five_sigma)
        #sm = np.mean(s[idx2])

        sm = np.mean(s)
        st = np.std(s)
        five_sigma = 5*st
        idx2 = np.where(np.abs(s-sm) < five_sigma)
        sm = np.mean(s[idx2]).astype(BIT_DEPTH)
        
        # G[idx] = G[idx] - sm
        G = G - sm

        print('Iteration: ', itera, five_sigma, '5*rms', sm, ' of ', iter)

    g = np.power(10, G/2**26) #/np.log(10,dtype='d')  + 0.5672334407 - 0.0018501610250685886#0.566 #exp(2.303)
    g[np.isneginf(g)] = 0
    g[np.isnan(g)] = 0
    g = (g*2**26).astype(BIT_DEPTH) 
    if verbose:
        plt.imshow(g, cmap='gray', vmin=0.95*2**26, vmax=1.05*2**26)
        plt.colorbar()
        plt.show()

    return g


def fdt_flat_fixpoint(observations, wavelength, npol, read_shits = 0, shifts = None, verbose = 1,
    expand = 1,thrd = 0,iter = 4, normalize = 1 , disp_method = 'Hough',
    inner_radius = 400, outer_radius = 800, steps = 20,shifts_file = False,imasize = [2048,2048]):

    image = []

    for i in observations:
       image.append(i.image[wavelength*4+npol])
    n_images = len(image)
    ys,xs = image[0].shape

    if read_shits == 1:
        try:
            print('... read user input shifts_file ...')
            centers = read_shifts(shifts_file+'_cnt_w'+str(wavelength)+'_n'+str(npol)+'.txt')
            radius  = read_shifts(shifts_file+'_rad_w'+str(wavelength)+'_n'+str(npol)+'.txt')
            for i in range(n_images):
                print('Image',i,'c: ',centers[i,0],',',centers[i,1],' rad: ', radius[i])
        except Exception:
            print("Unable to open fits file: {}",shifts_file+'_cnt_w'+str(wavelength)+'_n'+str(npol)+'.txt')        
    elif read_shits == 2:
        print('... shifts provided by user ...')
        centers = shifts[0]
        print(centers, '... read ')
        radius = shifts[1]
        print(radius, '... read ')
    else:
        print('... calculating shifts ...')

        if disp_method == 'Hough':

            centers, radius = do_hough(image, inner_radius, outer_radius, steps,verbose=False,threshold = 0.05)
        
            if shifts_file:
                _ = write_shifts(shifts_file+'_cnt_w'+str(wavelength)+'_n'+str(npol)+'.txt', centers)
                _ = write_shifts(shifts_file+'_rad_w'+str(wavelength)+'_n'+str(npol)+'.txt', radius )

        elif disp_method == 'FFT':
            print('TB checked. Input par "expand" should be negative number representing solar disk')
            image_dummy = np.zeros((n_images,ys,xs))
            for i in range(n_images):
                image_dummy[i,:,:] = image[i]
            s_y,s_x,_ = PHI_shifts_FFT(image_dummy,prec=5,verbose=False,norma=False,coarse_prec = 150)
            centers = np.zeros((n_images,2))
            radius = np.zeros((n_images))
            radius[:] = -expand
            expand = 5
            centers[:,0] = -s_x + xs//2
            centers[:,1] = -s_y + ys//2        
            for i in range(n_images):
                print('Image',i,'c: ',centers[i,0],',',centers[i,1],' rad: ', radius[i],'*')

        elif disp_method == 'circle':
            centers = np.zeros((n_images,2))
            radius = np.zeros((n_images))
            for i in range(n_images):
                centers[i,1],centers[i,0],radius[i] = find_center(image[i],sjump = 4,njumps = 100,threshold = 0.8)
                print, 'Image',i,'c: ',centers[i,0],',',centers[i,1],' rad: ', radius[i]
        else:
            pass
    #make sure we have integer numpy numbers in the centers 
    centers = np.array(centers).astype(int)
    mean_radii = np.mean(radius)

    for i in range(n_images):
        norma = np.mean(image[i][centers[i,1]-100:centers[i,1]+100,centers[i,0]-100:centers[i,0]+100])
        if normalize == 1:
            image[i] = image[i]/norma
            print('Normalization: ', norma)
        else:
            norma = 0
        pass

    if thrd != 0:
        gain = fdt_flat_gen_fixpoint(image, centers,iter=iter,thrd=thrd,verbose = verbose,imasize = imasize)
    else:
        gain = fdt_flat_gen_fixpoint(image, centers,iter=iter,radious=mean_radii,expand=expand,verbose = verbose, imasize = imasize)

        return gain, norma

def fdt_flat_testrun_fixpoint():
    '''
    Just for local test run in a folder at same level of SPGlib
    '''

    dir = '/Users/orozco/Dropbox_folder/SoPHI/PHI-COMMISSIONING/software-and-images/RSW1/Add-data/'
    files = ['solo_L0_phi-fdt-ilam_20200618T035946_V202007101227C_0066180100.fits',
      'solo_L0_phi-fdt-ilam_20200618T040546_V202007101223C_0066180125.fits',
      'solo_L0_phi-fdt-ilam_20200618T041146_V202007101231C_0066180150.fits',
      'solo_L0_phi-fdt-ilam_20200618T041746_V202009211029_0066180175_scorr.fits',
      'solo_L0_phi-fdt-ilam_20200618T042346_V202009211027_0066180200_scorr.fits',
      'solo_L0_phi-fdt-ilam_20200618T043004_V202011101020_0066180225_scorr.fits',
      'solo_L0_phi-fdt-ilam_20200618T043546_V202009211031_0066180250_scorr.fits',
      'solo_L0_phi-fdt-ilam_20200618T044146_V202009211031_0066180275_scorr.fits',
      'solo_L0_phi-fdt-ilam_20200618T044804_V202009291424C_0066180300.fits']
    files = [dir + s for s in files]

    dir = '/Users/orozco/Dropbox_folder/Python/VS-GitHub/SPGPylibs_tests/June-2020-flats/'
    files = list_fits(dir)

    dark_file = '/Users/orozco/Dropbox_folder/SoPHI/PHI-COMMISSIONING/software-and-images/RSW1/solo_L0_phi-fdt-ilam_20200618T000547_V202006221044C_0066181001_dark.fits'

    #create PHI image objects
    observations = []

    for i in files:
        observations.append(phidata(file=i))

    dark = phidata(file = dark_file)

    phidata.dark_offset = 0.98
    for i in observations:
       i.load()
       i.apply_dark(dark,verbose=True)

    wavelength = 0 
    npol = 0

    allgain = []
    norma = np.zeros((24))
    for wavelength in range(6):
      for npol in range(4):
        print(wavelength,npol,'................')
        gain, norma_out = fdt_flat_fixpoint(observations, wavelength, npol,read_shits = False, 
            shifts_file = 'shifts/shifts', correct_ghost = 0 , expand = 10, normalize = 0, 
            iter = 3,verbose=True)#,steps = -1)
        
        allgain.append(gain)
        norma[wavelength*4+npol] = norma_out

    with pyfits.open(files[0]) as hdu_list:
        hdu_list[0].data = allgain
        hdu_list.writeto('flat_fixpoint.fits', clobber=True)
