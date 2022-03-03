#=============================================================================
# Project: SoPHI
# File:    phifdt.py
# Author:  David Orozco Suárez (orozco@iaa.es)
# Contributors: Nestor Albelo (albelo@mps.mpg.de) and Hanna Streker (streckerh@iaa.es)
#-----------------------------------------------------------------------------
# Description: Class pipeline implementation of data reduction pipeline
#              Uses same modules as in the regular pipeline (check versions)
#-----------------------------------------------------------------------------

from math import nan
import json
import pprint
import numpy as np
from matplotlib import pyplot as plt
from .tools import printc,bcolors,timeit
from .phi_gen import * 
from .phi_utils import *
from .phi_fits import *
from .phi_reg import *
from .phifdt_pipe_modules import phi_correct_dark
# from .phifdt_flat import do_hough
from SPGPylibs.GENtools import *

#--------  GLOBALS  --------  
BIT_DEPTH: str = '>i4'   #"images are 32-bit big-endian integer"
#---------------------------

def fixpoint(imagen):
    if imagen.dtype != BIT_DEPTH:
        imagen = imagen.astype(BIT_DEPTH)       
    print('image bit depth ',imagen.dtype)
    return imagen

# def __Json__(self,json_file):
class Json():

    def __init__(self):

        self.json.file =            None 
        self.json.verbose =         None 
        self.json.input_data_dir =  None
        self.json.data_f =          None 
        self.json.shrink_mask =     None 
        self.json.center_method =   None 
        self.json.hough_params =    None 
        self.json.instrument =      None 
        self.json.flat_f =          None 
        self.json.dark_f =          None 
        self.json.dark_c =          None 
        self.json.flat_c =          None 
        self.json.flat_index =      None
        self.json.norm_f =          None 
        self.json.flat_scaling =    None 
        self.json.prefilter_fits =  None 
        self.json.prefilter =       None
        self.json.output_dir =      None 
        self.json.rte =             None 
        self.json.correct_fringes = None
        self.json.correct_ghost =   None
        self.json.putmediantozero = None 
        self.json.debug =           None 
        self.json.loopthis =        None 
        self.json.ItoQUV =          None 
        self.json.VtoQU =           None 
        self.json.realign =         None 
        self.json.ind_wave =        None 
        self.json.nlevel =          None

    def read_json(self,json_file):

        # =========================================================================== #
        # READING CONFIG FILE AND PRINTING
        # =========================================================================== #
        printc('--------------------------------------------------------------',bcolors.OKGREEN)
        printc(' Reading config json file '+json_file,bcolors.OKGREEN)
        with open(json_file) as j:
            CONFIG = json.load(j)

        self.json.file =            json_file
        self.json.verbose =         CONFIG['verbose']
        self.json.input_data_dir =  CONFIG['input_data_dir']
        self.json.data_f =          CONFIG['data_f']
        self.json.shrink_mask =     CONFIG['shrink_mask']
        self.json.center_method =   CONFIG['center_method'] 
        self.json.hough_params =    CONFIG['hough_params']
        self.json.instrument =      CONFIG['instrument']
        self.json.flat_f =          CONFIG['flat_f']
        self.json.dark_f =          CONFIG['dark_f']
        self.json.dark_c =          CONFIG['dark_c']
        self.json.flat_c =          CONFIG['flat_c']
        self.json.flat_index =      CONFIG['flat_index']
        self.json.norm_f =          CONFIG['norm_f']
        self.json.flat_scaling =    CONFIG['flat_scaling']
        self.json.prefilter_fits =  CONFIG['prefilter_fits']
        self.json.prefilter =       CONFIG['prefilter']
        self.json.output_dir =      CONFIG['output_dir']
        self.json.rte =             CONFIG['rte']
        self.json.correct_fringes = CONFIG['correct_fringes']
        self.json.correct_ghost =   CONFIG['correct_ghost']
        self.json.putmediantozero = CONFIG['putmediantozero']
        self.json.debug =           CONFIG['debug']
        self.json.loopthis =        CONFIG['loopthis']
        self.json.ItoQUV =          CONFIG['ItoQUV']
        self.json.VtoQU =           CONFIG['VtoQU']
        self.json.realign =         CONFIG['realign']
        self.json.ind_wave =        CONFIG['ind_wave']
        self.json.nlevel =          CONFIG['nlevel']

        # Prints the nicely formatted dictionary
        pprint.pprint(CONFIG)#, sort_dicts=False)

class phidata():

    dark_offset: float = 1   #

    import reprlib
    r = reprlib.Repr()
    r.maxlist = 4        # max elements displayed for lists
    r.maxstring = 100    # max characters displayed for strings

    def __init__(self, file: str = None):

        self.json         = Json()
        self.file         = file
        self.darkc: bool  = False
        self.flatc: bool  = False
        self.image        = None      #generic (same for flat, dark, etc)
        self.header       = None      #also generic (header info)
        self.imageSummary = None      
        self.imageSummary_head = None
        self.scaling = {"Present": [False,True], "scaling": [0,0], 'bit-depth': None}
        self.DID          = None

        #INIT variables should be the same as those that can be set with the JSON + internal stuff

    @staticmethod    
    def mprint(what,*args, **kw):
        print(phidata.r.repr(what),*args, **kw)
    
    def dscale(self):
        if self.image.dtype != phidata.BIT_DEPTH:
            self.image = self.image.astype(phidata.BIT_DEPTH)       
        printc('image bit depth ',self.image.dtype,color=bcolors.OKGREEN)

    def set_file(self,file):
        self.file = file

    def load_json(self,file):
        self.json.read_json(file)

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



def phifdt_test():
    '''
    Just for local test run in a folder at same level of SPGlib
    '''

    dark = phidata()

    phidata.dark_offset = 0.98
    for i in observations:
       i.load()
       i.apply_dark(dark,verbose=True)
