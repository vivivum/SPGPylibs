#=============================================================================
# Project: SoPHI
# File:    phi_fits.py
# Author:  David Orozco Suárez (orozco@iaa.es)
# Contributors: 
#-----------------------------------------------------------------------------
# Description: programs for accesing data and fits files
#-----------------------------------------------------------------------------

from astropy.io import fits as pyfits
from astropy.io.fits import getheader
import numpy as np
from .phi_utils import find_string

def fits_get(file,info = False,head = 0,scaling = False):
    '''helper function to load FITS data set
    if only file data is given, return data (16bits) + header: d,h = fits_get(file)
    if info = True prints fits info (headers structure)
    if head is provided, returns data and selected header
    if scaling = True, return the file scaling
        The scaling will be a list of two tuple elelemts
        scaling[0] ->  IMGFMT_16_0_S scaling 
            with scaling[0][0] = bool: True if present and scaling[0][1] the actual scaling 
        scaling[1] ->  IMGFMT_24_8 scaling
            with scaling[0][0] = bool: True if present and scaling[0][1] the actual scaling 
    '''
    if info == True:
        try:
            return pyfits.info(file)
        except Exception:
            print("Unable to open fits file: {}",file)        
            raise
    if scaling == True:
        index = 1
        scaling = {"Present": [False,True], "scaling": [0,0]}
        while True:
            try:
                dummy_head = getheader(file,index)
            except Exception:
                print("Unable to open fits file: {}",file)        
                raise
            if dummy_head['EXTNAME'] == 'PHI_FITS_imageSummary':
                with pyfits.open(file) as hdu_list:
                    header_data = hdu_list[index].data
                    #case 1 if that there is only ONE scaling (untouched data)
                    if len(header_data) == 1:
                        scaling["Present"][0] = False
                        scaling["Present"][1] = True
                        scaling["scaling"][0] = 0.
                        scaling["scaling"][1] = float(header_data[0][12])
                        return scaling
                    #case 2 if that there is more than TWO scaling data
                    if len(header_data) > 2:
                        #check the first one from below and if it is IMGFMT_16_0_S store it and continue
                        if header_data[-1][3] == 'IMGFMT_16_0_S':
                            scaling["Present"][0] = True
                            scaling["Present"][1] = True
                            scaling["scaling"][0] = float(header_data[-1][12])
                            scaling["scaling"][1] = float(header_data[-3][12])
                            return scaling
                        if header_data[-1][3] == 'IMGFMT_24_8':
                            scaling["Present"][0] = False
                            scaling["Present"][1] = True
                            scaling["scaling"][0] = 0.
                            scaling["scaling"][1] = float(header_data[-1][12])
                            return scaling
            index += 1
    with pyfits.open(file) as hdu_list:
        if head != 0:
            return hdu_list[head].data , hdu_list[head].header 
        else:
            return hdu_list[head].data.astype(np.dtype('float32')) , hdu_list[head].header

def fits_get_fpatimes(file,offset = None):
    '''
    for the moment provide the time corresponding to first exposure/wavelength 
    '''
    from datetime import datetime,timedelta
    fpa_head = 4
    # try:  
    dummy_head = getheader(file,0)
    NAXIS3 = dummy_head['NAXIS3']
    dummy_head = getheader(file,fpa_head)
    if dummy_head['EXTNAME'] != 'PHI_FITS_FPA_settings':
        print('ERROR .... fpa_head is not number 4')
        return
    adq_times = np.zeros((NAXIS3),dtype=float)
    with pyfits.open(file) as hdu_list:
        dat_fpa = hdu_list[fpa_head].data
        init = datetime.strptime(dat_fpa[0][1], '%Y-%m-%dT%H:%M:%S.%f')
        for i in range(NAXIS3):
            the_time = datetime.strptime(dat_fpa[i][1], '%Y-%m-%dT%H:%M:%S.%f')
            adq_times[i] = float((the_time - init).total_seconds()) 
        if offset:
            diff = init - offset
            print('diff: ', adq_times)
            print('diff: ', float(diff.total_seconds()) )
            adq_times = adq_times + float(diff.total_seconds()) 
    return adq_times,init
    # except Exception:
    #     print("Unable to open fits file: {}",file,' Other resons may happen (bug)')        
    #     return None

def list_fits(path = './', contain = None,remove_dir = False):
    '''Find all fits in directory.'''
    import os
    list_of_files = []
    for (dirpath, dirnames, filenames) in os.walk(path):
        for filename in filenames:
            if filename.endswith('.fits'):
                if contain != None:
                    _,exist = find_string(filename,contain)
                    if exist != -1:
                        print(exist)
                        if remove_dir:
                            list_of_files.append(filename)
                        else:
                            list_of_files.append(os.path.join(dirpath, filename)) 
                else:
                    list_of_files.append(os.path.join(dirpath, filename)) 
    return list_of_files

def fits_get_sampling(file,verbose = False):
    '''
    wave_axis,voltagesData,tunning_constant,cpos = fits_get_sampling(file)
    No S/C velocity corrected!!!
    cpos = 0 if continuum is at first wavelength and = 5 if continuum is at the end
    '''
    print('-- Obtaining voltages......')
    fg_head = 3
    try:
        with pyfits.open(file) as hdu_list:
            header = hdu_list[fg_head].data
            j = 0
            dummy = 0
            voltagesData = np.zeros((6))
            tunning_constant = 0.0
            ref_wavelength = 0.0
            for v in header:
                if (j < 6):
                    if tunning_constant == 0:
                        tunning_constant = float(v[4])/1e9
                    if ref_wavelength == 0:
                        ref_wavelength = float(v[5])/1e3
                    if np.abs(np.abs(v[2]) - np.abs(dummy)) > 5:
                        voltagesData[j] = float(v[2])
                        dummy = voltagesData[j] 
                        j += 1
        if verbose:
            print('     Voltages: ',voltagesData)
        d1 = voltagesData[0] - voltagesData[1]
        d2 = voltagesData[4] - voltagesData[5]
        if np.abs(d1) > np.abs(d2):
            cpos = 0
        else:
            cpos = 5
        if verbose:
            print('Continuum position at wave: ', cpos)
        wave_axis = voltagesData*tunning_constant + ref_wavelength  #6173.3356
        if verbose:
            print('     Data wave axis [mA]: ',wave_axis)

        return wave_axis,voltagesData,tunning_constant,cpos

    except Exception:
        print("Unable to open fits file: {}",file)        

def fits_get_part(file, wave, npol):
    '''helper function to load FITS data set
    wave and npol reref to wave and npol (wow!)
    '''
    try:
    #        with pyfits.open(file, do_not_scale_image_data=True, memmap=True, mode='denywrite') as hdu_list:
        with pyfits.open(file, mode='denywrite') as hdu_list:
            hdu_list.info()
            data = hdu_list[0].data.astype(np.dtype('float32')) 
            image = data[wave*4+npol,:,:]
        return image
    except Exception:
        print("Unable to open fits file: {}",file)        

        return None

def read_shifts(file):
    print('Reading ', file)
    content = np.loadtxt(file,dtype=np.int_)
    # with open(file, 'r') as reader:
    #     content = reader.readlines()
    return content

def write_shifts(file, content):
    print('writing ', file)
    np.savetxt(file, content, fmt='%5.0d')
    #     with open(file, 'w') as writer:
    #         formatted = [[format(v) for v in r] for r in content]
    #         writer.write(str(formatted))
    return 1

def set_level(file,what,towhat):
    _,exist = find_string(file,what)
    if exist != -1:
        return file.replace(what,towhat)        