from astropy.io import fits as pyfits
from astropy.io.fits import getheader
import numpy as np

def fits_get(file,info = False,head = 0,scaling = False):
    '''helper function to load FITS data set
    if only file data is given, return data (16bits) + header: d,h = fits_get(file)
    if info = True prints fits info (headers structure)
    if head is provided, returns data and selected header
    if scaling = True, return the file scaling
    '''
    try:
        if info == True:
            return pyfits.info(file)
        if scaling == True:
            index = 1
            while True:
                dummy_head = getheader(file,index)
                if dummy_head['EXTNAME'] == 'PHI_FITS_imageSummary':
                    with pyfits.open(file) as hdu_list:
                        header_data = hdu_list[index].data
                        if header_data[-1][3] == 'IMGFMT_16_0_S':
                            return float(header_data[-1][12])
                        else:
                            return float(header_data[-2][12])
                index += 1
        with pyfits.open(file) as hdu_list:
            if head != 0:
                return hdu_list[head].data , hdu_list[head].header 
            else:
                return hdu_list[head].data.astype(np.dtype('float32')) , hdu_list[head].header
    except Exception:
        print("Unable to open fits file: {}",file)        
        return None

def fits_get_sampling(file):
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
        print('     Voltages: ',voltagesData)
        d1 = voltagesData[0] - voltagesData[1]
        d2 = voltagesData[4] - voltagesData[5]
        if np.abs(d1) > np.abs(d2):
            cpos = 0
        else:
            cpos = 5
        print('Continuum position at wave: ', cpos)
        wave_axis = voltagesData*tunning_constant + ref_wavelength  #6173.3356
        print('     Data wave axis [mA]: ',wave_axis)

        return wave_axis,voltagesData,tunning_constant,cpos

    except Exception:
        print("Unable to open fits file: {}",file)        
        return None

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
