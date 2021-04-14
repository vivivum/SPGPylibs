from astropy.io import fits as pyfits
from astropy.io.fits import getheader

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


