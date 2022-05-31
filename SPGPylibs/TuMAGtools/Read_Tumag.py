"""
.. py:module:: TuMAGtools.Read_Image
.. module:: Read_Image
        :platform: Unix
        :synopsis: function for reading TuMAG images and headers
.. moduleauthor:: Pablo Santamarina  (SPG - IAA) <psanta@iaa.es>
"""

# ============================ IMPORTS ====================================== #
import matplotlib.pyplot as plt
# TuMAG Libs
import utils as ut
import numpy as np
import movie_tools as mt
# =========================================================================== #

def Read_Image(file):
    """
    This function reads TuMAG image and header information
    
    :param file: image file (can be a list)
    :type file: str

    :return: header, image.
    :rtype: Dict, np.ndarray(float)
        
    """
    
    Head = []
    Image = []

    for i in file:
        # Separate Header and image data 
        H, hl, I = ut.HeadernImageSeparator(file) # Header, header_lenth and image
        # Read header
        Head.append(ut.GetDatafromHeader(H))
        # Read Image
        Image.append(ut.read_image(file, Head, hl))
     
    return Head, Image

# =========================================================================== #

def Read_Image_example():
# Example of execution 

    Image_path = '2022_05_06_11_55_40_251_0_1.img'
    # Image_path = '2022_05_09_14_15_58_733_0_12801.img'
    H, I = Read_Image(Image_path)

    plt.figure()
    im = plt.imshow(I.T, cmap = 'inferno', vmin = 100, vmax = 150)
    plt.colorbar(im)
    plt.show()

def TuMAG_to_fits(Head,Image):

    '''
    convert to fits
    '''
    from astropy.io import fits as pyfits

    # default ones
    # SIMPLE  =                    T / conforms to FITS standard                      
    # BITPIX  =                    8 / array data type                                
    # NAXIS   =                    0 / number of array dimensions                     
    # EXTEND  =                    T                                                  

    hdu = pyfits.ImageHDU()
    hdu.header['BITPIX'] = 16
    hdu.header['NAXIS'] = 2
    Version = 'Fits TuMAG version 0.1'
    hdu.header.append(('Version', Version, 'Level version'))
    hdu.header.append(('NAXIS1', 2048, 'length of data axis 1'))
    hdu.header.append(('NAXIS2', 2048, 'length of data axis 2'))
    hdu.header.append(('BZERO', 65536, 'offset data range') )
    hdu.header.append(('TELESCOP', 'Sunrise3', 'Telescope')) 
    hdu.header.append(('INSTRUME', 'SUSI', 'Instrument') )
    hdu.header.append(('DATE_OBS', '2022-04-12T14:17:05.491', 'UTC Acquisition Time') )
    hdu.header.append(('TIMESYS', 'UTC', 'Time Standard') )
    hdu.header.append(('FILEORIG', 'original.img', 'Original instrument filename') )
    hdu.header.append(('BINNING', 2, 'binning (if thumbnails') )
    hdu.header.append(('BINNING', 2, 'binning (if thumbnails')) 
    hdu.header.append(('ROI_X1', 0, 'ROI START in X') )
    hdu.header.append(('ROI_Y1', 0, 'ROI START in Y') )
    hdu.header.append(('ROI_X2', 2048, 'ROI END in X') )
    hdu.header.append(('ROI_Y2', 2048, 'ROI END in Y') )
    hdu.header.append(('CDELT1', 0.005, 'pixel size in arcsec') )
    hdu.header.append(('CDELT2', 0.005, 'pixel size in arcsec') )
    # XCEN    =                 -nan / (arcsec)                                       
    # YCEN    =                 -nan / (arcsec)                                       
    # COSTHETA=          6.12323e-17 / cosine of heliocentric angle                   
    # GPS_TIME= '1969-12-31T23:59:59.000' / GPS Time                                  
    # GPS_LON =                 -nan / (deg) GPS Longitude of Balloon                 
    # GPS_LAT =                 -nan / (deg) GPS Latitude of Balloon                  
    # GPS_ALT =                 -nan / (m) GPS Altitude of Balloon                    
    # ELEV    =            90.000000 / (deg) Sun elevation                            
    # AZIMUTH =          -180.000000 / (deg) Sun azimuth                              
    # PARANGLE=                 -nan / (deg) Parallactic angle                        
    # P_ANGLE =           -14.315154 / (deg) P                                        
    # SOLAR_B0=            -6.325978 / (deg) B0                                       
    # SOLAR_L0=            72.971164 / (deg) L0                                       
    # EARTH_D =     147502272.446496 / (km) Sun Earth distance                        
    # SOLAR_R0=           973.267987 / (arcsec) Solar radius                          
    # CW_LOOP = 'open'               / CW Loop State                                  
    # PS_STATE= '15 arcsec'          / Pointing System state                          
    # FLAT_MOD= 'undefined'          / PS flatfield mode                              
    # AP_DOOR = 'switch error, opening' / Telescope Aperture Door State               
    # M2_XPOS =                   -1 / (micron) M2_XPos                               
    # M2_YPOS =                   -1 / (micron) M2_YPos                               
    # M2_ZPOS =                   -1 / (micron) M2_ZPos                               
    # M3_POS  =                   -1 / (micron) M3 Position Z                         
    # M4_POS  =                   -1 / (micron) M4 Position Z        

    hdu.data = np.zeros((2048,2048))
    hdu.writeto('hola.fits')





