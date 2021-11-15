#=============================================================================
# Project: SoPHI
# File:    phi_reg.py
# Author:  David Orozco Suárez (orozco@iaa.es)
# Contributors: 
#-----------------------------------------------------------------------------
# Description: 
#-----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from photutils import CircularAperture
from scipy.optimize import curve_fit
from .phi_gen import *
from .tools import * 
 
def sampling(N):
    """
    This function creates a grid of points with NxN dimensions for calling the
    Zernike polinomials.
    Output:
        X,Y: X and Y meshgrid of the detector
    """
    if N%2 != 0:
        print('Number of pixels must be an even integer!')
        return
    x=np.linspace(-N/2,N/2,N)
    y=np.copy(x)
    X,Y=np.meshgrid(x,y)
    return X,Y

def aperture(X,Y,N,R):
    """
    This function calculates a simple aperture function that is 1 within
    a circle of radius R, takes and intermediate value between 0
    and 1 in the edge and 0 otherwise. The values in the edges are calculated
    according to the percentage of area corresponding to the intersection of the
    physical aperture and the edge pixels.
    http://photutils.readthedocs.io/en/stable/aperture.html
    Input:
        X,Y: meshgrid with the coordinates of the detector ('sampling.py')
        R: radius (in pixel units) of the mask
    Output:
        A: 2D array with 0s and 1s
    """
    A=CircularAperture((N/2,N/2),r=R) #Circular mask (1s in and 0s out)
    A=A.to_mask(method='exact') #Mask with exact value in edge pixels
    A=A.to_image(shape=(N,N)) #Conversion from mask to image
    return A

@timeit
def PHI_shifts_FFT(data,norma=True,prec=100,verbose=False,coarse_prec = 1.5,sequential = False):

    """
    At least two images should be provided!
    s_y, s_x, simage = PHI_shifts_FFT(image_cropped,prec=500,verbose=True,norma=False)
    (row_shift, column_shift) deficed as  center = center + (y,x) 
    """

    #Normalization for each image
    sz,sy,sx = data.shape
    f=np.copy(data)
    if norma == True:
        norm=np.zeros(sz)
        for i in range(sz):
            norm[i]=np.mean(data[i,:,:])
            f[i,:,:]=data[i,:,:]/norm[i]

    #Frequency cut
    wvl=617.3e-9
    D = 0.14  #HRT
    foc = 4.125 #HRT
    fnum = foc / D
    fnum = 33.1 #FDT fnum
    nuc=1/(wvl*fnum) #Critical frequency (1/m)
    N=sx #Number of pixels per row/column (max. 2048)
    deltax = 10e-6 #Pixel size
    deltanu=1/(N*deltax)
    R=(1/2)*nuc/deltanu
    nuc=2*R#Max. frequency [pix]

    #Mask
    X,Y = sampling(N)
    mask = aperture(X,Y,N,R)

    #Fourier transform
    f0=f[0,:,:]
    #pf.movie(f0-f,'test.mp4',resol=1028,axis=0,fps=5,cbar='yes',cmap='seismic')
    F=np.fft.fft2(f0)

    #Masking
    F=np.fft.fftshift(F)
    F*=mask
    F=np.fft.ifftshift(F)

    #FJBM algorithm
    kappa=prec
    n_out=np.ceil(coarse_prec*2.*kappa)
    dftshift=np.fix(n_out/2)
    nr,nc=f0.shape
    Nr=np.fft.ifftshift(np.arange(-np.fix(nr/2),np.ceil(nr/2)))
    Nc=np.fft.ifftshift(np.arange(-np.fix(nc/2),np.ceil(nc/2)))
    kernc=np.exp((-1j*2*np.pi/(nc*kappa))*np.outer(\
    np.fft.ifftshift(np.arange(0,nc).T-np.floor(nc/2)),np.arange(0,n_out)-dftshift))
    kernr=np.exp((-1j*2*np.pi/(nr*kappa))*np.outer(\
    np.arange(0,n_out)-dftshift,np.fft.ifftshift(np.arange(0,nr).T-np.floor(nr/2))))

    row_shift=np.zeros(sz)
    col_shift=np.zeros(sz)
    shifted_image = np.zeros_like(data)

    if sequential == False:
        for i in np.arange(1,sz):
            g=f[i,:,:]
            G=np.fft.fft2(g)
            #Masking
            G=np.fft.fftshift(G)
            G*=mask
            G=np.fft.ifftshift(G)

            error,row_shift[i],col_shift[i],Gshift=dft_fjbm(F,G,kappa,dftshift,nr,\
            nr,Nr,Nc,kernr,kernc)
            shifted_image[i,:,:] = np.real(np.fft.ifft2(Gshift)) 
    if sequential == True:
        print('No fastidies')
        for i in np.arange(1,sz):
            g=f[i,:,:]
            G=np.fft.fft2(g)
            #Masking
            G=np.fft.fftshift(G)
            G*=mask
            G=np.fft.ifftshift(G)

            error,row_shift[i],col_shift[i],Gshift=dft_fjbm(F,G,kappa,dftshift,nr,\
            nr,Nr,Nc,kernr,kernc)
            shifted_image[i,:,:] = np.real(np.fft.ifft2(Gshift)) 
            F = np.copy(G) #Sequencial
            row_shift[i] = row_shift[i] + row_shift[i-1]
            col_shift[i] = col_shift[i] + col_shift[i-1]
 
    if verbose == True:
        plt.subplots()
        plt.plot(-row_shift,label='y shift',marker='o')
        plt.plot(-col_shift,label='x shift',marker='o')
        plt.ylabel('Shift [pix]')
        plt.legend()
        plt.show()

    return row_shift,col_shift,shifted_image

@timeit
def PHI_shifts_CC(data,norma=True,prec=100,verbose=False,sequential = False,ververbose=False):

    """
    cross-correlation of images
    """

    #Normalization for each image
    sz,sy,sx = data.shape
    apod_win = apod(sx,sy,0)
    
    f=np.copy(data)
    if norma == True:
        norm=np.zeros(sz,dtype=np.dtype('float32'))
        for i in range(sz):
            norm[i]=np.mean(data[i,:,:])
            f[i,:,:]=data[i,:,:]/norm[i]

    x_shift=np.zeros(sz,dtype=np.dtype('float32'))
    y_shift=np.zeros(sz,dtype=np.dtype('float32'))

    #Fourier transform reference image
    f0=f[0,:,:]
    F=fft.fft2(f0*apod_win)
    off_set = abs(F).max() * 1e-15

    #Masking
    # F=np.fft.fftshift(F)
    # F*=mask
    # F=np.fft.ifftshift(F)

    if sequential == False:
        for i in np.arange(1,sz):
            g=f[i,:,:]
            G=fft.fft2(g*apod_win)
            #f0, f1 = [fft.fft2(arr) for arr in (im0, im1)]

            #Masking
            # G=np.fft.fftshift(G)
            # G*=mask
            # G=np.fft.ifftshift(G)

            #cross-correlation
            cc = np.abs(fft.ifft2((F * G.conjugate()) / (1. + off_set))) #not notmalized!!!!
            cc_o = fft.fftshift(cc)
            
            max_position = np.argmax(cc_o)
            positions = np.array(list(np.unravel_index(max_position, cc_o.shape)))
            x_shift[i] = positions[0] - sx / 2
            y_shift[i] = positions[1] - sy / 2
            # print(positions)
            # print(x_shift[i],y_shift[i])
            # nxy = 3
            # ZZ = cc_o[positions[0]-nxy:positions[0]+nxy+1,positions[1]-nxy:positions[1]+nxy+1]
            # zoomed = ndii.zoom(ZZ, 10)
            # max_position = np.argmax(zoomed)
            # positions_zoomed = np.array(list(np.unravel_index(max_position, cc_o.shape)))
            # print(positions_zoomed)
            # x_shift[i] = positions_zoomed[0]/10.+positions[0]-nxy - sx / 2
            # y_shift[i] = positions_zoomed[1]/10.+positions[1]-nxy - sy / 2
            # print(x_shift[i],y_shift[i])
            #-----------------------------------

            if ververbose == True:
                plt.imshow(np.log10(np.abs(cc_o)))
                plt.show()

                #3D
                nxy = 3
                X,Y = np.meshgrid(np.arange(-nxy, nxy+1, 1), np.arange(-nxy, nxy+1, 1))
                ZZ = cc_o[positions[0]-nxy:positions[0]+nxy+1,positions[1]-nxy:positions[1]+nxy+1]

                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.plot_surface(X, Y, ZZ, cmap='plasma')
                ax.set_zlim(np.min(ZZ)-2,np.max(ZZ)+2)
                plt.show()
    if sequential == True:
        for i in np.arange(1,sz):
            g=f[i,:,:]
            G=fft.fft2(g*apod_win)
            #f0, f1 = [fft.fft2(arr) for arr in (im0, im1)]

            #Masking
            # G=np.fft.fftshift(G)
            # G*=mask
            # G=np.fft.ifftshift(G)

            #cross-correlation
            cc = np.abs(fft.ifft2((F * G.conjugate()) / (1. + off_set))) #not notmalized!!!!
            cc_o = fft.fftshift(cc)
            
            max_position = np.argmax(cc_o)
            positions = np.array(list(np.unravel_index(max_position, cc_o.shape)))
            x_shift[i] = positions[0] - sx / 2
            y_shift[i] = positions[1] - sy / 2

            F = np.copy(G) #Sequencial
            x_shift[i] = x_shift[i] + x_shift[i-1]
            y_shift[i] = y_shift[i] + y_shift[i-1]
            # print(positions)
            # print(x_shift[i],y_shift[i])
            # nxy = 3
            # ZZ = cc_o[positions[0]-nxy:positions[0]+nxy+1,positions[1]-nxy:positions[1]+nxy+1]
            # zoomed = ndii.zoom(ZZ, 10)
            # max_position = np.argmax(zoomed)
            # positions_zoomed = np.array(list(np.unravel_index(max_position, cc_o.shape)))
            # print(positions_zoomed)
            # x_shift[i] = positions_zoomed[0]/10.+positions[0]-nxy - sx / 2
            # y_shift[i] = positions_zoomed[1]/10.+positions[1]-nxy - sy / 2
            # print(x_shift[i],y_shift[i])
            #-----------------------------------

            if ververbose == True:
                plt.imshow(np.abs(cc_o))
                plt.colorbar()
                plt.show()

                #3D
                nxy = 3
                X,Y = np.meshgrid(np.arange(-nxy, nxy+1, 1), np.arange(-nxy, nxy+1, 1))
                ZZ = cc_o[positions[0]-nxy:positions[0]+nxy+1,positions[1]-nxy:positions[1]+nxy+1]

                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.plot_surface(X, Y, ZZ, cmap='plasma')
                ax.set_zlim(np.min(ZZ)-2,np.max(ZZ)+2)
                plt.show()

    if verbose == True:
        plt.plot(-x_shift,label='y shift',marker='o')
        plt.plot(-y_shift,label='x shift',marker='o')
        plt.ylabel('Shift [pix]')
        plt.xlabel('# of frame')
        plt.legend()
        plt.show()

    return x_shift,y_shift#,shifted_image

def dft_fjbm(F,G,kappa,dftshift,nr,nc,Nr,Nc,kernr,kernc):
    """
    Calculates the shift between a couple of images 'f' and 'g' with subpixel
    accuracy by calculating the IFT with the matrix multiplication tecnique.
    Shifts between images must be kept below 1.5 'dftshift' for the algorithm
    to work.
    Input: 
        F,G: ffts of images 'f' and 'g' without applying any fftshift
        kappa: inverse of subpixel precision (kappa=20 > 0.005 pixel precision)
    Output:
    """
    #DFT by matrix multiplication
    M=F*np.conj(G) #Cross-correlation
    CC=kernr @ M @ kernc
    CCabs=np.abs(CC)
    ind = np.unravel_index(np.argmax(CCabs, axis=None), CCabs.shape)
    CCmax=CC[ind]
    rloc,cloc=ind-dftshift
    row_shift=-rloc/kappa
    col_shift=-cloc/kappa
    rg00=np.sum(np.abs(F)**2)
    rf00=np.sum(np.abs(G)**2)
    error=np.sqrt(1-np.abs(CCmax)**2/(rg00*rf00))
    Nc,Nr=np.meshgrid(Nc,Nr)

    Gshift=G*np.exp(1j*2*np.pi*(-row_shift*Nr/nr-col_shift*Nc/nc)) 
    return error,row_shift,col_shift,Gshift

def shift_subp(im, shift=[0, 0], fill_value=0):
    '''define shift operator (subpixel)
        Input is y and x shifts (defined negative towards (0,0)
        new center = center + (x,y) 
        Note that image is defined as [sy,sx] so shifts = [sy (rows),sx (columns)]
    '''

    nr,nc = im.shape
    Nr = np.fft.ifftshift(np.arange(-np.fix(nr/2),np.ceil(nr/2)))
    Nc = np.fft.ifftshift(np.arange(-np.fix(nc/2),np.ceil(nc/2)))
    Nc,Nr = np.meshgrid(Nc,Nr)
    G = np.fft.fft2(im)
    Gshift = G*np.exp(1j*2*np.pi*(-shift[0]*Nr/nr-shift[1]*Nc/nc)) 
    Gshift = np.real(np.fft.ifft2(Gshift))
    return Gshift

def gaussian(coor ,height, x0, y0, sigma_x, sigma_y):
    x , y = coor
    f_gauss = height * np.exp(-((x-x0)**2/(2*sigma_x**2) + (y-y0)**2/(2*sigma_y**2)))
    return f_gauss.ravel()

def moments(data):
    """
    Usefull for 2D data
    Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def Gauss(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def Gauss2(x, a_1, x0_1, sigma_1,a_2, x0_2, sigma_2):
    return Gauss(x, a_1, x0_1, sigma_1)+Gauss(x, a_2, x0_2, sigma_2)

def fitgaussian(coor,data):
    print(type(coor))
    params = moments(data)
    popt, pcov = curve_fit(gaussian,coor, data, p0 = params)
    return popt,pcov

