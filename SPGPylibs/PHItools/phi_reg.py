import numpy as np
import matplotlib.pyplot as plt
from photutils import CircularAperture
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
        sys.exit()
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
