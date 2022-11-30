#=============================================================================
# Project: SoPHI
# File:    phi_gen.py
# Author:  David Orozco Suárez (orozco@iaa.es)
# Contributors: 
#-----------------------------------------------------------------------------
# Description: 
#-----------------------------------------------------------------------------

from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.signal import fftconvolve, savgol_filter
from scipy.signal.windows import tukey

import SPGPylibs.GENtools.plot_lib as plib


# __all__ = ['bar', 'baz']

def shift(matrix: np.ndarray, shift: list = None, fill_value: int = 0) -> np.ndarray:
    '''Shift operator
    Shift an image in 2D naively as in SOLO-PHI instrument.
    Faster and more efficient methods can be used in normal CPU.
    Input is a vector shift=[x,y] of x and y displacement
    +x -> positive; +y -> positive 
    fill_value = float. 
    This method does not have any boundary condition.
    '''
    try:
        dimy, dimx = matrix.shape
    except:
        raise ValueError("Input is not 2D matrix") 
    
    try:
        nx = shift[1]
        ny = shift[0]
    except:
        raise ValueError("Provided shift not in rigth format 'shift=[0, 0]' of not present") 

    e = np.empty_like(matrix)
    if nx > 0:
        e[:nx, :] = fill_value
        e[nx:, :] = matrix[:-nx, :]
    elif nx < 0:
        e[nx:, :] = fill_value
        e[:nx, :] = matrix[-nx:, :]
    else:
        e = matrix

    s = np.empty_like(matrix)
    if ny > 0:
        s[:, :ny] = fill_value
        s[:, ny:] = e[:, :-ny]
    elif ny < 0:
        s[:, ny:] = fill_value
        s[:, :ny] = e[:, -ny:]
    else:
        s = e

    return s

def generate_circular_mask(size, radius, r_width):
    """
    Create a circle mask of size = [dy,dx] with radius and r_width width
    """

    grids = np.mgrid[-size[0]//2:size[0]//2 + 1, -size[1]//2:size[1]//2+1]

    # [j][i] = r^2
    r2 = grids[0]**2 + grids[1]**2
    theta = np.arctan2(grids[1], grids[0])
    
    # get boolean value for inclusion in the circle
    outer_circle = r2 <= radius**2
    inner_circle = r2 < (radius - r_width)**2
    
    # back to integers
    outer_circle.dtype = inner_circle.dtype = np.int8
    annulus = outer_circle - inner_circle

    coords = np.where(annulus == 1)

    return annulus, coords

def gradient(xs, method='simply'):
    '''calculate gradient of real image using:
    Finite differences method="simply" (default)
    Prewitt method="prewitt" (central differences)
    sobel method="prewittsmooth"
         [-1 0 1]          [ 1  1  1]
    xdim [-1 0 1]    ydim  [ 0  0  0]
         [-1 0 1]          [-1 -1 -1]

         [-1 0 1]          [ 1  2  1]
    xdim [-2 0 2]    ydim  [ 0  0  0]
         [-1 0 1]          [-1 -2 -1]
    '''
    method = method.lower()
    if method == 'simply':
        print('Using simply for edges')
        sx = (xs - shift(xs, shift=[1, 0], fill_value=0))
        sy = (xs - shift(xs, shift=[0, 1], fill_value=0))
        return sx**2+sy**2
    elif method == 'prewitt':
        print('Using prewitt for edges')
        lxx = (shift(xs, shift=[-1, 0]) - shift(xs, shift=[+1, 0]))/2
        lyy = (shift(xs, shift=[0, -1]) - shift(xs, shift=[0, +1]))/2
        return lxx**2+lyy**2
    elif method == 'prewittsmooth':
        print('Using prewittsmooth for edges')
        lxx = np.zeros_like(xs)
        lyy = np.zeros_like(xs)
        for i in range(3):
            lxx = lxx + \
                (shift(xs, shift=[-1, -1+i]) - shift(xs, shift=[+1, -1+i]))
            lyy = lyy + \
                (shift(xs, shift=[-1+i, -1]) - shift(xs, shift=[-1+i, +1]))
        lxx = lxx/6
        lyy = lyy/6
        return lxx**2+lyy**2
    else:
        print('wrong method')
        return 0

def threshold_otsu(image, nbins=256):
    """Return threshold value based on Otsu's method.

    Parameters
    ----------
    image : (N, M) ndarray
        Grayscale input image.
    nbins : int, optional
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.

    Returns
    -------
    threshold : float
        Upper threshold value. All pixels with an intensity higher than
        this value are assumed to be foreground.

    Raises
    ------
    ValueError
         If `image` only contains a single grayscale value.

    References
    ----------
    .. [1] Wikipedia, http://en.wikipedia.org/wiki/Otsu's_Method

    Notes
    -----
    The input image must be grayscale.
    """
    if len(image.shape) > 2 and image.shape[-1] in (3, 4):
        raise ValueError("threshold_otsu is expected to work correctly only for "
                         "grayscale images; image shape {0} looks like an RGB image")

    # Check if the image is multi-colored or not
    if image.min() == image.max():
        raise ValueError("threshold_otsu is expected to work with images "
                         "having more than one color. The input image seems "
                         "to have just one color {0}.".format(image.min()))

    hist, bin_centers = histogram(image.ravel(), nbins)
    hist = hist.astype(float)

    # class probabilities for all possible thresholds
    weight1 = np.cumsum(hist)
    weight2 = np.cumsum(hist[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(hist * bin_centers) / weight1
    mean2 = (np.cumsum((hist * bin_centers)[::-1]) / weight2[::-1])[::-1]

    # Clip ends to align class 1 and class 2 variables:
    # The last value of `weight1`/`mean1` should pair with zero values in
    # `weight2`/`mean2`, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

    idx = np.argmax(variance12)
    return bin_centers[:-1][idx]

def histogram(image: np.ndarray, nbins: int=256):
    """Return histogram of image.

    Unlike `numpy.histogram`, this function returns the centers of bins and
    does not rebin integer arrays. For integer arrays, each integer value has
    its own bin, which improves speed and intensity-resolution.

    The histogram is computed on the flattened image: for color images, the
    function should be used separately on each channel to obtain a histogram
    for each color channel.

    Parameters
    ----------
    image : array
        Input image.
    nbins : int
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.

    Returns
    -------
    hist : array
        The values of the histogram.
    bin_centers : array
        The values at the center of the bins.

    See Also
    --------
    cumulative_distribution

    """
    sh = image.shape
    # For integer types, histogramming with bincount is more efficient.
    if np.issubdtype(image.dtype, np.integer):
        offset = 0
        image_min = np.min(image)
        if image_min < 0:
            offset = image_min
            image_range = np.max(image).astype(np.int64) - image_min
            # get smallest dtype that can hold both minimum and offset maximum
            offset_dtype = np.promote_types(np.min_scalar_type(image_range),
                                            np.min_scalar_type(image_min))
            if image.dtype != offset_dtype:
                # prevent overflow errors when offsetting
                image = image.astype(offset_dtype)
            image = image - offset
        hist = np.bincount(image.ravel())
        bin_centers = np.arange(len(hist)) + offset

        # clip histogram to start with a non-zero bin
        idx = np.nonzero(hist)[0][0]
        return hist[idx:], bin_centers[idx:]
    else:
        hist, bin_edges = np.histogram(image.flat, bins=nbins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
        return hist, bin_centers

def FindEdges(xs, threshold, method='simply', dthr=1, Otsu=None, verbose=False):
    """
    Use a method = [simply, prewitt, prewittsmooth] 
    for calculatin the image gradient and find the edges
    Otsu = True use Otsu to automatically threshold the image
    Threshold is input parameter for thresholding (factor wrt max)
    dthr = 2 or better to loop over different thresholds
    verbose = True for displaying
    Nosotros lo que hacíamos originalmente para binarizar era:

    -          Calcular derivada en X e Y usando el filtro prewit y la convolución
    -          Calcular el valor absoluto de cada derivada por separado
    -          Sumamos las dos imágenes
    -          Le restamos la original
    -          Aplicamos un threshold mayor o igual a 0
    """
    # prewitt filter gradient
    im_grad = gradient(xs.real, method=method)
    im_grad -= im_grad.min()

    # make binary image
    for _ in range(int(dthr)):
        image = np.copy(im_grad)
        image = image > image.max()*threshold
        image.dtype = np.int8
        nonzero = np.array(image.nonzero())
        density = float(nonzero[0].size)/image.size
        print("Signal density:", density*100.)
        if density > 0.15:
            if verbose == True:
                print("Too many points")
            threshold = threshold*2.
        elif density > 0.01 and density < 0.15:
            if verbose == True:
                print("Threshold OK - exit (between 0.15 and 0.01")
            break
        elif density < 0.01:
            if verbose == True:
                print("Too less points")
            threshold = threshold/2.

    if Otsu != None:
        print('Override current method. Using Otsu histogram')
        thresh = threshold_otsu(im_grad)
        image = np.copy(im_grad)
        image = image > image.max*thresh
        image.dtype = np.int8
    #show_one(image2)
    #np.save('test.npy', image)  # npzfile = np.load(outfile)
    #a = np.load('test.npy')
    #print(a)
    #raise SystemExit()
    # imgbin = image >= 1
    # plt.imshow(image)
    # plt.show()
    # show_one(image)
    # show_one(imgbin)
    if verbose == True:
        print('Stop in FindEdges. Close plot window to continue.')
        plib.show_one(image, title='FindEdges thresholded image')
    return image if dthr == 1 else (image, threshold)

def make_circles(radius, r_width):
    """
    Create a circle with radius and r_width width
    """
    grids = np.mgrid[-radius:radius +
                     1, -radius:radius+1]

    # [j][i] = r^2
    kernel_template = grids[0]**2 + grids[1]**2

    # get boolean value for inclusion in the circle
    outer_circle = kernel_template <= radius**2
    inner_circle = kernel_template < (radius - r_width)**2

    # back to integers
    outer_circle.dtype = inner_circle.dtype = np.int8
    return outer_circle - inner_circle

def find_Circles_ida(image, radii, r_width,verbose = False):
    """
    Perfrom a FFT Convolution over all the radii with the given annulus width.
    Smaller annulus width = more precise
    """
    acc = np.zeros((radii.size, image.shape[0], image.shape[1]))

    for i, r in enumerate(radii):
        C = make_circles(r, r_width)
        acc[i, :, :] = fftconvolve(image, C, 'same')
        print('Running fft convolution ', i, ' from a total of ', len(radii),end='\r')
        if verbose == True:
            print(image.shape,radii[i], r_width)
            plib.show_one(C,title='circle*')
            plib.show_one(acc[i], title='convolution')
    return acc

def votes(acc, radii):

    '''
    devielve: (circle_y, circle_x), radius, maxima, max_positions 
    c[0] = x
    c[1] = y   (The other way around of the definition!!!! which would be c[1] = x and c[0] = y) 
    '''
    maxima = []
    max_positions = []
    max_signal = 0
    print("calc: radius |  maxima  | max_position (x,y) | signal")

    #FIXME Hay que mejorarlo. Coger una caja de 5x5 y buscar el máximo en la caja
    for i, r in enumerate(radii):
        max_positions.append(np.unravel_index(acc[i].argmax(), acc[i].shape))
        maxima.append(acc[i].max())
        # use the radius to normalize
        signal = maxima[i]/np.sqrt(float(r))
    #        if signal > max_signal:
    #            max_signal = signal
        if maxima[i] > max_signal:
            max_signal = np.copy(maxima[i])
            (circle_y, circle_x) = max_positions[i]
            radius = np.copy(r)
        print("calc: %8.2f | %8.2f | %s | %8.2f" % (r, maxima[i], (max_positions[i]), signal))
    
    print(f"Last: {radius:8.2f} | {np.max(maxima):8.2f} | {(circle_x, circle_y)}")
    
    # Identify maximum. Note: the values come back as index, row, column
    #    max_index, circle_y, circle_x = np.unravel_index(acc.argmax(), acc.shape)

    return (circle_x, circle_y), radius, maxima, max_positions  # 

def bin_annulus(shape,radius, width, full = False):
    """
    This function creates a anulus mask of radius and width
    radius - width//2 < radius < radius + width//2 + 1
    """
    # try:
    #     rho
    # except NameError:

    rho = circle_grid(shape)  # FOR GENERATING CENTERS GLOVAL VARIABLE

    mask1 = rho <= radius + width//2 + 1
    mask2 = rho >= radius - width//2
    mask1.astype(np.int8)
    mask2.astype(np.int8)
    return mask1 if full ==  True else mask1 == mask2

def circle_grid(shape):
    """
    This function creates a grid of points with NxN dimensions.
    Output:
        X,Y: X and Y meshgrid of the detector
    """
    N = shape[0]
    if N % 2 != 0:
        print('Number of pixels must be an even integer!', N, N % 2)
        raise Exception
    x = np.linspace(-N/2, N/2, N)
    y = np.copy(x)
    X, Y = np.meshgrid(x, y)
    #globals()['rho'] = rho

    return np.sqrt(X**2 + Y**2)

def find_circle_hough(image,inner_radius, outer_radius, steps,method='prewitt',
            dhtr=10,normalize = False,verbose=False,Otsu = None,threshold = 0.15):
    '''
    Do Hough Transform
    '''
    imsize = image.shape

  ############################
  #Normalize images (using a box 100x100 in the central image)
  ############################

    if normalize == True:
        norma = np.mean(image[imsize[0]//2-100:imsize[0]//2 +
                                100, imsize[0]//2-100:imsize[0]//2+100])
        image = image/norma

  ############################
  #CALCULATE THE MASK GRADIENT FOR EACH IMAGE
  ############################

    binmask = []
    #threshold = 0.15
    #central image to determine threshold
    binmask, threshold = FindEdges(
        image, threshold, method=method, dthr=dhtr, verbose=verbose,Otsu=Otsu)

    #show_one(binmask)
    print(threshold)

  ############################
  #FIND CENTERS 
  ############################

    print('Analizing image........')
    r_width = (outer_radius - inner_radius)//steps * 2
    radii = np.linspace(inner_radius, outer_radius, steps)
    print('r_width',r_width,'radii',radii)
    acc_conv = find_Circles_ida(binmask, radii, r_width)
    #acc_conv = find_Circles(binmask, radii, r_width, verbose=verbose, full=True)
    center,radius,c,d = votes(acc_conv, radii)
    print('Found center [y,x]: ', center, ' and radius: ', radius)

    if verbose == True:
        fig = plt.figure(frameon=False)
        im1 = plt.imshow(binmask, cmap=plt.cm.gray, alpha=.5)
        circle_fit = bin_annulus(
            imsize, radius, 1, full=False).astype(float)
        dd = np.array(center)
        dx = dd[0] - imsize[0]//2
        dy = dd[1] - imsize[1]//2
        circle_fit = shift(circle_fit, shift=[dx,dy])
        im2 = plt.imshow(circle_fit, cmap=plt.cm.gray, alpha=.5)
        plt.show()

    return center, radius, threshold

def simple_shift(xs, shift=0, fill_value=0):
    '''define shift operator'''
    e = np.empty_like(xs)
    if shift > 0:
        e[:shift] = fill_value
        e[shift:] = xs[:-shift]
    elif shift < 0:
        e[shift:] = fill_value
        e[:shift] = xs[-shift:]
    else:
        e = xs
    return e

def find_center(im,sjump = 10,njumps = 50,threshold = 0.9):
    ys,xs = im.shape
    jumps = np.linspace(-sjump*njumps//2,sjump*njumps//2,njumps-1)
    rf = np.array([],dtype=float)
    xc = np.array([],dtype=float)
    yc = np.array([],dtype=float)
    for [i, j] in combinations(jumps, 2): # overall 36 combinations
        xi = xs//2 - int(i)#570#1024
        yi = ys//2 - int(j)#610#1024
        xcut = im[:,yi]
        ycut = im[xi,:]
        xcut = savgol_filter(xcut, 5, 3) # window size 51, polynomial order 3 
        ycut = savgol_filter(ycut, 5, 3) # window size 51, polynomial order 3 

    #calculate derivative
        xcut_d = (xcut - simple_shift(xcut, shift = 10))
        ycut_d = (ycut - simple_shift(ycut, shift = 10))
        xcut_d[:5] = 0
        ycut_d[:5] = 0
        xcut_d[-5:] = 0
        ycut_d[-5:] = 0

    #meter condicion de aumentar threshold si hay muchos o pocos/ninguno
        indices_x_max = np.asarray(np.where(xcut_d > xcut_d.max()*threshold)).flatten() 
        indices_x_min = np.asarray(np.where(xcut_d < xcut_d.min()*threshold)).flatten() 
        indices_y_max = np.asarray(np.where(ycut_d > ycut_d.max()*threshold)).flatten() 
        indices_y_min = np.asarray(np.where(ycut_d < ycut_d.min()*threshold)).flatten() 

        x1 = np.mean(indices_x_max*xcut_d[indices_x_max])/np.mean(xcut_d[indices_x_max]) - 5
        x2 = np.mean(indices_x_min*xcut_d[indices_x_min])/np.mean(xcut_d[indices_x_min]) - 5
        y1 = np.mean(indices_y_max*ycut_d[indices_y_max])/np.mean(ycut_d[indices_y_max]) - 5
        y2 = np.mean(indices_y_min*ycut_d[indices_y_min])/np.mean(ycut_d[indices_y_min]) - 5

        x0 = (x1+x2)/2 
        y0 = (y1+y2)/2 
        r1 = np.sqrt((y0-yi)*(y0-yi)+(x0-x1)*(x0-x1))
        r2 = np.sqrt((y0-y1)*(y0-y1)+(x0-xi)*(x0-xi))
        r = (r1+r2)/2

        rf = np.append(rf,r)
        xc = np.append(xc,x0)
        yc = np.append(yc,y0)

    rf_m = np.mean(rf[np.where((rf < (np.median(rf)+1)) & (rf > (np.median(rf)-1)))])
    xc_m = np.mean(xc[np.where((xc < (np.median(xc)+1)) & (xc > (np.median(xc)-1)))])
    yc_m = np.mean(yc[np.where((yc < (np.median(yc)+1)) & (yc > (np.median(yc)-1)))])

    print(xc_m,yc_m,rf_m)

    return xc_m,yc_m,rf_m

def FFTs(f, dir):
    """fft with shifting"""
    what = np.ascontiguousarray(f)

    if dir == 1:
        return np.fftshift(np.fftn(what))
    elif dir == -1:
        return np.fftshift(np.ifftn(np.ifftshift(what)))
    else:
        print('Select direction: 1 -> Direct FFT; -1 -> Inverse FFT')
        quit()
        return 0

def Laplacian(xs):
    "calculate gradient of real image using Laplacian filter"
    lxx = (shift(xs, shift=[1, 0]) - 2*xs + shift(xs, shift=[-1, 0]))/2
    lyy = (shift(xs, shift=[0, 1]) - 2*xs + shift(xs, shift=[0, -1]))/2
    return lxx**2+lyy**2

def rebin(arr, new_shape):
    """Rebin 2D array arr to shape new_shape by averaging."""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

def apod(nx,alpha):
    window = tukey(int(nx),alpha=alpha)
    return 1 - np.outer(window,window)

def expand_mask(mask,pixels = 1,step = 1):

    msk = shift(mask, shift=[-step, -step], fill_value=0) & shift(mask, shift=[0, step], fill_value=0) \
    & shift(mask, shift=[-step, 0], fill_value=0) & shift(mask, shift=[0, -step], fill_value=0) \
    & shift(mask, shift=[step, 0], fill_value=0) & shift(mask, shift=[1, step], fill_value=0) \
    & mask
    if pixels > 1:
        for _ in range(pixels):
            msk = shift(msk, shift=[-step, -step], fill_value=0) & shift(msk, shift=[0, step], fill_value=0) \
            & shift(msk, shift=[-step, 0], fill_value=0) & shift(msk, shift=[0, -step], fill_value=0)\
            & shift(msk, shift=[step, 0], fill_value=0) & shift(msk, shift=[1, step], fill_value=0) \
            & msk

    return msk

def sigma_mask(input,sigma=5,verbose=False,level=0.02):
    """ Input HAS to be a n,ny,nx array becouse we average in n"""
    shape = input.shape
    if len(shape) != 3:
        print('Input should be [n,ny,nx]')
        raise IndexError 
    md = np.mean(input,axis=0)
    md_std = np.std(input,axis=0)
    md_std = gaussian_filter(md_std,sigma=sigma) 
    if verbose:
        plt.imshow(md_std)
        plt.colorbar()
        plt.show()
    ms = np.ones_like(md)
    ms[md_std > level] = 0
    ms = ms.astype(int)
    if verbose:
        plt.imshow(ms.astype(int),interpolation=None,cmap='Greys')
        plt.colorbar()
        plt.show()
    return ms

def n_coord(x,y,d=1):
    w,h = x.shape
    assert w,h == y.shape
    if d ==1:
        print('Normalizing')
        return (2*x-w)/w, (2*y-h)/h
    if d ==-1:
        print('De-Normalizing')
        return (x+1)*w / 2 , (y+1)*h/2
    print('wrong direction')
    return 0,0
def dmodel(r,k1,k2):
    '''r is a map of the grid'''
    return 1 + k1*r + k2*r**2
def cmesh(w,h,**kwargs):
    X,Y = np.meshgrid(np.arange(w),np.arange(h))
    x,y = n_coord(X,Y,**kwargs)
    r = np.sqrt(x**2 + y**2) 
    return x,y,r

def test_distortion():
    # adjust k_1 and k_2 to achieve the required distortion
    k1 = 0.02
    k2 = 0.02
    margin = np.max([k1,k2])*4+1
    w,h = 2048,2048
    wn,hn = int(np.round(w*margin)),int(np.round(h*margin))
    print(w,h)
    x,y,r = cmesh(wn,hn)
    model = dmodel(r,k1,k2)
    #the model is applied like
    x_new = x * model
    y_new = y * model
    x_new2, y_new2 = n_coord(x_new,y_new,d=-1)

    # TODO incomplete program
    # grid = ndimage.map_coordinates(m, [y_new2.ravel(),x_new2.ravel()])

