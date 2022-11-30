"""
Description
============
Main software for calculating PHI/FDT flat fields.

Main program is called ``fdt_flat``

- File:    phifdt_flat.py

- Author:  David Orozco Suárez (orozco@iaa.es)

- Contributors: Nestor Albelo (albelo@mps.mpg.de)

- Description:

    Pipeline implementation for calculating the FDT flat field.
    It includes three different algorithms for flat calculations and to algorithms for
    finding the disc center coordinates: the circular Hough Transform in the frequency domain and
    a specially tailored simple circle fit program.

- Example:

.. literalinclude:: flat_processing_example/flat_processing_example.ipynb
   :language: ipynb

"""
from typing import Tuple

from SPGPylibs.GENtools import *
from .phi_fits import *
from .phi_reg import *
from .phi_utils import *
from .phifdt_pipe_modules import phi_correct_dark, phi_correct_ghost_single, phi_correct_prefilter, phi_correct_ghost

SIGMA_NORMA = 5


def Neckel(r: float) -> np.ndarray:
    """
    This functions provides de Neckel center to limb intensity variation *CLV* at 625 nm

    Provided r (solar radius) the program estimates the heliocentric angle :math:`{\\mu}`

    .. math::

        \mu = \sqrt{ 1 - r^2/r_\sun^2 } = \cos{\\theta}

    The *CLV* is given by:

    .. math::

        I = I_0  ( 1 - \sum_{k=1}^{n} a_k * \\mu^k)

        where

        \\mu = \cos{\\theta}

    :param r: the radius of the solar disk
    :type r: float
    :return: center to limb variation
    :rtype: np.ndarray
    """

    r_vector = np.arange(0, np.max(r), 1)
    mu = np.sqrt((1 - r_vector ** 2 / r ** 2))

    f1 = 0.31414418403067174
    f2 = 1.3540877312885482
    f3 = -1.827014432405401
    f4 = 2.355950448408068
    f5 = -1.6848471461910317
    f6 = 0.48767921486914473

    return f1 + f2 * mu + f3 * mu ** 2 + f4 * mu ** 3 + f5 * mu ** 4 + f6 * mu ** 5


def clv(r: float, r_val: float = None, coefficients: np.ndarray = False) -> np.ndarray:
    """
    This functions provides the center to limb intensity variation *CLV* model for PHI FDT data

    Provided r (solar radius) the program estimates the heliocentric angle :math:`{\\mu}`

    If the user does not provide coefficients, it takes default parameters for PHI FDT data
    coefficients = [1., 0.44984082, 0.11827295, 0.05089393]

    .. math::

        \mu = \sqrt{ 1 - r^2/r_\sun^2 } = \cos{\\theta}

    The *CLV* is given by:

    .. math::

        I = I_0  ( 1 - \sum_{k=1}^{n} a_k * \\mu^k)

        where

        \\mu = \cos{\\theta}

    :param r: the radius of the solar disk
    :type r: float
    :param r_val: if given, returns the value of the *CLV+ at the given solar radius
    :type r_val: float [in pixel units
    :param coefficients: the coefficients of the CLV function
    :type coefficients: float np.ndarray
    :return: center to limb variation
    :rtype: np.ndarray or float
    """
    # TODO homogenize this function with limb_darkening and newton.
    with contextlib.suppress(Exception):
        if not coefficients:
            coefficients = [1., 0.44984082, 0.11827295, 0.05089393]  # PHI values with n = 4 and solar radius of 789

    with contextlib.suppress(Exception):
        r_vector = np.arange(0, np.max(r), 1) if r_val is None else r_val

    mu = 1 - np.sqrt((1 - r_vector ** 2 / r ** 2))

    f = sum(coefficients[i + 1] * mu ** (i + 1) for i in range(len(coefficients) - 1))

    return coefficients[0] * (1 - f)


def azimuthal_mask(sx: int, sy: int, center: Tuple[int, int], r: float, coefficients: np.ndarray = None) -> np.ndarray:
    """
    This functions provides a CLV mask given the dimensions and the center and radius

    :param sx: x image dimensions
    :type sx: int
    :param sy: y image dimensions
    :type sy: int
    :param center: center (xc,yc) is a tuple with the center of the disk
    :type center: tuple(xc,yx)
    :param r: the radius of the solar disk
    :type r: int
    :param coefficients: the coefficients of the CLV function [optional]
    :type coefficients: float np.ndarray
    :return: a 2D image constructed from CLV parameters
    :rtype: np.ndarray

    *Example*

    .. code-block:: python

        mask = azimuthal_mask(2048,2048,[996,1245],786)

    """

    x_axis = np.arange(0, sx, 1)
    y_axis = np.arange(0, sy, 1)
    x_axis = x_axis - center[0] + 0.5  # if dimensions are odd add 0.5 for center
    y_axis = y_axis - center[1] + 0.5
    mask = np.zeros((sy, sx))
    for i, j in itertools.product(range(sx), range(sy)):
        pxy = np.sqrt(x_axis[i] ** 2 + y_axis[j] ** 2)
        if pxy < r:
            mask[int(j), int(i)] = clv(r, r_val=pxy, coefficients=coefficients)
    return mask


def fit_clv(data: np.ndarray, centers: tuple, verbose: bool = False):
    """
    Fits the center to limb o a 2D solar image using the Newton's method

    :param data: input image. Shall be a 2D image with the full Sun
    :type: np.ndarray
    :param centers: center of image [x,y]
    :type: tuple(integer)
    :param verbose: verbose
    :type: bool
    :return: an array with the coefficients of the CLV fit to use with ``clv.py``
    :rtype: np.ndarray
    """
    intensity, radius = azimutal_average(data, [centers[0], centers[1]])
    yd, xd = data.shape

    ints = np.zeros(int(np.sqrt(xd ** 2 + yd ** 2)))  # vector of sqrt(xd^^+yd^^)
    ints_rad = np.zeros(int(np.sqrt(xd ** 2 + yd ** 2)))

    ints[:len(intensity)] = intensity
    ints_rad[:len(intensity)] = radius

    # STEP --->>> FIT LIMB DATA
    der = (np.roll(intensity, 1) - intensity)  # simple derivative
    idx_max = np.where(der == np.max(der))  # look for peak position
    clv = intensity[:int(idx_max[0]) + 1]
    clv_rho = radius[:int(idx_max[0]) + 1]
    mu = np.sqrt((1 - clv_rho ** 2 / clv_rho[-1] ** 2))

    first_coefficient = 0.5
    i0 = 100
    locations = np.where(mu > 0.1)
    coefficients = newton(clv[locations], mu[locations], [i0, first_coefficient, 0.2, 0.2, 0.2], limb_darkening)

    if verbose:
        # TODO: clean ints_fit, etc... many unused parameters
        plt.plot(clv_rho, clv)
        plt.xlabel('Solar radius [pixel]')
        plt.ylabel('Intensity [DN]')
        plt.show()

        ints_fit = np.zeros(int(np.sqrt(xd ** 2 + yd ** 2)))
        ints_syn = np.zeros(int(np.sqrt(xd ** 2 + yd ** 2)))
        ints_fit_pars = np.zeros(5)

        fit, _ = limb_darkening(mu, coefficients)
        ints_fit[:len(fit)] = fit
        ints_fit_pars[:] = coefficients

        ints_syn = np.copy(ints)
        ints_syn[:len(fit)] = fit

        # STEP --->>> NORMALIZE
        ints_syn = ints_syn / ints_fit_pars[0]
        ints_fit = ints_fit / ints_fit_pars[0]
        ints = ints / ints_fit_pars[0]

        plt.plot(ints_fit, label='fitted clv')
        plt.plot(ints, '.', label='real clv')
        plt.plot(ints_syn, '--', label='synt clv')
        plt.xlabel('Heliocentric angle [' + r'$\theta$]')
        plt.ylabel('Intensity [DN]')
        plt.legend()
        plt.show()

    print('CLV coefficients:', coefficients)

    return coefficients


def centers_flat(n_images: int, inner_radius: int, outer_radius: int, steps: int, r_width: int, binmask: list,
                 imsize: tuple, verbose=None):
    """
    Main routine of Hough transform function

    :param verbose:
    :param n_images: number of images
    :type n_images: int
    :param inner_radius: smallest radius for votes
    :type inner_radius: int
    :param outer_radius: largest radius for votes
    :type outer_radius: int
    :param steps: number of subdivisions [+1]  between inner and outer radius
    :type steps: int
    :param r_width: width of the ring for cross-correlating with circle
    :type r_width: int
    :param binmask: input image
    :type binmask: list of images
    :param imsize: image size
    :type imsize: tuple(int)
    :return: the center and radius of the image calculated using Hough
    :rtype: tuple(float)
    """

    centers = []
    radius = []
    radii = np.linspace(inner_radius, outer_radius, steps + 1)
    printc('Analyzing ', n_images, ' images', color=bcolors.OKGREEN)

    for i in range(n_images):
        # acc_conv = find_Circles(
        #    binmask[i], radii_coarse, r_width_coarse, verbose=verbose, full=True)
        acc_conv = find_Circles_ida(binmask[i], radii, r_width)
        center, rad, c, d = votes(acc_conv, radii)

        centers.append(center)
        radius.append(rad)
        printc('Found center: ', centers[i], ' and radius: ', radius[i], color=bcolors.WARNING)
        if verbose:
            plt.figure(frameon=False)
            plt.imshow(binmask[i], cmap=plt.cm.gray, alpha=.5)
            circle_fit = bin_annulus(
                imsize, radius[i], 1, full=False).astype(float)
            dd = np.array(centers[i])
            dx = dd[0] - imsize[0] // 2
            dy = dd[1] - imsize[1] // 2
            circle_fit = shift(circle_fit, shift=[dx, dy])
            plt.imshow(circle_fit, cmap=plt.cm.gray, alpha=.5)
            plt.show()

    return centers, radius


@timeit
def do_hough(image, inner_radius, outer_radius, steps, org_centers=None, method='prewitt', save=False,
             dhtr=10, normalize=False, verbose=False, otsu=None, threshold=0.15):
    """
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
            (3): if step is a negative number then uses FM find model
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

    """

    imsize = image[0].shape
    n_images = len(image)
    if org_centers is None:
        org_centers = np.tile(np.array([0., 0.], dtype=np.int16), (n_images, 1))

    ############################
    # Normalize images (using a box 100x100 in the central image)
    ############################

    if normalize:
        norma = np.mean(image[0][imsize[0] // 2 - 100:imsize[0] // 2 + 100, imsize[0] // 2 - 100:imsize[0] // 2 + 100])
        if verbose:
            print('Normalization constant: ', norma, '[calculated with first image assumed to be central one]')

        for i in range(n_images):
            image[i] = image[i] / norma

    image_dummy, threshold = FindEdges(
        image[0], threshold, method=method, dthr=dhtr, verbose=verbose, Otsu=otsu)
    binmask = [image_dummy]
    for i in range(1, n_images):
        image_dummy = FindEdges(
            image[i], threshold, method=method, verbose=verbose, Otsu=otsu)
        binmask.append(image_dummy)

    #Mask to prevent edges
    circle_mask, _ = generate_circular_mask([imsize[0] - 1, imsize[1] - 1], 1020, 1020)
    for i in range(n_images):
        binmask[i] = binmask[i] * circle_mask
    ############################
    # FIND CENTERS - COARSE SEARCH
    ############################
    # Coarse and fine compressed in one call

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
    # FIND CENTERS - FINE SEARCH
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
        # FIND CENTERS - COARSE SEARCH
        #############################
        r_width = (outer_radius - inner_radius) // steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))

        printc('from: ', inner_radius, ' to: ', outer_radius, ' steps: ', steps, ' width: ', r_width,
               color=bcolors.OKGREEN)
        centers, radius = centers_flat(n_images, inner_radius, outer_radius, steps, r_width, binmask, imsize,
                                       verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            printc("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                   (i, org_centers[i, 0], org_centers[i, 1],
                    centers[i][1], centers[i][0], radius[i]), color=bcolors.FAIL)
        ###########################
        # FIND CENTERS - FINE SEARCH
        ###########################
        mean_r = np.int(np.mean(radius))
        inner_radius = mean_r - 32
        outer_radius = mean_r + 32
        steps = 16
        r_width = (outer_radius - inner_radius) // steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))

        printc('from: ', inner_radius, ' to: ', outer_radius, ' steps: ', steps, ' width: ', r_width,
               color=bcolors.OKGREEN)
        centers, radius = centers_flat(n_images, inner_radius, outer_radius, steps, r_width, binmask, imsize,
                                       verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            printc("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                   (i, org_centers[i, 0], org_centers[i, 1],
                    centers[i][1], centers[i][0], radius[i]), color=bcolors.FAIL)
        ################################
        # FIND CENTERS - VERY FINE SEARCH
        ################################
        mean_r = np.int(np.mean(radius))
        inner_radius = mean_r - 4
        outer_radius = mean_r + 4
        steps = 8
        r_width = (outer_radius - inner_radius) // steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))

        printc('from: ', inner_radius, ' to: ', outer_radius, ' steps: ', steps, ' width: ', r_width,
               color=bcolors.OKGREEN)
        centers, radius = centers_flat(n_images, inner_radius, outer_radius, steps, r_width, binmask, imsize,
                                       verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            printc("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                   (i, org_centers[i, 0], org_centers[i, 1],
                    centers[i][1], centers[i][0], radius[i]), color=bcolors.FAIL)
    elif steps < 0:
        ##################################
        # FIND CENTERS - FM SEARCH STRATEGY
        ##################################
        inner_radius = 128
        outer_radius = 1024
        steps = 32
        r_width = (outer_radius - inner_radius) // steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))
        printc('from: ', inner_radius, ' to: ', outer_radius, ' steps: ', steps, ' width: ', r_width,
               color=bcolors.OKGREEN)
        centers, radius = centers_flat(n_images, inner_radius, outer_radius, steps, r_width, binmask, imsize,
                                       verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            print("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                  (i, org_centers[i, 0], org_centers[i, 1],
                   centers[i][1], centers[i][0], radius[i]))
        mean_r = np.int(np.mean(radius))
        inner_radius = mean_r - 32
        outer_radius = mean_r + 32
        steps = 16
        r_width = (outer_radius - inner_radius) // steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))
        printc('from: ', inner_radius, ' to: ', outer_radius, ' steps: ', steps, ' width: ', r_width,
               color=bcolors.OKGREEN)
        centers, radius = centers_flat(n_images, inner_radius, outer_radius, steps, r_width, binmask, imsize,
                                       verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            print("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                  (i, org_centers[i, 0], org_centers[i, 1],
                   centers[i][1], centers[i][0], radius[i]))
        mean_r = np.int(np.mean(radius))
        inner_radius = mean_r - 8
        outer_radius = mean_r + 8
        steps = 8
        r_width = (outer_radius - inner_radius) // steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))
        printc('from: ', inner_radius, ' to: ', outer_radius, ' steps: ', steps, ' width: ', r_width,
               color=bcolors.OKGREEN)
        centers, radius = centers_flat(n_images, inner_radius, outer_radius, steps, r_width, binmask, imsize,
                                       verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            print("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                  (i, org_centers[i, 0], org_centers[i, 1],
                   centers[i][1], centers[i][0], radius[i]))
        mean_r = np.int(np.mean(radius))
        inner_radius = mean_r - 2
        outer_radius = mean_r + 2
        steps = 4
        r_width = (outer_radius - inner_radius) // steps * 2
        print(np.linspace(inner_radius, outer_radius, steps + 1))
        printc('from: ', inner_radius, ' to: ', outer_radius, ' steps: ', steps, ' width: ', r_width,
               color=bcolors.OKGREEN)
        centers, radius = centers_flat(n_images, inner_radius, outer_radius, steps, r_width, binmask, imsize,
                                       verbose=verbose)
        print('Image |   Original  |  Inferred   |   Radius')
        for i in range(n_images):
            print("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
                  (i, org_centers[i, 0], org_centers[i, 1],
                   centers[i][1], centers[i][0], radius[i]))
    else:
        print('NO HOUGH **** WRONG')

    if save:
        status = write_shifts('hough_centers.txt', (centers, radius))
        if status != 1:
            print('Error in write_shifts')

    return centers, radius


@timeit
def fdt_flat_gen(image, rel_centers, method, radious=0, thrd=0.05, iter=15,
                 bit_trun=0, verbose=0, expand=0, c_term=0, clv=0,
                 normalize=False):
    """
    fdt_flat_gen Khun-Lin-Lorantz algorithm

    _extended_summary_

    Args:
        image (list): [n_images][y,x]
        rel_centers (ndarrray): [n_images,2] where [:,0]=dx and [:,1]=dy. 
            Displacements are given with respect to image origin (0,0)
        method (string): Defaults to 'kll'.
            'kll': J. R. Kuhn, H. Lin, and D. Loranz,  PASP 103, 1097-1108, 1991    
            'chae': Adapted from: J. Chae, Solar Physics 221: 1–14, 2004
            'alter': James N. Caron, Marcos J. Montes, and Jerome L. Obermark,
                Review of Scientific Instruments 87, 063710 (2016); doi: 10.1063/1.4954730        
        radious:  radious of circular mask. In this case, the code uses the thrd to create the mask. Defaults to 0.
        thrd: threshold above which pixels are valid (assuming image is normalized to one). Defaults to 0.05.
        iter (int, optional): maximum number of iterations in the kll algorithm. Defaults to 15.
        verbose (int, optional): Verbosity level. Defaults to 0.
        expand: how much the circular mask is expanded (positive = schrinks the mask) in %. Defaults to 0.
        c_term (ndarray, optional): np.array(n_images) intensity factor correction for chae method  . Defaults to 0.
        clv (float, optional): _description_. Defaults to False.

    Returns:
        ndarray: flat
    """

    imsize = image[0].shape
    n_images = len(image)

    def do_normalization(array, n, type='zero'):
        """
        donorma _summary_

        _extended_summary_

        Args:
            array (_type_): _description_
            n (_type_): _description_
            type (str, optional): _description_. Defaults to 'zero'.

        Returns:
            _type_: _description_
        """
        idx = np.where(n > 0)
        s = array[idx]
        #  calculate average of gain table for normalization
        sm = np.mean(s)
        sm2 = np.mean(s ** 2)
        five_sigma = SIGMA_NORMA * np.sqrt(sm2 - sm * sm)
        idx2 = np.where(np.abs(s - sm) < five_sigma)
        sm = np.mean(s[idx2])
        if type == 'zero':
            array -= sm
        if type == 'one':
            array /= sm

        return array, sm

    ############################
    # set displacements of observed images (A) with respect Object image (centered) 
    ############################
    xyshifts = np.empty([n_images, 2], dtype=int)
    xyshifts[:, 0], xyshifts[:, 1] = rel_centers[:n_images, 0] - imsize[0] // 2, rel_centers[:n_images, 1] - imsize[
        1] // 2

    ############################
    # calculate masks
    ############################
    mask = np.zeros([n_images, imsize[0], imsize[1]], dtype=np.int8)

    if radious != 0:  # In case radius of solar disk is provided....
        maskn, dummy = generate_circular_mask([imsize[0] - 1, imsize[1] - 1], radious - expand, radious - expand)
        printc('Using circular mask', color=bcolors.OKGREEN)
        for i in range(n_images):
            mask[i] = shift(maskn, shift=[xyshifts[i, 0], xyshifts[i, 1]], fill_value=0)
    else:
        # find pixel coordinates with solar information (> thrd given by user, default = 0.05)
        # This step assumes input data has a mean value of one.
        printc('Setting mask depending n threshold level', thrd, color=bcolors.OKGREEN)
        for i in range(n_images):
            x, y = np.where(image[i] > thrd)
            mask[i][x, y] = 1

    circle_mask, _ = generate_circular_mask([imsize[0] - 1, imsize[1] - 1], 1020, 1020)
    for i in range(n_images):
        mask[i] = mask[i] * circle_mask

    ############################
    # NORMALIZATION Option 1 
    ############################
    norma = np.zeros(n_images)
    if normalize == 1:
        for i in range(n_images):
            image[i], norma[i] = do_normalization(image[i], mask[i], type='one')
            printc('Normalization image ', i, ': ', norma[i], color=bcolors.YELLOW)

    ############################
    # calculate a mask with the center to limb variation of the sun to flat the images
    ############################
    if clv:
        printc('Correcting CLV for masking spots', color=bcolors.OKGREEN)

        if clv < 0:
            printc('Reevaluating limb darkening', color=bcolors.OKBLUE)
            pars = fit_clv(image[0], rel_centers[0, :])
        else:
            pars = False
        clv = np.abs(clv)

        amask = azimuthal_mask(imsize[0], imsize[1], (imsize[0] // 2, imsize[1] // 2), radious, coefficients=pars)

        # find px outside threshols
        printc('... clv threshold ', clv, color=bcolors.OKBLUE)

        for i in range(n_images):
            image_dummy = image[i] / shift(amask, shift=[xyshifts[i, 0], xyshifts[i, 1]], fill_value=0)  # if any scale
            image_dummy[np.bitwise_not(np.isfinite(image_dummy))] = 0
            if normalize != 1:
                image_dummy, norma_ = do_normalization(image_dummy, mask[i], type='one')
            idx = np.where(image_dummy > (1 + clv))
            mask[i][idx] = 0
            idx = np.where(image_dummy < (1 - clv))
            mask[i][idx] = 0
            mask[i] = expand_mask(mask[i], pixels=1, step=2)

    ############################
    # DO LOG
    ############################

    D = np.log10(image)
    # replace NaNs and Infs by 0
    D[np.bitwise_not(np.isfinite(D))] = 0

    ############################
    # NORMALIZATION Option 2
    ############################
    if normalize == 2:
        tn = 0
        for i in range(n_images):
            norma = np.mean(D[i][np.where(mask[i] != 0)])
            print('Normalization to around zero in the log images: ', norma)
            D[i] = D[i] - norma
            tn += norma
        tn /= n_images

    if method == 'kll':

        ############################
        # CALCULATE CONSTANT
        ############################

        n = np.zeros([imsize[0], imsize[1]], dtype=int)
        sum_image = np.zeros([imsize[0], imsize[1]], dtype=float)

        # pairs = [] #for debugging purposes

        print('Rel centers: ', rel_centers)
        #  for [iq, ir] in itertools.combinations(range(n_images), 2): # overall 36 combinations
        for iq in range(1, n_images):  # loop in iq
            for ir in range(iq):  # loop in ir
                # shift of iq with respect ir
                dx = rel_centers[iq, 0] - rel_centers[ir, 0]
                dy = rel_centers[iq, 1] - rel_centers[ir, 1]
                if verbose == 2:
                    print('dx,dy', dx, dy, iq, ir)

                t_mask_1 = mask[ir] & shift(mask[iq], [-dx, -dy])
                t_mask_2 = mask[iq] & shift(mask[ir], [dx, dy])
                # t_mask = t_mask_1 & t_mask_2 # compound mask only used for mean (commented Sep 10-2022)

                t_image_1 = shift(D[iq], [-dx, -dy])
                t_image_2 = shift(D[ir], [dx, dy])
                aa = (D[iq] - t_image_2) * t_mask_2  # add _2
                bb = (D[ir] - t_image_1) * t_mask_1  # add _1
                image_pair = aa + bb
                sum_image += image_pair
                # pairs.append(image_pair) #for debugging purposes only

                # n += t_mask_1 # accumulate valid pixels first sumatorio
                # n += t_mask_2 # accumulate valid pixels second sumatorio
                t_mask = t_mask_1 + t_mask_2  # DOS  | vs + and /2.  Sep 10-2022
                n += t_mask  # DOS

        K = sum_image / n.astype(float)
        # replace NaNs and Infs by 0
        K[np.bitwise_not(np.isfinite(K))] = 0
        if verbose == 1:
            plt.imshow(K, cmap='gray')
            plt.clim(vmax=0.02, vmin=-0.02)
            plt.colorbar()
            plt.show()

        G = np.copy(K)
        if bit_trun == 1:
            K = np.int32(K * 256) / 256  # bit truncation
        k = np.power(10, K)
        if bit_trun == 1:
            k = np.int32(k * 256) / 256  # bit truncation

        for itera in range(iter):
            r_res = np.zeros(imsize, dtype=float)
            for iq in range(1, n_images):  # loop in iq
                for ir in range(iq):  # loop in ir
                    # shift of iq with respect ir
                    dx = rel_centers[iq, 0] - rel_centers[ir, 0]
                    dy = rel_centers[iq, 1] - rel_centers[ir, 1]
                    if verbose == 2:
                        print('dx,dy,sqrt(dx**2+dy**2)', dx, dy, np.sqrt(dx ** 2 + dy ** 2), iq, ir)
                    t_mask_1 = mask[ir] & shift(mask[iq], [-dx, -dy])
                    t_mask_2 = mask[iq] & shift(mask[ir], [dx, dy])

                    t_image_1 = shift(G, [-dx, -dy]) * t_mask_1
                    t_image_2 = shift(G, [dx, dy]) * t_mask_2
                    correction = (t_image_1 + t_image_2)
                    r_res += correction

                    # G = K + r_res / n.astype(np.float64)
            corr = r_res / n.astype(float)
            # corr = donorma(corr,n)
            G = K + corr
            # replace NaNs and Infs by 0
            G[np.bitwise_not(np.isfinite(G))] = 0

            idx = np.where(n > 0)
            s = G[idx]
            #  calculate average of gain table for normalization
            sm = np.mean(s)
            sm2 = np.mean(s ** 2)
            five_sigma = SIGMA_NORMA * np.sqrt(sm2 - sm * sm)
            idx2 = np.where(np.abs(s - sm) < five_sigma)
            sm = np.mean(s[idx2])
            # G[idx] = G[idx] - sm  
            G = G - sm
            # G,sm = donorma(G,n,dntype='zero')
            change = np.std(K - G)

            printc('Iteration: ', itera, ' of ', iter, 'Gain mean:', sm, 'using rms level', SIGMA_NORMA,
                   ' rms variation:', change, color=bcolors.OKBLUE)
            if verbose == 1:
                plt.imshow(G, cmap='gray')
                plt.clim(vmax=0.02, vmin=-0.02)
                plt.colorbar()
                plt.show()

                plt.imshow(K - G, cmap='gray')
                plt.clim(vmax=0.02, vmin=-0.02)
                plt.colorbar()
                plt.show()

        g = np.power(10, G)
        g[np.bitwise_not(np.isfinite(g))] = 0

        return g, norma

    elif method == 'chae2':

        tmask = np.sum(mask, axis=0)
        # mask_Ob,dummy = generate_circular_mask([imsize[0]-1, imsize[1]-1],radious,radious)  
        # mask_Ob,dummy = generate_circular_mask([imsize[0]-1, imsize[1]-1],radious - expand,radious - expand)  
        mask_Ob = shift(mask[0], shift=[-xyshifts[0, 0], -xyshifts[0, 1]], fill_value=0)

        # Constant term
        if c_term == 0:
            fit_c = 0
            c_term = np.log10(np.ones((n_images)))
        else:
            c_term = np.log10(np.ones((n_images)))

        flat = np.log10(np.ones_like(D[0]))
        Ob = np.zeros_like(D[0])

        for i in range(n_images):
            # shift input image to the center of the frame.
            Ob += shift(D[i], shift=-xyshifts[i, :])
        Ob = Ob / float(n_images)
        Ob[np.isneginf(Ob)] = 0.
        Ob[np.isnan(Ob)] = 0.

        idx = tmask >= 1

        for k in range(iter):

            numerator = np.zeros((imsize))
            for i in range(n_images):
                numerator += (
                        (c_term[i] + Ob - shift(D[i] - flat, shift=-xyshifts[i, :])) * mask_Ob)  # (i+xk,j+yk) Eq 8

            # Ob -= (numerator/mask_Ob/n_images)
            Ob -= (numerator / mask_Ob / tmask)
            Ob[np.isneginf(Ob)] = 0.
            Ob[np.isnan(Ob)] = 0.

            if verbose == 3:
                tshow = np.power(10, Ob)

                plt.imshow(tshow, cmap='gray', vmin=np.median(tshow[idx]) * 0.9, vmax=np.median(tshow[idx]) * 1.1)
                plt.show()

            numerator = np.zeros((imsize))
            for i in range(n_images):
                dummy = (c_term[i] + shift(Ob, shift=+xyshifts[i, :]) + flat - D[i]) * mask[i]
                dummy[np.bitwise_not(np.isfinite(dummy))] = 0
                numerator += dummy
                c_term[i] -= (np.sum(dummy[idx]) / np.sum(mask[i]))

            dummy = (numerator / tmask)
            flat -= dummy
            flat[np.isneginf(flat)] = 0.
            flat[np.isnan(flat)] = 0.

            if verbose == 3:
                plt.imshow(flat, cmap='gray', vmin=-0.02, vmax=0.02)
                plt.show()
            if verbose >= 1:
                print('Iter: ', k, ' STD: ', np.max(np.abs(dummy[idx])), np.exp(c_term))

            s = flat[idx]
            sm = np.mean(s)
            sm2 = np.mean(s ** 2)
            five_sigma = 5 * np.sqrt(sm2 - sm * sm)
            print('Iteration: ', k, five_sigma, '5*rms', sm, ' of ', k, ' STD: ', np.max(np.abs(dummy[idx])))

        flat, mf = do_normalization(flat, tmask)
        Ob = Ob + mf + np.mean(c_term)
        # flat = flat - np.mean(flat)
        # Ob = Ob + np.mean(flat) + np.mean(c_term)
        c_term = c_term - np.mean(c_term)

        flat = np.power(10, flat)
        flat[np.bitwise_not(np.isfinite(flat))] = 0

        if verbose >= 2:
            plt.imshow(flat, cmap='gray', vmin=0.95, vmax=1.05)  # vmin=)np.min(,vmax=0.05)
            plt.show()

        return flat, norma

    elif method == 'chae':

        tmask = np.sum(mask, axis=0)
        mask_Ob, dummy = generate_circular_mask([imsize[0] - 1, imsize[1] - 1], radious, radious)

        # Constant term
        if c_term == 0:
            fit_c = 0
            c_term = np.log10(np.ones(n_images))
        else:
            c_term = np.log10(np.ones(n_images))

        flat = np.log10(np.ones_like(D[0]))
        Ob = np.zeros_like(D[0])

        for i in range(n_images):
            # shift input image to the center of the frame.
            # 
            Ob += shift(D[i], shift=-xyshifts[i, :])
        Ob = Ob / float(n_images)
        idx = tmask >= 1

        for k in range(iter):

            numerator = np.zeros(imsize)
            for i in range(n_images):
                numerator += (
                        (c_term[i] + Ob - shift(D[i] - flat, shift=-xyshifts[i, :])) * mask_Ob)  # (i+xk,j+yk) Eq 8

            Ob -= (numerator / mask_Ob / n_images)
            Ob[np.isneginf(Ob)] = 0.
            Ob[np.isnan(Ob)] = 0.

            if verbose == 1:
                plt.imshow(Ob, cmap='gray', vmin=1, vmax=2)
                plt.show()

            numerator = np.zeros(imsize)
            for i in range(n_images):
                dummy = (c_term[i] + shift(Ob, shift=+xyshifts[i, :]) + flat - D[i]) * mask[i]
                dummy[np.bitwise_not(np.isfinite(dummy))] = 0
                numerator += dummy
                c_term[i] -= (np.sum(dummy) / np.sum(mask[i]))

            dummy = (numerator / tmask)
            flat -= dummy
            flat[np.isneginf(flat)] = 0.
            flat[np.isnan(flat)] = 0.

            if verbose == 1:
                plt.imshow(flat, cmap='gray', vmin=-0.02, vmax=0.02)
                plt.show()
            if verbose >= 1:
                print('Iter: ', k, ' STD: ', np.max(np.abs(dummy[idx])), np.exp(c_term))

            s = flat[idx]
            sm = np.mean(s)
            sm2 = np.mean(s ** 2)
            five_sigma = 5 * np.sqrt(sm2 - sm * sm)
            print('Iteration: ', k, five_sigma, '5*rms', sm, ' of ', k, ' STD: ', np.max(np.abs(dummy[idx])))

        # flat = flat - np.mean(flat)
        # Ob = Ob + np.mean(flat) + np.mean(c_term)
        # c_term = c_term - np.mean(c_term)
        flat, mf = do_normalization(flat, tmask)
        Ob = Ob + mf + np.mean(c_term)
        c_term = c_term - np.mean(c_term)

        flat = np.power(10, flat)
        flat[np.bitwise_not(np.isfinite(flat))] = 0

        if verbose >= 2:
            plt.imshow(flat, cmap='gray', vmin=0.95, vmax=1.05)  # vmin=)np.min(,vmax=0.05)
            plt.show()

        return flat, norma

    elif method == 'alter':

        # Extracting flat-field images from scene-based image sequences using phase
        # correlation
        # James N. Caron, Marcos J. Montes, and Jerome L. Obermark
        # Citation: Review of Scientific Instruments 87, 063710 (2016); doi: 10.1063/1.4954730
        # View online: https://doi.org/10.1063/1.4954730
        # View Table of Contents: http://aip.scitation.org/toc/rsi/87/6
        # Published by the American Institute of Physics

        Gf = np.zeros([imsize[0], imsize[1]], dtype=np.float64)
        n = np.zeros([imsize[0], imsize[1]], dtype=np.float64)
        sum_image = np.zeros([imsize[0], imsize[1]], dtype=np.float64)
        Gr = np.zeros([n_images, imsize[0], imsize[1]], dtype=np.float64)

        for itera in range(iter):
            ir = 0
            for iq in range(n_images):
                dx = rel_centers[iq, 0] - rel_centers[ir, 0]  # shift of iq with respect ir
                dy = rel_centers[iq, 1] - rel_centers[ir, 1]
                t_mask = mask[ir] * shift(mask[iq], [dx, dy])
                n += t_mask
                t_image = shift(D[iq] - Gf, [dx, dy])
                sum_image += (t_mask * t_image)

            K = sum_image / n
            K[np.isneginf(K)] = 0
            K[np.isnan(K)] = 0
            idx = np.where(K > 5)
            mask[0][idx] = 0
            idx2 = np.where(n == 0)
            K[idx2] = 0

            iq = 0
            for ir in range(n_images):
                dx = rel_centers[iq, 0] - rel_centers[ir, 0]  # shift of iq with respect ir
                dy = rel_centers[iq, 1] - rel_centers[ir, 1]
                Kn = shift(K * mask[0], [dx, dy])

                G = (D[ir] - Kn) * mask[ir]
                G[np.isneginf(G)] = 0
                G[np.isnan(G)] = 0
                Gr[ir, :, :] = G

            m = np.sum(mask, axis=0)
            Gf = np.sum(Gr, axis=0) / m
            Gf[np.isneginf(K)] = 0
            Gf[np.isnan(K)] = 0
            print('Iteration: ', itera)

        g = np.power(10, Gf)
        g[np.bitwise_not(np.isfinite(g))] = 0

        return g, norma
    else:
        return None, None


def fdt_flat(files, wavelength, npol, method='kll', dark=None, r_shifts=0, shifts=None, verbose=1,
             correct_ghost=0, expand=1, thrd=0., iter=4, normalize=False, disp_method='Hough', c_term=0,
             inner_radius=400, outer_radius=800, steps=20, shifts_file=False, imasize=None, single=False,
             clv=False,threshold=0.05,dhtr = 2):
    """
    fdt_flat _summary_

    _extended_summary_

    Args:
        files (_type_): _description_
        wavelength (_type_): _description_
        npol (_type_): _description_
        method (str, optional): _description_. Defaults to 'kll'.
        dark (nd.array, optional): Dark current. Dimensions should be same as input images. 
            Note that the correction is done directly substracting the dark from images. Defaults to None.
        r_shifts (int, optional): _description_. Defaults to 0.
        shifts (int, optional): Method for finding the image centers and Sun radious.
            0: find the disc center position and Sun radious using either the Hough transform of a circle fit method
                The method is given by [disp_method] keyword.
            1: read the shifts from files with a specific format and whose initial leadtext is given by [shifts_file]
            2: shifts are provided by the uses in [center] and [radious] keywords
            Defaults to None.
        verbose (int, optional): _description_. Defaults to 1.
        correct_ghost (int, optional): _description_. Defaults to 0.
        expand (int, optional): _description_. Defaults to 1.
        thrd (_type_, optional): _description_. Defaults to 0..
        iter (int, optional): _description_. Defaults to 4.
        normalize (bool, optional): _description_. Defaults to False.
        disp_method (str, optional): _description_. Defaults to 'Hough'.
        c_term (int, optional): _description_. Defaults to 0.
        inner_radius (int, optional): _description_. Defaults to 400.
        outer_radius (int, optional): _description_. Defaults to 800.
        steps (int, optional): _description_. Defaults to 20.
        shifts_file (bool, optional): _description_. Defaults to False.
        imasize (list, optional): _description_. Defaults to [2048,2048].
        single (bool, optional): This softwave is written explicitly for SOPHI since it uses SOPHI specially tailored read routines.
            However, one can use it with simulated data. In such a case, the keywords [pol_state] and wavelengt
            loss their meaning and one should set single = True. Defaults to False.
        clv (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
    """

    ############################
    # open 9 fits files and get one images corresponding to one polarization [pol_state] and one wavelength [wavelength]. 
    # This can be done all at once but I have this like that because I am lazy.
    ############################

    if imasize is None:
        imasize = [2048, 2048]
    if single:
        image = []
        header = []
        for i in files:
            img, hd = fits_get(i)
            image.append(img)
            header.append(hd)
    else:
        image = []
        header = []
        for i in files:
            im, hd = fits_get_part(i, wavelength, npol)
            image.append(im)
            header.append(hd)

    n_images = len(image)
    ys, xs = image[0].shape

    ############################
    # Correct Dark if not done in advance
    # dark current is provided in dark keyword
    ############################

    with contextlib.suppress(Exception):
        printc('...Correcting dark current...', color=bcolors.OKGREEN)
        for i in range(n_images):
            image[i] = np.abs(image[i] - dark)
            # REVIEW the dark correction is usually not O.K. in a 0.2-0.3%. I can remove the residual dark but 
            # this will affect the INTENSITY LEVEL. So, it is mandatory to correct the darks BEFOREHAND!!!
            # TODO PUT THIS AS AN ADDITIONAL INPUT OPTION
            # columns = np.mean(image[i][0:40,:],axis=0)
            # image[i] = image[i] - columns[np.newaxis,:]
            # image[i] = np.abs(image[i])
            # plt.imshow(image[i],vmax=100)
            # plt.show()

    ############################
    # Determine the image shifts
    # there are three different options: r_shifts = [0,1,2]
    # 0 find the disc center position and Sun radious using either the Hough transform of a circle fit method
    # 1 read the shifts from files with a specific format and whose initial leadtext is given by [shifts_file]
    # 2 shifts are provided by the uses in centers and radious keywords
    ############################

    if r_shifts == 1:
        try:
            printc('... read user input shifts_file ...', color=bcolors.OKGREEN)
            centers = read_shifts(shifts_file + '_cnt_w' + str(wavelength) + '_n' + str(npol) + '.txt')
            radius = read_shifts(shifts_file + '_rad_w' + str(wavelength) + '_n' + str(npol) + '.txt')
            printc("File: ", shifts_file + '_cnt_w' + str(wavelength) + '_n' + str(npol) + '.txt', color=bcolors.OKBLUE)
            for i in range(n_images):
                printc('Image', i, 'c: ', centers[i, 0], ',', centers[i, 1], ' rad: ', radius[i], color=bcolors.YELLOW)
        except Exception:
            printc("Unable to open file: {}", shifts_file + '_cnt_w' + str(wavelength) + '_n' + str(npol) + '.txt',
                   color=bcolors.FAIL)

    elif r_shifts == 2:
        printc('... shifts provided by user ...', color=bcolors.OKGREEN)
        centers = shifts[0]
        printc(centers, '... read ', color=bcolors.OKBLUE)
        radius = shifts[1]
        printc(radius, '... read ', color=bcolors.OKBLUE)

    elif r_shifts == 0:

        if disp_method == 'Hough':
            printc('... calculating shifts using Hough ...', color=bcolors.OKGREEN)

            centers, radius = do_hough(image, inner_radius, outer_radius, steps, verbose=verbose, threshold=threshold,dhtr=dhtr)

            if shifts_file:
                _ = write_shifts(shifts_file + '_cnt_w' + str(wavelength) + '_n' + str(npol) + '.txt', centers)
                _ = write_shifts(shifts_file + '_rad_w' + str(wavelength) + '_n' + str(npol) + '.txt', radius)

        elif disp_method == 'FFT':
            printc('... calculating shifts using FFT _NOT TESTED_ ...', color=bcolors.WARNING)
            printc('Input parameter "expand" should be negative number representing solar disk', color=bcolors.WARNING)
            image_dummy = np.zeros((n_images, ys, xs))
            for i in range(n_images):
                image_dummy[i, :, :] = image[i]
            s_y, s_x, _ = PHI_shifts_FFT(image_dummy, prec=5, verbose=False, norma=False, coarse_prec=150)
            centers = np.zeros((n_images, 2))
            radius = np.zeros(n_images)
            radius[:] = -expand
            expand = 5
            centers[:, 0] = -s_x + xs // 2
            centers[:, 1] = -s_y + ys // 2
            for i in range(n_images):
                printc('Image', i, 'c: ', centers[i, 0], ',', centers[i, 1], ' rad: ', radius[i], '*',
                       color=bcolors.YELLOW)

        elif disp_method == 'circle':
            printc('... calculating shifts using circle fit method ...', color=bcolors.OKGREEN)
            centers = np.zeros((n_images, 2))
            radius = np.zeros(n_images)
            for i in range(n_images):
                centers[i, 1], centers[i, 0], radius[i] = find_center(image[i], sjump=4, njumps=100, threshold=0.8)
                printc('Image', i, 'c: ', centers[i, 0], ',', centers[i, 1], ' rad: ', radius[i], color=bcolors.YELLOW)

            if shifts_file:
                _ = write_shifts(shifts_file + '_cnt_w' + str(wavelength) + '_n' + str(npol) + '.txt', centers)
                _ = write_shifts(shifts_file + '_rad_w' + str(wavelength) + '_n' + str(npol) + '.txt', radius)
        else:
            printc('Error in disp_method option. Given', disp_method, 'while options are [Hough,FFT,circle]',
                   color=bcolors.FAIL)
            raise Exception
    else:
        printc('Error in r_shift option. Given', r_shifts, 'while options are 0,1, or 2', color=bcolors.FAIL)
        raise Exception

    # make sure we have integer numpy numbers in the centers
    centers = np.array(centers).astype(int); mean_radii = np.mean(radius)

    loop = 0
    for i in range(n_images):

        if 'CRPIX1' in header[loop]:  # Check for existence
            header[loop]['CRPIX1'] = (round(centers[loop, 0], 2))
        else:
            header[loop].append('CRPIX1')
            header[loop]['CRPIX1'] = (round(centers[loop, 0], 2))

        if 'CRPIX2' in header[loop]:  # Check for existence
            header[loop]['CRPIX2'] = (round(centers[loop, 1], 2))
        else:
            header[loop].append('CRPIX2')
            header[loop]['CRPIX2'] = (round(centers[loop, 1], 2))
        loop += 1

    ############################
    # Correct ghost in images
    ############################

    if correct_ghost == 1:
        print('comprueba1 ', centers)
        for i in range(9):
            image[i], header[i] = phi_correct_ghost_single(image[i], header[i], radius[i],
                                                           verbose=verbose)  # ,center = centers[i,:])

    if thrd != 0:
        gain, norma = fdt_flat_gen(image, centers, method, iter=iter, verbose=verbose, c_term=c_term,
                                   normalize=normalize, thrd=thrd, clv=clv)
    else:
        gain, norma = fdt_flat_gen(image, centers, method, iter=iter, verbose=verbose, c_term=c_term,
                                   normalize=normalize, radious=mean_radii, expand=expand, clv=clv)
    return gain, norma


def fdt_flat_preprocessing(file: str = None, dark_f: str = None, verbose: bool = True, correct_ghost: int = 0,
                           correct_prefilter: bool = False,
                           prefilter_fits: str = '0000990710_noMeta.fits', version='01'):
    # TODO: This preprocessing should be done with the main phifdt_flat.py. So far it is here because of lack of time.

    if os.path.isfile(file):
        printc("Data file exist", bcolors.OKGREEN)
    else:
        printc("Data file do not exist in current location ", file, bcolors.FAIL)
        raise Exception("input files error")

    if os.path.isfile(dark_f):
        printc("Dark file exist", bcolors.OKGREEN)
    else:
        printc("Dark file do not exist in current location ", dark_f, bcolors.FAIL)
        raise Exception("input files error")

    if correct_prefilter:
        if os.path.isfile(prefilter_fits):
            printc("Prefilter file exist", bcolors.OKGREEN)
        else:
            printc("Prefilter file do not exist in current location ", prefilter_fits, bcolors.FAIL)
            raise Exception("input files error")

    # CHECK IF input is FITS OR FITS.GZ
    if file.endswith('.fits'):
        filetype = '.fits'
    elif file.endswith('.fits.gz'):
        filetype = '.fits.gz'
    else:
        raise ValueError("input data file type nor .fits neither .fits.gz")

    # -----------------
    # READ DATA
    # -----------------

    try:
        printc('-->>>>>>> Reading Data file: ' + file, color=bcolors.OKGREEN)
        data, header = fits_get(file)
    except:
        printc("ERROR, Unable to open fits file: {}", file, color=bcolors.FAIL)
        raise Exception("input files error")

    did = header['PHIDATID']
    acc = header['ACCACCUM']

    printc('-->>>>>>> data DID ' + did, color=bcolors.OKGREEN)

    printc('-->>>>>>> Reshaping data to [wave,Stokes,y-dim,x-dim] ', color=bcolors.OKGREEN)
    zd, yd, xd = data.shape
    data = np.reshape(data, (zd // 4, 4, yd, xd))
    data = np.ascontiguousarray(data)

    # -----------------
    # TAKE DATA DIMENSIONS
    # -----------------
    PXBEG1 = int(header['PXBEG1']) - 1
    PXEND1 = int(header['PXEND1']) - 1
    PXBEG2 = int(header['PXBEG2']) - 1
    PXEND2 = int(header['PXEND2']) - 1
    printc('Dimensions: ', PXBEG1, PXEND1, PXBEG2, PXEND2, color=bcolors.OKGREEN)

    if xd != (PXEND1 - PXBEG1 + 1) or yd != (PXEND2 - PXBEG2 + 1):
        printc('ERROR, Keyword dimensions and data array dimensions dont match ', color=bcolors.FAIL)
        return 0
    if xd < 2047:
        printc('         data cropped to: [', PXBEG1, ',', PXEND1, '],[', PXBEG2, ',', PXEND2, ']',
               color=bcolors.WARNING)

    data_scale = fits_get(file, get_scaling=True)

    # -----------------
    # READ AND CORRECT DARK FIELD
    # -----------------
    data, header = phi_correct_dark(dark_f, data, header, data_scale, verbose=verbose)

    # -----------------
    # GET INFO ABOUT VOLTAGES/WAVELENGTHS, determine continuum and new flat
    # -----------------
    printc('-->>>>>>> Obtaining voltages from data ', color=bcolors.OKGREEN)
    wave_axis, voltagesData, tunning_constant, cpos, ref_wavelength = fits_get_sampling(file)
    printc('          Data FG voltages: ', voltagesData, color=bcolors.OKBLUE)
    printc('          Continuum position at wave: ', cpos, color=bcolors.OKBLUE)
    printc('          Data ref_wavelength [mA]: ', ref_wavelength, color=bcolors.OKBLUE)
    printc('          Data wave axis [mA]: ', wave_axis, color=bcolors.OKBLUE)
    printc('          Data wave axis - axis[0] [mA]: ', wave_axis - wave_axis[0], color=bcolors.OKBLUE)
    dummy_1 = (voltagesData - np.roll(voltagesData, -1)) * (tunning_constant * 1000)
    dummy_2 = np.sort(np.abs(dummy_1))
    sampling = np.mean(dummy_2[0:-2])
    printc('          Data average sampling [mA]: ', sampling, ' using tunning constant: ', (tunning_constant * 1000),
           color=bcolors.OKBLUE)

    # -----------------
    # CORRECT PREFILTER 
    # -----------------

    if correct_prefilter:
        data, header = phi_correct_prefilter(prefilter_fits, header, data, voltagesData, verbose=verbose)

    # -----------------
    # FIND DATA CENTER 
    # -----------------
    # TODO: We can just save this information into the header and we are done!!!!

    # printc('-->>>>>>> finding the center of the solar disk (needed for masking) ',color=bcolors.OKGREEN)
    # try:
    #     if center_method == 'Hough':
    #         inner_radius,outer_radius,steps = hough_params
    #         c, radius,threshold = find_circle_hough(data[0,0,:,:],inner_radius,outer_radius,steps,threshold = 0.01,normalize=False,verbose=False)
    #         #c = np.roll(c,1)
    #         cx = c[0]
    #         cy = c[1]
    #         #TBE PUT IN CORRECT UNITS
    #     elif center_method  == 'circlefit':
    #         cy,cx,radius=find_center(data[0,0,:,:],sjump = 4,njumps = 100,threshold = 0.8)  #OJO Cy... Cx
    #         c = np.array([int(cx),int(cy)])   #El vector es [0,1,2,...] == [x,y,z,...] == [cx,cy,cz,...] Pero esto ultimo esta al reves
    #         radius = int(radius)
    #     elif center_method == None:
    #         #get from header
    #         cx = header['CRPIX1']
    #         cy = header['CRPIX2']
    #         c = np.array([int(cx),int(cy)])   #El vector es [0,1,2,...] == [x,y,z,...] == [cx,cy,cz,...] Pero esto ultimo esta al reves
    #         radius = header['RSUN_ARC']/header['CDELT1']
    #     else:  
    #         raise ValueError("ERROR in center determination method - check input 'circlefit','Hough',null/None") 
    # except ValueError as err:
    #     print(err.args)
    #     return 0

    # #Uptade header with new centers
    # if center_method == 'Hough' or center_method  == 'circlefit':
    #     printc('          Uptade header with new center:',color=bcolors.OKBLUE)
    #     printc('          OLD center:',color=bcolors.OKBLUE)
    #     printc('                  at: CRPIX1[x]=',header['CRPIX1'],' CRPIX2[y]=',header['CRPIX2'],' radius=',radius,color=bcolors.OKBLUE)
    #     header['history'] = ' CRPIX 1 and CRPIX2 uptated from ' + str(header['CRPIX1'])+ ' and ' + str(header['CRPIX2'])
    #     header['CRPIX1'] = (round(cx, 2))
    #     header['CRPIX2'] = (round(cy, 2))
    #     printc('          NEW center:',color=bcolors.OKBLUE)
    #     printc('                  at: CRPIX1[x]=',header['CRPIX1'],' CRPIX2[y]=',header['CRPIX2'],' radius=',radius,color=bcolors.OKBLUE)
    #     printc('ATTENTION: Keywords CRVAL1 and CRVAL2 are NOT updated but should be SET to zero!!!!',color=bcolors.FAIL)
    # else:
    #     printc('          Using header image center:',color=bcolors.OKBLUE)
    #     printc('                  at: CRPIX1[x]=',header['CRPIX1'],' CRPIX2[y]=',header['CRPIX2'],' radius=',radius,color=bcolors.OKBLUE)

    # -----------------
    # GHOST CORRECTION  
    # -----------------

    if correct_ghost:
        # here I am recalculating the correct center of the image. If I use the header it does not work
        cy, cx, radius = find_center(data[0, 0, :, :], sjump=4, njumps=100, threshold=0.8)  # OJO Cy... Cx
        c = np.array([int(cx),
                      int(cy)])
        # El vector es [0,1,2,...] == [x,y,z,...] == [cx,cy,cz,...] Pero esto ultimo esta al reves
        radius = int(radius)
        data, header = phi_correct_ghost(data, header, radius, verbose=verbose)

    printc('---------------------------------------------------------', color=bcolors.OKGREEN)

    # -----------------
    # SAVE DATA
    # -----------------
    # basically replace L1 by L1.5
    outfile = set_level(file, 'L1', 'L1.5')
    outfile = set_level(outfile, 'ilam', 'ilam_offpoint')
    outfile = append_id(outfile, filetype, version, did)

    printc(' Saving data to:', outfile)

    data = data.astype(float)

    with pyfits.open(file) as hdu_list:
        hdu_list[0].data = data
        hdu_list[0].header = header
        hdu_list.writeto(outfile, overwrite=True)

    return
