################################################################################
################################################################################
# KLL LIBRARY 
################################################################################
################################################################################

import numpy as np
import matplotlib.pyplot as plt
from .tools import *
from .phi_gen import *
from .phi_utils import *
from .phi_fits import *
from .phi_reg import *
from SPGPylibs.GENtools import *

@timeit
def do_hough(image,inner_radius, outer_radius, steps, org_centers=None,method='prewitt',save=False,
            dhtr=10,normalize = False,verbose=False,Otsu = None,threshold = 0.15):
    '''
    Do Hough Transform
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

    centers = []
    radius = []
    r_width_coarse = (outer_radius - inner_radius)//steps * 2
    radii_coarse = np.linspace(inner_radius, outer_radius, steps)
    print('Analizing ', n_images, ' images (coarse search)')

    for i in range(n_images):
    #acc_conv = find_Circles(
    #    binmask[i], radii_coarse, r_width_coarse, verbose=verbose, full=True)
        acc_conv = find_Circles_ida(binmask[i], radii_coarse, r_width_coarse)
        center,rad,c,d = votes(acc_conv, radii_coarse)

        centers.append(center)
        radius.append(rad)
        print('Found center: ', centers[i], ' and radius: ', radius[i])
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

    print('Image |   Original  |  Inferred   |   Radius')
    for i in range(n_images):
        print("  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
            (i, org_centers[i, 0], org_centers[i, 1],
            centers[i][0], centers[i][1], radius[i]))

    ############################
    #FIND CENTERS - FINE SEARCH
    ############################

    centers_fine = []
    radius_fine = []

    mean_r = np.mean(radius)
    print('pp',mean_r)
    inner_radius = mean_r-20
    outer_radius = mean_r+20
    steps = 20
    r_width_fine = 5
    radii_fine = np.linspace(inner_radius, outer_radius, steps)
    print('Analizing ', n_images, ' images (fine case)')

    for i in range(n_images):
        acc_conv = find_Circles_ida(binmask[i], radii_fine, r_width_fine,verbose=False)
        center,rad,c,d = votes(acc_conv, radii_fine)
        centers_fine.append(center)
        radius_fine.append(rad)
        print('Found center: ', centers_fine[i],
            ' and radius: ', radius_fine[i])
        if verbose == True:
            fig = plt.figure(frameon=False)
            im1 = plt.imshow(binmask[i], cmap=plt.cm.gray, alpha=.5)
            circle_fit = bin_annulus(
                imsize, radius_fine[i], 1, full=False).astype(float)
            dd = np.array(center)
            dx = dd[0] - imsize[0]//2
            dy = dd[1] - imsize[1]//2
            circle_fit = shift(circle_fit, shift=[dx,dy])
            im2 = plt.imshow(circle_fit, cmap=plt.cm.gray, alpha=.5)
            plt.show()

    print('Method  |  Image |   Original  |  Inferred   |   Radius')
    for i in range(n_images):
        print(" Coarse  %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
            (i, org_centers[i, 0], org_centers[i, 1],
            centers[i][0], centers[i][1], radius[i]))
        print(" Fine    %2.0f  | (%4.0f,%4.0f) | (%4.0f,%4.0f) |  %6.2f" %
            (i, org_centers[i, 0], org_centers[i, 1],
            centers_fine[i][0], centers_fine[i][1], radius_fine[i]))
    
    if save == True:
        status = write_shifts('hough_centers.txt', (centers_fine,radii_fine))
        if status != 1:
            print('Error in write_shifts')
    return centers_fine, radius_fine

@timeit  
def fdt_flat_gen(image, rel_centers, method, radious = 0, thrd=0.05, iter=15, \
    bit_trun = 0,verbose = 0, expand=0, c_term = 0):
    '''
    Khun-Lin-Lorantz algorithm ()
    Input: 
    image -> [n_images][y,x]
    rel_centers -> [n_images,2] where [:,0]=dx and [:,1]=dy 
        Displacements are given with respect to image origin (0,0)
    radious -> radious of circular mask. Default = 0. In this case, the code uses the thrd to create the mask
    thrd -> threshold above which pixels are valid. Default = 0.05 (assuming image is normalized to one)
    iter -> maximum number of iterations in the kll algorithm.  Default = 15
    expand -> how much the circular mask is expanded (positive = schrinks the mask) in %.
    verbose -> 0 = do nothing.
    method = 'kll'
        'kll': 
        'chae': 
        'alter': 
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

    D = np.log10(image) 
    # replace NaNs and Infs by 0
    D[np.isneginf(D)] = 0
    D[np.isnan(D)] = 0

    if method == 'kll':

        ############################
        # CALCULATE CONSTANT
        ############################

        n = np.zeros([imsize[0], imsize[1]], dtype=np.float64)
        sum_image = np.zeros([imsize[0], imsize[1]], dtype=np.float64)

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
        
        K = sum_image / n.astype(np.float64) 
        # replace NaNs and Infs by 0
        K[np.isneginf(K)] = 0
        K[np.isnan(K)] = 0
        if verbose == 1:
            plt.imshow(K,cmap='gray')
            plt.clim(vmax=0.02,vmin=-0.02)
            plt.colorbar()
            plt.show()

        G = np.copy(K)
        if bit_trun == 1:
            K = np.int32(K * 256) / 256  # bit truncation
        k = np.power(10, K)
        if bit_trun == 1:
            k = np.int32(k * 256) / 256  # bit truncation

        for itera in range(iter):
            r_res = np.zeros([2048, 2048], dtype=np.float64)
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

            G = K + r_res / n.astype(np.float64) 
            # replace NaNs and Infs by 0
            G[np.isneginf(G)] = 0
            G[np.isnan(G)] = 0

            idx = np.where(n > 0)
            s = G[idx]
            #  calculate average of gain table for normalization
            sm = np.mean(s)
            sm2 = np.mean(s**2)
            five_sigma = 5*np.sqrt(sm2-sm*sm)
            idx2 = np.where(np.abs(s-sm) < five_sigma)
            sm = np.mean(s[idx2])
            G[idx] = G[idx] - sm

            print('Iteration: ', itera, five_sigma, '5*rms', sm, ' of ', iter)
            if verbose == 1:
                plt.imshow(G, cmap='gray', vmin=-0.05, vmax=0.05)
                plt.colorbar()
                plt.show()

        g = np.power(10, G,dtype='d')#/np.log(10,dtype='d')  + 0.5672334407 - 0.0018501610250685886#0.566 #exp(2.303)
        g[np.isneginf(g)] = 0
        g[np.isnan(g)] = 0

        return g

    elif method == 'chae':

        tmask = np.sum(mask,axis=0)
        mask_Ob,dummy = generate_circular_mask([imsize[0]-1, imsize[1]-1],radious,radious)  

        #Constant term
        if c_term == 0:
            fit_c = 0
            c_term = np.log10(np.ones((n_images)))
        else:
            c_term = np.log10(np.ones((n_images)))
        
        flat = np.log10(np.ones_like(D[0]))
        Ob = np.zeros_like(D[0])

        for i in range(n_images):
            # shift input image to the center of the frame.
            # 
            Ob += shift(D[i], shift = -xyshifts[i,:])
        Ob = Ob / float(n_images)
        idx = tmask >= 1

        for k in range(iter):

            numerator = np.zeros((imsize))
            for i in range(n_images):
                numerator += ((c_term[i] + Ob - shift(D[i] - flat, shift = -xyshifts[i,:]))*mask_Ob) #(i+xk,j+yk) Eq 8

            Ob -= (numerator/mask_Ob/n_images)
            Ob[np.isneginf(Ob)] = 0.
            Ob[np.isnan(Ob)] = 0.

            if verbose == 3:
                plt.imshow(Ob,cmap='gray',vmin=1,vmax=2)
                plt.show()
            
            numerator = np.zeros((imsize))
            for i in range(n_images):
                dummy = (c_term[i] + shift(Ob, shift = +xyshifts[i,:]) + flat - D[i])*mask[i]
                numerator += dummy
                c_term[i] -= ( np.sum(dummy) / np.sum(mask[i]) )
            
            dummy = (numerator/tmask)
            flat -= dummy
            flat[np.isneginf(flat)] = 0.
            flat[np.isnan(flat)] = 0.

            if verbose == 3:
                plt.imshow(flat,cmap='gray',vmin=-0.02,vmax=0.02)
                plt.show()
            if verbose >= 1:
                print('Iter: ',k, ' STD: ',np.max(np.abs(dummy[idx])),np.exp(c_term))

            s = flat[idx]
            sm = np.mean(s)
            sm2 = np.mean(s**2)
            five_sigma = 5*np.sqrt(sm2-sm*sm)
            print('Iteration: ', k, five_sigma, '5*rms', sm, ' of ', k, ' STD: ',np.max(np.abs(dummy[idx])))

        flat = flat - np.mean(flat)
        Ob = Ob + np.mean(flat) + np.mean(c_term)
        c_term = c_term - np.mean(c_term)

        flat = np.power(10, flat,dtype='d')
        flat[np.isneginf(flat)] = 0
        flat[np.isnan(flat)] = 0

        if verbose >= 2:
            plt.imshow(flat,cmap='gray',vmin=0.95,vmax=1.05)#vmin=)np.min(,vmax=0.05)
            plt.show()

        return flat
    
    elif method == 'alter':

        #Extracting flat-field images from scene-based image sequences using phase
        #correlation
        #James N. Caron, Marcos J. Montes, and Jerome L. Obermark
        #Citation: Review of Scientific Instruments 87, 063710 (2016); doi: 10.1063/1.4954730
        #View online: https://doi.org/10.1063/1.4954730
        #View Table of Contents: http://aip.scitation.org/toc/rsi/87/6
        #Published by the American Institute of Physics

        Gf = np.zeros([imsize[0], imsize[1]],dtype=np.float64)
        n = np.zeros([imsize[0], imsize[1]],dtype=np.float64)
        sum_image = np.zeros([imsize[0], imsize[1]],dtype=np.float64)
        Gr = np.zeros([n_images,imsize[0], imsize[1]],dtype=np.float64)

        iter = 1
        for itera in range(iter):
            ir = 0
            for iq in range(n_images):
                dx = rel_centers[iq,0] - rel_centers[ir,0] #shift of iq with respect ir
                dy = rel_centers[iq,1] - rel_centers[ir,1]
                t_mask = mask[ir] * shift(mask[iq], [dx,dy]) 
                n += t_mask
                t_image = shift(D[iq] - Gf, [dx,dy]) 
                sum_image += t_mask * t_image

            K = sum_image / n
            K[np.isneginf(K)] = 0
            K[np.isnan(K)] = 0
            idx = np.where(K > 5)
            mask[0][idx] = 0
            idx2 = np.where(n == 0)
            K[idx2] = 0

            iq = 0
            for ir in range(n_images):
                dx = rel_centers[iq,0] - rel_centers[ir,0] #shift of iq with respect ir
                dy = rel_centers[iq,1] - rel_centers[ir,1]
                Kn = shift(K * mask[0], [dx,dy])

                G = (D[ir] - Kn) * mask[ir]
                G[np.isneginf(G)] = 0
                G[np.isnan(G)] = 0
                Gr[ir,:,:] = G

            m = np.sum(mask,axis=0)
            Gf = np.sum(Gr,axis=0) / m
            Gf[np.isneginf(K)] = 0
            Gf[np.isnan(K)] = 0
            print('Iteration: ', itera)

        g = np.power(10, Gf)
        g[np.isneginf(g)] = 0
        g[np.isnan(g)] = 0

        return g
    else:
        return None

def fdt_flat(files, wavelength, npol, method, dark = None, read_shits = 0, shifts = None, verbose = 1,
    correct_ghost = 0,expand = 1,thrd = 0,iter = 4, normalize = 1 , disp_method = 'Hough', c_term = 0,
    inner_radius = 400, outer_radius = 800, steps = 20):
    '''
    The Dark, if provided, should have the same scaling as the data and same size!!!!!!!
    This program does not take care of sizes. For that go to fdt_pipeline
    '''

    ############################
    # open the 9 FITS images and get one pol and wave. 
    # This can be done all at once but I have this like that because I am lazy.
    ############################
    
    image = [fits_get_part(i,wavelength,npol) for i in files]
    n_images = len(image)
    ys,xs = image[0].shape
    ############################
    # Correct Dark if not done in advance
    ############################

    try: 
        print('...Dark correction...')
        for i in range(n_images):
            image[i] = image[i] - dark
    except:
        pass

    #TODO To be implemented (detailes dark correction)
    # try:
    #     for i in range(n_images):
    #         row = image[i][1800:,:].sum(axis=0) / (2048.-1800.)
    #         dark2 = np.zeros_like(dark)
    #         dark2 = dark2 + row[:np.newaxis]
    #         image[i] = np.abs(image[i] - dark2)
    # except:
    #     pass
    
    if read_shits == 1:
        try:
            print('... read user input shifts_file ...')
            centers = read_shifts(shifts+'_cnt_w'+str(wavelength)+'_n'+str(npol)+'.txt')
            radius  = read_shifts(shifts+'_rad_w'+str(wavelength)+'_n'+str(npol)+'.txt')
            for i in range(n_images):
                print('Image',i,'c: ',centers[i,0],',',centers[i,1],' rad: ', radius[i])
        except Exception:
            print("Unable to open fits file: {}",shifts+'_cnt_w'+str(wavelength)+'_n'+str(npol)+'.txt')        
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
            _ = write_shifts(shifts+'_cnt_w'+str(wavelength)+'_n'+str(npol)+'.txt', centers)
            _ = write_shifts(shifts+'_rad_w'+str(wavelength)+'_n'+str(npol)+'.txt', radius )

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
            centers = np.array([n_images,2])
            radius = np.array([n_images])
            for i in range(n_images):
                centers[i,1],centers[i,0],radius[i] = find_center(image[i],sjump = 4,njumps = 100,threshold = 0.8)
                print, 'Image',i,'c: ',centers[i,0],',',centers[i,1],' rad: ', radius[i]
        else:
            pass
    #make sure we have integer numpy numbers in the centers 
    centers = np.array(centers).astype(int)
    mean_radii = np.mean(radius)

    if correct_ghost == 1:
        coef = [-1.98787669,1945.28944245]
        print(' Ghost corrrection...')
        poly1d_fn = np.poly1d(coef)
        sh = poly1d_fn(centers[4,:]).astype(int) #np.array([ -1.99350209*centers[4,0] + 1948.44866543,-1.98963222*centers[4,1] + 1949.61650596]).astype(int)
        reflection = image[4] - shift(image[4], shift=sh) * 0.004
        reflection = shift(reflection, shift=[-centers[4,1]+1024,-centers[4,0]+1024])
        for i in range(9):
            sh = poly1d_fn(centers[i,:]).astype(int) 
            image[i] = image[i] - shift(reflection, shift=[sh[1]+centers[i,1]-1024,sh[0]+centers[i,0]-1024]) * 0.004

    #PROBAR CON TODA LA IMAGES TODO 1
    for i in range(n_images):
        norma = np.mean(image[i][centers[i,1]-100:centers[i,1]+100,centers[i,0]-100:centers[i,0]+100])
        if normalize == 1:
            image[i] = image[i]/norma
            print('Normalization: ', norma)
        else:
            norma = 0
        pass

    if thrd != 0:
        gain = fdt_flat_gen(image, centers,method,iter=iter,thrd=thrd,verbose = verbose, c_term = c_term)
    else:
        gain = fdt_flat_gen(image, centers,method,iter=iter,radious=mean_radii,expand=expand,verbose = verbose, c_term = c_term)

        return gain, norma


import sys
sys.path.append('../SPGPylibs/')
import SPGPylibs as spg
def fdt_flat_testrun():
    '''
    Just for local test
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

    dark_file = '/Users/orozco/Dropbox_folder/SoPHI/PHI-COMMISSIONING/software-and-images/RSW1/solo_L0_phi-fdt-ilam_20200618T000547_V202006221044C_0066181001_dark.fits'
    dark, _ = fits_get(dark_file)
    
    scaling_dark = fits_get(dark_file,scaling = True)
    scaling_flat = fits_get(files[0],scaling = True)
    dark = dark * scaling_flat / scaling_dark

    wavelength = 0
    npol = 0

    allgain = []
    norma = np.zeros((24))
    # for wavelength in range(6):
      # for npol in range(4):
    for wavelength in range(1):
      for npol in range(1):
        print(wavelength,npol,'................')
        gain, norma_out = fdt_flat(files, wavelength, npol, 'kll',dark = dark,read_shits = False, 
            shifts_file = 'shifts/shifts', correct_ghost = 0 , expand = 10, normalize = 0, 
            iter = 3, method = 'kll',verbose=0)
        # gain, norma_out = fdt_flat(files,wavelength,npol,'kll',dark=dark,read_shits = False, shifts_file = ' '
        #     correct_ghost=0 , expand = 10, normalize=0,thrd = 0.2, iter = 3, method = 'kll',verbose=0)
        allgain.append(gain)
        norma[wavelength*4+npol] = norma_out

    np.savetxt('flats_norma.txt', norma)
    
    header = getheader(files[0]) 

    hdu = pyfits.PrimaryHDU(allgain)
    hdu = pyfits.PrimaryHDU(header=header)
    hdu.writeto('flats.fits', clobber=True, output_verify='fix')
