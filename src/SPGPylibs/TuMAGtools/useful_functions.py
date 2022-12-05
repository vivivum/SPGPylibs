# ============================ IMPORTS ====================================== #

import os
import numpy as np
import glob
import matplotlib.pyplot as plt

from alive_progress import alive_bar
import datetime as dt

from astropy.io import fits

# Own Libs
import utils as ut

# ============================= CONFIG ====================================== #

plt.style.use('default')

# =========================================================================== #

def Read_Image(path, PlotFlag = False, printHeader_flag = False, vmin = 0, vmax = 4000):
    """
    PlotFlag: imshow de la imagen si True
    printHeader_flag: saca el Header por terminal si True
    vmin, vmax: límites del imshow si PlotFlag es True
    """

    # Separate Header and image data
    H, hl, I = ut.HeadernImageSeparator(path) # Header, header_lenth and image

    # Read header
    Head = ut.GetDatafromHeader(H)

    # Read Image
    Im  = ut.read_image(path, Head['Roi_x_size'], Head['Roi_y_size'], hl)

    if PlotFlag:
        plt.figure()
        im = plt.imshow(Im.T, cmap = 'inferno', vmin = vmin, vmax = vmax)
        plt.colorbar(im)

    if printHeader_flag:
        for key in Head:
            print(key, ' : ', Head[key])

    return Head, Im

def Read_Header(path, printHeader_flag = False):
    """
    Lo mismo que la función anterior pero sólo lee la cabecera.
    """

    # Separate Header and image data
    H, hl, I = ut.HeadernImageSeparator(path) # Header, header_lenth and image

    # Read header
    Head = ut.GetDatafromHeader(H)

    # Read Image
    if printHeader_flag:
        for key in Head:
            print(key, ' : ', Head[key])

    return Head

def folder_image_reader(path):
    """
    Si le pasas el path de una carpeta, lee todas las imágenes
    y las cabeceras y los pone en un array. El primer índice es el
    índice del número de imágenes
    """

    Images_paths = sorted(glob.glob(os.path.join(path, '*.img')))

    IMAGES = []
    DATA = []

    print('Reading headers in folder:', path)
    with alive_bar(len(Images_paths)) as bar:
        for img in Images_paths:
            H , I = Read_Image(img)
            IMAGES.append(I)
            DATA.append(H)
            bar()

    return DATA, np.array(IMAGES)


def folder_header_reader(path):
    """
    Lo mismo que la función anterior pero sólo con la cabecera.
    """

    Images_paths = sorted(glob.glob(os.path.join(path, '*.img')))

    DATA = {}

    print('Reading headers in folder:', path)
    with alive_bar(len(Images_paths)) as bar:
        for img in Images_paths:
            image_name = os.path.basename(img)
            H = Read_Header(img)
            DATA[image_name] = H
            bar()

    return DATA

def voltage_vs_time(path, PlotFlag = True):
    """
    Al pasarle el path de una carpeta, pinta el voltaje comandado y leído
    frente al tiempo. También lo devuelve en 3 arrays.
    """

    Comm = []
    Time = []
    Read = []

    Images = sorted(glob.glob(os.path.join(path, '*')))

    Nimag = len(Images)
    print('Reading headers in folder:', path)
    with alive_bar(Nimag) as bar:
        for img in Images:

            H  = Read_Header(img)

            C = (-1) ** (int(H['EtalonSign']) + 1) \
                        * (int(H['EtalonDN']) * 4999 / (2 ** 12 - 1))

            R = 5000 * ( ( ( 2 * H['EtalonVoltsReading'] ) / 4095 ) - 1 )

            col = os.path.basename(img).split('_')
            T = dt.datetime(year = int(col[0]),
                                      month = int(col[1]),
                                      day = int(col[2]),
                                      hour = int(col[3]),
                                      minute = int(col[4]),
                                      second = int(col[5]),
                                      microsecond = 1000 * int(col[6]))

            Comm.append(C)
            Read.append(R)
            Time.append(T)

            bar()

    if PlotFlag:
        t0 = Time[0]
        TT = [(x -t0).total_seconds() for x in Time]

        plt.figure()
        plt.plot(TT, Comm, color = 'deeppink', lw = 5, label = 'Commanded')
        plt.plot(TT, Read, color = 'k', lw = 2, label = 'Read', zorder = 10)
        plt.legend()
        plt.xlabel('Seconds since first image')
        plt.ylabel('Volts [V]')
        plt.grid(True)
        plt.show()

    return Comm, Read, Time #Comandado, leído y tiempo


def voltage_scan(path):
    """
    Si le pasas el path de una carpeta, ordena las imágenes por voltaje
    y promedia todas las imágenes para cada voltaje. Devuelve dos listas.
    """

    Images = sorted(glob.glob(os.path.join(path, '*')))

    Nimag = len(Images)
    print('Reading headers in folder:', path)

    Volts     = []
    Intensity = []

    with alive_bar(Nimag) as bar:
        for img in Images:

            H, I = Read_Image(img)

            size = H['Roi_x_size']
            mid = int(size / 2)
            deltax=250

            Intensity.append(np.mean(I[mid - deltax : mid + deltax, mid\
             - deltax : mid + deltax]))

            C = (-1) ** (int(H['EtalonSign']) + 1) \
                        * (int(H['EtalonDN']) * 4999 / (2 ** 12 - 1))

            Volts.append(C)

            bar()

    return Volts, Intensity


def Dark_current(path):
    """
    Promedia todas las imágenes para cada cámara. Sirve para construir
    darks y flats. Si en la carpeta hay sólo imágenes de una cámara,
    saca una única variable.
    """

    Images_0 = glob.glob(os.path.join(path, '*_0_*'))
    Images_1 = glob.glob(os.path.join(path, '*_1_*'))

    print('Found images -> Cam 1 :', len(Images_0), 'Cam 2 :', len(Images_1))

    H = Read_Header(Images_0[0])
    size = H['Roi_x_size']

    dc0 = np.zeros((size, size))
    dc1 = np.zeros((size, size))

    if len(Images_0) > 0:

        print('Computing Darks for cam 1 :')
        with alive_bar(len(Images_0)) as bar:
            for img in Images_0:

                H, I = Read_Image(img)

                dc0 += I

                bar()
        dc0 /= len(Images_0)

    if len(Images_1) > 0:
        print('Computing Darks for cam 2 :')
        with alive_bar(len(Images_1)) as bar:
            for img in Images_1:

                H, I = Read_Image(img)

                dc1 += I

                bar()

        dc1 /= len(Images_1)

    if len(Images_0) > 0 and len(Images_1) > 0:
        return dc0, dc1

    elif len(Images_0) > 0:
        return dc0

    else:
        return dc1

# =========================================================================== #

pathfolder='D:/20220504/'#'C:/Users/Francisco J. Bailén/Dropbox (IdAdA)/Trabajo/TuMag/End to end tests/Datos e2e/Kiruna/'
pathsubfolder='Test5/517nm'
Folder = pathfolder+pathsubfolder

#C, R, T = voltage_vs_time(Folder)
V, I = voltage_scan(Folder)

plt.plot(V,I)
plt.savefig(Folder+'.png',dpi=480)
plt.show()
quit()

H, D = folder_image_reader(Folder)
#print(H[5])

print(np.mean(D[0,:,:]))
plt.imshow(D[0,:,:],cmap='gray')
plt.show()
quit()

#plt.figure()
#plt.plot(V, I, 'deeppink', lw = 3)


dc1 = Dark_current(Folder)
