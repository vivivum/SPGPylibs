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

from mpl_toolkits.axes_grid1 import make_axes_locatable

def show_row(im,title,svmin=0,svmax=0,xlabel='Pixel',ylabel='Pixel',zoom=1,block=True):

    ni = len(im)
    fig, maps = plt.subplots(1,ni,figsize=(18*zoom,6*zoom))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    for i in range(ni):
        plot = maps[i].imshow(im[i], cmap='viridis',vmin = svmin,vmax=svmax)
        maps[i].set_title(title[i])
        maps[i].set_ylabel(ylabel)
        maps[i].set_xlabel(xlabel)
        colorbar(plot)
    plt.show(block = block)
    return

#Function for color bar
def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

def Read_Image(path, PlotFlag = False, printHeader_flag = False, vmin = 0, vmax = 4000):
    """
    This function reads TuMAG image and header information
    
    :param Path: Path to the image file (including image name)
    :type Path: str
    :param PlotFlag: If True the image gets represented. Default value = False
    :type PlotFlag: Boolean
    :param printHeader_flag: If True, Header info will be printed out. Default = False
    :type printHeader_flag: Boolean

    :return: header, image.
    :rtype: Dict, np.ndarray(float)
        
    """
     
    # Separate Header and image data 
    H, hl, I = ut.HeadernImageSeparator(path) # Header, header_lenth and image
     
    # Read header
    Head = ut.GetDatafromHeader(H) 
    # Read Image
    Im  = ut.read_image(path, Head, hl) 
     
    if PlotFlag:
        plt.imshow(Im, cmap = 'Greys', vmin = vmin, vmax = vmax)
        plt.colorbar()
        plt.show()
         
    if printHeader_flag:
        for key in Head:
             print(key, ' : ', Head[key])
     
    return Head, Im

# =========================================================================== #

def Read_Image_example():
# Example of execution 

    Image_path = '2022_05_06_11_55_40_251_0_1.img'
    # Image_path = '2022_05_09_14_15_58_733_0_12801.img'
    H, I = Read_Image(Image_path,PlotFlag=True,vmin=120,vmax=140)
    # H  = ut.Read_Header(Image_path)
    Im, H = ut.Read_TuMAG(Image_path)

    plt.figure()
    im = plt.imshow(I.T, cmap = 'inferno', vmin = 100, vmax = 150)
    plt.colorbar(im)
    plt.show()
        
def Read_scan():
# Example of execution 

    Image_path = '/Volumes/New Volume/TuMAG/FW26_517/'
    import glob, os
    os.chdir(Image_path)
    files = glob.glob("*.img")
    im = []
    for f in files:
        print(Image_path+f)
        H, I = Read_Image(Image_path+f,PlotFlag=True)#,vmin=120,vmax=140)
        im.append(I)

    curva = np.zeros((len(im)))
    for j in range(len(im)):
        curva[j] = np.sum(im[j])

    plt.figure()
    plt.plot(curva)
    plt.show()

def test_light():
    # import Image_reader as tm
    import glob, os
    paths = ['/Volumes/New Volume/TuMAG/Telescope_heights/p0/',
    '/Volumes/New Volume/TuMAG/Telescope_heights/p1/',
    '/Volumes/New Volume/TuMAG/Telescope_heights/p2/',
    '/Volumes/New Volume/TuMAG/Telescope_heights/p3/',
    '/Volumes/New Volume/TuMAG/Telescope_heights/p4/',
    '/Volumes/New Volume/TuMAG/Telescope_heights/p5/']
    c = []
    images = []
    for i in paths:
        Image_path = i
        os.chdir(Image_path)
        files = sorted(glob.glob("*.img"))

        curva = []
        l = 0
        for f in files:
            l += 1
            print(Image_path+f)
            H, I = Read_Image(Image_path+f)#,PlotFlag=True,vmin=120,vmax=1000)
            #images.append(I)
            curva.append(np.sum(I[1024-128:1024+128,1024-128:1024+128]))
            if l == 10:
                images.append(I)
        c.append(curva)

    plt.plot(c[0],label='0º')
    plt.plot(c[1],label='9.18º')
    plt.plot(c[2],label='17.93º')
    plt.plot(c[3],label='26.98º')
    plt.plot(c[4],label='36.25º')
    plt.plot(c[5],label='45.13º')
    plt.xlabel('position')
    plt.ylabel('sum')
    plt.legend()
    plt.show()

    level = np.mean(images[0])
    # plt.imshow((images[0]-images[1])/level*100,vmin=-3,vmax=3,label='0º - 9.10º')
    # plt.colorbar()
    # plt.show()
    # plt.imshow((images[0]-images[2])/level*100,vmin=-30,vmax=30,label='0º - 17.93º')
    # plt.colorbar()
    # plt.show()
    # plt.imshow((images[0]-images[3])/level*100,vmin=-30,vmax=30,label='0º - 26.98º')
    # plt.colorbar()
    # plt.show()
    # plt.imshow((images[0]-images[4])/level*100,vmin=-30,vmax=30,label='0º - 36.25º')
    # plt.colorbar()
    # plt.show()
    # plt.imshow((images[0]-images[5])/level*100,vmin=-30,vmax=30,label='0º - 45.13º')
    # plt.colorbar()
    # plt.show()

    show_row(((images[0]-images[1])/level*100,(images[0]-images[2])/level*100,(images[0]-images[3])/level*100
        ,(images[0]-images[4])/level*100,(images[0]-images[5])/level*100)
        ,('0º - 9.10º','0º - 17.93º','0º - 26.98º','0º - 36.25º','0º - 45.13º'),svmin=-4,svmax=4)

    return curva

def Read_grid_data():
# Example of execution 

    Image_path = '/Volumes/New Volume/TuMAG/Sun_obs/step6/'
    import glob, os
    os.chdir(Image_path)
    files = sorted(glob.glob("*.img"))
    print(len(files))
    Mg = 0
    Fe = -1
    index = 0
    files = files[0:(10 + 8) * 1]
    n_ciclos = int(len(files)/(10 + 8))
    ImMg = np.zeros((2048,2048,10,n_ciclos))
    ImFe = np.zeros((2048,2048,10,n_ciclos))
    for f in files:
        print(f,len(f),Mg,Fe,index)
        H, I = Read_Image(Image_path+f)#,PlotFlag=True)#,vmin=120,vmax=140)
        print(H)
        if Mg >= 0:
            #Mg case
            ImMg[:,:,Mg,index//(10+8)] += I
            Mg += 1
            if Mg > 9:
                Mg = -1
                Fe = 0
        if Fe >= 0:
            #Mg case
            ImFe[:,:,Fe,index//(10+8)] += I
            Fe += 1
            if Fe > 7:
                Mg = 0
                Fe = -1
        index += 1
        print(index//(10+8))


Read_grid_data()
#test_light()
# Read_scan()
# Read_Image_example()



