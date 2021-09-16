from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

plt.style.use('seaborn-white')
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (5, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'10',
         'ytick.labelsize':'10'}
pylab.rcParams.update(params)
PLT_RNG = 3

#Function for color bar
def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

def show_one(img,vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Image no title',cbarlabel='Some units',save=False,cmap='gray'):

    plt.figure(figsize=(6, 6))
    ax = plt.gca()
    if vmin == None and vmax == None:
        im = ax.imshow(img, cmap=cmap,vmin=img.mean() - PLT_RNG * img.std(),
           vmax=img.mean() + PLT_RNG * img.std(), interpolation='none')
    elif vmin == None:
        im = ax.imshow(img, cmap=cmap,vmin=img.mean() - PLT_RNG * img.std(),
           vmax=vmax, interpolation='none')
    elif vmax == None:
        im = ax.imshow(img, cmap=cmap,vmin=vmin,
           vmax=img.mean() + PLT_RNG * img.std(), interpolation='none')
    else:
        im = ax.imshow(img, cmap=cmap,vmin=vmin,
           vmax=vmax, interpolation='none')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(cbarlabel)
    if save:
        plt.savefig(save,dpi=300)
        plt.close()
    else:
        plt.show()

    return

def show_two(im1,im2,vmin=[None,None],vmax=[None,None],block=True,pause=0.1,title=['',''],xlabel='Pixel',ylabel='Pixel'):
    
    fig, maps = plt.subplots(1,2,figsize=(8,8))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)

    if vmin[0] == None and vmax[0] == None:
        vmin[0]=im1.mean() - PLT_RNG * im1.std()
        vmax[0]=im1.mean() + PLT_RNG * im1.std()
    if vmin[1] == None and vmax[1] == None:
        vmin[1]=im2.mean() - PLT_RNG * im2.std()
        vmax[1]=im2.mean() + PLT_RNG * im2.std()

    im = maps[0].imshow(im1, cmap='gray',vmin=vmin[0],vmax=vmax[0])
    maps[0].set_title(title[0])
    maps[0].set_xlabel(xlabel)
    maps[0].set_ylabel(ylabel)
    colorbar(im)

    im = maps[1].imshow(im2, cmap='gray',vmin=vmin[1],vmax=vmax[1])
    maps[1].set_title(title[1])
    maps[1].set_xlabel(xlabel)
    colorbar(im)

    plt.show(block=block)
    plt.pause(pause)
    plt.close()
    return

def show_three(im1,im2,im3,vmin=[None,None,None],vmax=[None,None,None],block=True,pause=0.1,title=['','',''],xlabel='Pixel',ylabel='Pixel',save=False,cmap='gray'):

    #vmin = [None,None,None]
    #vmax = [None,None,None]

    fig, maps = plt.subplots(1,3,figsize=(12,5))
    #fig, maps = plt.subplots(1,3, sharex=True, sharey=True,figsize=(15,5))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)

    if vmin[0] == None and vmax[0] == None:
        vmin[0]=im1.mean() - PLT_RNG * im1.std()
        vmax[0]=im1.mean() + PLT_RNG * im1.std()
    if vmin[1] == None and vmax[1] == None:
        vmin[1]=im2.mean() - PLT_RNG * im2.std()
        vmax[1]=im2.mean() + PLT_RNG * im2.std()
    if vmin[2] == None and vmax[2] == None:
        vmin[2]=im3.mean() - PLT_RNG * im3.std()
        vmax[2]=im3.mean() + PLT_RNG * im3.std()

    im = maps[0].imshow(im1,vmin=vmin[0],vmax=vmax[0], interpolation='none',cmap=cmap)
    maps[0].set_title(title[0])
    maps[0].set_xlabel(xlabel)
    maps[0].set_ylabel(ylabel)
    colorbar(im)

    im = maps[1].imshow(im2,vmin=vmin[1],vmax=vmax[1], interpolation='none',cmap=cmap)
    maps[1].set_title(title[1])
    maps[1].set_xlabel(xlabel)
    colorbar(im)

    im = maps[2].imshow(im3,vmin=vmin[2],vmax=vmax[2], interpolation='none',cmap=cmap)
    maps[2].set_title(title[2])
    maps[2].set_xlabel(xlabel)
    colorbar(im)

    if save:
        plt.savefig(save,dpi=300)
        plt.close()
    else:
        plt.show(block=block)
        plt.pause(pause)
        plt.close()
    return

def show_four_row(im1,im2,im3,im4,svmin=0,svmax=0,title=['','','',''],xlabel='Pixel',ylabel='Pixel',save=False,zoom=1,block=True):

    fig, maps = plt.subplots(1,4,figsize=(12*zoom,6*zoom))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    yd,xd = im1.shape
    for i in range(4):
        if i == 0:
            dummy = np.copy(im1)
        if i == 1:
            dummy = np.copy(im2)
        if i == 2:
            dummy = np.copy(im3)
        if i == 3:
            dummy = np.copy(im4)
        try:
            if len(svmin) == 1 or len(svmax) == 1:
                vmin = dummy.mean() - PLT_RNG * dummy.std()
                vmax = dummy.mean() + PLT_RNG * dummy.std()
            else:
                vmin = svmin[i]
                vmax = svmax[i]
        except:
                vmin = np.min(dummy[yd//2-50:yd//2+50,xd//2-50:xd//2+50])
                vmax = np.max(dummy[yd//2-50:yd//2+50,xd//2-50:xd//2+50])
                lim = np.max([np.abs(vmin),vmax])
                vmin = -lim
                vmax = lim

        dummy[dummy<vmin] = vmin
        dummy[dummy>vmax] = vmax
        im = maps[i].imshow(dummy, cmap='gray',vmin = vmin,vmax=vmax)
        maps[i].set_title(title[i])
        if i == 0:
            maps[i].set_ylabel(ylabel)
        maps[i].set_xlabel(xlabel)
        colorbar(im)

    if save != False:
        plt.savefig(save)
        plt.clf()
        return
    plt.show(block = block)

    return

def show_hist(x,bins=40,title=' ',xlabel='',ylabel='',leave='close',color='green'):

    n, bins, patches = plt.hist(x, bins, facecolor=color, alpha=0.75)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    if leave == 'close':
        plt.show()

def show_six_row(image,plrt = 0,title='n/a',vmin=None,vmax=None):

    fig, maps = plt.subplots(1,6,figsize=(24,4))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    for i in range(6):
        el_mean = image[i][:,:].mean() 
        el_std = image[i][:,:].std()
        el_min = el_mean - (PLT_RNG + plrt) * el_std 
        el_max = el_mean + (PLT_RNG + plrt)  * el_std
        try:
            im = maps[i].imshow(image[i][:,:],vmin=vmin,vmax=vmax)
        except:
            im = maps[i].imshow(image[i][:,:],vmin=el_min,vmax=el_max)
        maps[i].set_title(title+'{: d}'.format(i))

        colorbar(im)
    plt.show()
    return
def squar(n):
    column = np.int(np.sqrt(n))
    # if column <= n:
    #     row = column + 1
    #     if row*column <= n:
    #         column = column + 1
    # else:
    row = column
    return row,column

def show_all(image,save=False):
    ishape = image.shape
    row,column = squar(ishape[2])
    print(row,column,ishape)
    fig, maps = plt.subplots(row,column, sharex='col', sharey='row',figsize=(12,12))
    plt.subplots_adjust(top=0.92)
    for i in range(ishape[2]):
        #themedian = median(sp[:,:,i])
        #image = image_histogram_equalization(sp[:,:,i])[0]
        im = maps[divmod(i, column)].imshow(image[:,:,i],\
        cmap='gray',vmin=image[:,:,i].mean() - PLT_RNG * image[:,:,i].std(),
               vmax=image[:,:,i].mean() + PLT_RNG * image[:,:,i].std(),\
                interpolation='none')
        colorbar(im)
    if save != False:
        plt.savefig(save)
        plt.clf()
        return
    plt.show()
    return

def doplots(im1,im2,im3):

    font = {'family' : 'serif','weight' : 'normal','size'   : 6}
    plt.rc('font', **font)
    plt.rc('figure',figsize=[8,8])
    fig, maps = plt.subplots(2, 2, sharex='col', sharey='row')
    plt.subplots_adjust(hspace=0.3, wspace=0.3, top=0.92)
    im = maps[0,0].imshow(im1, cmap='gray')
    maps[0,0].set_title('Restored image')
    colorbar(im)
    im = maps[0,1].imshow(im2, cmap='gray')
    maps[0,1].set_title('Original image')
    colorbar(im)
    im = maps[1,0].imshow(im1-im2, cmap='gray',vmax=0.05,vmin=-0.05)
    maps[1,0].set_title('|Original-restored|')
    colorbar(im)
    im = maps[1,1].imshow(im3, cmap='gray',vmax=0.05,vmin=-0.05)
    maps[1,1].set_title('Defocused')
    colorbar(im)
#    plt.savefig('resultados'+str(i)+'.png')
    plt.show()
    plt.close()
    return

def twoimages(im1,im2,vmin=[None,None],vmax=[None,None],block=True,pause=0.1):

    fig, maps = plt.subplots(1,2)
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    im = maps[0].imshow(im1, cmap='gray',vmin=vmin[0],vmax=vmax[0])
    maps[0].set_title('im1')
    colorbar(im)
    im = maps[1].imshow(im2, cmap='gray',vmin=vmin[1],vmax=vmax[1])
    maps[1].set_title('im2')
    colorbar(im)
#    plt.savefig('resultados'+str(i)+'.png')
    plt.show(block=block)
    plt.pause(pause)
    plt.close()
    return

def fourimages(im1,im2,im3,im4,vmin=[None,None,None,None],vmax=[None,None,None,None],block=True,pause=0.1):

    fig, maps = plt.subplots(2,2)
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    im = maps[0,0].imshow(im1, cmap='gray',vmin=vmin[0],vmax=vmax[0])
    maps[0,0].set_title('im1')
    colorbar(im)
    im = maps[1,0].imshow(im2, cmap='gray',vmin=vmin[1],vmax=vmax[1])
    maps[1,0].set_title('im2')
    colorbar(im)
    im = maps[0,1].imshow(im3, cmap='gray',vmin=vmin[2],vmax=vmax[2])
    maps[0,1].set_title('im3')
    colorbar(im)
    im = maps[1,1].imshow(im4, cmap='gray',vmin=vmin[3],vmax=vmax[3])
    maps[1,1].set_title('im4')
    colorbar(im)
#    plt.savefig('resultados'+str(i)+'.png')
    plt.show(block=block)
    plt.pause(pause)
    plt.close()
    return

def images3x3(image):
    fig, maps = plt.subplots(3, 3, sharex='col', sharey='row',figsize=(8,8))
    plt.subplots_adjust(top=0.92)
    for i in range(9):
        #themedian = median(sp[:,:,i])
        #image = image_histogram_equalization(sp[:,:,i])[0]
        im = maps[divmod(i, 3)].imshow(image[:,:,i],cmap='gray')
        colorbar(im)
    plt.show()
    plt.close()
    return

def two_plot(x,d1,d2):
    ###def array(*args, **kwargs):

    plt.plot(x,d1)
    plt.plot(x,d2)
    plt.show()

    return

def i_h_e(image, number_bins=256):
    # from http://www.janeriksolem.net/2009/06/histogram-equalization-with-python-and.html

    # get image histogram
    image_histogram, bins = np.histogram(image.flatten(), number_bins, normed=True)
    cdf = image_histogram.cumsum() # cumulative distribution function
    cdf = 255 * cdf / cdf[-1] # normalize

    # use linear interpolation of cdf to find new pixel values
    image_equalized = np.interp(image.flatten(), bins[:-1], cdf)

    return image_equalized.reshape(image.shape), cdf

