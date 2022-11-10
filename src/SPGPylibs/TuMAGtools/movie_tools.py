import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FFMpegWriter, FuncAnimation
plt.rcParams['animation.ffmpeg_path']='/Users/orozco/Downloads/ffmpeg'

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

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

def make_movie(data,filename,figsize=[9,9],fps=15,cbar='no',cmap='gray',vmin=0,vmax=0,title='title',xtitle='Pixels',ytitle='Pixels'):

        import matplotlib.animation as nanim 

        FFMpegWriter = nanim.writers['ffmpeg']
        metadata = dict(title='Movie', artist='dos',comment='Movie support!')
        writer = FFMpegWriter(fps = fps, metadata=metadata,bitrate=1800)

        plt.rc('lines', linewidth=2, color='r')
        plt.rc('figure', figsize=figsize)

        font = {'family' : 'serif',
                'weight' : 'normal',
                'size'   : 10}
        plt.rc('font', **font)
        fig  = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_title(title)
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)

        if (vmin != 0) and (vmax != 0):
                im = ax.imshow(data[:,:,0],clim=(vmin,vmax), interpolation = 'none',cmap=cmap)
        else:
                im = ax.imshow(data[:,:,0], interpolation = 'none',cmap=cmap)
        if cbar:
                plt.colorbar(im)
        yd,xd,n = data.shape
        with writer.saving(fig, filename, 100):
                for i in range(n-1):
                        im.set_data(data[:,:,i+1],)
                        if (vmin != 0) and (vmax != 0):
                                im.set_clim = ([vmin,vmax])
                        writer.grab_frame()
                        plt.close()

def usage():
    make_movie(d,'blos.mp4',figsize=[9,9],fps=2,cbar=True,cmap='Greys',vmin=-40,vmax=40,title='blos')