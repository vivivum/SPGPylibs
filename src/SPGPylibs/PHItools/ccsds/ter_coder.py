# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import os
import sys
import numpy as np
from astropy.io import fits as pyfits
from dataclasses import dataclass

_print = print
def print(*args, **kw):
    _print('--> ',  *args, **kw)

# %%
@dataclass
class ter:
   ''' data type is fixed to uint16 or int16 '''
   TER_coder_sw_path: str = 'TER_iaa_coder.jar'
   TER_decoder_sw_path: str = 'TER_iaa_decoder.jar'
   input_folder_path: str = './'
   code_ext: str = '.rec'
   imageName_input: str = 'inputImage'
   imageName_output_coded: str = 'outputImage_coded'
   imageName_output_decoded: str = 'outputImage_decoded'
   image_width: int = 2048
   image_height: int = 2048
   endianess: str = 'b'          # 'b' = big endian; 'l' = little endian
   target_bpp: int = 0           # 0 = lossless
   dc_stop: int = 0
   bitplaneStop: int = 0 #
   #### --customWtFlag
   # - bitplaneStop = 0 --> lossless
   # - bitplaneStop = 23 --> super lossy (hasta 8)
   stage_stop: int = 4 #number of bit planes 0-4
   output_folder_path: str = './'
   decode_ext: str = '.raw'
   scale_value: float = (0,0)
   io_format: str = '<i2'   #'<i2' #depends on signedPixels value and endianess
   #signedPixels: int = 1         # 1 = signed; 0 = unsigned
   _signedPixels: int = 1
   verbose: int = 1
   
# Para recorrer de mejor a peor calidad de imagen descomprimida, las pruebas tendrían que ir en un bucle anidado tal que así:
# for bitplaneStop=0:1:23
#     for stage_stop=4:-1:0
#         do_compression(bitplane, stage_stop);
#     end for;
# end for;

#    print(ter.__dict__.keys()) 
#    print(ter.__dict__.values()) 

   @property
   def signedPixels(self) -> int:
       return self._signedPixels
   @signedPixels.setter
   def signedPixels(self, value):
       if value == 1:
           fmt = 'i2'
       if value == 0:
           fmt = 'u2'
       if self.endianess == 'b':
           fmt = '>'+fmt
       if self.endianess == 'l':
           fmt = '<'+fmt
       self.io_format = fmt   #'<i2' #depends on signedPixels value and endianess
       self.pinfo('setting i/o: ',fmt)
       self._signedPixels = int(value)


     #'<u2') #uint16
     #'<i2') #int16

   def pinfo(self,*args, **kw):
       if self.verbose:
         print(*args, **kw)
   def set_format(self):
       ## condition to check whether var is suitable or not
       if self.signedPixels == 1:
         fmt = 'i2'
       if self.signedPixels == 0:
         fmt = 'u2'
       if self.endianess == 'b':
         fmt = '>'+fmt
       if self.endianess == 'l':
         fmt = '<'+fmt
       self.io_format = fmt
       return  
     
   def save_raw(self,data):
       self.set_format()
       self.pinfo('saving... '+self.input_folder_path+self.imageName_input+'.raw')
       data.astype(self.io_format).tofile(self.input_folder_path+self.imageName_input+'.raw')
       return

   def read_output(self):
       #formatis = '<i2'
       self.pinfo('reading... ',self.output_folder_path+self.imageName_output_decoded+'.raw')
       #return np.fromfile(self.output_folder_path+self.imageName_output_decoded+'.raw',dtype=formatis).reshape(self.image_width, self.image_height)
       return np.fromfile(self.output_folder_path+self.imageName_output_decoded+'.raw',dtype=self.io_format).reshape(self.image_width, self.image_height)

   def read_input(self):
       self.pinfo('reading... ',self.input_folder_path+self.imageName_input+'.raw')
       return np.fromfile(self.input_folder_path+self.imageName_input+'.raw',dtype=self.io_format).reshape(self.image_width, self.image_height)

   def code(self) -> None:
       ''' assumes there is already a file with the proper format'''

       # input_imageName_fullPath
       # get last character of input_folder_path
       last = self.input_folder_path[-1]
       if last != '/':
           input_imageName_fullPath = self.input_folder_path + '/'
       else:
           input_imageName_fullPath = self.input_folder_path

       input_imageName_fullPath = input_imageName_fullPath + self.imageName_input+'.raw'

       # output image name in full path format
       last = self.output_folder_path[-1]
       if last != '/':
           output_imageName_fullPath = self.output_folder_path + '/'
       else:
           output_imageName_fullPath = self.output_folder_path

       output_imageName_fullPath = self.output_folder_path + self.imageName_output_coded

       # TER compression parameters
       code_call = 'java -jar ' + self.TER_coder_sw_path + ' -i '  + input_imageName_fullPath

       # numOfBlocks is always strip mode
       numOfblocks = self.image_width//8

       # -cl Indicated the coded word length for each segment. 
       #   Valid values are:
       #     0 - 8-bit word
       #     1 - 16-bit word
       #     2 - 24-bit word
       #     3 - 32-bit word
       parameters = ' -cl 1 -bs ' + str(numOfblocks)

       # segByteLimit
       # segment_size_blocks = image_width/8;
       # segment_size_coeff  = segment_size_blocks*64 = image_width/8*64 = 
       #                          image_width*8		
       # segByteLimit    = uint16(segment_size_coeff*target_bpp_vector/8) = 
       #                       uint16(image_width*8*target_bpp_vector/8) = 
       #                       uint16(image_width*target_bpp_vector)
       segByteLimit = int(self.image_width*self.target_bpp)

       if segByteLimit > 0:
           parameters += ' -bl ' + str(segByteLimit)

       # dc_stop, bitplaneStop, stage_stop in SOPHI ar set as follow
       parameters += ' -dc ' + str(self.dc_stop)
       parameters += ' -bp ' + str(self.bitplaneStop)
       parameters += ' -ss ' + str(self.stage_stop)

       # % (-g) Geometry of raw image data. Parameters are:
       # % 1- zSize (number of image components)
       # % 2- ySize (image height)
       # % 3- xSize (image width)
       # % 4- data type. Possible values are:
       # %   0- boolean (1 byte)
       # %   1- unsigned int (1 byte)
       # %   2- unsigned int (2 bytes)
       # %   3- signed int (2 bytes)
       # %   4- signed int (4 bytes)
       # %   5- signed int (8 bytes)
       # %   6- float (4 bytes)
       # %   7- double (8 bytes)
       # % 5- Byte order (0 if BIG ENDIAN, 1 if LITTLE ENDIAN)    
       # % 6- 1 if 3 first components are RGB, 0 otherwise.
       parameters += ' -g 1 ' + str(self.image_height) + ' ' + str(self.image_width)

       # signedPixels
       if self.signedPixels == 0:
           parameters += ' 2' # uint16
       else:
           parameters += ' 3'  # int16

       # endianess
       if self.endianess == 'b':
       	parameters += ' 0 0' # big endian
       elif self.endianess == 'l':
           parameters += ' 1 0'  # little endian
       else:
           self.pinfo('ednianess error')

       code_call += parameters + ' -o ' + output_imageName_fullPath
       self.pinfo(self.io_format)

       self.pinfo(code_call)
       stream = os.popen(code_call)
       output = stream.read()
       self.pinfo(output)

       return 

   def decode(self) -> None:
       ''' assumes there is already a file with the proper format'''

       # input_imageName_fullPath
       # get last character of input_folder_path
       last = self.input_folder_path[-1]
       if last != '/':
           input_imageName_fullPath = self.input_folder_path + '/'
       else:
           input_imageName_fullPath = self.input_folder_path

       input_imageName_fullPath = input_imageName_fullPath + self.imageName_output_coded+'.rec'

       # output image name in full path format
       last = self.output_folder_path[-1]
       if last != '/':
           output_imageName_fullPath = self.output_folder_path + '/'
       else:
           output_imageName_fullPath = self.output_folder_path

       output_imageName_fullPath = self.output_folder_path + self.imageName_output_decoded+'.raw'

       # TER compression parameters
       parameters = ' -g 1 ' + str(self.image_height) + ' ' + str(self.image_width)

       # signedPixels
       if self.signedPixels == 0:
           parameters += ' 2' # uint16
       else:
           parameters += ' 3'  # int16

       # endianess
       if self.endianess == 'b':
       	parameters += ' 0 0' # big endian
       elif self.endianess == 'l':
           parameters += ' 1 0'  # little endian
       else:
           self.pinfo('ednianess error')

       #parameters += ' -vp 1 -vm 1 '
       code_call = 'java -jar ' + self.TER_decoder_sw_path + ' -i '  + input_imageName_fullPath + parameters + ' -o ' + output_imageName_fullPath

       self.pinfo(code_call)
       stream = os.popen(code_call)
       output = stream.read()
       self.pinfo(output)

       return 

   def convert(self, img, type_limits, direction):
       '''    convert(img, 0, 255, np.uint8) '''
       if direction == 1:
           imin = img.min()
           imax = img.max()
           #print('check 1',imin,imax)
           self.scale_value = (imax,imin)
           #print('check self',self.scale_value)
           a = (type_limits[1] - type_limits[0]) / (imax - imin)
           b = type_limits[1] - a * imax
           new_img = (a * img + b).astype(self.io_format)
       if direction == -1:
           imax,imin = self.scale_value
           #print('check -1',imin,imax)
           #print('check self',self.scale_value)

           a = (type_limits[1] - type_limits[0]) / (imax - imin)
           b = type_limits[1] - a * imax
           new_img = (img.astype(float) - b)/a
       #print('a,b= ',a,b)
       return new_img

   def scale_image(self,data,direction = 1, dynamic = 0):
       # signedPixels
       self.set_format()
       info = np.iinfo(self.io_format)
       type_limits = [info.min,info.max]

       if dynamic != 0:
           img =  self.convert(data,type_limits,direction) #data
       else:
        if direction == 1:
            if self.signedPixels == 0:
                scale = np.max(data)
            if self.signedPixels == 1:
                scale = np.max((np.abs(np.min(data)),np.max(data)))            
            self.scale_value = (scale,0)
            data = data.astype(np.float64) / self.scale_value[0] * type_limits[1] # normalize the data to [ -1 -- 1 ]
            img =  data.astype(self.io_format)       
        if direction == -1:
            img = data.astype(float) / type_limits[1] * self.scale_value[0]

       self.pinfo('Scaling...',self.scale_value)
       return img
   def check_options(self):

       #options = ter.__dict__
       #ter = ter()
       options = vars(self)
       print('...............................................')
       print('Available options:    [default values]')
       print('···············································')
       for item in options.items():
           extension = len(item[0]) 
           spaces = 20 - extension
           string_val = " " * spaces
           print("%s: %s %s" % (item[0],string_val,item[1]))
       print('···············································')
       #signedPixels = 1,     # 1 = signed; 0 = unsigned
       #endianess = 'b'         # 'b' = big endian; 'l' = little endian
       #target_bpp = 5     # 0 = lossless
       
       #print(self.__dict__.keys())
       #print(self.__dict__.values())
       
   def readfits(self,file,info=False,head=0):
   
       if info == True:
           return pyfits.info(file)
       '''helper function to load FITS data set'''
       hdu_list = pyfits.open(file)
       if head != 0:
           return hdu_list[head].data , hdu_list[head].header
       else:
           return hdu_list[head].data.astype(np.dtype('float32')) , hdu_list[head].header

   def PSNR(self,original: float, compressed: float) -> float: 
       info = np.iinfo(self.io_format)
       type_limits = [info.min,info.max]
       mse = np.mean((original - compressed) ** 2) 
       if(mse == 0):  # MSE is zero means no noise is present in the signal . 
                   # Therefore PSNR have no importance. 
           return 100
       max_pixel = type_limits[1]
       psnr = 20 * np.log10(max_pixel / np.sqrt(mse)) 
       return psnr 


''' Drawing defs'''
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def show_one(img: float,vmax=None,vmin=None,xlabel='pixel',ylabel='pixel',title='Image no title',cbarlabel='Some units',save=None):

   plt.figure(figsize=(6, 6))
   ax = plt.gca()
   if vmin == None and vmax == None:
       im = ax.imshow(img, cmap='gray',vmin=img.mean() - PLT_RNG * img.std(),
          vmax=img.mean() + PLT_RNG * img.std(), interpolation='none')
   elif vmin == None:
       im = ax.imshow(img, cmap='gray',vmin=img.mean() - PLT_RNG * img.std(),
          vmax=vmax, interpolation='none')
   elif vmax == None:
       im = ax.imshow(img, cmap='gray',vmin=vmin,
          vmax=img.mean() + PLT_RNG * img.std(), interpolation='none')
   else:
       im = ax.imshow(img, cmap='gray',vmin=vmin,
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

def show_four_row(im1,im2,im3,im4,svmin=[None,None,None,None],svmax=[None,None,None,None],block=True,pause=0.1,title=['','','',''],xlabel='Pixel',ylabel='Pixel',cagonento=0,save=False):

    fig, maps = plt.subplots(1,4,figsize=(24,6))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    if cagonento == 0:
        svmin=[None,None,None,None]
        svmax=[None,None,None,None]
    for i in range(4):
        if i == 0:
            dummy = np.copy(im1)
        if i == 1:
            dummy = np.copy(im2)
        if i == 2:
            dummy = np.copy(im3)
        if i == 3:
            dummy = np.copy(im4)
        if svmin[i] == None and svmax[i] == None:
            svmin[i]=dummy.mean() - PLT_RNG * dummy.std()
            svmax[i]=dummy.mean() + PLT_RNG * dummy.std()

        dummy[dummy<svmin[i]] = svmin[i]
        dummy[dummy>svmax[i]] = svmax[i]
        im = maps[i].imshow(dummy, cmap='gray',vmin = svmin[i],vmax=svmax[i],interpolation='none')
        maps[i].set_title(title[i])
        if i == 0:
            maps[i].set_ylabel(ylabel)
        maps[i].set_xlabel(xlabel)
        colorbar(im)

    if save != False:
        plt.savefig(save)
    plt.show()

    return

# ter = ter()
# ter.check_options()
# ter.image_height = 288
# ter.image_width = 288

# files = ['mhddata_noi05_small.fits','mhddata_noi1_small.fits','mhddata_noi15_small.fits']

# #STOKES I
# ter.signedPixels=0
# bpp = np.arange(3,10,0.1)
# rms = np.zeros((3,len(bpp)))
# rms_scaled = np. zeros((3,len(bpp)))
# i = 0
# for s in files:
#     j = 0
#     data, header = ter.readfits(s)
#     datau = data[:,:,0,0]
#     #datas = data[:,:,1,0]
#     datau_scaled = ter.scale_image(datau,1,1)
#     ter.save_raw(datau_scaled)
#     norma_o = np.mean(datau)
#     norma_d = np.mean(datau_scaled)
#     for bpp_i in bpp:
#         ter.target_bpp = bpp_i.round(decimals=2)
#         print(ter.target_bpp)
#         ter.code()
#         ter.decode()
#         out_data = ter.read_output()
#         out_data_descaled = ter.scale_image(out_data,-1,1)
#         rms1 = np.std(datau - out_data_descaled)
#         rms2 = np.std(datau_scaled - out_data)
#         rms[i,j] = rms1
#         rms_scaled[i,j] = rms2
#         j += 1
#     i += 1
# plt.plot(bpp,rms[0,:]/norma_o)
# plt.plot(bpp,rms[1,:]/norma_o)
# plt.plot(bpp,rms[2,:]/norma_o)
# plt.plot(bpp,rms_scaled[0,:]/norma_d)
# plt.plot(bpp,rms_scaled[1,:]/norma_d)
# plt.plot(bpp,rms_scaled[2,:]/norma_d)
# plt.show()

# quit()

# ###CHECK CODE
# # %%
# ter = ter()
# ter.check_options()
# ter.image_height = 288
# ter.image_width = 288

# data, header = ter.readfits('mhddata_noi05_small.fits')

# datau = data[:,:,0,0]
# datas = data[:,:,1,0]
# sy,sx = datau.shape
# ter.image_height = sy
# ter.image_width = sx


# # %%
# '''un-signed image'''
# ter.check_options()
# ter.signedPixels=0
# datau_scaled = ter.scale_image(datau,1,1)
# datau_back = ter.scale_image(datau_scaled,-1,1)
# #show_one(datau_scaled,title='datau_scaled',vmax=2**16,vmin=-0)
# #show_one(datau_back-datau,title='datau_scaled - original')#,vmax=100,vmin=-100)

# ter.save_raw(datau_scaled)
# data_read = ter.read_input()
# ter.check_options()
# #show_one(data_read,title='datau_scaled',vmax=2**16,vmin=-0)
# #show_one(data_read-datau_scaled,title='datau_scaled - original')

# #show_four_row(datau_scaled,datau_back-datau,data_read,data_read-datau_scaled,svmin=[0,None,0,None],svmax=[2**16,None,2**16,None],block=True,pause=0.1,\
# #    title=['datau_scaled','datau_scaled - original','datau_scaled_read_disk','datau_scaled_read_disk - datau_scaled'],xlabel='Pixel',ylabel='Pixel')

# ter.target_bpp = 3.2
# ter.check_options()
# ter.code()
# ter.check_options()
# ter.decode()

# in_data = ter.read_input()
# out_data = ter.read_output()
# out_data_descaled = ter.scale_image(out_data,-1,1)
# # show_four_row(in_data,out_data,datau,out_data_descaled,svmin=[0,None,0,None],svmax=[2**16,None,2**16,None],block=True,pause=0.1,\
# #     title=['datau_scaled_read_disk','datau_read_decompress','original','datau_read_decompress_descaled'],xlabel='Pixel',ylabel='Pixel')
# # show_four_row(in_data-out_data,datau,datau-out_data_descaled,out_data_descaled,svmin=[0,None,0,None],svmax=[2**16,None,2**16,None],block=True,pause=0.1,\
# #     title=['datau_scaled_read_disk-datau_read_decompress','datau_scaled_read_disk','original-datau_read_decompress_descaled','datau_read_decompress_descaled'],xlabel='Pixel',ylabel='Pixel')

# bpp = np.arange(3,10,0.1)
# rms = np. array([])
# norma = np.mean(datau)
# for i in bpp:
#     ter.target_bpp = i.round(decimals=2)
#     print(ter.target_bpp)
#     ter.code()
#     ter.decode()
#     out_data = ter.read_output()
#     out_data_descaled = ter.scale_image(out_data,-1,1)
#     new_rms = rms,np.std(datau - out_data_descaled)
#     rms = np.append(rms,new_rms)
# plt.plot(bpp,rms/norma)
# plt.show()
# quit()


# # %%

# quit()
# # %%
# '''signed image'''
# ter.signedPixels=1

# datas_scaled = ter.scale_image(datas,1)
# datas_back = ter.scale_image(datas_scaled,-1)
# show_one(datas_back-data[:,:,1,0],vmax=100,vmin=-100)

# ter.save_raw(datas_scaled)
# data_read = ter.read_input()
# print(np.max(data_read),np.min(data_read))
# show_one(data_read)
# show_one(data_read-datas_scaled)
# show_one(data[:,:,1,0])

# ter.check_options()


# # %%
# ter.code()


# # %%
# ter.decode()


# # %%
# in_data = ter.read_input()
# out_data = ter.read_output()
# print(in_data.shape,out_data.shape)
# show_one(in_data)
# show_one(out_data[0:50,0:50])
# show_one(datas[0:50,0:50])
# show_one(ter.scale_image(out_data[0:50,0:50],-1))


# # %%
# print(in_data.shape,out_data.shape)


# # %%



# # %%



# # %%



