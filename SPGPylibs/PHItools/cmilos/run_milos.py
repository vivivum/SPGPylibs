import pymilos
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import subprocess
import time

DTYPE_INT = np.intc
DTYPE_DOUBLE = np.float_

from scipy.io import readsav
s = readsav('cmilos_testdata/test_data.save')
datos = s.data
wave_axis = s.wave_axis

datos = np.tile(datos, 4) 

print('shape1',datos.shape)                                                                                                                                                                                                                      
datos = np.einsum('ijkl->klij',datos)
y,x,p,l = datos.shape
print('shape2',datos.shape)                                                                                                                                                                                                                      
print('shape2',y,x,p,l)                                                                                                                                                                                                                      

options = np.zeros((4))#,dtype=DTYPE_INT)
options[0] = len(wave_axis) #NLAMBDA
options[1] = 15 #MAX_ITER
options[2] = 0 #CLASSICAL_ESTIMATES [0,1]
options[3] = 0 #RFS [0,1,2]
# options[4] = 0 #FWHM(in A)
# options[5] = 0 #DELTA(in A)
# options[6] = 0 #NPOINTS

npx = x
npy = y

print('----------- pmilos---------')
start = time.process_time()
out = pymilos.pmilos(options,datos[0:npy,0:npx,:,:],wave_axis)
print(time.process_time() - start)
print(out.shape)

out = np.reshape(out,(npy,npx,12))                                                                                                                                                                                        

#use milos 

sdata = datos[0:npy,0:npx,:,:]
y,x,p,l = sdata.shape

print('----------- milos---------')
start = time.process_time()

filename = 'dummy_in.txt'
with open(filename,"w") as f:
    for i in range(x):
        for j in range(y):
            for k in range(l):
                f.write('%e %e %e %e %e \n' % (wave_axis[k],sdata[j,i,0,k],sdata[j,i,1,k],sdata[j,i,2,k],sdata[j,i,3,k]))
del sdata

rte_on = subprocess.call("time ./milos.MacBook-Pro-de-David 6 15 0 0 dummy_in.txt  >  dummy_out.txt",shell=True)
print(rte_on)
print(time.process_time() - start)


# cmd = subprocess.Popen('ls -l', shell=True, stdout=PIPE)
# for line in cmd.stdout.readlines():
#     print line


res = np.loadtxt('dummy_out.txt')
npixels = res.shape[0]/12.
result = np.zeros((y*x,12)).astype(float)
for i in range(y*x):
    result[i,:] = res[i*12:(i+1)*12]
result = np.einsum('ijk->jik',result.reshape(y,x,12))

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(121)
im1 = ax1.imshow(result[:,:,8], interpolation='None')

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im1, cax=cax, orientation='vertical')

ax2 = fig.add_subplot(122)
im2 = ax2.imshow(out[:,:,8], interpolation='None')

divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im2, cax=cax, orientation='vertical')

plt.show()

			# outputdata[cnt_model] = contador;
			# outputdata[cnt_model+1] = iter;
			# outputdata[cnt_model+2] = initModel.B;
			# outputdata[cnt_model+3] = initModel.gm;
			# outputdata[cnt_model+4] = initModel.az;
			# outputdata[cnt_model+5] = initModel.eta0;
			# outputdata[cnt_model+6] = initModel.dopp;
			# outputdata[cnt_model+7] = initModel.aa;
			# outputdata[cnt_model+8] = initModel.vlos;
			# outputdata[cnt_model+9] = initModel.S0;
			# outputdata[cnt_model+10] = initModel.S1;
			# outputdata[cnt_model+11] = chisqrf;

import sys
sys.exit()
#prepare input
input_data = datos[:,:,0:1,0:1].astype(DTYPE_DOUBLE)  #get data and set dtype
nlon,npol,ny,nx = input_data.shape
input_data = input_data.flatten(order='C')#.copy(order='c')

wave_axis = wave_axis.astype(DTYPE_DOUBLE)            #get wave_axis dtype

length = len(input_data)
output_data = np.zeros((length//len(wave_axis)//4 * 12),dtype=DTYPE_DOUBLE)


print(options,input_data,wave_axis)
print(options.flags,input_data.flags,wave_axis.flags,output_data.flags)

pymilos.py_milos(options,input_data,wave_axis,output_data)


    # options = options.astype(DTYPE_INT)
    # options = options.copy(order='c').astype(DTYPE_INT)
    # wave_axis = wave_axis.copy(order='C') 
    # input_data = input_data.flatten().copy(order='C') #input_data.flatten()#(order='C')
    # length = len(input_data)
    # output_data = np.zeros((length//options[0]//4,12), dtype=DTYPE_DOUBLE,order='C')

    # print(options.shape,input_data.shape,output_data.shape)
    # print(options.flags,input_data.flags,output_data.flags)

    # x = np.array(list(range(10)), '>i4') # big endian
    # newx = x.byteswap().newbyteorder() # force native byteorder

    # py_milos_wrap(input_data, wave_axis, output_data)

    # return output_data

#out = pymilos.py_milos(options,wave_axis,input_data)

