import numpy as np
from scipy.integrate import simps

def COG(stokes,wave,wave_axis,lande_factor=0,cpos = False):

    '''
    cog: data input data = [wave,stokes,x,y], wave_axis = 'wavelength axis in mA centered in line', rx,ry bondaries

    wavelength = 6173.3356
    wave_axis = np.array([-160,-80,0,80,160,300])/1000.+wavelength
    '''
    
    # check dimensions of input 
    # posibilitis are: 

    # if blos then lande_factor != 0


    #multiple profiles#
    l,s,sy,sx = stokes.shape
    try:
    lpos = int(cpos)
    except:
    lpos = l

    Ic = stokes[lpos-1,0]
    t1 = wave_axis[:,np.newaxis,np.newaxis] * ( Ic[np.newaxis,:,:] - stokes[:,0,:,:] ) 
    tc = ( Ic[np.newaxis,:,:] - stokes[:,0,:,:] )
    Itn = simps(t1, x=wave_axis,axis=0) / simps(tc, x=wave_axis,axis=0)

    vlos = -(wave - Itn ) * 2.99792458e+5 / wave 

    t_plus = stokes[:,0,:,:] + stokes[:,3,:,:]
    t_minus = stokes[:,0,:,:] - stokes[:,3,:,:]
    t1 = wave_axis[:,np.newaxis,np.newaxis] * ( Ic[np.newaxis,:,:] - t_plus ) #Ic*0.5 ??
    tc1 = Ic[np.newaxis,:,:] - t_plus 
    t2 = wave_axis[:,np.newaxis,np.newaxis] * ( Ic[np.newaxis,:,:] - t_minus ) #Ic*0.5 ??
    tc2 = Ic[np.newaxis,:,:] - t_minus
    l_plus = simps(t1, x=wave_axis,axis=0) / simps(tc1, x=wave_axis,axis=0)
    l_minus = simps(t2, x=wave_axis,axis=0) / simps(tc2, x=wave_axis,axis=0)

    blos = (l_plus - l_minus) / 2 / 4.67e-13/ (lande_factor * wave**2)

    return vlos,blos
