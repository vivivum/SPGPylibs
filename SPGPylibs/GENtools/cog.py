import numpy as np
from scipy.integrate import simps

def cog(input_data,wave,wave_axis,lande_factor=0,cpos = False):

    '''
    Calculates the velocity [and longitudinal field] of a profile or set of profiles using CoG technique.
    input_data: the input data may have dimmensions of:
        [wave] = 1D the program assumes Stokes I profile (one profile) - returns Vlos
        [wave, Stokes] = 2D the program assumes Stokes I,Q,U,V [length 4] (one profile) - return Vlos + Blos
        [wave, X , Y] = 3D the program assumes Stokes I only in an X and Y image - return Vlos 
        [wave, Stokes, X , Y] = 4D the program assumes Stokes I,Q,U,V in an X and Y image - return Vlos + Blos

    wave: line central wavelength in Angstrom.
    wave_axis = wavenelgth axis in Angstrom.
    cpos: position where the continuum of the line is located. Default is last position.
    In case of Stokes parameters are given, the user has to provide the lande_factor of the transition
    '''
    
    # check dimensions of input 
    # posibilitis are: 
    ndim = input_data.ndim
    if ndim == 1:
        l = input_data.shape
        print('Single Stokes I profile')
        #check continuum position
        try:
            lpos = int(cpos)
        except:
            lpos = l

        Ic = input_data[lpos-1]
        t1 = wave_axis * ( Ic - input_data ) 
        tc = ( Ic - input_data )
        Itn = simps(t1, x=wave_axis) / simps(tc, x=wave_axis)
        vlos = -(wave - Itn ) * 2.99792458e+5 / wave 

        return vlos

    elif ndim == 2:
        l,s = input_data.shape
        print('Single Full Stokes profile')
        if lande_factor==0:
            print('Warning: lande factor is zero pelotero')
        #check continuum position
        try:
            lpos = int(cpos)
        except:
            lpos = l

        Ic = input_data[lpos-1,0]
        t1 = wave_axis * ( Ic - input_data[:,0] ) 
        tc = ( Ic - input_data[:,0] )
        Itn = simps(t1, x=wave_axis) / simps(tc, x=wave_axis)
        vlos = -(wave - Itn ) * 2.99792458e+5 / wave 

        t_plus = input_data[:,0] + input_data[:,3]
        t_minus = input_data[:,0] - input_data[:,3]
        t1 = wave_axis * ( Ic - t_plus ) #Ic*0.5 ??
        tc1 = Ic - t_plus 
        t2 = wave_axis * ( Ic - t_minus ) #Ic*0.5 ??
        tc2 = Ic - t_minus
        l_plus = simps(t1, x=wave_axis) / simps(tc1, x=wave_axis)
        l_minus = simps(t2, x=wave_axis) / simps(tc2, x=wave_axis)
        blos = (l_plus - l_minus) / 2 / 4.67e-13/ (lande_factor * wave**2)

        return vlos,blos

    elif ndim == 3:
        l,sy,sx = input_data.shape
        print('Multiple Stokes I profiles [l,x,y]: ', l, sx, sy )
        #check continuum position
        try:
            lpos = int(cpos)
        except:
            lpos = l

        Ic = input_data[lpos-1,:,:]
        t1 = wave_axis[:,np.newaxis,np.newaxis] * ( Ic[np.newaxis,:,:] - input_data ) 
        tc = ( Ic[np.newaxis,:,:] - input_data )
        Itn = simps(t1, x=wave_axis,axis=0) / simps(tc, x=wave_axis,axis=0)
        vlos = -(wave - Itn ) * 2.99792458e+5 / wave 

        return vlos

    elif ndim == 4:
        l,s,sy,sx = input_data.shape
        print('Multiple Full Stokes profiles [l,stokes,x,y]: ', l, s, sx, sy )
        if lande_factor==0:
            print('Warning: lande factor is zero pelotero')
        #check continuum position
        try:
            lpos = int(cpos)
        except:
            lpos = l

        Ic = input_data[lpos-1,0,:,:]
        t1 = wave_axis[:,np.newaxis,np.newaxis] * ( Ic[np.newaxis,:,:] - input_data[:,0,:,:] ) 
        tc = ( Ic[np.newaxis,:,:] - input_data[:,0,:,:] )
        Itn = simps(t1, x=wave_axis,axis=0) / simps(tc, x=wave_axis,axis=0)
        vlos = -(wave - Itn ) * 2.99792458e+5 / wave 

        t_plus = input_data[:,0,:,:] + input_data[:,3,:,:]
        t_minus = input_data[:,0,:,:] - input_data[:,3,:,:]
        t1 = wave_axis[:,np.newaxis,np.newaxis] * ( Ic[np.newaxis,:,:] - t_plus ) #Ic*0.5 ??
        tc1 = Ic[np.newaxis,:,:] - t_plus 
        t2 = wave_axis[:,np.newaxis,np.newaxis] * ( Ic[np.newaxis,:,:] - t_minus ) #Ic*0.5 ??
        tc2 = Ic[np.newaxis,:,:] - t_minus
        l_plus = simps(t1, x=wave_axis,axis=0) / simps(tc1, x=wave_axis,axis=0)
        l_minus = simps(t2, x=wave_axis,axis=0) / simps(tc2, x=wave_axis,axis=0)
        blos = (l_plus - l_minus) / 2 / 4.67e-13/ (lande_factor * wave**2)

        return vlos,blos

    else:
        print('No input data or wrong dimentions',ndim)  
        return     


