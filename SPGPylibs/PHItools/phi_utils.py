#hough

from datetime import timedelta
import matplotlib.dates as mdates
import numpy as np
from astropy import units as u 
import spiceypy
import warnings
warnings.filterwarnings('ignore')

def dot(K,L):
    return sum(i*j for i,j in zip(K, L)) if len(K) == len(L) else 0

def angle(K,L):
    dot_prod = dot(K,L)
    n_K = np.linalg.norm(K, axis=0)
    n_L = np.linalg.norm(L, axis=0)
    c = dot_prod/n_K/n_L
    angle = np.arccos(np.clip(c, -1, 1))
    return angle*180/np.pi

def sphere2cart(r, phi, theta):
    x = r*np.cos(theta)*np.cos(phi)
    y = r*np.cos(theta)*np.sin(phi)
    z = r*np.sin(theta)
    return (x, y, z) 

def cart2sphere(x,y,z):
    r = np.sqrt(x**2+ y**2 + z**2)            # r
    theta = np.arctan2(z,np.sqrt(x**2+ y**2))     # theta
    phi = np.arctan2(y,x)                        # phi
    return (r, theta, phi)

def phi_orbit(init_date, object, end_date = None, frame = 'ECLIPJ2000', resolution = 1, kernel = None):
    ''' get basic orbitar parameters using Spice kernels
    Need to have spicepy installed.
    '''

    spice_files = ['Orbits-data/kernels/lsk/naif0012.tls',\
        'Orbits-data/kernels/spk/de421.bsp',\
        'Orbits-data/kernels/pck/pck00010.tpc',\
        'Orbits-data/kernels/spk/solo_ANC_soc-orbit_20200210-20301120_L004_V1_00062_V01.bsp']

    #Load the kernals 
    spiceypy.furnsh(spice_files)

    # Calculate the timerange
    # times -->> TDB seconds past the J2000 epoch.
   
    if end_date == None:
        sc_time = init_date
    else:
        sc_time = []
        while init_date < end_date:
            sc_time.append(init_date)
            init_date += timedelta(days=resolution)

    sc_time_spice = [spiceypy.str2et(t.strftime('%Y-%m-%d %H:%M')) for t in sc_time]

    solo, lT = spiceypy.spkezr('Solar Orbiter', sc_time_spice,frame , 'NONE', 'Sun')
    pos_solo = np.array(solo)[:, :3] * u.km
    vel_solo = np.array(solo)[:, 3:] * u.km / u.s

    sc_x = pos_solo[:, 0].to(u.au)
    sc_y = pos_solo[:, 1].to(u.au)
    sc_z = pos_solo[:, 2].to(u.au)
    sc_vx = vel_solo[:, 0]
    sc_vy = vel_solo[:, 1]
    sc_vz = vel_solo[:, 2]

    sc_sp = np.sqrt(sc_vx**2 + sc_vy**2 + sc_vz**2)

    sc_r = np.sqrt(sc_x**2 + sc_y**2 + sc_z**2)
    elevation = np.rad2deg(np.arcsin(sc_z / sc_r))
    angle = np.rad2deg(np.arcsin((sc_x.to_value()**2 + sc_y.to_value()**2 + sc_z.to_value()**2 ) / sc_r.to_value()))

    sc_r, sc_lat, sc_lon = cart2sphere(sc_x,sc_y,sc_z)

    rad_sun = 1391016 # solar radius in km
    dit_sun = 149597870 #Sun-Earth distance 1 AU
    asiz = np.arctan(rad_sun/dit_sun) #Sun angular size
    sun_arc = asiz*180/np.pi*60*60  #Sun angular size in arcsec

    s_size = sun_arc/sc_r #,'solar size in arcsec from solo'

    sc_vr = (sc_x*sc_vx + sc_y*sc_vy + sc_z*sc_vz)/np.sqrt(sc_x**2 + sc_y**2 + sc_z**2)
    sc_vt = (sc_x*sc_vx + sc_y*sc_vy + sc_z*sc_vz)/np.sqrt(sc_x**2 + sc_y**2 + sc_z**2)
    
    screc=np.rec.array([mdates.date2num(sc_time),sc_r,sc_lon,sc_lat,sc_x,sc_y,sc_z,sc_vx,sc_vy,sc_vz,\
        sc_sp, sc_vr, elevation, angle, s_size],\
        dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),\
        ('y','f8'),('z','f8'),('vx','f8'),('vy','f8'),('vz','f8'),('sp','f8'),\
        ('vr','f8'),('elevation','f8'),('angle','f8'),('s_size','f8')])
    
    return screc

def phi_orbit_test():

    from datetime import datetime
    import matplotlib.pyplot as plt

    starttime = datetime(2020, 11, 30)
    endtime = datetime(2022, 9, 30)
    res_in_days = 0.1
#kernel = None

    solo_info = phi_orbit(starttime, 'Solar Orbiter', end_date = endtime,resolution = res_in_days)

    plt.plot_date(solo_info.time,solo_info.vr ,'-')
    plt.xlabel('Time')
    plt.ylabel('Velocity [km/s]')
    plt.show()

    when = [datetime(2020, 5, 15,19,0,0),datetime(2020, 5, 21,14,00,0)]

    solo_heeq = phi_orbit(when, 'Solar Orbiter')

    for i in range(len(when)):
        print('___________________________________________')
        print('Time: ',when[i])
        print('   Solo Sun center [HEEQ,degree]: ',solo_heeq.lat[i]*180/np.pi)
        print('   Longitud (wrt Earth): ',solo_heeq.lon[i]*180/np.pi)
        print('   Distance [AU]: ',solo_heeq.r[i])
        print('   Solar size in arcsec from solo: ',solo_heeq[i].s_size)
        print('   FDT solar diameter in pixels from solo: ',solo_heeq[i].s_size/3.61)
        print('-------------------------------------------')
        

    return None