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

def phi_orbit(init_date, object, end_date = None, frame = 'ECLIPJ2000', resolution = 1, kernel = None, kernel_dir = None):
    ''' get basic orbitar parameters using Spice kernels
    Need to have spicepy installed.
    The program goes to the kernel dir for loading the mk (meta-kernel) file that takes into account all stuff.

    Default frame is Mean ecliptic and equinox of J2000 (ECLIPJ2000)
    ECLIPJ2000  -- >> Ecliptic coordinates based upon the J2000 frame.

    Kernels provided by SOLO (see readme in the kernels)
      SPICE Frame Name            Common names and other designators
      --------------------        --------------------------------------
      SUN_ARIES_ECL               HAE, Solar Ecliptic (SE)
      SUN_EARTH_CEQU              HEEQ, Stonyhurst Heliographic
      SUN_INERTIAL                HCI, Heliographic Inertial (HGI)
      EARTH_SUN_ECL               GSE of Date, Hapgood
      SUN_EARTH_ECL               HEE of Date
      EARTH_MECL_MEQX             Mean Ecliptic of Date (ECLIPDATE)
    
     The Heliocentric Earth Ecliptic frame (HEE) is defined as follows (from [3]):
       -  X-Y plane is defined by the Earth Mean Ecliptic plane of date,
           therefore, the +Z axis is the primary vector,and it defined as
            the normal vector to the Ecliptic plane that points toward the
            north pole of date;
       -  +X axis is the component of the Sun-Earth vector that is
            orthogonal to the +Z axis;
       -  +Y axis completes the right-handed system;
       -  the origin of this frame is the Sun's center of mass.

     The Heliocentric Inertial Frame (HCI) is defined as follows (from [3]):

       -  X-Y plane is defined by the Sun's equator of epoch J2000: the +Z
          axis, primary vector, is parallel to the Sun's rotation axis of
          epoch J2000, pointing toward the Sun's north pole;
       -  +X axis is defined by the ascending node of the Sun's equatorial
          plane on the ecliptic plane of J2000;
       -  +Y completes the right-handed frame;
       -  the origin of this frame is the Sun's center of mass.

    The Heliocentric Earth Equatorial (HEEQ) frame is defined as follows (from [3]
      and [4]):

       -  X-Y plane is the solar equator of date, therefore, the +Z axis 
            is the primary vector and it is aligned to the Sun's north pole
            of date;
       -  +X axis is defined by the intersection between the Sun equatorial
            plane and the solar central meridian of date as seen from the Earth.
            The solar central meridian of date is defined as the meridian of the
            Sun that is turned toward the Earth. Therefore, +X axis is the
            component of the Sun-Earth vector that is orthogonal to the +Z axis;
       -  +Y axis completes the right-handed system;
       -  the origin of this frame is the Sun's center of mass.

    SOLO mission specific generic frames are (see solo_ANC_soc-sci-fk_V07.tf):

      SOLO_SUN_RTN                Sun Solar Orbiter Radial-Tangential-Normal
      SOLO_SOLAR_MHP              S/C-centred mirror helioprojective
      SOLO_IAU_SUN_2009           Sun Body-Fixed based on IAU 2009 report
      SOLO_IAU_SUN_2003           Sun Body-Fixed based on IAU 2003 report
      SOLO_GAE                    Geocentric Aries Ecliptic at J2000 (GAE)
      SOLO_GSE                    Geocentric Solar Ecliptic at J2000 (GSE)
      SOLO_HEE                    Heliocentric Earth Ecliptic at J2000 (HEE)
      SOLO_VSO                    Venus-centric Solar Orbital (VSO)

   Heliospheric Coordinate Frames developed for the NASA STEREO mission:

      SOLO_ECLIPDATE              Mean Ecliptic of Date Frame
      SOLO_HCI                    Heliocentric Inertial Frame
      SOLO_HEE_NASA               Heliocentric Earth Ecliptic Frame
      SOLO_HEEQ                   Heliocentric Earth Equatorial Frame ***
      SOLO_GEORTN                 Geocentric Radial Tangential Normal Frame

   Heliocentric Generic Frames(*):

      SUN_ARIES_ECL               Heliocentric Aries Ecliptic   (HAE)
      SUN_EARTH_CEQU              Heliocentric Earth Equatorial (HEEQ)
      SUN_EARTH_ECL               Heliocentric Earth Ecliptic   (HEE)
      SUN_INERTIAL                Heliocentric Inertial         (HCI)

    '''

    import os 

    KERNELS_FULL_PATH_DIRECTORY = \
        '/Users/orozco/Dropbox_folder/Python/VS-GitHub/SPGPylibs/SPGPylibs/PHItools/Orbits-data/kernels/'
    try:
        if kernel_dir != None:
            KERNELS_FULL_PATH_DIRECTORY = kernel_dir
    except:
        pass

    REQUIRED_KERNELS = ['mk/solo_ANC_soc-flown-mk_v107_20210412_001.tm']
    try:
        if kernel != None:
            REQUIRED_KERNELS = kernel
    except:
        pass    

    #Here the function moves to KERNELS_FULL_PATH_DIRECTORY folder to furnsh the corresponding kernels
    
    #Get urrent folder
    cd = os.getcwd()
    try:
        os.chdir(KERNELS_FULL_PATH_DIRECTORY)
    except Exception:
        print("Unable to go to kernel folder location: {}",KERNELS_FULL_PATH_DIRECTORY)        
        return None
    try:
        #Load the kernals 
        spiceypy.furnsh(REQUIRED_KERNELS)
    except:
        raise IOError("Error reading kernel files",KERNELS_FULL_PATH_DIRECTORY+REQUIRED_KERNELS)

    #back to working directory
    os.chdir(cd)

    # Calculate the timerange
    # times -->> TDB seconds past the J2000 epoch.
   
    if end_date == None:
        sc_time = init_date
    else:
        sc_time = []
        while init_date < end_date:
            sc_time.append(init_date)
            init_date += timedelta(days=resolution)

    # Convert UTC to ephemeris time
    sc_time_spice = [spiceypy.str2et(t.strftime('%Y-%m-%d %H:%M')) for t in sc_time]


    solo, lT = spiceypy.spkezr(object, sc_time_spice, frame , 'NONE', 'Sun')
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

    solo_info = spg.phi_orbit(starttime, 'Solar Orbiter', end_date = endtime,resolution = res_in_days, frame = 'ECLIPJ2000')

    plt.plot_date(solo_info.time,solo_info.vr ,'-')
    plt.xlabel('Time')
    plt.ylabel('Velocity [km/s]')
    plt.show()

    solo_info = spg.phi_orbit(starttime, 'Solar Orbiter', end_date = endtime,resolution = res_in_days, frame = 'SOLO_HEEQ')

    plt.plot_date(solo_info.time,solo_info.vr ,'-')
    plt.xlabel('Time')
    plt.ylabel('Velocity [km/s]')
    plt.show()

    when = [datetime(2020, 5, 15,19,0,0),datetime(2020, 5, 21,14,00,0)]

    solo_j2000 = spg.phi_orbit(when, 'Solar Orbiter', frame = 'ECLIPJ2000')
    solo_heeq = spg.phi_orbit(when, 'Solar Orbiter', frame = 'SOLO_HEEQ')

    for i in range(len(when)):
        print('____________________J2000 - HEEQ_______________________')
        print('Time: ',when[i])
        print('   Solo Sun center [degree]: ',solo_j2000.lat[i]*180/np.pi,solo_heeq.lat[i]*180/np.pi)
        print('   Longitud (wrt Earth): ',solo_j2000.lon[i]*180/np.pi,solo_heeq.lon[i]*180/np.pi)
        print('   Distance [AU]: ',solo_j2000.r[i],solo_heeq.r[i])
        print('   Solar size in arcsec from solo: ',solo_j2000[i].s_size,solo_heeq[i].s_size)
        print('   FDT solar diameter in pixels from solo: ',solo_j2000[i].s_size/3.61,solo_heeq[i].s_size/3.61)
        print('-------------------------------------------')

        
def phi_orbit_conjuntion():

    from datetime import datetime
    import matplotlib.pyplot as plt

    when = [datetime(2021, 2, 1,12,30,3),\
        datetime(2021, 2, 2 ,12,30,2),\
        datetime(2021, 2, 3 ,12,30,2),\
        datetime(2021, 2, 4 ,12,30,2),\
        datetime(2021, 2, 5 ,0,30,2),\
        datetime(2021, 2, 5 ,6,30,2),\
        datetime(2021, 2, 5 ,14,30,2),\
        datetime(2021, 2, 5 ,22,30,2),\
        datetime(2021, 2, 6 ,12,30,2),\
        datetime(2021, 2, 7 ,12,30,2),\
        datetime(2021, 2, 9 ,12,30,3),\
        datetime(2021, 2, 10,12,30,2),\
        datetime(2021, 2, 11,12,30,2),\
        datetime(2021, 2, 12,12,30,2),\
        datetime(2021, 2, 21,6,00,2)]

    # 2021-02-01 12:30:03
    # 2021-02-02 12:30:02
    # 2021-02-03 12:30:02
    # 2021-02-04 12:30:02
    # 2021-02-05 00:30:02
    # 2021-02-05 06:30:02
    # 2021-02-05 14:30:02
    # 2021-02-05 22:30:02
    # 2021-02-06 12:00:02
    # 2021-02-07 12:30:02
    # 2021-02-09 12:30:03
    # 2021-02-10 12:30:02
    # 2021-02-11 12:30:02
    # 2021-02-12 12:30:02
    # 2021-02-21 06:00:02FF

    solo = phi_orbit(when, 'Solar Orbiter')

    for i in range(len(when)):
        print('___________________________________________')
        print('Time: ',when[i])
        print('   Solo Sun center [HEEQ,degree]: ',solo.lat[i]*180/np.pi)
        print('   Solo radial velocity [km/s]: ',solo.vr[i])
        print('   Longitud (wrt Earth): ',solo.lon[i]*180/np.pi)
        print('   Distance [AU]: ',solo.r[i])
        print('   Solar size in arcsec from solo: ',solo[i].s_size)
        print('   FDT solar diameter in pixels from solo: ',solo[i].s_size/3.61)
        print('-------------------------------------------')
    
    for i in range(len(when)):
        print('   Time, Solo radial velocity [km/s]: ',when[i],solo.vr[i])

    # t = spiceypy.str2et('2021-2-1 12:30:03') 
    # solo, lT = spiceypy.spkezr('Solar Orbiter', t,'SOLO_HEEQ', 'NONE', 'Sun') 


    return None