from tqdm import tqdm
import numpy as np
from matplotlib import pyplot as plt
import time

pardata = {'delta1': np.array([225, 225, 315, 315]),
            'delta2' : np.array([234, 125.26, 54.74, 305.26]),
            'theta1' : 0,
            'theta2' : 45,
            'pol_angle' : 0,
            'rot_inst' : 0}

def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            t = (te - ts) * 1000
            if t > 1e3:
                print('%r  %2.2f seconds' %
                      (method.__name__, t/1000.))
            elif t > 1e6:
                print('%r  %2.2f mim' %
                      (method.__name__, t/1000./60.))
            else:
                print('%r  %2.2f ms' %
                      (method.__name__, t))

        return result
    return timed 

def svd_solve(A, b=None, w_cut=1e-10):
    """
    This function solves the system of equations Ax=b by calculating the
    inverse of A using the SVD method: x=A^(-1)*b; A^(-1)=V*S^(-1)*U'
    Inputs:
        A: 2D array of dimensions nxm (n>=m)
        b: 1D array of dimensions n
        w_cut: cut-off frequency for singular values (fraction of the maximum).
        Diagonal elements S^(-1) are zero for the positions of S where its
        value is less than w_cut.
    """
    # Ab = np.abs(A)
    # A = np.where(Ab < np.min(Ab)*10.,0,A)
    # #Ac = np.zeros_like(b)
    # #Ac[:,0] = A[:,0]
    # #b = np.where(Ac == 0.0, 0 ,b)
    # bb = np.abs(b)
    # b = np.where(bb < np.min(bb)*10.,0,b)

    U, S, Vt = np.linalg.svd(A)
    sigma = w_cut*np.max(S)
    Sinv = np.where(S < sigma, 0, (1/S))
    # print(S)
    # print(sigma)
    zeros_Sinv = np.argwhere(Sinv == 0)
    Ainv = np.dot(np.transpose(Vt)*Sinv, np.transpose(U))
    if b is None:
        return Ainv
    else:
        delta_a = np.dot(Ainv, b)
        return delta_a

def cal_unit(alpha, delta, theta, angle_rot=None):
    '''
    #TODO implement derivatives (it is in lib_pmp_v1.4)
    #theta = off set of the retarder
    #delta = retardance
    #alpha = off set of the polarizer(angle)
    '''
    theta_r = theta * np.pi/180.
    delta_r = delta * np.pi/180.
    alpha_r = alpha * np.pi/180.
    c2 = np.cos(2*theta_r)
    s2 = np.sin(2*theta_r)

    cl = 0.5*np.matrix( [1.,
                        (c2**2.+s2**2*np.cos(delta_r))*np.cos(2*alpha_r) +
                        c2*s2*(1-np.cos(delta_r))*np.sin(2*alpha_r),
                        (s2**2.+c2**2*np.cos(delta_r))*np.sin(2*alpha_r) + 
                        c2*s2*(1-np.cos(delta_r))*np.cos(2*alpha_r),
                        s2*np.sin(delta_r)*np.cos(2*alpha_r)-c2*np.sin(delta_r)*np.sin(2*alpha_r)] ).T
                        #transpose just to have a column matrix
    if angle_rot is None:
        return cl
    else:
        return np.matmul(rotation_matrix(angle_rot),cl)

# theta_r = THETA * !dpi/180d0
# delta_r = DELTA * !dpi/180d0
# alpha_r = ALPHA * !dpi/180d0
# c2 = cos(2d0*theta_r)
# s2 = sin(2d0*theta_r)

# CL = 0.5*[1d0,$
#       (c2^2d0+s2^2d0*cos(delta_r))*cos(2d0*alpha_r) + c2*s2*(1d0-cos(delta_r))*sin(2d0*alpha_r),$
#       (s2^2d0+c2^2d0*cos(delta_r))*sin(2d0*alpha_r) + c2*s2*(1d0-cos(delta_r))*cos(2d0*alpha_r),$
#       s2*sin(delta_r)*cos(2d0*alpha_r)-c2*sin(delta_r)*sin(2d0*alpha_r)]

def test_cal_unit():
    theta = np.arange(0, 360, 1)
    out = np.zeros((4,360))
    for i in theta:
        out[:,i] = cal_unit(0.,75.,i,angle_rot=0).flatten()
    plt.plot(theta, out[0, :], label = 'I1')
    plt.plot(theta, out[1, :], label = 'I2')
    plt.plot(theta, out[2, :], label = 'I3')
    plt.plot(theta, out[3, :], label = 'I4')
    plt.title('rot = 0 deg')
    plt.legend()
    plt.show()
    for i in theta:
        out[:, i] = cal_unit(0., 75., i, angle_rot=60).flatten()
    plt.plot(theta, out[0, :], label='I1')
    plt.plot(theta, out[1, :], label='I2')
    plt.plot(theta, out[2, :], label='I3')
    plt.plot(theta, out[3, :], label='I4')
    plt.title('rot = 60 deg')
    plt.legend()
    plt.show()

def rotation_matrix(angle_rot):
    c, s = np.cos(2*angle_rot*np.pi/180), np.sin(2*angle_rot*np.pi/180)
    return np.matrix([[1, 0, 0, 0], [0, c, s, 0], [0, -s, c, 0], [0, 0, 0, 1]])

def test_rotation_matrix():
    return


def pol_lin(alpha):
    #OJO ALPHA IS = 0 WHEN THE POLARIZER IS HORIZONTAL!!!! Vertical = > THETA = 90
    #this is a right handed system
    # y
    # |
    # |
    # |
    # --------- -> X

    angle = alpha * np.float64(np.pi/180.)
    c2 = np.cos(2*angle, dtype=np.float64)
    s2 = np.sin(2*angle, dtype=np.float64)

    pl = 0.5*np.matrix([[1 , c2    , s2   , 0],
                        [c2, c2**2 , c2*s2, 0],
                        [s2, c2*s2 , s2**2, 0],
                        [0 , 0     , 0    , 0]])

    return pl


def rotating_wp(alpha,retardance):
    '''
    a retarder
    '''
    #OJO ALPHA IS = 0 fast axis horizontal  
    #this is a right handed system
    # y
    # |
    # |
    # |
    # --------- -> X

    c2 = np.cos(2*alpha * np.pi/180.)
    s2 = np.sin(2*alpha * np.pi/180.)
    cr = np.cos(retardance * np.pi/180.)
    sr = np.sin(retardance * np.pi/180.)

    rt = np.matrix([[1,              0,              0,      0],
                    [0, c2**2+s2**2*sr, c2*s2*(1-cr)  , -s2*sr],
                    [0, c2*s2*(1-cr)  , c2**2+c2**2*sr, c2*sr ],
                    [0, s2*sr         , c2*sr         , cr    ]])

    return rt


def test_rotating_wp(retardance=150.):

    theta = np.arange(0, 360, 1)
    polange = np.arange(0, 90, 0.25)
    out = np.zeros((360, 4))
    for i in range(len(theta)):
        out[i, :] = np.matmul(pol_lin(polange[0]), rotating_wp(
            theta[i], retardance))[0, :]
    plt.plot(theta, out[:, 0], label='I1')
    plt.plot(theta, out[:, 1], label='I2')
    plt.plot(theta, out[:, 2], label='I3')
    plt.plot(theta, out[:, 3], label='I4')
    #TODO ESTABA POR AQUI
    I2 = 0.5*(np.cos(2*theta*np.pi/180.)**2+np.sin(2*theta*np.pi /
                                               180.)**2*np.sin(retardance*np.pi/180))
    I3 = 0.5*np.cos(2*theta*np.pi/180.)*np.sin(2*theta*np.pi /
                                               180.)*(1-np.cos(retardance*np.pi/180))
    I4 = -0.5*np.sin(2*theta*np.pi/180.)*np.sin(retardance*np.pi/180)
    plt.plot(theta, I2, label='I2b', linestyle='--')
    plt.plot(theta, I3, label='I3b', linestyle='--')
    plt.plot(theta, I4, label='I4b', linestyle='--')
    plt.title('rot = '+str(retardance)+' deg')
    plt.legend()
    plt.show()

def lcvr(theta,delta):

    '''
    lcvr is just a retarder
    '''
    #quarter-wave plate, fast axis horizontal  THETA = 0, DELTA = 90
    #quarter-wave plate, fast axis vertical  THETA = 90, DELTA = 90
    #theta is fast axis, delta retardance
    #angles are measured positive counterclockwise
    #z axis positive direction is along the direction of light propagation
    #this is a right handed system
    # y
    # |
    # |
    # |
    # --------- -> X
    # angle in degree (changed internally into rad)
    # matrix [row,columns]

    deltar = np.deg2rad(delta)
    thetar = np.deg2rad(theta)
    c2 = np.cos(2*thetar)
    s2 = np.sin(2*thetar)

    return np.matrix([[1 , 0                         , 0                          , 0                  ],
                      [0 , c2**2+s2**2*np.cos(deltar), c2*s2*(1-np.cos(deltar))   , -s2*np.sin(deltar) ],
                      [0 , c2*s2*(1-np.cos(deltar))  , s2**2+c2**2*np.cos(deltar) , c2*np.sin(deltar)  ],
                      [0 , s2*np.sin(deltar)         , -c2*np.sin(deltar)         , np.cos(deltar)    ]])

def instrument_model(pardata=pardata):

    # get variables from global if not given
    delta1 = pardata['delta1']
    delta2 = pardata['delta2']
    theta1 = pardata['theta1']
    theta2 = pardata['theta2']
    pol_angle = pardata['pol_angle']
    rot_inst = pardata['rot_inst']

    mod_matrix = np.zeros((4,4),dtype=np.float64)
    for i in range(4):
        LC1 = lcvr(theta1,delta1[i])
        LC2 = lcvr(theta2,delta2[i])
        PL = pol_lin(pol_angle)
        MR = rotation_matrix(rot_inst)
        dummy = np.matmul(np.matmul(PL, np.matmul(LC2, LC1)),MR)
        mod_matrix[i, :] = dummy[0,:]
    return mod_matrix

#   ; ventana = Depo + retarder
#   ; lenses =  Nothing
#   ;  m123 = mirror(45,0.97)##mirror(22,0.97)##mirror(45,0.97)
#   ; mirror = function mirror,theta,T ; tres a 45 grados -> espejo,t,angle
#   ; etalon = Fran equation
#   ;  MR=rotacion(angrot)
#   ;  TELESCOPE = LCVR(10d0,angrot)


def pol_cal_model(alpha=0, delta=75, theta=0,pardata=pardata, angle_rot=None,plot=None):

    if not isinstance(theta, np.ndarray):
        theta=np.arange(0, 360, 1)
    '''
    #TODO implement derivatives (it is in lib_pmp_v1.4)
    #theta = off set of the retarder
    #delta = retardance
    #alpha = off set of the polarizer(angle)
    '''

#    alpha = 0
#    delta = 75
#    theta = np.arange(0, 360, 1)
#    polange = np.arange(0, 90, 0.25)
    out = np.zeros((len(theta),4))
    pm = 2.*instrument_model(pardata=pardata)
    ipm = svd_solve(pm)
    for i in range(len(theta)):
        dummy = np.matmul(pm, cal_unit(
            alpha, delta, theta[i], angle_rot=angle_rot))
        out[i, :] = dummy.flatten()
    if plot is None:
        pass
    else:
        plt.plot(theta, out[:, 0], label='I1')
        plt.plot(theta, out[:, 1], label='I2')
        plt.plot(theta, out[:, 2], label='I3')
        plt.plot(theta, out[:, 3], label='I4')
        plt.legend()
        plt.show()
    return pm, out


@timeit
def random_model(sample=10000, theta=0):
    if len(theta) == 0:
        theta = np.arange(0, 360, 1)
    arr = np.random.rand(12, sample)*360
    modeli = np.zeros((sample, len(theta), 4))
    modelo = np.zeros(( sample, 4, 4))

    for i in tqdm(range(sample)):
        pardata['delta1'] = arr[0:4, i]
        pardata['delta2'] = arr[4:8, i]
        pardata['rot_inst'] = arr[8, i]
        pardata['theta1'] = 0#arr[9, i]/360*1-0.5
        pardata['theta2'] = 45#arr[10, i]/360*1+44.5
        pardata['pol_angle'] = 90#arr[11, i]
        a, b = pol_cal_model(pardata=pardata, theta=theta)

        modeli[i,:,:] = b
        modelo[i,:,:] = a

#np.random.shuffle(arr)
    return modeli, modelo


def mirror(theta,t):
    '''
    Capitani et al, 1989, Sol. Phys., 120, 173
    '''
    a = theta * np.pi/180.

    return [[ (t**2. + 1) , (t**2. - 1) ,     0      ,     0      ],
        [  (t**2. - 1) , (t**2. + 1) ,     0      ,     0      ],
        [       0      ,      0      , 2*t*np.cos(a) , 2*t*np.sin(a) ],
        [       0      ,      0      , -2*t*np.sin(a), 2*t*np.cos(a) ]] / 2.

def smirror(Phase_shift,R, Angle_incidence):
    '''
    '''
    #reflectance = r
    #Phase_shift
    alpha = Phase_shift * np.pi/180. #small angle app
    #angle of incidence
    beta = Angle_incidence * np.pi/180.

    a = (R+1.)/2.
    b = (R-1.)/2.

    return [[ a          , b*np.cos(2.*beta)                                                 ,    b*np.sin(2.*beta)                                              ,     0                                     ],
    [  b*np.cos(2.*beta) , a*np.cos(2.*beta)**2.-np.sqrt(R)*np.cos(alpha)*np.sin(2.*beta)**2 , (a+np.sqrt(R)*np.cos(alpha))*np.sin(4.*beta)/2.                   ,  np.sqrt(R)*np.sin(alpha)*np.sin(2.*beta) ],
    [  b*np.cos(2.*beta) , (a+np.sqrt(R)*np.cos(alpha))*np.sin(4.*beta)/2.                   , a*np.sin(2.*beta)^2.-np.sqrt(R)*np.cos(alpha)*np.cos(2.*beta)**2. , -np.sqrt(R)*np.sin(alpha)*np.cos(2.*beta) ],
    [       0            , -np.sqrt(R)*np.sin(alpha)*np.sin(2.*beta)                         , np.sqrt(R)*np.sin(alpha)*np.cos(2.*beta)                          , -np.sqrt(R)*np.cos(alpha)                 ]]
