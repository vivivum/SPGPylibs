import numpy as np
from .tools import *

@timeit
def shift(matrix, shift=[0, 0], fill_value=0):
    '''Shift operator
    Shift an image in 2D naively as in SOLO-PHI instrument.
    Faster and more efficient methods can be used in normal CPU.
    Input is a vector shift=[x,y] of x and y displacement
    +x -> positive; +y -> positive 
    fill_value = float. 
    This method does not have any boundary condition.
    '''
    try:
        dimy, dimx = matrix.shape
    except:
        raise ValueError("Input is not 2D matrix") 
    
    try:
        nx = shift[1]
        ny = shift[0]
    except:
        raise ValueError("Provided shift not in rigth format 'shift=[0, 0]' of not present") 

    e = np.empty_like(matrix)
    if nx > 0:
        e[:nx, :] = fill_value
        e[nx:, :] = matrix[:-nx, :]
    elif nx < 0:
        e[nx:, :] = fill_value
        e[:nx, :] = matrix[-nx:, :]
    else:
        e = matrix

    s = np.empty_like(matrix)
    if ny > 0:
        s[:, :ny] = fill_value
        s[:, ny:] = e[:, :-ny]
    elif ny < 0:
        s[:, ny:] = fill_value
        s[:, :ny] = e[:, -ny:]
    else:
        s = e

    return s

