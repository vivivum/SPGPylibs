import numpy as np
from itertools import combinations

from typing import List,Tuple

def bit_mask_serie(image_list: List[np.ndarray], blos: bool = False,method: str = None) -> Tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Generators have a ``Yields`` section instead of a ``Returns`` section.

    Args:
        n (int): The upper limit of the range to generate, from 0 to `n` - 1.

    Returns:
        The return value. True for success, False otherwise.


    Examples:
        Examples should be written in doctest format, and should illustrate how
        to use the function.

        >>> print([i for i in example_generator(4)])
        [0, 1, 2, 3]

    """
    y,x = image_list[0].shape
    n = len(image_list)
    mask = [np.zeros([y,x], dtype=np.int8) for i in image_list]
    cmask = np.ones((y,x), dtype=np.int8)
    mean = np.zeros((y,x), dtype=np.float32)
    #determine rms of image central area
    rg = 100 
    if blos == True:
        rms = 10
    else:
        rms = image_list[0][int(y/2-rg):int(y/2+rg),int(x/2-rg):int(x/2+rg)].std()

    print(rms)
    #determine mask of each image
    for i in range(len(image_list)): 
        mask[i][np.where(np.abs(image_list[i]) > rms)] = 1

    for iter in combinations(mask,2):
        cmask *= iter[0] & iter[1]

    #cmask contains only the common stuff

    if method == 'JBR':
        for i in image_list: 
            mean += np.abs(i)
        mean /= np.mean(mean)
        cmask = cmask * mean
        cmask[np.where(cmask) == 0] = 1
    else:
        for i in image_list: 
            mean += i
        mean /= len(image_list)
        flat = mean * cmask

    return mask,cmask,flat

def _bit_mask_serie_test(image_list,blos = False):

    pass
    #mask,cmask,flat = bit_mask_serie(bb,blos=True)

if __name__ == "__main__":

    _bit_mask_serie_test()
