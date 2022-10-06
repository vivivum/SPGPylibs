from scipy import ndimage
import numpy as np


def idl_smooth(image: np.ndarray, size: int = 3) -> np.ndarray:
    """
    :param image: nd.array of two dimensions
    :param size: integer indicating the number of pixels to smooth
    :return: nd.ndarray of two dimensions with the smooth image
    """

    return ndimage.uniform_filter(image, size=size)
