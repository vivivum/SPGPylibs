def idlsmooth(image,size=3):
    '''
    IDL smooth function smooth(imahge,size,/edge_truncate)
    '''
    from scipy import ndimage, misc
    return ndimage.uniform_filter(image, size=size)
