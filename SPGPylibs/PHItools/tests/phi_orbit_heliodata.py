import sys
sys.path.append('../SPGPylibs/')
import SPGPylibs as spg
import matplotlib.pyplot as plt

def phi_orbit_heliodata(): 

    from datetime import datetime
    import matplotlib.pyplot as plt

    '''Find all fits in directory.'''
    adonde = '/Users/orozco/Dropbox_folder/SoPHI/PHI-COMMISSIONING/software-and-images/STP-142/helio/'
    list_of_files = spg.list_fits(path = adonde)
    list_of_files = sorted(list_of_files)

    lostimes = []
    for a in list_of_files:
        lostimes.append(spg.fits_get_fpatimes(a))

    solo_heeq = spg.phi_orbit(lostimes, 'Solar Orbiter')#, frame = 'SOLO_HEEQ')
    print(solo_heeq.vr)

if __name__ == "__main__":

        phi_orbit_heliodata()
