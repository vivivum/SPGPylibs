import sys
sys.path.append('../SPGPylibs/')
import SPGPylibs as spg
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits

def pipe_test():

    '''
    This is an example use of the program phifdt_pipe.py
    First look the info in file phifdt_pipe.py for IN/OUT keywords and options.


    '''
    wavelength = 6173.3356
    prefilter = '/Users/orozco/Dropbox_folder/SoPHI/PHI-COMMISSIONING/software-and-images/RSW1/0000990710_noMeta.fits'


    # data_f = '/Users/orozco/Dropbox_folder/SoPHI/PHI-COMMISSIONING/software-and-images/Nov-2020-STP122/solo_L0_phi-fdt-ilam_0658980568_V202012012104C_0051180402.fits'
    dir = 'nov-data-wcs/'
    # with pyfits.open(data_f) as hdu_list:
    #     datos = hdu_list[0].data
    # with pyfits.open(dir+'wcs.fits') as hdu_list:
    #     hdu_list[0].data = datos
    #     hdu_list.writeto(dir+'l0_wcs.fits', clobber=True)
    # return
    data_f = dir + 'l0_wcs.fits'
    dark_f = '/Users/orozco/Dropbox_folder/SoPHI/PHI-COMMISSIONING/software-and-images/RSW1/solo_L0_phi-fdt-ilam_20200618T000547_V202006221044C_0066181001_dark.fits'
    flat_f = '/Users/orozco/Dropbox_folder/Python/VS-GitHub/SPGPylibs_tests/flats_kll_23April2021_June24_spg2.fits'
    # hdu_list.writeto('flats_kll_23April2021_June24_spg2.fits',overwrite=True) #kll + expand 20
    # hdu_list.writeto('flats_kll_23April2021_June24_spg3.fits',overwrite=True) #chae + expand 20  NO BUENO
    # hdu_list.writeto('flats_kll_23April2021_June24_spg4.fits',overwrite=True) #kll + circle + expand 20
    # hdu_list.writeto('flats_kll_23April2021_June24_spg5.fits',overwrite=True) #kll + circle + expand 20 + ghost
    # hdu_list.writeto('flats_kll_23April2021_June24_spg6.fits',overwrite=True) #kll + circle + expand 20 + normalize NO BUENO

    flat_f = '/Users/orozco/Dropbox_folder/Python/VS-GitHub/SPGPylibs_tests/flats_kll_23April2021_June24_spg2.fits'
    outfile = 'solo_L1_phi-fdt-ilam_0658980568_V202012012104C_0051180402' #test CI
    spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=False,directory = dir,
        rte = False,correct_ghost=False,putmediantozero = True,fringes = False,
        individualwavelengths = 0,vqu = 0, do2d = 0,verbose = False, debug = False,instrument = 'FDT40',
        prefilter_fits = prefilter,mask_margin = 1,realign = 0,flat_c=True)
    # spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=True,directory = dir,
    #     rte = 'CE+RTE',correct_ghost=True,putmediantozero = True,fringes = False,
    #     individualwavelengths = 0,vqu = 0, do2d = 0,verbose = False, debug = False,instrument = 'FDT40',
    #     prefilter_fits = prefilter,mask_margin = 1,realign = 0,flat_c=True)

    return
    flat_f = '/Users/orozco/Dropbox_folder/Python/VS-GitHub/SPGPylibs_tests/flats_kll_23April2021_June24_spg3.fits'
    outfile = 'Nov_map2' #test CI
    dir = './'
    spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=False,directory = dir,
        rte = 'CE',correct_ghost=True,putmediantozero = 0,fringes = False,
        individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
        prefilter_fits = prefilter,mask_margin = 6,realign = 0,flat_c=True)
    flat_f = '/Users/orozco/Dropbox_folder/Python/VS-GitHub/SPGPylibs_tests/flats_kll_23April2021_June24_spg4.fits'
    outfile = 'Nov_map3' #test CI
    dir = './'
    spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=False,directory = dir,
        rte = 'CE',correct_ghost=True,putmediantozero = 0,fringes = False,
        individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
        prefilter_fits = prefilter,mask_margin = 6,realign = 0,flat_c=True)
    flat_f = '/Users/orozco/Dropbox_folder/Python/VS-GitHub/SPGPylibs_tests/flats_kll_23April2021_June24_spg5.fits'
    outfile = 'Nov_map4' #test CI
    dir = './'
    spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=False,directory = dir,
        rte = 'CE',correct_ghost=True,putmediantozero = 0,fringes = False,
        individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
        prefilter_fits = prefilter,mask_margin = 6,realign = 0,flat_c=True)
    flat_f = '/Users/orozco/Dropbox_folder/Python/VS-GitHub/SPGPylibs_tests/flats_kll_23April2021_June24_spg6.fits'
    outfile = 'Nov_map5' #test CI
    dir = './'
    spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=False,directory = dir,
        rte = 'CE',correct_ghost=True,putmediantozero = 0,fringes = False,
        individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
        prefilter_fits = prefilter,mask_margin = 6,realign = 0,flat_c=True)

    return
    outfile = 'Nov_map2' #realineados pesos 20,20,10
    outfile = 'Nov_map3' #realineados pesos 10,10,4
    outfile = 'Nov_map4' #change cmilos

    # OJO         
    # shift_w =  wave_axis[3] - wavelength
    # wave_axis = wave_axis - shift_w

    dir = 'spg6/'

    # outfile = 'Nov_map_test1' #Only dark + flat
    # spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=False,directory = dir,
    #     rte = 'RTE',correct_ghost=False,putmediantozero = 0,fringes = False,
    #     individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
    #     prefilter_fits = prefilter,mask_margin = 6,realign = 0,flat_c=True)
    # move_file = subprocess.call("mv dummy_out.txt "+dir+outfile+".txt",shell=True)
    # print(move_file)
    # outfile = 'Nov_map_test2' #Only dark + flat + fringes
    # spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=False,directory = dir,
    #     rte = 'RTE',correct_ghost=False,putmediantozero = 0,fringes = True,
    #     individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
    #     prefilter_fits = prefilter,mask_margin = 6,realign = 0,flat_c=True)
    # move_file = subprocess.call("mv dummy_out.txt "+dir+outfile+".txt",shell=True)
    # print(move_file)
    # outfile = 'Nov_map_test3' #Only dark + flat + fringes + ghost
    # spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=False,directory = dir,
    #     rte = 'RTE',correct_ghost=True,putmediantozero = 0,fringes = True,
    #     individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
    #     prefilter_fits = prefilter,mask_margin = 6,realign = 0,flat_c=True)
    # move_file = subprocess.call("mv dummy_out.txt "+dir+outfile+".txt",shell=True)
    # print(move_file)
    # outfile = 'Nov_map_test4' #Only dark + flat + fringes + ghost + medianzero
    # spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=False,directory = dir,
    #     rte = 'RTE',correct_ghost=True,putmediantozero = 1,fringes = True,
    #     individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
    #     prefilter_fits = prefilter,mask_margin = 6,realign = 0,flat_c=True)
    # move_file = subprocess.call("mv dummy_out.txt "+dir+outfile+".txt",shell=True)
    # print(move_file)
    # outfile = 'Nov_map_test5' #Only dark + flat + fringes + ghost + medianzero + prefilter
    # spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=True,directory = dir,
    #     rte = 'RTE',correct_ghost=True,putmediantozero = 1,fringes = True,
    #     individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
    #     prefilter_fits = prefilter,mask_margin = 6,realign = 0,flat_c=True)
    # move_file = subprocess.call("mv dummy_out.txt "+dir+outfile+".txt",shell=True)
    # print(move_file)
    # outfile = 'Nov_map_test6' #Only dark + flat + fringes + ghost + medianzero + prefilter + individual waves
    # spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=True,directory = dir,
    #     rte = 'RTE',correct_ghost=True,putmediantozero = 1,fringes = True,
    #     individualwavelengths = 1,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
    #     prefilter_fits = prefilter,mask_margin = 6,realign = 0,flat_c=True)
    # move_file = subprocess.call("mv dummy_out.txt "+dir+outfile+".txt",shell=True)
    # print(move_file)
    # outfile = 'Nov_map_test7' #Only dark + flat + fringes + ghost + medianzero + prefilter + individual waves + realign
    # spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=True,directory = dir,
    #     rte = 'RTE',correct_ghost=True,putmediantozero = 1,fringes = True,
    #     individualwavelengths = 1,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
    #     prefilter_fits = prefilter,mask_margin = 6,realign = 1,flat_c=True)
    # move_file = subprocess.call("mv dummy_out.txt "+dir+outfile+".txt",shell=True)
    # print(move_file)
    # outfile = 'Nov_map_test8' #Only dark + flat + fringes + ghost + medianzero + prefilter + realign
    # spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=True,directory = dir,
    #     rte = 'RTE',correct_ghost=True,putmediantozero = 1,fringes = True,
    #     individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
    #     prefilter_fits = prefilter,mask_margin = 6,realign = 1,flat_c=True)
    # move_file = subprocess.call("mv dummy_out.txt "+dir+outfile+".txt",shell=True)
    # print(move_file)

    return

    file1 = dir+outfile+'_blos_rte.fits'
    with pyfits.open(file1) as hdu_list:
        b_los = hdu_list[0].data = b_los
    fileo = 'solo_L1_phi-fdt-blos_20201118T021009_V202105171449C_0051180402.fits'
    with pyfits.open(fileo) as hdu_list:
        hdu_list[0].data = b_los
        hdu_list.writeto(directory+outfile+'_wcs.fits', clobber=True)

    loop = [0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.0,1.01,1.02,1.03,1.04,1.05,1.05,1.07,1.08,1.09,1.1]
    loop =[ 0]
    for i in loop:
        spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=True,
            rte = 'RTE',correct_ghost=True,putmediantozero = 1,fringes = True,
            individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 0, debug = False,instrument = 'FDT40',
            prefilter_fits = prefilter,mask_margin = 6,realign = 1,loopthis=i,flat_c=True)

    return
    #RSW1 SCIENCE SCAN
    data_f = '../RSW1/solo_L0_phi-fdt-ilam_20200618T110546_V202006221044C_0046180001_image.fits'
    dark_f = '../RSW1/solo_L0_phi-fdt-ilam_20200618T000547_V202006221044C_0066181001_dark.fits'
    flat_f = '/Users/orozco/Dropbox_folder/Python/VS-GitHub/SPGPylibs_tests/flats_kll_23April2021_June24_spg.fits'
    outfile = 'sqw1'

    spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=True,
        rte = 'RTE',correct_ghost=True,putmediantozero = 0,fringes = True,
        individualwavelengths = 0,vqu = 0, do2d = 0,verbose = True, debug = False,instrument = 'FDT40',
        prefilter_fits = prefilter,mask_margin = 3,realign = 0,flat_c=True)

    return
    

    #July data (One image)
    data_f = '/Users/orozco/Dropbox_folder/SoPHI/PHI-COMMISSIONING/software-and-images/RSW1/solo_L0_phi-fdt-ilam_20200618T110546_V202006221044C_0046180001_image.fits'
    dark_f = '/Users/orozco/Dropbox_folder/SoPHI/PHI-COMMISSIONING/software-and-images/RSW1/solo_L0_phi-fdt-ilam_20200618T000547_V202006221044C_0066181001_dark.fits'
    flat_f = '/Users/orozco/Dropbox_folder/Python/VS-GitHub/SPGPylibs_tests/flats_kll_23April2021_June24_spg4.fits'
    outfile = 'test2.fits'

    vl = spg.phifdt_pipe(data_f,dark_f,flat_f,outfile=outfile,normalize = 0,prefilter=0,
        inner_radius = 500, outer_radius = 550, steps = 25,rte = True,correct_ghost=0,putmediantozero = 1,
        individualwavelengths = 0,vqu = 0, do2d = 0,verbose = 1, debug = False,instrument = 'FDT40',
        prefilter_fits = prefilter,mask_margin = 10)

    quit()


    #November Spot data (me sirve par aver la rotación)
    data_f = '../Nov-2020-STP122/solo_L0_phi-fdt-ilam_0658980568_V202012012104C_0051180402.fits'
    dark_f = '../RSW1/solo_L0_phi-fdt-ilam_20200618T000547_V202006221044C_0066181001_dark.fits'
    flat_f = 'flatfile.fits'
    vl = run_main(data_f,dark_f,flat_f,outfile='nov-data-preliminar.fits',normalize = 0,prefilter=1,verbose=1,rte=True)#,do2d = 2 )#,cavity = cavity)
    quit()


    quit()

    cavity,h = phi.fits_read('HRT_cavity_map_IP5.fits')
    cavity = cavity * 0.3513e-3/6173.*300000. #A/V 
    cavity = cavity - np.median(cavity[800:1200,800:1200])

    if txt == '-19May':
        print(txt)
        data[19,:,:] = data[18,:,:] 
        #flat[4,:,:] = flat[6,:,:]
    if txt == '-21May':
        print(txt)
        #data[19,:,:] = data[18,:,:] 
        #flat[4,:,:] = flat[6,:,:]
    if txt == '-18Jun':
        print(txt)
        #data[19,:,:] = data[18,:,:] 
        #flat[4,:,:] = flat[6,:,:]


    vl = run_main(data_f,dark_f,flat_f,label=str(0),outfile='rsw1P0flat.fits',index = [5,0,1,2,3,4])#,cavity = cavity)


    dir = '/Users/orozco/Dropbox_folder/SoPHI/PHI-COMMISSIONING/software-and-images/PHI-FDT/'

    files_ip8 = ['solo_L0_phi-fdt-ilam_20200519T124648_V202006061716C_0065250200.fits',
    'solo_L0_phi-fdt-ilam_20200519T125248_V202006061710C_0065250225.fits',
    'solo_L0_phi-fdt-ilam_20200519T125848_V202006061709C_0065250250.fits',
    'solo_L0_phi-fdt-ilam_20200519T130449_V202006061658C_0065250275.fits',
    'solo_L0_phi-fdt-ilam_20200519T131048_V202006061720C_0065250300.fits',
    'solo_L0_phi-fdt-ilam_20200519T131648_V202006061716C_0065250325.fits',
    'solo_L0_phi-fdt-ilam_20200519T132248_V202006061717C_0065250350.fits',
    'solo_L0_phi-fdt-ilam_20200519T132848_V202006061707C_0065250375.fits',
    'solo_L0_phi-fdt-ilam_20200519T133448_V202006061700C_0065250400.fits']
    files_ip8 = [dir + s for s in files_ip8]

if __name__ == "__main__":

        pipe_test()
