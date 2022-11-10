
import json
def json_generator():

    json_file = './HRT.json'
    dict_HRT = { 

        "__comment1__":"data and calibration data goes here",
        "data_f" : "solo_L1_phi-hrt-ilam_20201117T170209_V202108301639C_0051170001.fits",
        "flat_f" : "/data/slam/home/sinjan/fits_files/april_avgd_2020_flat.fits",
        "dark_f" : "../fits_files/solo_L0_phi-fdt-ilam_20200228T155100_V202002281636_0022210004_000.fits",
        
        "__comment2__":"input/output type + scaling goes here",
        'L1_input' : True, 
        'L1_8_generate': False, 
        'scale_data' : True,  
        'accum_scaling' : True, 
        'bit_conversion' : True, 
        
        "__comment3__":"Data reduction options go here",
        'dark_c' : True,
        'flat_c' : True, 
        'norm_f' : True, 
        'clean_f' : True, 
        'sigma' : 59, 
        'clean_mode' : "V", 
        'flat_states' : 24, 
        'prefilter_f': None,
        'fs_c' : True, 
        'demod' : True, 
        'norm_stokes' : True, 
        'ItoQUV' : False,
        'ctalk_params' : None, 
        'rte' : False, 
        'p_milos' : False, 
        
        "__comment4__":"Output options go here",
        'output_dir' : './',  
        'out_demod_file' : False,  
        'out_demod_filename' : None, 
        'out_rte_filename' : None,  
        'config_file' : True 
    }

    with open(json_file, 'w', encoding="utf-8", newline='\r\n') as outfile:
        json.dump(dict_HRT, outfile, indent=4, ensure_ascii=False)



        #     prefilter = True, prefilter_fits = '0000990710_noMeta.fits',
        #     realign = False, verbose = True, outfile=None, mask_margin = 2, correct_fringes = False,
        #     individualwavelengths = False,correct_ghost = False,putmediantozero=True,directory = './',
        #     rte = False, debug = False,nlevel = 0.3,loopthis=0,
        #     cross_talk_IQUV = False, cross_talk_VQU = False, do2d = 0,output_dir='./'):

    # data_f -> input data (single file for FDT - )
    # dark_f -> dark file (scaling calculated inside pipe) 
    # flat_f -> flat file (no scaling - just normalized to one)
    # dark_c -> correct dark or not (default true)
    # flat_c -> correct flat or not (default true)
    # instrument = 'FDT40' 
    # hough_params = [...,...,...] -> three values array
    #   inner_radius = 250, outer_radius = 600, steps = 100 : initial values for finding sun center
    #   if not present or none uses data header information
    #   default is set to [250,600,100] 
    # center_method = ['circlefit','hough']
    #   Default is 'circlefit'. If set to 'hough' uses the given hough_params parameters.
    #   If find_center is set to None then uses header information, in any case.

    json_file = './FDT.json'
    fits_version = '01'
    dict_FDT = { 

        "__comment1__":"data and calibration data goes here",
        "data_f" : "solo_L1_phi-hrt-ilam_20201117T170209_V202108301639C_0051170001.fits",
        "flat_f" : "/data/slam/home/sinjan/fits_files/april_avgd_2020_flat.fits",
        "dark_f" : "../fits_files/solo_L0_phi-fdt-ilam_20200228T155100_V202002281636_0022210004_000.fits",
        'input_data_dir' : './',  
        'verbose' : False,
        "__comment2__":"Data reduction options go here",
        'dark_c' : True,
        'flat_c' : True,
        'instrument' : 'FDT40',
        'hough_params' : None,
        'center_method' : 'circlefit',
        'shrink_mask' : 2,
        'norm_f' : False, 
        'flat_scaling' : 1,
        'flat_index' : False,
        'prefilter': False,
        'prefilter_fits' : '0000990710_noMeta.fits',
        'rte' : False, 
        'correct_fringes' : 'manual',
        'correct_ghost' : False,
        'putmediantozero' : True,
        'ItoQUV' : False,
        'VtoQU' : False,
        'realign' : False,
        'ind_wave' : False,
        'nlevel' : 0.3,
        'debug' : False,
        'vers': fits_version,        #desired version number, only 2 characters 01 -> 99, if not specified, '01' default
        'RTE_code': 'cmilos',        #specify which code will be used in the inversion. Default is cmilos (uses .txt ASCII input files)
        "__comment3__":"Output options go here",
        'output_dir' : './',  

    }

    with open(json_file, 'w', encoding="utf-8", newline='\r\n') as outfile:
        json.dump(dict_FDT, outfile, indent=4, ensure_ascii=False)
