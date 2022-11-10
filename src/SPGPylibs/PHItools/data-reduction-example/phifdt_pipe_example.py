try:
    import SPGPylibs as spg
    print('Using local installation')
except:
    import sys
    sys.path.append('../SPGPylibs/')
    import SPGPylibs as spg
    print('Using local distribution')

def pipe_test():

    '''
    This is an example use of the program phifdt_pipe.py
    First look the info in file phifdt_pipe.py for IN/OUT keywords and options.
    '''

    import pathlib
    #lib_path = pathlib.Path().resolve()
    lib_path = pathlib.Path(__file__).parent.resolve()
    lib_path = str(lib_path)+'/'
    
    dir = '/Users/orozco/Dropbox_folder/Python/VS-GitHub/SPGPylibs_run/Pipeline-run-tests/data-reduction-example'
    input_data_dir = dir + '/data/'
    list_of_files = sorted(spg.list_fits(inpath = input_data_dir,contain='L1',remove_dir=True))
    for i in list_of_files:
        print(i)

    dark_f = dir + '/cal_files/solo_L0_phi-hrt-ilam_20210428T130238_V202106300924C_0164281001_dark.fits'
    flat_f = dir + '/cal_files/flats_kll_23April2021_June24_spg5.fits'
    prefilter = dir + '/cal_files/0000990710_noMeta.fits'
    output_dir = dir + '/red/'

    print(lib_path)
    spg.phifdt_pipe(json_input = lib_path+'FDT_test.json')
    return
    #only one
    spg.phifdt_pipe(data_f=list_of_files[0],input_data_dir = input_data_dir, dark_f = dark_f,flat_f = flat_f,prefilter=False,
        instrument = 'FDT40', prefilter_fits = prefilter,shrink_mask = 1,correct_ghost=False, verbose = False,
        correct_fringes=False,ItoQUV=True,center_method = None,output_dir = output_dir,rte = 'RTE',RTE_code="cmilos",
        vers = '01')

    # for i in list_of_files:
    #     spg.phifdt_pipe(data_f=i,input_data_dir = input_data_dir, dark_f = dark_f,flat_f = flat_f,prefilter=False,
    #         instrument = 'FDT40', prefilter_fits = prefilter,shrink_mask = 1,correct_ghost=False, verbose = False,
    #         correct_fringes=False,ItoQUV=True,center_method = None,output_dir = lib_path + '/red/',rte = 'RTE',cmilos="pymilos")
    #     return
                
    #test with json

if __name__ == "__main__":

        pipe_test()
