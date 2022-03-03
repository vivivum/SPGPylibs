import sys
sys.path.append('../SPGPylibs/')
import SPGPylibs as spg

def pipe_test():

    '''
    This is an example use of the program phifdt_pipe.py
    First look the info in file phifdt_pipe.py for IN/OUT keywords and options.
    '''

    import pathlib
    #lib_path = pathlib.Path().resolve()
    lib_path = pathlib.Path(__file__).parent.resolve()
    lib_path = str(lib_path)

    
    input_data_dir = lib_path + '/data/'
    list_of_files = sorted(spg.list_fits(inpath = input_data_dir,contain='L1',remove_dir=True))
    for i in list_of_files:
        print(i)

    dark_f = lib_path + '/cal_files/solo_L0_phi-hrt-ilam_20210428T130238_V202106300924C_0164281001_dark.fits'
    flat_f = lib_path + '/cal_files/flats_kll_23April2021_June24_spg5.fits'
    prefilter = lib_path + '/cal_files/0000990710_noMeta.fits'

    for i in list_of_files:
        spg.phifdt_pipe(data_f=i,input_data_dir = input_data_dir, dark_f = dark_f,flat_f = flat_f,prefilter=True,
            instrument = 'FDT40', prefilter_fits = prefilter,shrink_mask = 1,correct_ghost=False, verbose = False,
            correct_fringes='manual',ItoQUV=True,center_method = None,output_dir = lib_path + '/red/',rte = 'RTE')
                

if __name__ == "__main__":

        pipe_test()
