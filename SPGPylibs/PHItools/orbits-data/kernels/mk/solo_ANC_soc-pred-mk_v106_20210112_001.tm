KPL/MK

Meta-kernel for Solar Orbiter Dataset v106 -- Predicted 20210112_001
============================================================================

   This meta-kernel lists the Solar_Orbiter Predicted SPICE kernels
   that provide information for the Predicted scenario.

   The kernels listed in this meta-kernel and the order in which
   they are listed are picked to provide the best data available and
   the most complete coverage for the Solar_Orbiter Predicted scenario.

   This meta-kernel was generated with the Auxiliary Data Conversion
   System version: ADCSng v2.5.2.


Usage of the Meta-kernel
---------------------------------------------------------------------------

   This file is used by the SPICE system as follows: programs that make use
   of this kernel must "load" the kernel normally during program
   initialization. Loading the kernel associates the data items with
   their names in a data structure called the "kernel pool". The SPICELIB
   routine FURNSH loads a kernel into the pool.

   The kernels listed below can be obtained from the ESA SPICE FTP server:

      ftp://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels/


Implementation Notes
---------------------------------------------------------------------------

   It is recommended that users make a local copy of this file and
   modify the value of the PATH_VALUES keyword to point to the actual
   location of the Solar_Orbiter SPICE data set's ``data'' directory on
   their system. Replacing ``/'' with ``\'' and converting line
   terminators to the format native to the user's system may also be
   required if this meta-kernel is to be used on a non-UNIX workstation.


-------------------

   This file was created on January 12, 2021 by Marc Costa Sitja ESA/ESAC.
   The original name of this file was solo_ANC_soc-pred-mk_V106_20210112_001.tm.


   \begindata

     PATH_VALUES       = ( '..' )

     PATH_SYMBOLS      = ( 'KERNELS' )

     KERNELS_TO_LOAD   = (

                           '$KERNELS/ck/solo_ANC_soc-sc-iboom-ck_20180930-21000101_V01.bc'
                           '$KERNELS/ck/solo_ANC_soc-sc-oboom-ck_20180930-21000101_V01.bc'
                           '$KERNELS/ck/solo_ANC_soc-sc-fof-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-eui-fsi-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-eui-hri-euv-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-eui-hri-lya-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-metis-euv-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-metis-vis-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-metis-m0-tel-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-phi-fdt-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-phi-hrt-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-solohi-ck_20180930-21000101_V03.bc'
                           '$KERNELS/ck/solo_ANC_soc-spice-sw-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-spice-lw-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-stix-ck_20180930-21000101_V02.bc'
                           '$KERNELS/ck/solo_ANC_soc-default-att_20200210-20301120_L004_V1_00062_V01.bc'

                           '$KERNELS/fk/solo_ANC_soc-sc-fk_V08.tf'
                           '$KERNELS/fk/solo_ANC_soc-ops-fk_V02.tf'
                           '$KERNELS/fk/solo_ANC_soc-sci-fk_V07.tf'
                           '$KERNELS/fk/earth_topo_050714.tf'
                           '$KERNELS/fk/estrack_v03.tf'

                           '$KERNELS/ik/solo_ANC_soc-epd-ik_V02.ti'
                           '$KERNELS/ik/solo_ANC_soc-eui-ik_V01.ti'
                           '$KERNELS/ik/solo_ANC_soc-metis-ik_V03.ti'
                           '$KERNELS/ik/solo_ANC_soc-phi-ik_V01.ti'
                           '$KERNELS/ik/solo_ANC_soc-solohi-ik_V01.ti'
                           '$KERNELS/ik/solo_ANC_soc-spice-ik_V02.ti'
                           '$KERNELS/ik/solo_ANC_soc-stix-ik_V02.ti'
                           '$KERNELS/ik/solo_ANC_soc-swa-ik_V03.ti'

                           '$KERNELS/lsk/naif0012.tls'

                           '$KERNELS/pck/pck00010.tpc'

                           '$KERNELS/pck/earth_070425_370426_predict.bpc'

                           '$KERNELS/sclk/solo_ANC_soc-sclk-fict_20000101_V01.tsc'

                           '$KERNELS/spk/solo_ANC_soc-orbit_20200210-20301120_L004_V1_00062_V01.bsp'
                           '$KERNELS/spk/de421.bsp'
                           '$KERNELS/spk/earthstns_itrf93_050714.bsp'
                           '$KERNELS/spk/estrack_v03.bsp'

                         )

   \begintext


SPICE Kernel Dataset Version
--------------------------------------------------------------------------

   The SPICE Kernel Dataset version of the kernels present in this
   meta-kernel is provided by the following keyword (please note that
   this might not be the last version of the SPICE Kernel Dataset):

   \begindata

      SKD_VERSION = 'v106_20210112_001'

   \begintext


Contact Information
--------------------------------------------------------------------------

   If you have any questions regarding this file contact the
   ESA SPICE Service (ESS) at ESAC:

           Marc Costa Sitja
           (+34) 91-8131-457
           esa_spice@sciops.esa.int, marc.costa@esa.int,


End of MK file.