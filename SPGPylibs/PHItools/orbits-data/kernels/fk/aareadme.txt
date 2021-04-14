Solar Orbiter FK files
===========================================================================

   This ``aareadme.txt'' file describes the contents of the kernels/fk
   directory of the Solar Orbiter SPICE Kernel Dataset.

   It was last modified on July 4, 2019 by Marc Costa Sitja, ESAC/ESA.


Brief Summary
--------------------------------------------------------

   This directory contains the SPICE Frames Definition Kernel files for the
   SOLAR ORBITER spacecraft, its structures, and science instruments, as
   well as for the Ground Earth Stations. Moreover, it contains frame
   definitions for science operations and for data analysis and scientific
   investigations.


File naming conventions and Current FK Kernels Set
--------------------------------------------------------

   The following table provides information on which file contains which
   frames, their file naming conventions and some particular details:

      File                         Contents
     ---------------------------  -----------------------------------------
      earthfixediau.tf             Makes the IAU_EARTH coincide with
                                   the Earth fixed reference frame.

      estrack_vNN.tf               ESA Ground Stations topocentric frames.

      new_norcia_topo.tf           ESA New Norcia topocentric frame.

      earth_topo_050714.tf         NASA DSN Ground Stations topocentric
                                   frames.

      solo_ANC_soc-sci-fk_VNN.tf   SOLAR ORBITER relevant scientific
                                   reference frames.

      solo_ANC_soc-ops-fk_VNN.tf   SOLAR ORBITER  science operations relevant
                                   frames.

      solo_ANC_soc-sc-orb-norm-fk_VNN.tf

                                   Mapping/overwrite of SOLO_PRF to
                                   SOLO Orbit-Normal Pointing.

      solo_ANC_soc-sc-eclip-norm-fk_VNN.tf

                                   Mapping/overwrite of SOLO_PRF to
                                   SOLO Ecliptical-Normal Pointing.

      solo_ANC_soc-sc-equat-norm-fk_VNN.tf

                                   Mapping/overwrite of SOLO_PRF to
                                   SOLO Equatorial-Normal Pointing.

      solo_ANC_soc-sc-fk_VNN.tf    SOLAR ORBITER spacecraft, structures,
                                   instruments and sensors frames.

   where

           NN       version number -- two digits (required)

                    If multiple versions of a Frames Kernel file are
                    provided, always use the latest version (unless
                    earlier version is needed for some special reasons.)


Other directory contents
--------------------------------------------------------

    aareadme.txt         This file.


Particulars
--------------------------------------------------------

    Important information about the Mapping/overwrite frame definition
    kernels (solo_ANC_soc-sc-orb-norm-fk_VNN.tf,
    solo_ANC_soc-sc-eclip-norm-fk_VNN.tf and
    solo_ANC_soc-sc-equat-norm-fk_VNN.tf):

       These frames have been implemented to overwrite the SOLO_PRF
       frame definition provided in the Solar Orbiter Frames Definitions
       kernel (  solo_ANC_soc-sc-fk_VNN.tf ) and map it to the
       SOLO_ORBIT_NORM, SOLO_ECLIP_NORM and SOLO_EQUAT_NORM frames defined
       in the Solar Orbiter Science Operations Frames Definitions kernel
       (solo_ANC_soc-ops-fk_VNN.tf).

       In order to make use of these frames' kernel, ONLY ONE of these files
       MUST BE LOADED AFTER the SOLAR ORBITER frames definition kernel and
       the Solar Orbiter Science Operations Frames Definition kernel.

       NOTE THAT BY USING ANY OF THESE KERNEL, THE SOLO_PRF FRAME
       WILL BE  MAPPED TO THE CORRESPONDING SCIENCE OPERATIONS FRAME, AND
       ANY CK PROVIDING ORIENTATION FOR THE SOLO_PRF FRAME WILL
       NOT BE USED BY THE APPLICATION SOFTWARE, EVEN IF IT IS LOADED IN
       THE KERNEL POOL.


Kernel File Details
--------------------------------------------------------

   The most detailed description of the data in an FK file is provided in
   metadata included inside the descriptive text areas of the file. This
   information can be viewed using any text editor.


Contact Information
--------------------------------------------------------

   If you have any questions regarding this file contact the
   ESA SPICE Service (ESS) at ESAC:

           Marc Costa Sitja
           (+34) 91-8131-457
           esa_spice@sciops.esa.int, marc.costa@esa.int

   or the Solar Orbiter Science Operations Center at ESAC:

           sol_soc@esa.int


References and required readings
--------------------------------------------------------

   1. ``Frames Required Reading'', NAIF Document

   2. ``Kernel Pool Required Reading'', NAIF

   3. ``C-Kernel Required Reading'', NAIF

   4. ``SOC-Provided Ancillary Data for Solar Orbiter'',
      SOL-SGS-TN-0017, A. Walsh, Issue 0, Revision 3, 22/09/2018


End of aareadme file.