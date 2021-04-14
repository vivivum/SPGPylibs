Solar Orbiter SPK files
===========================================================================

   This ``aareadme.txt'' file describes the contents of the kernels/spk
   directory of the Solar Orbiter SPICE Kernel Dataset.

   It was last modified on July 4, 2019 by Marc Costa Sitja, ESAC/ESA.


Brief Summary
--------------------------------------------------------

   This directory contains the SPICE SP-Kernel files for the SOLAR ORBITER
   mission, including mission analysis, nominal and operational spacecraft
   trajectory SPKs, generic planetary and satellite ephemeris SPKs.


File naming conventions
--------------------------------------------------------

   Naming Scheme for Solar Orbiter spacecraft SPKs:

     The naming scheme for the Solar Orbiter spacecraft trajectory SPKs is:

           solo_ANC_soc-orbit_YYYYMMDD-YYYYMMDD_[VOEM]_VNN.bsp

     where

           YYYYMMDD  coverage start and stop times in TDB (required);

           VOEM      reference to the source OEM file version (optional);

           NN        version number, starting from 01 (required; e.g. 01)
                     from the source VOEM.


   Naming Scheme for Generic Planetary Ephemeris SPKs

     The naming scheme for generic planetary SPKs is:

           deNNN.bsp

     where

           NNN       DE version (required; e.g. 421);


   Naming Scheme for ESA ESTRACK ground stations SPKs

     The naming scheme for ESA ESTRACK ground stations SPKs is:

           estrack_vNN.bsp

     where

           NN       version (required; e.g. 01);


   Naming Scheme for ESA New Norcia ground station SPKs

     The naming scheme for ESA New Norcia ground stations SPKs is:

           new_norcia.bsp


Current SPK Kernels Set
--------------------------------------------------------

   solo_ANC_soc-orbit_YYYYMMDD-YYYYMMDD_V01.bsp

      SPICE SPK file containing CrEMA candidate trajectories from Mission
      Analyis for the complete mission.


   de432s.bsp

      SPICE SPK file containing JPL planetary ephemerides version DE432,
      created by NAIF, JPL.


   mar097_YYYYMMDD_YYYYMMDD.bsp

      SPICE SPK file containing JPL Martian satellite ephemerides version
      MAR097, created by NAIF, JPL.


   estrack_vNN.bsp

      SPICE SPK file that defines the position for each of the ESA ESTRACK
      ground stations, as seen from the center of the Earth, in the
      topocentric frame defined in the corresponding ESTRACK FK.


   new_norcia.bsp

      SPICE SPK file that defines the position of the ESA New Norcia
      ground station, as seen from the center of the Earth, in the
      topocentric frame defined in the corresponding New Norcia FK.


Other directory contents
--------------------------------------------------------

   aareadme.txt         This file.


Particulars
--------------------------------------------------------

     Nothing to report.


Kernel File Details
--------------------------------------------------------

     The most detailed description of the data in a binary SPK file is
     provided in metadata included inside the comment area of the file.
     This information can be viewed using the utility programs COMMNT and
     SPACIT included in the NAIF Toolkit.


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

     1. ``SPK Required Reading'', NAIF Document

     2. ``SOC-Provided Ancillary Data for Solar Orbiter'',
        SOL-SGS-TN-0017, A. Walsh, Issue 0, Revision 3, 22/09/2018


End of aareadme file.