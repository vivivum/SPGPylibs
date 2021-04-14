Solar Orbiter SCLK files
===========================================================================

   This ``aareadme.txt'' file describes the contents of the kernels/sclk
   directory of the Solar Orbiter SPICE Kernel Dataset.

   It was last modified on July 4, 2019 by Marc Costa Sitja, ESAC/ESA.


Brief Summary
--------------------------------------------------------

   This directory contains the SPICE Spacecraft Clock-Kernel files for the
   SOLAR ORBITER spacecraft on-board clock.


File naming conventions
--------------------------------------------------------

   The naming scheme for SOLAR ORBITER SCLKs is:

        solo_ANC_soc-sclk-[fict]_YYYYMMDD_V01.tsc

   where

        TYPE     SCLK data type contained within the file:

                 fict:   ``fictional'' data (optional);

        YYYYMMDD start time of applicability, i.e. either release date for
                 fictional SCLKs or date of the last TCP used to generate
                 the SCLK in step mode.


Current SCLK Kernels Set
--------------------------------------------------------

   solo_ANC_soc-sclk_YYYYMMDD_V01.tsc

      SPICE SCLK file containing time correlation data for the SOLAR ORBITER
      on-board clock. Created by the ESA SPICE Service (ESS).


  solo_ANC_soc-sclk-fict_20000101_V01.tsc

      SPICE SCLK file assuming an absolutely constant clock rate for
      SOLAR ORBITER. This file is used for working with planning data.
      Created by the ESA SPICE Service (ESS).


Other directory contents
--------------------------------------------------------

   aareadme.txt         This file.


Particulars
--------------------------------------------------------

   Nothing to report.


Kernel File Details
--------------------------------------------------------

   The most detailed description of the data in an SCLK file is provided
   in metadata included inside the descriptive text areas of the file.
   This information can be viewed using any text editor.


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

   1.   ``SCLK required Reading'', NAIF Document

   2.   ``TIME Required Reading'', NAIF Document

   3.   ``KERNEL Pool Required Reading'', NAIF Document

   4.   ``SOC-Provided Ancillary Data for Solar Orbiter'',
        SOL-SGS-TN-0017, A. Walsh, Issue 0, Revision 3, 22/09/2018


End of aareadme file.