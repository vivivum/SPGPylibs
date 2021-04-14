Solar Orbiter IK files
===========================================================================

   This ``aareadme.txt'' file describes the contents of the kernels/ik
   directory of the Solar Orbiter SPICE Kernel Dataset.

   It was last modified on July 4, 2019 by Marc Costa Sitja, ESAC/ESA.


Brief Summary
--------------------------------------------------------

   This directory contains the SPICE Instrument Kernel files for the
   SOLAR ORBITER remote sensing instruments.


File naming conventions and Current FK Kernels Set
--------------------------------------------------------

   The naming scheme for SOLAR ORBITER IKs is:

         solo_ANC_soc-SENSOR-ik_vNN.ti

   where

          SENSOR   Solar Orbiter instrument identifier(required):

                       edp:      Energetic Particle Detector (EPD)
                       eui:      Extreme Ultraviolet Imager (EUI)
                       metis:    Multi Element Telescope for Imaging and
                                 Spectroscopy (Metis)
                       phi:      Polarimetric and Helioseismic Imager (PHI)
                       solohi:   SOLO Heliospheric Imager (SOLOHI)
                       spice:    Spectral Imaging of the Coronal Environment
                                 (SPICE)
                       stix:     Spectrometer Telescope for Imaging X rays
                                 (STIX)
                       str:      Star Trackers (STR)
                       swa:      Solar Wind Analyzer (SWA)

         NN       version number (required);

                  If multiple versions of an IK file are provided for an
                  instrument, always use the latest version (unless
                  earlier version is needed for some special reasons.)


Other directory contents
--------------------------------------------------------

   aareadme.txt         This file.


Particulars
--------------------------------------------------------

   Nothing to report.


Kernel File Details
--------------------------------------------------------

   The most detailed description of the data in an IK file is provided in
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

   1. ``Kernel Pool Required Reading'', NAIF

   2. Header of GETFOV SPICE API, latest version

   3. ``SOC-Provided Ancillary Data for Solar Orbiter'',
      SOL-SGS-TN-0017, A. Walsh, Issue 0, Revision 3, 22/09/2018


End of aareadme file.