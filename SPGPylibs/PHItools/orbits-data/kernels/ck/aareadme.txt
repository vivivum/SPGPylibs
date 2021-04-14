Solar Orbiter CK files
===========================================================================

     This ``aareadme.txt'' file describes the contents of the kernels/ck
     directory of the Solar Orbiter SPICE Kernel Dataset.

     It was last modified on January 4, 2020 by Marc Costa Sitja, ESAC/ESA.


Brief Summary
--------------------------------------------------------

     CK (Camera-Matrix Kernel) files contain time varying orientation for
     the SOLAR ORBITER spacecraft, its structures, and science instruments.


File naming conventions
--------------------------------------------------------

   SOLAR ORBITER spacecraft CK:

      The naming scheme for the as-planned default and as-planned predicted
      S/C Roll SOLAR ORBITER S/C CKs are:

          solo_ANC_soc-default-att_YYYYMMDD-YYYYMMDD_VOEM_VNN.bc
          solo_ANC_soc-pred-roll-att_YYYYMMDD-YYYYMMDD_VOEM_VFECS_VNN.bc

      where

          YYYYMMDD  coverage start and stop times in TDB (required);

          VOEM      reference to the source OEM file version (required);

          VFECS     reference to the source FECS file version (required);

          NN        version number, starting from 01 (required; e.g. 01)
                     from the source VOEM.


      The naming scheme for the as-flown SOLAR ORBITER spacecraft CKs is:

          solo_ANC_soc-flown-att_YYYYMMDD[-YYYYMMDD]_VNN.bc

      where

          YYYYMMDD  coverage start time or day in TDB (required);

          YYYYMMDD  coverage stop time or day in TDB (optional);

          NN        version number, starting from 01 (required; e.g. 01)
                     from the source VOEM.


   SOLAR ORBITER structures CK:

      The naming scheme for the SOLAR ORBITER sensors line of sight CKs are:

          solo_ANC_soc-sc-STRUCT-ck_YYYYMMDD-YYYYMMDD_VNN.bc

      where

          STRUCT    Solar Orbiter Line of Sight sensors:

                       boom:  Instrument Boom
                       iboom: Instrument Boom Inboard Segment
                       oboom: Instrument Boom Outboard Segment
                       sa:    Solar Arrays
                       hga:   High Gain Antenna
                       mga:   Medium Gain Antenna
                       fof:   Flight Optical Frame

          YYYYMMDD  coverage start and stop times in TDB (required);

          VOEM      reference to the source OEM file version (optional);

          NN        version number, starting from 01 (required; e.g. 01)
                    from the source VOEM.

   SOLAR ORBITER sensors line of sight CK:

      The naming scheme for the SOLAR ORBITER sensors line of sight CKs are:

          solo_ANC_soc-SENSOR-ck_YYYYMMDD-YYYYMMDD_VNN.bc

      where

          SENSOR    Solar Orbiter Line of Sight sensors (required):

                       eui-fsi:      Extreme Ultraviolet Imager (EUI) FSI
                       eui-hri-euv:  Extreme Ultraviolet Imager (EUI) EUV
                       eui-hri-lya:  Extreme Ultraviolet Imager (EUI) LYA
                       metis-euv:    Multi Element Telescope for Imaging and
                                     Spectroscopy (Metis)
                       metis-m0-tel: Multi Element Telescope for Imaging and
                                     Spectroscopy (Metis) IEO-M0 Boom
                       metis-vis:    Multi Element Telescope for Imaging and
                                     Spectroscopy (Metis) VIS
                       phi-fdt:      Polarimetric and Helioseismic Imager
                                     (PHI) FDT
                       phi-hrt:      Polarimetric and Helioseismic Imager
                                     (PHI) HRT
                       solohi:       SOLO Heliospheric Imager (SOLOHI)
                       spice-lw:     Spectral Imaging of the Coronal
                                     Environment (SPICE) LW
                       spice-sw:     Spectral Imaging of the Coronal
                                     Environment (SPICE) LW
                       stix:         Spectrometer Telescope for Imaging X rays
                                     (STIX)

          YYYYMMDD  coverage start and stop times in TDB (required);

          NN        version number, starting from 01 (required; e.g. 01)
                    from the source VOEM.


   If multiple versions of a C-Kernel file are provided, always use the latest
   version (unless earlier version is needed for some special reasons.)


Current CK Kernels Set
--------------------------------------------------------

    The current kernels correspond mostly to test data.


Other directory contents
--------------------------------------------------------

     aareadme.txt         This file.


Particulars
--------------------------------------------------------

     Nothing to report.


Kernel File Details
--------------------------------------------------------

     The most detailed description of the data in a binary CK file is
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

     1. ``Frames Required Reading'', NAIF Document

     2. ``C-Kernel Required Reading'', NAIF Document

     3. ``SCLK Required Reading'', NAIF Document

     4. ``TIME Required Reading'', NAIF Document

     5. ``SOC-Provided Ancillary Data for Solar Orbiter'',
        SOL-SGS-TN-0017, A. Walsh, Issue 0, Revision 3, 22/09/201


End of aareadme file.