SOLAR-ORBITER SPICE Kernel Dataset
==================================

SPICE is an essential tool for scientists and engineers alike in the
planetary science field for Solar System Geometry. Please visit the
[ESA SPICE Service (ESS) web site](https://www.cosmos.esa.int/web/spice/).



The Solar Orbiter SPICE Kernel Dataset (SKD) consists on a complete set of SPICE
Kernels that cover the whole mission lifespan including reconstructed and long
term predicted trajectory and orientation. The SKD has the following directory
structure:

* root directory (named after the mission: SOLAR-ORBITER)
  * kernels
    * ck: Kernels that contain orientation for the spacecraft and some of its structures, (solar arrays, for instance).
    * fk: Kernels that define reference frames needed for the Mission.
    * ik: Kernels for the instruments on board the spacecraft.
    * lsk: Leapseconds kernel.
    * mk: Meta-kernel files (a.k.a "furnsh" files) that provide lists of kernels suitable for a given mission period.
    * pck: Kernels that define planetary constants.
    * sclk: Spacecraft clock coefficients kernels.
    * spk: Orbit kernels, for the spacecraft and other solar system bodies.
  *	misc
    * release notes: Includes the Release Notes of the SKD.
    * other: different directories of interest.

The SKD also include text files that help to understand the contents of every directory. Those files are named aareadme.txt.


Installation
------------

In order to use Git to obtain the operational subset of the SKD, the user needs
to have Git installed. After installint Git clone the repository with

    git clone https://repos.cosmos.esa.int/socci/projects/SPICE_KERNELS/repos/SOLAR-ORBITER

In order to run the SKD in SPICE outside of the mk directory of the Git
repository the user needs to modify the following the PATH_VALUE variable
of the meta-kernel:

    PATH_VALUES       = ( '..' )

It is recommended for users to make a local copy of this file and
modify the value of the PATH_VALUES keyword to point to the actual
location of the Solar Orbiter SPICE data set's 'data' directory on
their system. Replacing '/' with '\' and converting line
   terminators to the format native to the user's system may also be
   required if this meta-kernel is to be used on a non-UNIX workstation.

Finally please note that in order to use the SPICE kernels you need to install
your prefered SPICE Library. You can find instructions on how to install
SPICE [here](https://naif.jpl.nasa.gov/naif/toolkit.html)



Function and Purpose
--------------------

ESA has a number of science missions under development and in operation that are dedicated to the study of our Solar System (i.e. MEX, Rosetta, ExoMars, BepiColombo, Solar Orbiter and JUICE). The Science Operations Centres for these missions, located at the European Space Astronomy Centre (ESAC) in Spain, are responsible for all science operations planning, data processing and archiving tasks, being the essential interface between the science instruments and the spacecraft, and with the scientific community.

From the concept study phase to the day-to-day science operations, these missions produce and use auxiliary data (spacecraft orbital state information, attitude, event information and relevant spacecraft housekeeping data) to assist science planning, data processing, analysis and archiving.

Within the Solar System scientific community, the SPICE information system is the 'de facto' standard for auxiliary data handling and geometry computations, and has been adopted and is extensively used in ESA missions. SPICE is developed and maintained by the Navigation and Ancillary Information Facility (NAIF) of the Jet Propulsion Laboratory (JPL).

The SPICE format includes definitions for orbit, attitude and event information, and the data files describing these categories are called ‘SPICE kernels’. These kernels contain sufficient information to compute derived geometrical or positional values using the existing NAIF SPICE toolkit.

For ESA missions, there are three different sources of SPICE kernels:
* The spacecraft orbit (SPK), attitude (CK) and clock correlation (SCLK) kernels, produced regularly for each mission from Mission Analysis and Flight Dynamics products and spacecraft telemetry.
* Generic kernels, such as planetary constants (PCK), leap seconds (LSK) and ephemerides (SPK), mostly provided by NAIF
* Instrument kernels (IK) and the frames specification kernel (FK), created in collaboration with the Instruments Teams

The Solar Orbiter SKD is pre-operational and therefore it is updated with the natural development of instrument and frames kernels and the updated trajectories provided by mission analysis.


The SKD contain two different types of kernels:
* Setup kernels (STK) FK, IK, PCK, LSK: they are typically text files and they are developed by ESS and are reviewed and iterated with the SGS and with the Instrument Teams when need be during the whole duration of the mission. The STKs include the following information:
  *	Set of Reference Frames of interest for geometry computations.
  *	FoV and boresight modeling for science payload.
  *	Study trajectory default orientation for S/C.
  *	Physical models for natural bodies of the mission.
* Time-varying kernels (TVK) SPK, CK, SCLK, MK: are either text or binary files and are generated with the Auxiliary Data Conversion System (ADCSng) and the source data is provided by the Flight Dynamics in terms of OEMs, AEMs and Housekeeping TM data. The TVKs include the following information:
  *	Predicted attitude and predicted/reconstructed trajectory.
  *	OBT to UTC/CAL time conversion.
  *	Reconstructed trajectory and measured orientation for S/C.
  *	Orientation of Solar Arrays and HGA (if applicable).
  *	Position of scans or turn-tables or articulations of payload.


Meta-kernel management
----------------------

Meta-kernel (MK) file management is a critical element of any SPICE dataset.
For ESA Planetary Missions the MKs of the operational SKDs
are also used to control the version of the SPICE dataset in the
ESA SPICE FTP server. In general missions in operations will have an
"operational" MK. This MK will contain the latest SPICE
kernels available for the mission including those kernels generated
by ADCSng in a periodical basis. Since the periodicity of this
updates ranges from hours to weeks, the MK in the Meta-Kernel
directory does not have a version number and is regularly
overwritten. Nevertheless in order to keep traceability and
backwards compatibility with former MKs, a duplicate of every
generated operational MK is stored. This duplicated MK filename is
appended with a time tag and a version that corresponds to the
extended version of the SKD.


Environmental Considerations
----------------------------

 At least a basic knowledge of the SPICE system is needed in order to
 use these kernels. The SPICE toolkit provides versions in Fortran
 (SPICELIB), C (CSPICE), IDL (icy), Matlab (mice) and Python (SpiceyPy),
 and the user can choose any one that suits him/her.

 The SPICELIB routine FURNSH, CSPICE function furnsh_c, icy routine
 cspice_furnsh and SpiceyPy routine furnsh load a kernel file into the
 kernel pool as shown below.

    CALL FURNSH  ( 'kernel_name' )
    furnsh_c     ( "kernel_name" );
    cspice_furnsh, 'kernel_name'
    cspice_funsh( 'kernel_name' );
    spiceypy.furnsh( 'kernel_name' )


Documentation
-------------

The HTML documentation for the Solar Orbiter SKD is available [here](ftp://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/misc/SOLAR-ORBITER.html).


How to Help
-----------

Please note that the Solar Orbiter SKD might contain bugs and is being constantly updated.
Feedback is always welcome. If you discover a bug please do not hesitate to contact the [ESA SPICE Service](spice@sciops.esa.int).
