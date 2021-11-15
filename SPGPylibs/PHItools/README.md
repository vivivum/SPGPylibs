<div style="width:800px">

<img src="../SPGLOGO-LR.png" align="right" width=100px />

## `SPGPylibs/PHItools`
--------------------------

These PHItools functions have been developed during the Solar Orbiter PHI commissioning and cruise phase activities. They provide tools for performing data analysis and calibration tasks. These libraries are not optimized for performance or anything. The user of these libraries should take care of their proper use. This software will be in continuous development until end of 2022. 

Functions present in SPGPylibs.PHItools are listed below. There are also two folders containing a beta version of the c-milos code (cmilos) and the SOLO kernels (orbits-data). 

The c-milos code is not optimized for performance and should be compiled in each machine (cmilos/lib directory). It uses gcc standard compiler. There is a version for Python using Cython and can be compiled in the main cmilos folder. There is no guarantee it will compile succesfully (under development). 

The kernels should be updated to the last version by the user at <https://www.cosmos.esa.int/web/spice/solar_orbiter>

Main functions present in SPGPylibs.PHItools are listed below.

-------------------------- 
</div>

#### *phifdt_flat.py*  
- `do_hough()`             <span style="float:right; width:45em;">Find the center of full disk images using the Hough transform</span> 
- `fdt_flat_gen()`         <span style="float:right; width:45em;">FDT kll flatfield determination core code.</span> 
- `fdt_flat()`             <span style="float:right; width:45em;">Determine image displacements and call fdt_flat_gen for flat determination</span>      
- `fdt_flat_testrun()`     <span style="float:right; width:45em;">FDT flat generation test program</span>

#### *phi_utils.py*      
- `dot()`             <span style="float:right; width:45em;">dot product</span> 
- `angle()`             <span style="float:right; width:45em;">angle between two 2D-vectors</span>      
- `angle2()`             <span style="float:right; width:45em;">angle between two 3D-vectors</span>      
- `crossp()`             <span style="float:right; width:45em;">Cross product between two 3D vectors</span>      
- `sphere2cart()`         <span style="float:right; width:45em;">Coordinate transformation</span> 
- `cart2sphere()`     <span style="float:right; width:45em;">Coordinate transformation</span>
- `cart2polar()`      <span style="float:right; width:45em;">Coordinate transformation</span>
- `polar2cart()`       <span style="float:right; width:45em;">Coordinate transformation</span>
- `phi_orbit()`     <span style="float:right; width:45em;">Wrapper to scipypy for orbit calculations</span>
- `phi_orbit_test()`     <span style="float:right; width:45em;">tests program</span>
- `allen_clv()`     <span style="float:right; width:45em;">Allen CLV coefficients</span>
- `azimutal_average()`     <span style="float:right; width:45em;">Azimuthally averages the intensity of solar disk</span>
- `limb_darkening()`     <span style="float:right; width:45em;">Limb darkening function</span>
- `newton()`     <span style="float:right; width:45em;">Newton LS</span>
- `get_time()`    <span style="float:right; width:45em;">get time (h) from datetime format</span>
- `genera_2d()`    <span style="float:right; width:45em;">generates a 2D image from a 1D cut of dim/2</span>
- `running_mean()`    <span style="float:right; width:45em;">Calculated the running mean over a 1D vector</span>
- `find_string()`    <span style="float:right; width:45em;">find a string in an string</span>
- `rotate_grid()`    <span style="float:right; width:45em;">Rotate a grid or a point</span>

#### *phi_gen.py*      
- `shift()`             <span style="float:right; width:45em;">Shift 2D images at pixel level</span> 
- `generate_circular_mask()`             <span style="float:right; width:45em;">Circular mask generator</span>      
- `gradient()`         <span style="float:right; width:45em;">Calculate gradient of image using different methods</span> 
- `threshold_otsu()`     <span style="float:right; width:45em;">Return threshold value based on Otsu's method</span>
- `histogram()`     <span style="float:right; width:45em;">Return histogram of image</span>
- `FindEdges()`     <span style="float:right; width:45em;">Find edges of image using simply, prewitt or prewittsmooth methods</span>
- `make_circles()`     <span style="float:right; width:45em;">Create a circle (or mask)</span>
- `find_Circles_ida()`     <span style="float:right; width:45em;">Hough transform circle finder</span>
- `votes()`     <span style="float:right; width:45em;">Hough transform votes program</span>
- `bin_annulus()`     <span style="float:right; width:45em;">creates a anulus mask</span>
- `circle_grid()`     <span style="float:right; width:45em;">creates a grid of points with NxN dimensions</span>
- `find_circle_hough()`     <span style="float:right; width:45em;">Do Hough Transform</span>
- `simple_shift()`     <span style="float:right; width:45em;">shifts elements in a vector</span>
- `find_center()`     <span style="float:right; width:45em;">find center of Sun using vertical and horizontal cuts</span>
- `FFTs()`     <span style="float:right; width:45em;">Fourier transform wrapper</span>
- `Laplacian()`     <span style="float:right; width:45em;">calculate gradient of real image using Laplacian filter</span>
- `rebin()`     <span style="float:right; width:45em;">Rebin 2D array arr to shape new_shape by averaging</span>
- `apod()`     <span style="float:right; width:45em;"> Turkey apodization mask</span>

#### *phi_reg.py*      
- `sampling()`             <span style="float:right; width:45em;">creates a grid of points with NxN dimensions</span> 
- `aperture()`             <span style="float:right; width:45em;">calculates a simple aperture function</span>      
- `PHI_shifts_FFT()`         <span style="float:right; width:45em;">calculate relative shifts between two images</span> 
- `PHI_shifts_CC()`         <span style="float:right; width:45em;">Image cross-correlation</span> 
- `shift_subp()`         <span style="float:right; width:45em;">Shift (subpixel) an image</span> 
- `gaussian()`         <span style="float:right; width:45em;">2D gaussian function</span> 
- `moments()`         <span style="float:right; width:45em;">Moments of 2D gaussian function</span> 
- `Gauss()`         <span style="float:right; width:45em;">1D gauss function</span> 
- `Gauss2()`         <span style="float:right; width:45em;">Two 1D gauss functions</span> 
- `fitgaussian()`         <span style="float:right; width:45em;">fit 2D gaussian fucntion to image</span> 

#### *phi_fits.py*      
- `fits_get()`             <span style="float:right; width:45em;">Function to load PHI files</span> 
- `fits_get_fpatimes()`         <span style="float:right; width:45em;">Get FPA times</span>
- `list_fits()`         <span style="float:right; width:45em;">List available fits in a folder</span>
- `fits_get_sampling()`             <span style="float:right; width:45em;">get wavelength / voltage sampling in file</span>      
- `fits_get_part()`         <span style="float:right; width:45em;">Get single npol and wave from PHI data</span> 
- `read_shifts()`         <span style="float:right; width:45em;">For KLL centers</span> 
- `write_shifts()`         <span style="float:right; width:45em;">For KLL centers</span> 

#### *phifdt_pipe.py*      
- `phifdt_pipe()`             <span style="float:right; width:45em;">FDT pipeline - See file header for detailed usage</span> 

#### *phihrt_pipe.py*      
- `----()`             <span style="float:right; width:45em;">--</span> 

#### *phifdt_pipe_modules.py*      
- `phi_correct_dark()`             <span style="float:right; width:45em;">Read and correct dark field</span> 
- `interpolateImages()`             <span style="float:right; width:45em;">interpolate 2D images</span> 
- `phi_correct_prefilter()`             <span style="float:right; width:45em;">Correct prefilter</span> 
- `applyPrefilter_dos()`             <span style="float:right; width:45em;">PHI prefilter. Modified version from K. Albert.</span> 
- `phi_apply_demodulation()`             <span style="float:right; width:45em;">Demodulate data</span> 
- `crosstalk_ItoQUV()`             <span style="float:right; width:45em;">Evaluate and correct crosstalk from Stokes I to QUV</span> 
- `cross_talk_QUV()`             <span style="float:right; width:45em;">Evaluate and correct crosstalk from Stokes V to QU</span> 
- `crosstalk_ItoQUV2d()`             <span style="float:right; width:45em;">Evaluate and correct crosstalk from Stokes I to QUV in 2D</span> 
- `phi_correct_ghost()`             <span style="float:right; width:45em;">Ghost correction (Under development)</span> 
- `phi_correct_fringes()`             <span style="float:right; width:45em;">Finge correction (Under development)</span> 

#### *tools.py*      
- `printc()`            <span style="float:right; width:45em;">--</span> 
- `countcall()`            <span style="float:right; width:45em;">--</span> 
- `timeit()`            <span style="float:right; width:45em;">--</span> 
- `fix_path()`             <span style="float:right; width:45em;">--</span> 

#### *phi_rte.py*      
- `phi_rte()`            <span style="float:right; width:45em;">RTE wrapper (in development)</span> 

#### *inputs_json folder*      
- `json_gnerator()`            <span style="float:right; width:45em;">Generate a generic FDT json input file (HRT Not implemented)</span> 

#### *tests folder*      
- `Examples`            <span style="float:right; width:45em;">...</span> 

#### *polcal_py folder*      
- `polcal_lib_v2()`            <span style="float:right; width:45em;">Polarization calibration software</span> 

#### *cmilos folder*      
- `cmilos_testdata`        <span style="float:right; width:45em;">RTE test data</span> 
- `pymilos.pyx`            <span style="float:right; width:45em;">Milos Cython wrapper </span> 
