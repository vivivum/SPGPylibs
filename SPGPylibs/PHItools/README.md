<div style="width:800px">

<img src="../SPGLOGO-LR.png" align="right" width=100px />

## `SPGPylibs/PHItools`
--------------------------

These PHItools functions have been developed during the Solar Orbiter PHI commissioning and cruise phase activities. They provide tools for performing data analysis and calibration tasks. These libraries are not optimized for performance or anything. The user of these libraries should take care of their proper use. 

Functions present in SPGPylibs.PHItools are listed below.

-------------------------- 
</div>

#### *phifdt_flat.py*      
- `do_hough()`             <span style="float:right; width:45em;">n_images flat gen test run</span> 
- `fdt_flat()`             <span style="float:right; width:45em;">KLL</span>      
- `fdt_flat_gen()`         <span style="float:right; width:45em;">FDT kll flatfield determination wrapper.</span> 
- `fdt_flat_testrun()`     <span style="float:right; width:45em;">FDT flat gen test run</span>

#### *phi_utils.py*      
- `dot()`             <span style="float:right; width:45em;">dot product</span> 
- `angle()`             <span style="float:right; width:45em;">angle between two vectors</span>      
- `sphere2cart()`         <span style="float:right; width:45em;">Coordinate transformation</span> 
- `cart2sphere()`     <span style="float:right; width:45em;">Coordinate transformation</span>
- `cart2polar()`      <span style="float:right; width:45em;">Coordinate transformation</span>
- `polar2cart()`       <span style="float:right; width:45em;">Coordinate transformation</span>
- `phi_orbit()`     <span style="float:right; width:45em;">Wrapper to scipypy for orbit calculations</span>
- `phi_orbit_test()`     <span style="float:right; width:45em;">tests program</span>
- `allen_clv()`     <span style="float:right; width:45em;">Allen CLV coefficients</span>
- `azimutal_average()`     <span style="float:right; width:45em;">Azimutically averages the intensity of solar disk</span>
- `limb_darkening()`     <span style="float:right; width:45em;">Limb darkening function</span>
- `newton()`     <span style="float:right; width:45em;">Newton LS</span>
- `get_time()`    <span style="float:right; width:45em;">get time (h) from datetime format</span>

#### *phi_gen.py*      
- `shift()`             <span style="float:right; width:45em;">Shift function</span> 
- `generate_circular_mask()`             <span style="float:right; width:45em;">Circular mask</span>      
- `gradient()`         <span style="float:right; width:45em;">--</span> 
- `threshold_otsu()`     <span style="float:right; width:45em;">--</span>
- `histogram()`     <span style="float:right; width:45em;">--</span>
- `FindEdges()`     <span style="float:right; width:45em;">--</span>
- `make_circles()`     <span style="float:right; width:45em;">--</span>
- `find_Circles_ida()`     <span style="float:right; width:45em;">--</span>
- `votes()`     <span style="float:right; width:45em;">tests program</span>
- `bin_annulus()`     <span style="float:right; width:45em;">tests program</span>
- `circle_grid()`     <span style="float:right; width:45em;">tests program</span>
- `find_circle_hough()`     <span style="float:right; width:45em;">tests program</span>
- `FFTs()`     <span style="float:right; width:45em;">tests program</span>
- `Laplacian()`     <span style="float:right; width:45em;">tests program</span>
- `rebin()`     <span style="float:right; width:45em;">tests program</span>

#### *phi_reg.py*      
- `sampling()`             <span style="float:right; width:45em;">--</span> 
- `aperture()`             <span style="float:right; width:45em;">--</span>      
- `PHI_shifts_FFT()`         <span style="float:right; width:45em;">--</span> 
- `PHI_shifts_CC()`         <span style="float:right; width:45em;">--</span> 
- `shift_subp()`         <span style="float:right; width:45em;">--</span> 
- `gaussian()`         <span style="float:right; width:45em;">--</span> 
- `moments()`         <span style="float:right; width:45em;">--</span> 
- `fitgaussian()`         <span style="float:right; width:45em;">--</span> 

#### *phi_fits.py*      
- `fits_get()`             <span style="float:right; width:45em;">--</span> 
- `fits_get_sampling()`             <span style="float:right; width:45em;">--</span>      
- `fits_get_part()`         <span style="float:right; width:45em;">--</span> 
- `read_shifts()`         <span style="float:right; width:45em;">--</span> 
- `write_shifts()`         <span style="float:right; width:45em;">--</span> 
- `fits_get_fpatimes()`         <span style="float:right; width:45em;">--</span>

#### *phifdt_pipe.py*      
- `----()`             <span style="float:right; width:45em;">--</span> 

#### *phihrt_pipe.py*      
- `----()`             <span style="float:right; width:45em;">--</span> 

#### *tools.py*      
- `countcall()`            <span style="float:right; width:45em;">--</span> 
- `timeit()`             <span style="float:right; width:45em;">--</span> 

