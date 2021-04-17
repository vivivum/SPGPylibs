<div style="width:800px">

<img src="SPGLOGO-LR.png" align="right" width=100px />

## `SPGPylibs/PHItools`
--------------------------


These PHItools functions have been developed during the Solar Orbiter PHI commissioning and cruise phase activities. They provide tools for performing data analysis and calibration tasks. These libraries are not optimized for performance or anything. The user of these libraries should take care of their proper use. 

Functions present in SPGPylibs.PHItools are listed below.

-------------------------- 
</div>

| <div style="width:220px">*phihrt_flat.py* </div> |  <div style="width:400px"><font size="2"> description</font> </div>    |
|-------------------------------|-----------------------------|
|` `                      | <font size="2">  </font> |


| <div style="width:220px">*phifdt_flat.py* </div> |  <div style="width:400px"><font size="2"> description</font> </div>    |
|-------------------------------|-----------------------------|
|`do_hough()`     | <font size="2"> n_images flat gen test run </font>               |
|`fdt_flat()`                      | <font size="2"> KLL </font> |
|`fdt_flat_gen()`     | <font size="2"> FDT kll flatfield determination wrapper. </font>               |
|`fdt_flat_testrun()`     | <font size="2"> FDT flat gen test run </font>               |


| <div style="width:220px">*phi_utils.py* </div> |  <div style="width:400px"><font size="2"> description</font> </div>    |
|-------------------------------|-----------------------------|
|`dot()`                      | <font size="2"> dot product </font> |
|`angle()`     | <font size="2"> angle between two vectors </font>               |
|`sphere2cart()`     | <font size="2"> Coordinate transformation </font>               |
|`cart2sphere()`     | <font size="2"> Coordinate transformation </font>               |
|`phi_orbit()`     | <font size="2"> Wrapper to scipypy for orbit calculations </font>               |
|`phi_orbit_test()`     | <font size="2"> tests program </font>               |


|  <div style="width:220px">*phi_gen.py*</div>  |  <div style="width:400px"><font size="2"> description</font> </div>    |
|-------------------------------|-----------------------------|
|`shift()`                      | <font size="2"> Shift function</font> |
|`generate_circular_mask()`     | <font size="2"> Circular mask</font>               |
|`gradient()`     | <font size="2">  - </font>               |
|`threshold_otsu()`     | <font size="2"> - </font>               |
|`histogram()`     | <font size="2"> Circular mask</font>               |
|`FindEdges()`     | <font size="2"> Circular mask</font>               |
|`make_circles()`     | <font size="2"> Circular mask</font>               |
|`find_Circles_ida()`     | <font size="2"> Circular mask</font>               |
|`votes()`     | <font size="2"> Circular mask</font>               |
|`bin_annulus()`     | <font size="2"> Circular mask</font>               |
|`circle_grid()`     | <font size="2"> Circular mask</font>               |
|`find_circle_hough()`     | <font size="2"> Circular mask</font>               |


|  <div style="width:220px">*phi_reg.py*</div>  |  <div style="width:400px"><font size="2"> Image registering tools</font> </div>    |
|-------------------------------|-----------------------------|
|`sampling()`                      | <font size="2">  </font> |
|`aperture()`                      | <font size="2">  </font> |
|`PHI_shifts_FFT()`                      | <font size="2">  </font> |


|  <div style="width:220px">*phi_fits.py*</div>  |  <div style="width:400px"><font size="2"> PHI I/O tools </font> </div>    |
|-------------------------------|-----------------------------|
|`fits_get()`                      | <font size="2">  </font> |
|`fits_get_sampling()`     | <font size="2">  </font>               |
|`fits_get_part()`     | <font size="2">  </font>               |
|`read_shifts()`     | <font size="2">  </font>               |
|`write_shifts()`     | <font size="2">  </font>               |


|  <div style="width:220px">*phifdt_pipe.py*</div>  |  <div style="width:400px"><font size="2"> FDT pipeline</font></div>     |
|-------------------------------|-----------------------------|
|` `                      | <font size="2">  </font> |


|  <div style="width:220px">*phihrt_pipe.py*</div>  |  <div style="width:400px"><font size="2"> HRT pipeline</font></div>     |
|-------------------------------|-----------------------------|
|` `                      | <font size="2">  </font>          |


|  <div style="width:220px"> *tools.py* </div>             |  <div style="width:400px"> <font size="2"> Tools </font> </div>     |
|-------------------------------|-----------------------------|
|`countcall()`                  | <font size="2">  </font> |
|`timeit()`                     | <font size="2">  </font>               |

