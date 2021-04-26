<div style="width:800px">

<img src="SPGLOGO-LR.png" align="right" width=100px />

## `SPGPylibs`
--------------------------

This repository contains libraries and functions developed by the Solar Physics Group of the [Instituto de Astrofísica de Andalucía](https://www.iaa.es). Most of them are for the analysis and pre-procesing of data coming from the Polarimetric and Helioseismic Imager of the Solar Orbiter ESA mission and from the Tunable Magnetograph (TuMAG) of the Sunrise balloom borne mission. This repository also contains tools for fully modeling the behabiour of Fabry-Pérot etalons and for Phase-Diversity  image reconstruction. More details can be found in the readme file located in the different package folder.

For more information about the members of the group visit <http://spg.iaa.es/>

Inversions tools of the SPG can be found in the main web and in the following repositories:
- <https://github.com/IAA-InvCodes/P-MILOS>
- <Desire>
- <SIR>

</div>

Installation
------------

SPGPylibs is based in Python 3+ and can be installed in your distribution using pip (**not active yet**) or directly downloading it from here.

```shell
pip install SPGPylibs
```

Getting started
===============

```python
import SPGPylibs as spg
```

Available packages
------------

- `FPtools`             <span style="float:right; width:45em;">Fabry-Pérot modeling tools</span> 
- `GENtools`             <span style="float:right; width:45em;">Generic tools</span> 
- `PHItools`             <span style="float:right; width:45em;">Solar Orbiter PHI related software</span> 
- `PDtools`             <span style="float:right; width:45em;">Phase-Diversity reconstruction tools</span> 
- `TuMAGtools`             <span style="float:right; width:45em;">Sunrise TuMAG related software</span> 

Contributors
------------

	-D. Orozco Suárez
	-P. Santamarina
	-J. Blanco
    -A. Sui
    -A. Dorantes
    -H. Strecker
    -A. Moreno Vacas

----
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)