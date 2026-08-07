## PIGW 
This is the distribution of the PIGW code.
This code originated from a module in KSSOLV (see 10.1021/acs.jpca.1c03762) for reference, and is now well-separated from KSSOLV.

## Installation
Quick installation instructions for the impatient:
`./configure [options]`
` make all`
"make" alone prints a list of acceptable targets. Binaries go in bin/.

## Want to know more?
For more information, see the Yambo [main web-site](https://www.yambo-code.eu/)
Yambo is also a flagship code of the MaX Centre of Excellence [MaX web-site](https://www.max-centre.eu)

For specific documentation visit the [educational web-site](https://www.yambo-code.eu/wiki/) and related subsections
* [Getting started](https://www.yambo-code.eu/wiki/index.php?title=Tutorials)
* [Download](https://www.yambo-code.eu/wiki/index.php?title=Download)
* [Install](https://www.yambo-code.eu/wiki/index.php?title=Installation)

For support please refer to the Yambo [forum web-site](https://www.yambo-code.eu/forum) 

## License
All the material included in this distribution is free software; you can redistribute it and/or modify it under the terms of the BSD 3-Clause License.

These programs are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the BSD 3-Clause License for more details (see LICENSE file).

## AUTHORS
Please refer to the AUTHORS file

## ACKNOWLEDGEMENTS

In all source files the developers are included with their initials.
<!-- Zhengbang Zhou is the main developer currently responsible for the code development. -->

For acknowledging this work please refer to the following article:

- Z. Zhou, H. Ma, W. Wu, W. Gao, J. Yang, M. Shao, and W. Hu,
  "A fast low-rank inversion algorithm of dielectric matrix in GW approximation",
  arXiv:2403.12340 (2024).
  https://arxiv.org/abs/2403.12340

For more info please refer to the AUTHORS file


## Known Issues
Please refer to the ISSUES file

## Code structure
Yambo is composed of the following components:

* driver
* services
* apps
* controllers
* plugins

## Services

In the following the different libraries are grouped on the basis of the dependence level. Each group
depends only on lower level groups.

Each group is listed as folder and strings that label routines and modules. All paths refer to services.

### Level 0
* core: <none> 
* stderr: STDERR 
* strings: STRINGS 
* cloud: CLOUD 
* shell_operations: SH
* units: UNITS 
* numerics: NUM 
* vectors_and_matrices: V, M
* parallel/core: PARALLEL, PAR

### Level 1
* parsing_and_init: PARSER, it

### Level 2
* openmp: OPENMP
* bosons: bosons 
* timing: TIMING
* communication: COM, MSG

### Level 3
* memory: memory, MEM 
* debug: debug 
* parallel: PARALLEL (no modules, just operations)
* output: OF 
* GPUs: GPU
* io: IO_srv,DESC
* linear_algebra: LA, SLK, PAR_MATRIX, linear_algebra
* lattices: LAT, KPT, DL, RL
* service/FFT: FFT 
* service/xc_functions: XC
* service/pseudo_potentials: PP, pseudo
* service/wave_functions: WF
* service/electrons: EL, electrons [**io to complete**] 
* service/Kleinman-Bylander: PP_KB 
* service/interpolate: INTERPOLATION
* service/observables: OBSERVABLE

### Level 4
* io_more: **temporary, to be removed**
* controllers: SERVICES

## PACKAGEs

## APPs

All paths below refer to apps folder

### Level 0 (in services)
* frequencies: FREQUENCIES
* cloud: CLOUD_apps
* RIM: RIM [**io_RIM and io_RIM_W** to check]
* dipoles: DIPOLES [**io folder to rename**]

### Level \infty (engines) ###
* acfdt: ACFDT 

## Controllers
They are the closing components. At the end of the compilation they connect all services/apps/packages.

## Tested compilers 
* gcc/13.2.0/openmpi-5.0.5
* gcc/13.2.0/mpich-4.2.2

