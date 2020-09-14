3D Filaments
-------------

MATLAB scripts and functions for the simulation of flexible, inextensible
untwistable filaments in three dimensions.

Installation and use
-------------

Clone or download the contents of the repository, and run 

main

in MATLAB from the directory. For best performance, first compile
dz_free_space.m into a MEX file using MATLAB Coder, and edit main.m to call
this MEX file on line 110-111. This can be achieved by running the included
script

codegen_dz_free_space

from within MATLAB, providing a supported compiler has been set up (see 'help
codegen' in MATLAB). Sample MEX files for Mac OS are provided, but are not
guaranteed to work on your system.

Extensions
-------------
To include the effects of internally generated moments, specify them in
integrated_internal_moments.m (integrated from s to 1) following the example
sinusoidal moments. Similarly, specify non-zero intrinsic curvatures in
intrinsic_curvature.m

License
-------------

This work is available under a CC BY 4.0 Attribution license. The authors,
Benjamin J. Walker, Kenta Ishimoto, and Eamonn A. Gaffney, along with the
published article "A new basis for filament simulation in three dimensions" by
the same authors should be appropriately cited in all use.