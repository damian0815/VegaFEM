VegaFEM
=======

Unofficial mirror of Vega-FEM (http://run.usc.edu/vega/) with a fix or two.

* fixes a compile error with friend definition in mat3d
* fixes some issues using SPOOLES solver (see also https://github.com/damiannz/spooles)

On the branch `spooles`, the SPOOLES solver is used on OSX, which is faster (it's expected to be in `../spooles.2.2` relative to the VegaFEM root, see `Makefile-headers/Makefile-header.osx`).

INSTALL.txt
-----------

Installation and compilation instructions are available in the PDF User's Manual, 
located in the "documentation" subfolder.

Quick start, main functionality (Linux, Mac OS X):
==================================================

1. Go to the subfolder "Makefile-headers", and launch either "./selectLinux" or "./selectMacOSX".
2. Go to the subfolder "utilities/interactiveDeformableSimulator" and type "make".
3. Go to the subfolder "examples/turtle" and launch "./run_turtle". LMB = apply force on model, RMB = camera movement. You can also try "examples/asianDragon" and other examples of course.

You can perform steps 1., 2., by running "./build".

Quick start, model reduction (Linux, Mac OS X):
===============================================

1. Build main functionality (above).
2. Edit "Makefile-headers/Makefile-header": fill in the sections "needed for model reduction", "needed for reducedDynamicSolver-rt", and optionally, "needed for LargeModalDeformationFactory"
3. Launch "./buildModelReduction"
4. Go to the subfolder "utilities/reducedDynamicSolver-rt" and launch "./run_A". LMB = apply force on model, RMB = camera movement

--------
Jernej Barbic, USC, 2013


