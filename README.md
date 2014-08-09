VegaFEM
=======

Unofficial mirror of Vega-FEM (http://run.usc.edu/vega/) with a fix or two.

* fixes a compile error with friend definition in mat3d
* fixes some issues using SPOOLES solver (see also https://github.com/damiannz/spooles)

On the branch `spooles`, the SPOOLES solver is used on OSX, which is faster (it's expected to be in `../spooles.2.2` relative to the VegaFEM root, see `Makefile-headers/Makefile-header.osx`).
