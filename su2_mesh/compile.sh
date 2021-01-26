#!/bin/bash

gfortran -O3 -o su2mesh.x class_ifile.f90 mod_class_mcs.f90 mod_hermite_functions.f90 mod_vector_search.f90 mod_class_path2d.f90 mod_boundary.f90 mod_points_distribution.f90 main.f90

rm *.mod
