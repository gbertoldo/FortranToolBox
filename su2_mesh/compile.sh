#!/bin/bash

gfortran -O3 -o su2mesh.x  mod_boundary.f90 class_ifile.f90 main.f90

rm *.mod
