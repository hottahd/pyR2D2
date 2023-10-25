# R2D2_py
Python code for analyzing  the result of R2D2

# Upgrade resolution

if you want to use fortran utility such as `d_x`, `d_y`, `d_z`, `interp`, and `spherical2cartesian`, compile fortran programs at `R2D2_py/R2D2/fortran_src`

```
make
```

We need `gfortran` compiler but you can edit `R2D2_py/R2D2/fortran_src/Makefile` so that the other compilers are used.