# TriTime2d Integration into Jupyter Notebooks

Here are a few notebooks with TriTime2d integrated into Python workflows. The examples 
are based on the Murphy and Herrero Seismica submission. Before running the demos, the 
Python wrapper needs to be compiled first — instructions are provided below.

The `Demo` notebooks provide examples on how to construct meshes and calculate 
traveltimes using _TriTime2d_ for a variety of different settings. Once `_tritime2d.so` 
is compiled, these notebooks will run. The `Comparison_TwoLayer` notebook provides 
comparisons with other traveltime solvers, namely Mark Noble's Eik2d solver and Podvin 
and Lecomte's traveltime solver. As a result, additional libraries are required:

- Eik2d and instructions on how to call it from Python can be found [here](https://github.com/Mark-Noble/FTeik-Eikonal-Solver).
- Podvin and Lecomte's C subroutine is also required — [link or instructions needed]

## How to Compile the Wrapper

To compile, run the makefile in this directory:
```make -f makefile_tritime2d```

This assumes that the Fortran files are located in the same directory structure as in 
this repository and that the compiler is `gfortran`. It is also possible to specify the 
Intel compiler `ifx`, as well as the location of the Fortran files if you choose to 
store them elsewhere:
```make -f makefile_tritime2d COMPILER=ifx SRC=../..```

On successful compilation, a `_tritime2d.so` file will be generated in the current 
directory. This is read by `tritime2d.py`, which is in turn imported into the notebooks using `import tritime2d as tt`.


# References
Podvin, P. and Lecomte, I., (1991). Finite difference computation of traveltimes in very contrasted velocity models: a massively parallel approach and its associated tools, Geophysical Journal International,105(1), 271–284

Noble, M., Gesret A. and Belayouni N., (2014). Accurate 3-D finite difference computation of traveltimes in strongly heterogeneous media, Geophys.J.Int.,199,(3),1572-158.
