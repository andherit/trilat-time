# TriTime2d Integration into Jupyter Notebooks

Here are a few notebooks with TriTime2d integrated into Python workflows. The examples 
are based on the Murphy and Herrero Seismica submission. Before running the demos, the 
Python wrapper needs to be compiled first — instructions are provided below.

The `Diffraction`,`Ramp`,`Velocity_Gradient` notebooks provide examples on how to construct meshes and calculate 
traveltimes using _TriTime2d_ for a variety of different settings. A `requirements.txt` file is provided containing a list of python libraries. To install run:
```pip install -r requirements.txt  ```
Once `_tritime2d.so` is compiled and python libraries installed, these notebooks will run. The `TwoLayer` notebook provides comparisons with other traveltime solvers, namely Mark Noble's Eik2d solver and Podvin and Lecomte's traveltime solver. As a result, additional libraries are required:

- Eik2d and instructions on how to call it from Python can be found [here](https://github.com/Mark-Noble/FTeik-Eikonal-Solver).
- Podvin and Lecomte's C subroutine is used, instructions on how to install it are provied below. 

## How to Compile the Wrapper

### Prerequisites per platform

The compile command below (`make -f makefile_tritime2d`) is identical on every platform — the 
makefile fixes the output filename to `_tritime2d.so`, and Python's `ctypes.CDLL` loads a shared 
library by that exact path regardless of its internal format or extension. What differs between 
platforms is simply whether `make` and a Fortran compiler are available beforehand:

- **Linux**: `make` and `gfortran` are typically available already, or installable via your 
  package manager (e.g. `apt install make gfortran`).
- **macOS**: install Xcode Command Line Tools (`xcode-select --install`) for `make`, and 
  `gfortran` via Homebrew (`brew install gcc`, which includes `gfortran`).
- **Windows**: `make` and `gfortran` are not installed by default. Install them through 
  [MSYS2](https://www.msys2.org/) or [WSL](https://learn.microsoft.com/en-us/windows/wsl/) (which 
  gives you a Linux environment), then run the same command from that shell.

To compile, run the makefile in this directory:
```make -f makefile_tritime2d```

This assumes that the Fortran files are located in the same directory structure as in 
this repository and that the compiler is `gfortran`. It is also possible to specify the 
Intel compiler `ifx`, as well as the location of the Fortran files if you choose to 
store them elsewhere:
```make -f makefile_tritime2d COMPILER=ifx SRC=../..```

On successful compilation, a `_tritime2d.so` file will be generated in the current 
directory. This is read by `tritime2d.py`, which is in turn imported into the notebooks using `import tritime2d as tt`.

## Installing Podvin and Lecomte's subroutine

The `PL` directory contains `Time_2d.c` and `time_2d.h`, Podvin and Lecomte's finite-difference 
traveltime solver. Unlike `tritime2d`, this is called directly from `Comparison_TwoLayer.ipynb` as 
a C shared library via Python's `ctypes` module, so it must be compiled separately (not through 
`makefile_tritime2d`).

From the `PL` directory, compile with `gcc`:
```
cd PL
gcc -shared -fPIC -o time_2d.so Time_2d.c
```

On Windows (e.g. with MinGW), build a `.dll` instead:
```
gcc -shared -o time_2d.dll Time_2d.c
```

Once compiled, set `so_file` in `Comparison_TwoLayer.ipynb` to the path of your compiled library, 
e.g.:
```python
so_file = r"PL/time_2d.so"
time_so = ctypes.CDLL(so_file)
time_2d = time_so.time_2d
```

# References
Podvin, P. and Lecomte, I., (1991). Finite difference computation of traveltimes in very contrasted velocity models: a massively parallel approach and its associated tools, Geophysical Journal International,105(1), 271–284

Noble, M., Gesret A. and Belayouni N., (2014). Accurate 3-D finite difference computation of traveltimes in strongly heterogeneous media, Geophys.J.Int.,199,(3),1572-158.
