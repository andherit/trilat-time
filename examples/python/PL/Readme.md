
# Podvin and Lecomte's traveltime solver

The `Comparison_TwoLayer` notebook benchmarks TriTime2d against Podvin and Lecomte's finite-difference traveltime solver (Podvin & Lecomte, 1991; see References). This is not our code and if used should be referenced correctly - see comment at start of file. To compile into a shared library: 
```
gcc -shared -fPIC -o time_2d.so Time_2d.c

```
Then inside `Comparison_TwoLayer.ipynb` point to the location of the shared library:

```python
so_file = r"path/to/time_2d.so"
time_so = ctypes.CDLL(so_file)
time_2d = time_so.time_2d
```
(On Windows, build a `.dll` instead: `gcc -shared -o time_2d.dll Time_2d.c`.)

# References
Podvin, P. and Lecomte, I., (1991). Finite difference computation of traveltimes in very contrasted velocity models: a massively parallel approach and its associated tools, Geophysical Journal International,105(1), 271–284