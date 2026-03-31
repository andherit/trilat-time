# Examples

This directory contains usage examples of **Trilat-time**, illustrating different workflows and use cases.

---

## Structure

```
examples/
  ├── fortran/
  └── python/
```

* `fortran/` : standalone Fortran examples (recommended for performance and production use)
* `python/`  : Jupyter notebooks using the Python wrapper (for prototyping and visualization)

---

## Fortran examples

Each example is organized in its own directory and provides:

* a complete workflow (mesh → traveltime computation)
* a `Makefile` for compilation
* input/output files when required
These examples demonstrate how to use the solver directly in Fortran and are the reference for performance-oriented usage.

---

## Python examples

Python examples are provided as Jupyter notebooks and illustrate:

* mesh handling (e.g. with Gmsh)
* velocity model construction
* traveltime computation via the wrapper
* visualization

Before running the notebooks, the Python wrapper must be compiled.

## Notes

* Fortran examples provide the most direct and efficient usage of the solver
* Python examples are intended for convenience, testing, and visualization
* Both approaches rely on the same core solver implementation

