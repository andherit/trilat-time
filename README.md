# Trilat-time

## First-arrival traveltime computation on unstructured triangular meshes

### Why?

Rapidly estimating traveltimes of seismic or tsunami waves in heterogeneous media is a central problem in geophysics. Trilat-time provides a hybrid method based on triangular meshes that is particularly well-suited for complex geometries (curved interfaces, irregular boundaries).

The method is described in:

**Murphy & Herrero (2026), submitted to Seismica**

---

## Overview

Trilat-time computes **first-arrival traveltimes** on 2D unstructured triangular meshes using a trilateration-based propagation scheme.

The project is centered on a **Fortran core solver**, with optional Python tools for preprocessing and visualization.

Typical workflows:

1. **Pure Fortran** (recommended for performance)
2. **Hybrid** (Python for mesh generation, Fortran for computation)
3. **Python wrapper** (ease of use, prototyping)

---

## Core solver (Fortran)

### Entry point

```
call timeonevsall2d(amesh, velocity, time, nton, adiff)
```

### Arguments

* `amesh` : mesh structure (nodes + triangular connectivity)
* `velocity` : velocity defined per cell
* `time` : nodal traveltime array (input/output)
* `nton` : node-to-node connectivity structure (must be precomputed)
* `adiff` : diffraction structure

### Required preprocessing

```
call compnton(amesh, nton)
```

### Memory management

```
call free_nton(nton)
```

---

## Python interface (optional)

A Python wrapper is provided for convenience and rapid testing.

```
traveltime = tritime2d(x, y, z, cells, cell_vel, initial_time, diff_nodes=None, fast=False)
```

### Notes

* The wrapper reconstructs internal Fortran structures at each call
* This introduces overhead compared to pure Fortran usage
* It is therefore best suited for:

  * prototyping
  * testing
  * mesh preprocessing pipelines

### Diffraction mode

* `fast = True`  → fast mode (reduced diffraction modeling)
* `fast = False` → full mode (more accurate, slower)

---

## Mesh and velocity

Meshes can be:

* generated externally (e.g. Gmsh)
* or constructed directly in Fortran (examples provided)

Velocity is defined **per cell**.

---

## Examples

See `examples/`:

* Two-layer model
* Diffraction example
* Ramp geometry
* Velocity gradient

These illustrate both Fortran and Python workflows.

---

## Build

### Fortran

Each example directory provides a `Makefile`.

### Python wrapper

```
cd examples/python
make
```

---

## Notes for advanced users

* The solver computes **first-arrival times only**
* Performance is optimal when used directly in Fortran
* Reusing precomputed structures (e.g. `nton`) can significantly reduce runtime in repeated simulations

---

## Summary

Trilat-time provides:

* A robust solver for complex geometries
* A high-performance Fortran implementation
* Flexible workflows depending on user needs

---

## License / Citation

(To be completed upon publication.)

