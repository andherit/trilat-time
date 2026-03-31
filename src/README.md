# Trilat-time – `src/` Technical Documentation

This document describes the internal architecture, data structures, and API of the Fortran core.

---

## Modules overview

Source files in `src/`:

* `generic.f90`        : numerical kinds, constants, utilities
* `lists.f90`          : linked-list containers used for sparse graph structures
* `LAT_mesh.f90`       : mesh definition and topology (node ↔ node, node ↔ cell)
* `LAT_time.f90`       : diffraction and auxiliary runtime structures
* `time.f90`           : main solver (trilateration propagation)
* `iomod.f90`          : I/O utilities (VTK export)
* `wrapper.f90`        : C/Python interface layer

---

## Numerical conventions

Defined in `generic`:

* `pr`  : working real precision (double)
* `pin` : integer precision
* `infinity` : large sentinel value used for uninitialized traveltime

All computations assume **double precision** (`real(pr)`).

---

## Core data structures

### Mesh (`LAT_mesh`)

```fortran
type mesh
  integer(pin) :: Nnodes, Ncells
  real(pr), allocatable :: px(:), py(:), pz(:)
  integer(pin), allocatable :: cell(:,:)
end type
```

* Nodes are indexed from **1 (Fortran convention)**
* `cell(i,:)` gives the 3 node indices of triangle `i`

---

### Node-to-node connectivity (`nton`)

Constructed with:

```fortran
call compnton(amesh, nton)
```

`nton` is an array of linked lists:

* For each node `i`, `nton(i)` stores all neighboring nodes
* Each edge record contains:

  * neighbor node id
  * edge length (`donedge`)
  * adjacent cells (`cellonedge(2)`)

This structure is central to the solver and acts as a **sparse graph of the mesh**.

---

### Memory management (critical)

`nton` uses pointer-based linked lists and must be freed explicitly:

```fortran
call free_nton(nton)
```

Failure to call this results in memory leaks.

---

### Diffraction structure (`LAT_time`)

```fortran
type diff
  logical :: fast
  integer(pin) :: NdiffNodes
  integer(pin), allocatable :: nodes(:)
end type
```

* `fast = .true.`  → single-pass (no secondary diffraction)
* `fast = .false.` → multi-pass with secondary sources

---

## Solver API

### Main routine

```fortran
call timeonevsall2d(amesh, velocity, time, nton, adiff)
```

### Inputs

* `amesh`     : mesh structure
* `velocity`  : array of size `Ncells`
* `nton`      : node-to-node connectivity (precomputed)
* `adiff`     : diffraction control

### In/Out

* `time` : nodal traveltime array

  * input: initialized values (sources finite, others = `infinity`)
  * output: first-arrival times

---

## Initialization helpers

### Single-node source

```fortran
call pre_timeonevsall2d_onvertex(amesh, k, time, nton)
```

* initializes `time` to `infinity`
* sets `time(k) = 0`
* builds `nton`

---

## Internal algorithm (high-level)

The solver implements a **trilateration-based fast marching scheme**:

1. Initialize active node list (`ntodo`) from finite times
2. Iteratively:

   * extract node with minimum time
   * propagate to neighbors using multiple operators:

     * edge propagation
     * planar propagation
     * curved (circular) propagation
     * head waves
3. Update neighboring nodes if a smaller time is found
4. Optionally reloop for diffraction sources

The algorithm uses a **linked-list priority approximation** instead of a heap.

---

## Operator hierarchy

At each update, candidate times are computed and ranked:

* curved (face)
* head wave
* edge
* planar

Mode selection ensures consistent physics and avoids invalid updates.

---

## Diffraction handling

If `fast = .false.`:

* additional loops are performed
* secondary sources are activated
* traveltime field is updated iteratively

If `fast = .true.`:

* only the first propagation is performed
* computation stops after first pass

---

## Wrapper interface (`wrapper.f90`)

Provides a C-compatible entry point:

```fortran
subroutine tritime2d(...) bind(C)
```

Responsibilities:

* reconstruct `mesh` from flat arrays
* build `nton`
* call solver
* free `nton`

This wrapper is intended for Python (ctypes) usage.

---

## Performance notes

* Main cost components:

  * `compnton` (topology construction)
  * iterative propagation loop
* Reusing `nton` across runs can significantly reduce runtime
* Pure Fortran usage avoids wrapper overhead

---

## Debugging

Controlled via:

```fortran
integer(pin) :: verbose
```

Levels:

* `0` : silent
* `1` : workflow
* `2` : detailed
* `3` : debug (stores intermediate fields)

---

## Summary

The `src/` directory contains a modular Fortran implementation of a
trilateration-based first-arrival solver built around:

* explicit mesh topology (`nton`)
* linked-list propagation
* multiple physical operators
* optional diffraction reprocessing

Understanding the lifecycle of `nton` and the initialization of `time`
is essential for correct and efficient usage.

