# Two-layer example

This example is a **didactical introduction** to the Fortran workflow of Trilat-time. It shows how to compute first-arrival traveltimes on a simple two-layer model, export the result and visualize the wavefronts.

![Two-layer model and resulting wavefronts](two_layers.png)

The figure shows a source located in the upper layer, a horizontal interface separating two constant-velocity media, and the resulting first-arrival field. In this example the upper layer is slower than the lower layer, so the solution contains both direct arrivals and refracted arrivals across the interface. The model is loaded from an ASCII VTK file containing the mesh, the cellwise velocity, and a theoretical traveltime field used for comparison.

---

## Purpose of this example

This example is meant to illustrate the **minimal Fortran workflow** around the core solver:

1. load a mesh and the associated velocity model,
2. define the source,
3. initialize the solver,
4. compute the traveltime field,
5. compare with a reference solution,
6. write the results to a VTK file.

The mesh-loading routine is only a support utility here. The important part is the sequence of calls around the solver.

---

## Files

* `two_layers.f90` : main Fortran program
* `Makefile` : build instructions for the example
* `input_velocity.vtk` : input mesh, velocity, and theoretical traveltime field
* `two_layers.png` : illustration of the model and the computed wavefronts
* `result.vtk` : output written by the example after execution

---

## How to run

Compile the example with:

```bash
make
```

Then run:

```bash
./two_layers
```

The program writes the numerical traveltime, the theoretical traveltime, and the relative error into `result.vtk`.

---

## Step-by-step explanation

### 1. Build the physical setup

The program starts by activating the fast mode for diffraction:

```fortran
adiff%fast = .true.
```

In this example diffraction is not used, so the diffraction structure only serves to tell the solver to perform a single fast pass. In the solver, `adiff%fast = .true.` exits after the first traveltime computation rather than performing secondary diffraction loops.

### 2. Load the model

The call

```fortran
call load_model(amesh, velocity, theo_time)
```

reads the input VTK file and fills three objects:

* `amesh` : the triangular mesh,
* `velocity` : the velocity value for each cell,
* `theo_time` : the theoretical first-arrival time at each node.

For this example, the mesh-loading details are not the main point. What matters is that after this call the program has a complete mesh and a cellwise velocity model ready for the solver.

### 3. Choose the source node

The source is prescribed by its node index:

```fortran
k = 7
```

The comments in the program indicate that this corresponds to a source at `(5000,1500)` in the upper layer. Since the code is written in Fortran, indexing is **1-based**.

### 4. Initialize the traveltime problem

The example then calls:

```fortran
call pre_timeonevsall2d_onvertex(amesh, k, traveltime, nton)
```

This helper routine performs two key operations:

* it allocates and initializes the traveltime array,
* it constructs the node-to-node connectivity structure `nton`.

Inside `pre_timeonevsall2d_onvertex`, the traveltime array is initialized to `infinity` everywhere and set to zero at the source node `k`. The same routine also allocates `nton` and computes it with `compnton`. This is exactly the standard initialization path for a source located on a mesh vertex. The helper is defined in `time.f90`, where it allocates `time` and `nton`, calls `compnton`, initializes `time` to `infinity`, and sets `time(k)=0`.

### 5. Compute the traveltime field

The core computation is performed by:

```fortran
call timeonevsall2d(amesh, velocity, traveltime, nton, adiff)
```

This is the main solver entry point. It takes the mesh, the cellwise velocity, the initialized nodal traveltime array, the precomputed node-to-node connectivity, and the diffraction control structure. In the source code, `timeonevsall2d` is documented as the routine that computes the first-arrival field from all nodes already carrying a finite initial time.

Internally, the solver propagates time through the triangular mesh using several local operators and updates the nodal first-arrival field until no better estimate is found.

### 6. Free the connectivity structure

After the solver call, it is safer to release the `nton` structure explicitly:

```fortran
call free_nton(nton)
```

This step is important because `nton` is built from pointer-allocated linked lists. These internal list nodes are **not** released automatically when the array goes out of scope. The dedicated cleanup routine `free_nton` is provided in `LAT_mesh.f90` exactly for this reason.

For clarity, the central part of the example should therefore be:

```fortran
call cpu_time(start_time)
call pre_timeonevsall2d_onvertex(amesh, k, traveltime, nton)
call timeonevsall2d(amesh, velocity, traveltime, nton, adiff)
call free_nton(nton)
call cpu_time(finish_time)
```

### 7. Export and visualize the result

After the traveltime field is computed, the program writes the results to a VTK file (`result.vtk`). The output typically includes:

- mesh geometry,
- cellwise velocity,
- numerical traveltime.

You can inspect the result using **ParaView**:

1. Open `result.vtk` in ParaView
2. Apply a color map to the `time` field
3. Add contour filters to visualize wavefronts
4. Optionally overlay the mesh to inspect resolution and interface geometry

This is the recommended way to explore and validate the computed traveltime field.
---

## What this example teaches

This example is the reference starting point if you want to understand the Fortran usage of Trilat-time.

It shows the essential solver pattern:

```fortran
call load_model(...)
call pre_timeonevsall2d_onvertex(...)
call timeonevsall2d(...)
call free_nton(...)
```

Once this sequence is clear, more advanced examples mainly differ by:

* how the mesh is built,
* how the initial condition is prescribed,
* whether diffraction is used,
* how the outputs are post-processed.

---

## Notes

* This example uses a **source on a mesh node**.
* It uses **fast mode** (`adiff%fast = .true.`), so no secondary diffraction loop is performed.
* The theoretical reference field is read from the input VTK file and used only for validation.
* The mesh-loading routine is intentionally not the focus here; it is just a convenient way to provide a ready-to-run didactical case.
* The mesh has been realized with `create_two_layers.py` which uses the gmsh api for python. 
