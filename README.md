<h1 align=center>First arrival time through trilaterations</h1>

# Why ?

Rapidly estimating how long seismic or tsunami waves take to travel through complex environments
is important for many problems in the geosciences. We present a new hybrid method that uses triangular
meshes to calculate wave travel times for a given velocity model. The triangular approach performs
better when the boundaries of the model or the interfaces between materials are irregular or curved. 

From Murphy and Herrero (2026) - Submitted to Seismica.

# What is it exactly ?

The main kernel is a set of Fortran modules. With different examples, we
show how to use the code in different context, for instance using only
fortran or only python using the associated wrapper or also with a
mixed approach pre and post processing python for the meshes and 
a fortran approach for the main computation

Keep in mind that by construction, today the code is slightly faster if it
is used outside the wrapper.

# the main module

The main module is `time.f90` and the main routine is called `timeonevsall`

In the directory `examples`, we present different models which allows you to understand how
to call and initialize the variables. in each case, a `Makefile` is present
to show how to compile the code.

On the fortran side, the call is

```
call timeonevsall2d(amesh,velocity,time,nton,adiff)
```
where

`amesh` is the mesh structure (in)

`velocity` is the velocity field defined on each cell of the mesh (in)

`time` is the first arrival time defined on each node of the mesh (inout)

`nton` is the edge array and their attributes to be computed before (in)

`adiff` is the diffraction structure defining is the diffraction mode is on and where are the diffractors.

On the python side, with the wrapper, we trasform the use in a function
```
traveltime = tritime2d(x, y, z, cells, cell_vel, initial_travel, diff_nodes, fast=True)
```
Because we cannot use a fortran structure like `amesh` in the python code, the
wrapper splits the structures for the call where

`x, y, z` are the positions of the nodes in an array format

`cells` are the definition of the cells through their node indexes

`cell_vel` is the velocity array

`initial_travel` is the initial time array

`diff_nodes` is the array of the diffraction nodes 

`fast=True`indicates if the code take into account the diffraction
