<h1 align=center>First arrival time through trilaterations</h1>

# Why ?

Rapidly estimating how long seismic or tsunami waves take to travel through complex environments
is important for many problems in the geosciences. We present a new hybrid method that uses triangular
meshes to calculate wave travel times for a given velocity model. Tthe triangular approach performs
better when the boundaries of the model or the interfaces between materials are irregular or curved. 

From Murphy and Herrero (2026) - Submitted to Seismica.

# What is it exactly ?

It is a set of Fortran modules. The main modules is `time.f90` and the main routine is called
`timeonevsall`

In the directory `examples`, we present four different models which allows you to understand how
to call and initialize the variables.

```
call timeonevsall2d(amesh,velocity,time,nton,adiff)
```
where

`amesh` is the mesh structure (in)

`velocity` is the velocity field defined on each cell of the mesh (in)

`time` is the first arrival time defined on each node of the mesh (inout)

`nton` is the edge array and their attributes to be computed before (in)

`adiff` is the diffraction structure defining is the diffraction mode is on and where are the diffractors.
