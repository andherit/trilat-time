program mesh_gradient

   use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t
   use gmsh
   use generic,  only : pin, pr
   use LAT_mesh, only : mesh
   use iomod,    only : dumpmeshvtk, dumpcellattributevtk, dumpnodeattributevtk

   implicit none

   type(gmsh_t) :: gm
   type(mesh)   :: amesh

   integer(c_int) :: p1, p2, p3, p4
   integer(c_int) :: l1, l2, l3, l4
   integer(c_int) :: cloop, surf
   integer(c_int) :: field_id

   integer(c_size_t), allocatable :: node_tags(:)
   integer(c_size_t), allocatable :: tri_tags(:)
   integer(c_size_t), allocatable :: tri_node_tags(:)

   real(c_double), allocatable :: coord(:)
   real(c_double), allocatable :: param(:)

   integer(pin), allocatable :: tag2local(:)

   real(pr), allocatable :: velocity(:)
   real(pr), allocatable :: theory(:)

   integer(pin) :: i, j
   integer(pin) :: max_node_tag, local_tag
   integer(pin) :: source_node
   integer      :: vtk

   real(pr), parameter :: xmin = 0._pr
   real(pr), parameter :: xmax = 20000._pr
   real(pr), parameter :: ymin = 0._pr
   real(pr), parameter :: ymax = 10000._pr

   real(pr), parameter :: xs = 0._pr
   real(pr), parameter :: ys = 0._pr

   real(pr), parameter :: v0   = 500._pr
   real(pr), parameter :: grad = 0.4_pr

   real(pr), parameter :: htop = 20._pr
   real(pr), parameter :: hbot = 400._pr

   real(pr) :: x, y, yc
   real(pr) :: vsrc, vrec, r2, arg

!-------------------------------------------------------------------------------
! 1. Build the rectangular geometry and generate the Gmsh mesh
!-------------------------------------------------------------------------------

   call gm%initialize()
   call gm%option%setNumber("General.Terminal", 1.0d0)

   call gm%model%add("vertical velocity gradient example")

   p1 = gm%model%geo%addPoint(real(xmin,c_double), real(ymin,c_double), 0.0d0, &
                              real(htop,c_double), 1)
   p2 = gm%model%geo%addPoint(real(xmax,c_double), real(ymin,c_double), 0.0d0, &
                              real(htop,c_double), 2)
   p3 = gm%model%geo%addPoint(real(xmax,c_double), real(ymax,c_double), 0.0d0, &
                              real(hbot,c_double), 3)
   p4 = gm%model%geo%addPoint(real(xmin,c_double), real(ymax,c_double), 0.0d0, &
                              real(hbot,c_double), 4)

   l1 = gm%model%geo%addLine(p1, p2, 1)
   l2 = gm%model%geo%addLine(p2, p3, 2)
   l3 = gm%model%geo%addLine(p3, p4, 3)
   l4 = gm%model%geo%addLine(p4, p1, 4)

   cloop = gm%model%geo%addCurveLoop([l1, l2, l3, l4], 1)
   surf  = gm%model%geo%addPlaneSurface([cloop], 1)

   call gm%model%geo%synchronize()

!  Background target size:
!
!     h(y) = 20 + 0.038 y
!
!  With y in meters. The global min/max options clip it to [20, 400].

   field_id = gm%model%mesh%field%add("MathEval", 1)
   call gm%model%mesh%field%setString(field_id, "F", "20 + 0.038*y")
   call gm%model%mesh%field%setAsBackgroundMesh(field_id)

   call gm%option%setNumber("Mesh.MeshSizeFromPoints", 0.0d0)
   call gm%option%setNumber("Mesh.MeshSizeFromCurvature", 0.0d0)
   call gm%option%setNumber("Mesh.MeshSizeExtendFromBoundary", 0.0d0)
   call gm%option%setNumber("Mesh.MeshSizeMin", real(htop,c_double))
   call gm%option%setNumber("Mesh.MeshSizeMax", real(hbot,c_double))
   call gm%option%setNumber("Mesh.ElementOrder", 1.0d0)
   call gm%option%setNumber("Mesh.Algorithm", 5.0d0)

   call gm%model%mesh%generate(2)

!-------------------------------------------------------------------------------
! 2. Extract Gmsh nodes and linear triangular elements
!-------------------------------------------------------------------------------

   call gm%model%mesh%getNodes(node_tags, coord, param, returnParametricCoord=.false.)
   call gm%model%mesh%getElementsByType(2, tri_tags, tri_node_tags)

   call gm%finalize()

   if (size(node_tags) == 0) stop "Gmsh returned no nodes."
   if (size(tri_tags) == 0)  stop "Gmsh returned no triangular cells."

   if (size(tri_node_tags) /= 3 * size(tri_tags)) then
      stop "Unexpected triangle connectivity size from Gmsh."
   endif

!-------------------------------------------------------------------------------
! 3. Fill the internal trilat_time mesh structure
!-------------------------------------------------------------------------------

   amesh%Nnodes = int(size(node_tags), pin)
   amesh%Ncells = int(size(tri_tags), pin)

   allocate(amesh%px(amesh%Nnodes))
   allocate(amesh%py(amesh%Nnodes))
   allocate(amesh%pz(amesh%Nnodes))
   allocate(amesh%cell(amesh%Ncells,3))

   max_node_tag = int(maxval(node_tags), pin)
   allocate(tag2local(max_node_tag))
   tag2local = 0

   do i = 1, amesh%Nnodes
      local_tag = int(node_tags(i), pin)
      tag2local(local_tag) = i

      amesh%px(i) = real(coord(3*i-2), pr)
      amesh%py(i) = real(coord(3*i-1), pr)
      amesh%pz(i) = real(coord(3*i  ), pr)
   enddo

   do i = 1, amesh%Ncells
      do j = 1, 3
         local_tag = int(tri_node_tags(3*(i-1)+j), pin)

         if (local_tag < 1 .or. local_tag > max_node_tag) then
            stop "Triangle references a node tag outside the tag map."
         endif

         if (tag2local(local_tag) == 0) then
            stop "Triangle references an unknown node tag."
         endif

         amesh%cell(i,j) = tag2local(local_tag)
      enddo
   enddo

!-------------------------------------------------------------------------------
! 4. Locate the source node
!-------------------------------------------------------------------------------

   source_node = tag2local(1)

!-------------------------------------------------------------------------------
! 5. Compute cell-centered velocity
!-------------------------------------------------------------------------------

   allocate(velocity(amesh%Ncells))

   do i = 1, amesh%Ncells
      yc = ( amesh%py(amesh%cell(i,1)) &
           + amesh%py(amesh%cell(i,2)) &
           + amesh%py(amesh%cell(i,3)) ) / 3._pr

      velocity(i) = v0 + grad * yc
   enddo

!-------------------------------------------------------------------------------
! 6. Compute node-centered analytical time for linear velocity gradient
!
!    v(y) = v0 + grad*y
!
!    T = acosh( 1 + grad^2*r^2/(2*vs*vr) ) / grad
!
!    The acosh is evaluated as log(a + sqrt(a*a - 1)) to avoid depending
!    on compiler support for the acosh intrinsic.
!-------------------------------------------------------------------------------

   allocate(theory(amesh%Nnodes))

   vsrc = v0 + grad * ys

   do i = 1, amesh%Nnodes
      x = amesh%px(i) - xs
      y = amesh%py(i) - ys

      vrec = v0 + grad * amesh%py(i)
      r2 = x*x + y*y

      arg = 1._pr + grad * grad * r2 / (2._pr * vsrc * vrec)
      arg = max(arg, 1._pr)

      theory(i) = log(arg + sqrt(max(arg*arg - 1._pr, 0._pr))) / grad
   enddo

   theory(source_node) = 0._pr

!-------------------------------------------------------------------------------
! 7. Write input.vtk in ASCII legacy POLYDATA format
!-------------------------------------------------------------------------------

   open(newunit=vtk, file="input.vtk", status="replace", action="write", &
        form="formatted")

   call dumpmeshvtk(vtk, amesh, "vertical velocity gradient example")
   call dumpcellattributevtk(vtk, amesh, velocity, "velocity", .true.)
   call dumpnodeattributevtk(vtk, amesh, theory, "time", .true.)

   close(vtk)

!-------------------------------------------------------------------------------
! 8. Diagnostics
!-------------------------------------------------------------------------------

   write(*,'(a)') "Wrote input.vtk"
   write(*,'(a,i0)') "Number of nodes: ", amesh%Nnodes
   write(*,'(a,i0)') "Number of cells: ", amesh%Ncells
   write(*,'(a,i0)') "Source node, Fortran index: ", source_node
   write(*,'(a,i0)') "Source node, VTK index:     ", source_node - 1
   write(*,'(a,1x,es12.5,1x,es12.5)') "Velocity min/max:", &
                                      minval(velocity), maxval(velocity)
   write(*,'(a,1x,es12.5,1x,es12.5)') "Theory time min/max:", &
                                      minval(theory), maxval(theory)

end program mesh_gradient
