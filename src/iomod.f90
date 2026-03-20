module iomod
use generic
use LAT_mesh
implicit none
!! support module for managing the input and the output of the mesh
!! in different formats
!!   - Ascii legacy VTK
!!     - dumpmeshvtk
!!     - dumpnodeattributevtk
!  - gmsh type
!  - Standard Rupture Format

contains

!###############################################################################
subroutine dumpmeshvtk(dev,amesh)
  type(mesh) :: amesh
  integer(pin) :: dev

  integer(pin) :: i,j

  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a8)') 'distance'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  if( pr == 4) then
     write(dev,'(a7,i10,a6)') 'POINTS ',amesh%Nnodes,' float'
  else
     write(dev,'(a7,i10,a7)') 'POINTS ',amesh%Nnodes,' double'
  endif

  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a9,2i10)') 'POLYGONS ',amesh%Ncells,amesh%Ncells*4
  do i=1,amesh%Ncells
     write(dev,'(a2,$)') '3 '
     write(dev,*) (amesh%cell(i,j)-1,j=1,3)
  enddo
  return
end subroutine dumpmeshvtk
!###############################################################################
subroutine dumpnodeattributevtk(dev,amesh,field,attname,init)
  type(mesh) :: amesh
  real(pr), dimension(amesh%Nnodes) :: field
  integer(pin) :: dev
  character*(*) :: attname
  logical :: init

  integer(pin) :: i

  if (init) write(dev,'(a11,i10)') 'POINT_DATA ',amesh%Nnodes
  if( pr == 4) then
    write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' float 1'
  else
     write(dev,'(a8,a,a10)') 'SCALARS ',trim(attname),' double 1'
  endif
  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,amesh%Nnodes
     write(dev,*) field(i)
  enddo
  return
end subroutine dumpnodeattributevtk
!###############################################################################
subroutine dumpcellattributevtk(dev,amesh,field,attname,init)
  type(mesh) :: amesh
  real(pr), dimension(amesh%Ncells) :: field
  integer(pin) :: dev
  character*(*) :: attname
  logical :: init

  integer(pin) :: i

  if (init) write(dev,'(a11,i10)') 'CELL_DATA ',amesh%Ncells
  if( pr == 4) then
    write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' float 1'
  else
     write(dev,'(a8,a,a10)') 'SCALARS ',trim(attname),' double 1'
  endif
  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,amesh%Ncells
     write(dev,*) field(i)
  enddo
  return
end subroutine dumpcellattributevtk
end module iomod
