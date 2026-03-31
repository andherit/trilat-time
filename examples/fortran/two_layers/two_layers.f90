program two_layers
!! test of the one vs all time routine on a 2D mesh with two layers of different velocity.
!! We use the fast option (no diffraction)
use generic
use lists
use LAT_mesh
use LAT_time
use time
use iomod
implicit none

   type(mesh) :: amesh
   real(pr), dimension(:), allocatable :: traveltime,velocity,theo_time,error_time
   real(pr) :: start_time,finish_time
   type(diff) :: adiff
   integer(pin) :: i,k
   type(containern), dimension(:), allocatable :: nton
   
! fast option on
! the disffraction array is not used in this case.
   adiff%fast=.true.
! load the two-layer model in ascii vtk
   call load_model(amesh,velocity,theo_time)
! source is at (5000,1500) and the fortran code is 1-based.
   k=7
! pre-computation of :
! - the node-to-node edge array nton
! - the initialisation of the time array
! - the initialisation of the time at the source node to zero 
!   call pre_timeonevsall2d_onvertex(amesh,k,traveltime,nton)
! computation
   write(0,*) 'starting time computing'
   call cpu_time(start_time)
   call pre_timeonevsall2d_onvertex(amesh,k,traveltime,nton)
   call timeonevsall2d(amesh,velocity,traveltime,nton,adiff)
   call cpu_time(finish_time)
   write(*,*) 'cpu time for timeonevsall2d : ',finish_time-start_time
! relative error in percent with respect to the theoretical solution
   allocate(error_time(amesh%Nnodes))
   do i=1,amesh%Nnodes
      if (i == k) then
         error_time(i)=0._pr
      elseif (abs(theo_time(i)) <= water_level(1._pr)) then
         error_time(i)=0._pr
      else
         error_time(i)=100._pr*(traveltime(i)-theo_time(i))/theo_time(i)
      endif
   enddo
! freeing the memory of the nton array
   call free_nton(nton)
! dumping the result in ascii VTK
   write(*,*) 'dumping time in result.vtk'
   open(10,file='result.vtk',form = 'formatted')
   call dumpmeshvtk(10,amesh)
   call dumpcellattributevtk(10,amesh,velocity,'velocity',.true.)
   call dumpnodeattributevtk(10,amesh,traveltime,'time',.true.)
   call dumpnodeattributevtk(10,amesh,theo_time,'theo_time',.false.)
   call dumpnodeattributevtk(10,amesh,error_time,'error_time',.false.)
   close(10)

contains
!#########################################################################
subroutine load_model(amesh,velocity,theo_time)
! load the two-layer vtk file called input_velocity.vtk
! 
use generic
use LAT_mesh
implicit none
type(mesh) :: amesh
real(pr), dimension(:), allocatable :: velocity
real(pr), dimension(:), allocatable :: theo_time

integer :: i,j,idum,cell_list_size
character(len=32) :: keyword, datatype

! read VTK header and get POINTS count
  open(12,file="input_velocity.vtk")
  do i=1,4
     read(12,*)
  enddo
  read(12,*) keyword, amesh%Nnodes, datatype
  allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
  do i=1,amesh%Nnodes
     read(12,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
! read POLYGONS header and get cell count
  read(12,*) keyword, amesh%Ncells, cell_list_size
  allocate(amesh%cell(amesh%Ncells,3))
! first number is the cell type, we ignore it, then the 3 nodes of the triangle
  do i=1,amesh%Ncells
     read(12,*) idum,(amesh%cell(i,j),j=1,3)
  enddo
! vtk files are 0-based, we need to add 1 to the cell connectivity
  amesh%cell=amesh%cell+1
! reading the velocity attribute.
! the file has a void line and the CELL_DATA header is 3 lines long. thus we need to read 4 lines before the velocity data
  do i=1,4
     read(12,*)
  enddo
! allocate the velocity array and read the velocity value for each cell
  allocate(velocity(amesh%Ncells))
   do i=1,amesh%Ncells
       read(12,*) velocity(i)
   enddo
! reading the theoretical time node attribute.
! the file has a void line and the POINT_DATA header is 3 lines long.
  do i=1,4
     read(12,*)
  enddo
! allocate the theoretical time array and read the value for each node
  allocate(theo_time(amesh%Nnodes))
   do i=1,amesh%Nnodes
       read(12,*) theo_time(i)
   enddo
  close(12)
  return
end subroutine load_model
!###############################################################################
end program two_layers
