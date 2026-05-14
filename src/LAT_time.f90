module LAT_time
  use generic
  use LAT_mesh
  implicit none

  type diff
    logical :: fast
    integer(pin) :: NdiffNodes
    integer(pin),allocatable,dimension(:) :: nodes
  end type diff

  type interface   ! replace reflect 
    logical :: check = .false.  ! can we do this (preset value?)
    integer(pin) :: nnodes 
    integer(pin),allocatable,dimension(:) :: nodes
  end type interface 

  type reflect
    logical :: present
    integer(pin) :: Nnodes
    integer(pin),allocatable,dimension(:) :: nodes
  end type reflect

  real(pr), allocatable, dimension(:) :: kappa
  integer(pin), allocatable, dimension(:) :: mode
  

end module LAT_time
