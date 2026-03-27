subroutine tritime2d(Nnodes, Ncells, px, py, pz, cell_flat, &
                               velocity, NdiffNodes, diff_nodes, fast_int, &
                               traveltime) bind(C, name='tritime2d')
  use iso_c_binding
  use generic,  only: pr, pin, infinity
  use LAT_mesh, only: mesh,compnton,free_nton
  use LAT_time, only: diff
  use lists,    only: containern
  use time,     only: timeonevsall2d
  implicit none

  ! --- inputs from Python ---
  integer(c_int), value, intent(in) :: Nnodes, Ncells, NdiffNodes, fast_int
  real(c_double),        intent(in) :: px(Nnodes), py(Nnodes), pz(Nnodes)
  real(c_double),        intent(in) :: velocity(Ncells)
  integer(c_int),        intent(in) :: cell_flat(Ncells*3)   ! row-major from numpy
  integer(c_int),        intent(in) :: diff_nodes(NdiffNodes)
  ! --- in/out: Python pre-sets source nodes to 0, rest to infinity ---
  real(c_double),        intent(inout) :: traveltime(Nnodes)

  type(mesh) :: amesh
  type(diff) :: adiff
  type(containern), dimension(:), allocatable :: nton
  real(pr), dimension(Nnodes) :: tt
  real(pr), dimension(Ncells) :: vel

  ! build mesh type
  amesh%Nnodes = int(Nnodes, pin)
  amesh%Ncells = int(Ncells, pin)
  allocate(amesh%px(Nnodes), amesh%py(Nnodes), amesh%pz(Nnodes))
  amesh%px = real(px, pr);  amesh%py = real(py, pr);  amesh%pz = real(pz, pr)
  allocate(amesh%cell(Ncells, 3))
  ! numpy sends row-major [Ncells,3]; reshape reads column-major, so transpose
  amesh%cell = transpose(reshape(cell_flat, [3, Ncells]))

  ! build diff type
  adiff%fast = (fast_int /= 0)
  adiff%NdiffNodes = int(NdiffNodes, pin)
  allocate(adiff%nodes(max(1, NdiffNodes)))
  if (NdiffNodes > 0) adiff%nodes(1:NdiffNodes) = int(diff_nodes, pin)


  ! copy to local pr-kind arrays (c_double == pr in practice, but avoids creating temporary arrays)
  vel = velocity
  tt  = traveltime

  ! computing the node-to-node edge array nton
  allocate(nton(amesh%Nnodes))
  call compnton(amesh,nton)
  call timeonevsall2d(amesh, vel, tt, nton, adiff)

  ! copy result back to c_double for Python
  traveltime = tt

  ! free pointer-allocated linked-list nodes (not freed automatically on return)
  call free_nton(nton)

end subroutine tritime2d
