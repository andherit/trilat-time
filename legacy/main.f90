program main_generic
use LAT_time
use LAT_mesh_util
use calc_time

implicit none
   type(diff) :: adiff
   type(mesh) :: amesh
   real(pr), dimension(:), allocatable :: delta
!   real(pr), dimension(:), allocatable :: time,velocity,
   integer(pin) :: i,j,imx,ptisup
   real(pr) :: mx,perc,vup,vdown,ac,di,x,t,tmin,x1,x2
   integer(pin) :: nxp,nyp
   integer(pin) :: edge_no,n1,n2,cell_no,other_cell !,find_face_cell
   logical :: minnotf
   real(pr) :: vel,a,b,y,io
   real(pr) :: finish,start,dx
   real(pr), parameter :: pi=acos(-1._pr)
   real(pr),allocatable, dimension(:) :: theo_time,res
   real(pr),allocatable, dimension(:) :: res_sface,res_shead,res_splane,res_sedge,res_scplane

   ! read in mesh file
   call read_mesh_file(amesh,adiff)
   write(0,*) 'starting lateration time computing'
   call cpu_time(start)

   call timeonevsall2dV2(amesh,adiff)

   call cpu_time(finish)

   call calc_theory(amesh,theo_time)

   allocate(res(amesh%Nnodes))
   allocate(res_sface(amesh%Nnodes))
   allocate(res_shead(amesh%Nnodes))
   allocate(res_splane(amesh%Nnodes))
   allocate(res_sedge(amesh%Nnodes))
   allocate(res_scplane(amesh%Nnodes))

   ! do i = 1,amesh%Nnodes 
   !    if (theo_time(i) > 0._pr) then 
   !       res(i) = (time(i) - theo_time(i))/theo_time(i)
   !    else 
   !       res(i) = 0._pr
   !    endif
   !    res_sface(i) = sface(i) - theo_time(i)
   !    res_shead(i) = shead(i) - theo_time(i)
   !    res_splane(i) = splane(i) - theo_time(i)
   !    res_sedge(i) = sedge(i) - theo_time(i)
   !    res_scplane(i) = scplane(i) - theo_time(i)
   ! enddo 

!================================================================
!======= write to file =======
   write(*,*) 'output results'
   open(10,file='prova.vtk',form = 'formatted')
   call dumpmeshvtk_jup(10,amesh)
   call dumpcellattributevtk(10,amesh,velocity,'velocity',.true.)
   call dumpnodeattributevtk(10,amesh,time,'time',.true.)

   call dumpnodeattributevtk(10,amesh,mode,'mode',.false.)
   ! call dumpnodeattributevtk(10,amesh,kappa,'kappa',.false.)
   ! call dumpnodeattributevtk(10,amesh,theo_time,'Theoretical_times',.false.)
   ! call dumpnodeattributevtk(10,amesh,res,'Residuals',.false.)		!uncomment for spherical residual check
   ! call dumpnodeattributevtk(10,amesh,res_sface,'Res_Face',.false.)		!uncomment for spherical residual check
   ! call dumpnodeattributevtk(10,amesh,res_shead,'Res_Head',.false.)		!uncomment for spherical residual check
   ! call dumpnodeattributevtk(10,amesh,res_splane,'Res_Plane',.false.)		!uncomment for spherical residual check
   ! call dumpnodeattributevtk(10,amesh,res_sedge,'Res_Edge',.false.)		!uncomment for spherical residual check
   ! call dumpnodeattributevtk(10,amesh,res_scplane,'Res_Cplane',.false.)		!uncomment for spherical residual check
   ! call dumpnodeattributevtk(10,amesh,origidnode,'id_from',.false.)		!uncomment for spherical residual check

  close(10)

end program main_generic
!###############################################################################
