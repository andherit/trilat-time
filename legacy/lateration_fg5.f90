program lateration_fg5
use lateration_wk3
implicit none

   type(mesh) :: amesh
   real, dimension(:), allocatable :: time,delta,velocity
   real :: vup,vdown,ac,di,t,tmin,x1,x2,x,finish,start
   integer :: i,j,nu
   logical :: minnotf

   
   di=1.
   call loadmesh(amesh)
   vup=1000.
   vdown=2000.
   ac=asin(vup/vdown)
   allocate(velocity(amesh%Ncells))
! velocity field
   do i=1,amesh%Ncells
      nu=0
      do j=1,3
         if (amesh%py(amesh%cell(i,j)) <= 700.) nu=nu+1
      enddo
      if (nu == 3) then
         velocity(i)=vup
      else
         velocity(i)=vdown
      endif
   enddo
! time initialization
   allocate(time(amesh%Nnodes))
   time=infinity
   time(7)=0.
   write(*,*) 'source position :',amesh%px(7),amesh%py(7)
! lateration call
   write(0,*) 'time computing'
   call cpu_time(start)
   call timeonevsall2d(amesh,time,velocity)
   call cpu_time(finish)
   write(*,*) ' cpu time : ',finish-start
! error computation
   write(0,*) 'starting error estimation'
   allocate(delta(amesh%Nnodes))
   do i=1,amesh%Nnodes
! direct waves
   if (amesh%py(i) .le. 700.) then
      delta(i)=sqrt((amesh%px(i)-amesh%px(7))**2.+(amesh%py(i)-amesh%py(7))**2.)/vup
! conic waves
      if (abs(amesh%px(i)-amesh%px(7)) > 150.*tan(ac)+(700.-amesh%py(i))*tan(ac)) then
         x1=150.*tan(ac)
         x2=(700.-amesh%py(i))*tan(ac)
         delta(i)=min((150.+700.-amesh%py(i))/cos(ac)/vup+(abs(amesh%px(i)-amesh%px(7))-x1-x2)/vdown,delta(i))
      endif
      if (delta(i) /= 0.) delta(i)=(delta(i)-time(i))/delta(i)*100.
   else
! refracted waves
         minnotf=.true.
         tmin=infinity
         x=0.
         tmin=sqrt((amesh%px(i)-amesh%px(7))**2.+(amesh%py(i)-700.)**2.)/vdown+150./vup
         do while (minnotf)
            x=x+di
   t=sqrt((abs(amesh%px(i)-amesh%px(7))-x)**2.+(amesh%py(i)-700.)**2.)/vdown+sqrt(x**2.+150.**2.)/vup
            if (t < tmin) then
               tmin=t
            else
               minnotf=.false.
            endif
         enddo
         delta(i)=tmin
      delta(i)=(delta(i)-time(i))/delta(i)*100.
      endif
   enddo
   write(*,*) 'dumping time'
   open(10,file='prova.vtk',form = 'formatted')
   call dumpmeshvtk(10,amesh)
   call dumpcellattributevtk(10,amesh,velocity,'velocity',.true.)
   call dumpnodeattributevtk(10,amesh,time,'time',.true.)
   call dumpnodeattributevtk(10,amesh,delta,'delta',.false.)
  close(10)
end program lateration_fg5
!###############################################################################
subroutine loadmesh(amesh)
use lateration_wk3
implicit none

   type(mesh) :: amesh

   integer :: i,j,idum

  amesh%Nnodes=13401
  allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
  open(12,file="noblefig5.vtk")
  do i=1,5
     read(12,*)
  enddo
  do i=1,amesh%Nnodes
     read(12,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  read(12,*)
  amesh%Ncells=26400
  allocate(amesh%cell(amesh%Ncells,3))
  do i=1,amesh%Ncells
     read(12,*) idum,(amesh%cell(i,j),j=1,3)
  enddo
  amesh%cell=amesh%cell+1
  return
end subroutine loadmesh
!###############################################################################
