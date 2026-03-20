module LAT_velocity
use generic
use LAT_time

implicit none

contains


!########################################################
subroutine assign_velocity(amesh,cell_id,adiff)
use LAT_mesh
use LAT_time

type(mesh) :: amesh
type(diff) :: adiff
real(pr) :: c_x,c_y,c_z,r

integer(pin),allocatable,dimension(:) :: cell_id
real(pr),allocatable, dimension(:) :: vlayer,v_grad,depth0
real(pr),allocatable, dimension(:) :: radius,v_profile!,velocity
real(pr) :: bary_depth
integer(pin) :: i,j,nlines,is,ibc,id
integer(pin),allocatable,dimension(:) :: s_nodes,bc_nodes
real(pr) :: x,y,z,surf_radius,bc_radius


write(*,*) 'read velocity file'
call read_velocity_jup(vlayer,v_grad,depth0)
!call read_velocity(nlines,radius,v_profile)
write(*,*) 'finished reading  velocity file'

! apply velocity to layer
do i = 1,amesh%Ncells
     id = cell_id(i)
     if (abs(v_grad(id)) < epsilon(vlayer(1))) then
           velocity(i) = vlayer(id)
     else
           bary_depth = sum(amesh%py(amesh%cell(i,:)))/3._pr ! depth of barycentre of cell
           if (bary_depth < depth0(id)) then
              write(*,*) 'ERROR: velocity gradient incorrect'
              stop
           endif
           velocity(i) = vlayer(id) + v_grad(id)*(bary_depth-depth0(id))
     endif
enddo


!
! allocate(s_nodes(amesh%Nnodes),bc_nodes(amesh%Nnodes))
!
! s_nodes = 0
! bc_nodes = 0
!
! surf_radius = maxval(radius(:)) !
! !write(*,*) 'radius of sphere....', surf_radius
!
! bc_radius =  5961._pr   ! boundary radius that could cause diffraction
! !bc_radius =  3482._pr
! is = 0
! ibc = 0
! do i = 1,amesh%Ncells
! ! get baricentre of cell radius
!     c_x = sum(amesh%px(amesh%cell(i,:)))/3.
!     c_y = sum(amesh%py(amesh%cell(i,:)))/3.
!     c_z = sum(amesh%pz(amesh%cell(i,:)))/3.
!     r = sqrt(c_x**2+c_y**2+c_z**2)
! ! calc  velocity for cell
!     velocity(i) = find_velocity(nlines,v_profile,radius,r)
! ! find layer boundaries
!     do j = 1,3
!        x = amesh%px(amesh%cell(i,j))
!        y = amesh%py(amesh%cell(i,j))
!        z = amesh%pz(amesh%cell(i,j))
!        r = sqrt(x**2+y**2+z**2)
!        if (abs(r-surf_radius)<=  1._pr) then
!            is = is+1
!            ! log surface node
!            s_nodes(is) = amesh%cell(i,j)
!        elseif (abs(r-bc_radius)<= 1._pr) then
!            ibc = ibc+1
!            ! log boundary nodes
!            bc_nodes(ibc) = amesh%cell(i,j)
!        endif
!     enddo
! enddo
!
!

!write(*,*)'number of surface nodes......',is
!time%Nsurf = is
!allocate(amesh%s_nodes(amesh%Nsurf))
!time%s_nodes = s_nodes(1:is)


adiff%fast = .true.  ! turn off diffraction
if (.not.adiff%fast) then
  write(*,*)'number of diffraction nodes.......',ibc
  adiff%Nnodes = ibc
  allocate(adiff%nodes(adiff%Nnodes))
  adiff%nodes = bc_nodes(1:ibc)
else
  adiff%Nnodes = 1
  allocate(adiff%nodes(adiff%Nnodes))
  adiff%nodes = 1
endif
!time%Ndiff = ibc
!allocate(amesh%diff_nodes(amesh%Ndiff))
!time%diff_nodes(1:ibc) = bc_nodes(1:ibc)


!write(*,*)'number of diffraction nodes.......',ibc
!time%Ndiff = ibc
!allocate(amesh%diff_nodes(amesh%Ndiff))
!time%diff_nodes(1:ibc) = bc_nodes(1:ibc)

! if you want to check all nodes, uncomment next lines
! amesh%Ndiff = amesh%Nnodes
! write(*,*)'number of diffraction nodes.......',amesh%Ndiff
! allocate(amesh%diff_nodes(amesh%Ndiff))
! do i =1,amesh%Nnodes
!    amesh%diff_nodes(i) = i
! enddo


return
end subroutine assign_velocity
!==============================================================================
function find_velocity(nlines,v_profile,radius,z)
  real(pr) :: z,find_velocity,m,c
  real(pr),dimension(nlines) :: v_profile,radius
  integer :: nlines,i

  find_velocity = 0._pr
  ! find correct depth
  ! assuming radius starts from outside and goes to core
  do i = 1,nlines
    if (z >= radius(i)) then
       exit
    endif
  enddo

  if (i >= nlines) then
! point is inside core
     find_velocity = v_profile(nlines)
  elseif(maxval(radius) - z <= 35.) then
     find_velocity = 6.5
  else
!          find_velocity = v_profile(i-1)
!     if (abs(radius(i+1)-radius(i))<= epsilon(1._pr)) then
!        find_velocity = v_profile(i)
!     else
! linear interpolate between nearest two velocity points
        m = (z-radius(i-1))/(radius(i)-radius(i-1))
        find_velocity = v_profile(i-1)*(1._pr-m)+v_profile(i)*m
!     endif
  endif
  ! write(*,*) 'z, radius used', z,radius(i),radius(i-1)
  ! write(*,*) 'inner most value', v_profile(nlines)
  ! write(*,*) 'outer most value', v_profile(1)
  ! write(*,*)'value choosen',find_velocity
  ! pause
end function find_velocity
!==============================================================================
subroutine read_velocity_jup(v,dv,d0)
  real(pr),allocatable,dimension(:) :: v,dv,d0
  integer(pin) :: i,id,ipos_end
  integer(pin) :: nlayers
  character(len=60) :: line,fname

  fname = 'velocity.txt'
  call check_file(fname)
  open(13,file='velocity.txt',form = 'formatted', action = 'read')


  read(13,'(a)') line
  ipos_end = scan(line,":",back=.false.)
  read (line(ipos_end+1:60),*) nlayers
  write(*,*) nlayers

  read(13,*)  ! header
  allocate(v(nlayers),dv(nlayers),d0(nlayers))
  do i = 1,nlayers
      read(13,*) id, v(i), dv(i),d0(i)
!      write(*,*) id, v(i), dv(i),d0(i)
  enddo
  close(13)


  return
end subroutine read_velocity_jup
!==============================================================================
subroutine check_file(fname)
  character(60),intent(in) :: fname
  logical :: file_exists

  inquire(file=trim(adjustl(fname)), exist=file_exists)
  if (.not.file_exists) then
      write(*,*) 'ERROR: an input file was not found'
      write(*,*) 'Programme stopping '
      stop
  endif

 return
end subroutine
!==============================================================================
subroutine read_velocity(nlines,radius,v_profile)

real(pr),allocatable,dimension(:) :: depth,radius,vp,vs,v_profile
integer :: i,ierr, nlines
character :: line*100
logical :: file_exists
ierr = 0


inquire(file="IASP91.csv", exist=file_exists)
if (.not.file_exists) then
    write(*,*) 'ERROR: no velocity file found'
    write(*,*) 'Programme stopping '
    stop
endif

open(10,file='IASP91.csv')
 nlines= 0
do while (ierr == 0)
   nlines = nlines+1
   read(10,'(a)',iostat=ierr) line
enddo
rewind(10)
nlines = nlines-1
allocate(depth(nlines),radius(nlines),vp(nlines),vs(nlines))
write(*,*) 'number of lines   ',nlines
do i = 1,nlines
 read(10,*) depth(i),radius(i),vp(i),vs(i)
! write(*,120) a,b,c,d
enddo
close(10)
!write(*,*) 'last values   ',depth(nlines),radius(nlines),vp(nlines),vs(nlines)
!120 format (f4.2, 1x , f7.2, 1x, f6.4, 1x, f6.4)
!120 format (f4.2, f7.2, f6.4, f6.4)
allocate(v_profile(nlines))
v_profile = vp


return
end subroutine read_velocity
!==============================================================================

end module LAT_velocity
