!###############################################################################
! Need a comment - making a test
! 20210305AH Check Required - This routine does not seem to be used
!            To be moved later on the submodule dumpvtk
!###############################################################################
!this is another test
subroutine dumpcellattributevtk(dev,amesh,field,attname,init)
 type(mesh) :: amesh
  real(pr), dimension(amesh%Ncells) :: field
  integer(pin) :: dev
  character*(*) :: attname
  logical :: init

  integer(pin) :: i

  if (init) write(dev,'(a10,i10)') 'CELL_DATA ',amesh%Ncells
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
