module LAT_mesh_util_legacy

contains

!==============================================================================
subroutine mesh_connectivity(amesh)
  !  Create a matrix that defines neighbouring cells
  !  Dimension :  (No of cells , 3)
  !  Row relates to cells
  !  Columns give the cell numbers of the 3 cells with adjoining faces to current cell
  !  If a cell has less than 3 neighbours the index of the current cell is placed in the column
  type(mesh),intent (inout) :: amesh
  integer(pin), dimension(3) :: node_id,neighbour_element
  integer(pin),allocatable,dimension(:) :: common_node
  integer(pin) :: i,j,k,Edge_no,noedges
  integer(pin), dimension(2) :: test_face
  integer(pin), allocatable, dimension(:,:) ::   Face_array,dummy_Edge
  logical :: new_face
  Edge_no = 0
  allocate(Face_array(3*amesh%Ncells,2))
  allocate(dummy_Edge(3*amesh%Ncells,4))
  dummy_Edge = 0

  allocate(amesh%EToE(amesh%Ncells,3))
  ! find boundary of mesh
  do i = 1,amesh%Ncells
       amesh%EToE(i,:) = i   ! initial assumption : there is no overlap
       node_id = amesh%cell(i,:)
       noedges = 0
       neighbour_element = 0
  ! check node id against rest of elements
       do k = 1,amesh%Ncells
         if (i /= k) then
             if (any(node_id(1).eq.amesh%cell(k,:)).or.      &
                 any(node_id(2).eq.amesh%cell(k,:)).or.      &
                 any(node_id(3).eq.amesh%cell(k,:))) then
  ! compare all the nodes of the two elements
                call intersect(node_id,amesh%cell(k,:),common_node)
                if (size(common_node)==2) then
            ! two common nodes = face
                      noedges = noedges+1
                      neighbour_element(noedges) = k
                endif
              endif
         endif

       enddo
       if (noedges > 0) then
         do k = 1,noedges
            amesh%EToE(i,k) = neighbour_element(k)
         enddo
       endif
  ! if any of k index values in EToE = i value => there are less than 3 faces
  ! 3 faces = internal element
  ! 2 faces = boundary element
  ! 1 face = edge element
  !loop inside cells & find faces
         do j = 1,3
          new_face = .true.
          if (j < 3) then
           test_face = amesh%cell(i, (/ j, j+1 /))
          else
           test_face = amesh%cell(i ,(/ j, 1 /))
          endif
  ! check if face is already in list
          if (Edge_no > 0) then
            do k = 1,Edge_no
              if (test_face(1) == Face_array(k,1).and.test_face(2) == Face_array(k,2)          &
              .or.test_face(2) == Face_array(k,1).and.test_face(1) == Face_array(k,2)) then
                 new_face =.false.
                 dummy_Edge(k,2) = i
               endif
           enddo
          endif
          if(new_face) then
              Edge_no = Edge_no+1
              Face_array(Edge_no,:) = test_face
              dummy_Edge(Edge_no,1) = i
              dummy_Edge(Edge_no,3:4) = test_face
           endif
      enddo
  enddo
  amesh%Nedges = Edge_no
  allocate(amesh%Edge(Edge_no,4))
  amesh%Edge = 0
  amesh%Edge = dummy_Edge(1:Edge_no,:)
  deallocate(dummy_Edge)
  write(*,*) 'No of faces.......', amesh%Nedges

  return
end subroutine mesh_connectivity
!==============================================================================
subroutine intersect(A,B,v)
  !  Find the common values to both A and B
  !
  !
  integer(pin),dimension(:),intent(in) :: A, B
  integer(pin), allocatable,dimension(:),intent(out) :: v
  integer(pin), allocatable,dimension(:) :: dummy
  integer(pin),allocatable,dimension(:) :: haveit

  integer(pin) :: i, inx

  allocate(dummy(size(A)))
  inx = 0
  do  i = 1,size(A)
         if (allocated(haveit)) deallocate(haveit)

         call find_1dI(B,'==',A(i),haveit)
         if ( haveit(1) > 0 ) then
                inx = inx +1
                dummy(inx) = A(i)
         endif
  enddo

  allocate(v(inx))
  v = dummy(1:inx)
  deallocate(dummy,haveit)

  return
  end subroutine intersect
!==============================================================================
subroutine find_1dI(array,condt,target_value,idx)
  implicit none
  integer(pin),dimension(:),intent(in) :: array
  integer(pin),allocatable,dimension(:), intent(out) :: idx
  integer(pin),allocatable,dimension(:) :: loc
  integer(pin) :: i,n,target_value,inx
  character(2),intent(in) :: condt

  n = size(array,1)
  allocate(loc(n))
  loc = 0
  inx = 0

  if (condt == '==') then
         do i = 1,n
                if(array(i) == target_value) then
                       inx = inx+1
                       loc(inx) = i
                endif
         enddo
  else
         do i = 1,n
                if(array(i) /= target_value) then
                       inx = inx+1
                       loc(inx) = i
                endif
         enddo
  endif


  if (inx > 0) then  ! target_value is in array
    allocate(idx(inx))
    idx = loc(1:inx)
  else   ! target_value is not in array
    allocate(idx(1))
    idx = -1 !loc(1)

  endif
  deallocate(loc)

  return

  end subroutine find_1dI
!========================================================
!==============================================================================
subroutine read_vtk_mesh_matlab(omesh,mesh_file)
  type(mesh) :: omesh
  character(30),intent(in) :: mesh_file
  integer,allocatable,dimension(:) :: layer_bc
  integer :: i,c1,c2,c3,element_type
  integer :: max_cells,t_cells,t_start,nnodes
  integer :: ipos_start,ipos_end
  logical :: tri_cells_start
  character(len=60) :: line

!  integer :: a,b,c
!  integer,allocatable,dimension(:) :: counter,idx
  open(13,file=trim(adjustl(mesh_file)),form = 'formatted', action = 'read')
  do i = 1,4
         read(13,*)
  enddo
  read(13,'(a)') line
  ipos_start = scan(line,"S",back=.true.)
  ipos_end = scan(line,"d",back=.true.)
  read (line(1+ipos_start:ipos_end),*) nnodes !no unused note at zero
  omesh%Nnodes = nnodes
  write(*,*) 'No. of nodes....',omesh%Nnodes
  allocate(omesh%px(omesh%Nnodes),omesh%py(omesh%Nnodes),omesh%pz(omesh%Nnodes))
  omesh%px = 0._pr
  omesh%py = 0._pr
  omesh%pz = 0._pr
  do i = 1,omesh%Nnodes
     read(13,*) omesh%px(i),omesh%py(i),omesh%pz(i)
!     print*,amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  read(13,*)
  read(13,'(a)')  line
!  write(*,'(a)') line
  ipos_start = scan(line,"CELLS",back=.true.)
!  print*,'ipos ........',ipos_start
  read (line(1+ipos_start:),*) max_cells

  omesh%Ncells = max_cells
  allocate(omesh%cell(omesh%Ncells,3))
  write(*,*)  'no. of cells........ ',omesh%Ncells
  do i = 1,omesh%Ncells
   read(13,*) c1,c2,c3
   omesh%cell(i,:) = (/ c1, c2,c3 /)

  enddo
  close(13)

end subroutine read_vtk_mesh_matlab
!==============================================================================

end module LAT_mesh_util_legacy
