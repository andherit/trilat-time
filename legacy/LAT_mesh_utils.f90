module LAT_mesh_util
use generic
use LAT_mesh
implicit none

interface dumpnodeattributevtk
    module procedure idumpnodeattributevtk
    module procedure rdumpnodeattributevtk
end interface 


contains

!###############################################################################
subroutine read_mesh_file(amesh,adiff)
  use LAT_time
  use LAT_velocity
  integer(pin),allocatable,dimension(:) :: cell_id
  type(mesh) :: amesh
  type(diff) :: adiff
  character(30) :: mesh_file

  ! read unstructured mesh from file
  write(*,*) 'enter read_vtk_mesh'
  !  mesh_file = 'geoid_cross_section.vtk'
  mesh_file = 'test2.vtk'
!   call read_vtk_mesh_jup(amesh,cell_id,mesh_file)
   call read_vtk_mesh_jupv2(amesh,mesh_file,adiff)

   write(*,*) 'exit read_vtk_mesh'

! initiate time & velocity arrays
!  call initiate_time(amesh)

  ! assign velocities to cells and define possible diffraction points
!  call assign_velocity(amesh,cell_id,adiff)

  call mesh_connectivity(amesh)

  return
end subroutine read_mesh_file
!###############################################################################
subroutine initiate_time(amesh)
 use LAT_time
 type(mesh) :: amesh
 integer(pin)  :: nodetodump
 real(pr) :: stime

! set up time array
 allocate(time(amesh%Nnodes))
  time = infinity




  ! select a node as source for travel time
  !nodetodump = 3367 !3367 away from any of the ellipse arcs - GMSH with mshspc=100
  !time(nodetodump) = 0._pr

  ! read source file
  open(13,file='source.txt',form = 'formatted', action = 'read')
  read(13,*)    ! header
  read(13,*)   nodetodump, stime
  close(13)

  time(nodetodump) = stime

! initiate velocity field
  allocate(velocity(amesh%Ncells))
  velocity = 0._pr

  return
end subroutine initiate_time
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

  integer(pin),dimension(3,2) :: cell_side
  integer(pin),dimension(2) :: face
  integer :: Faces,NFaces, id
  integer(pin), allocatable, dimension(:,:) ::  FtoNode,FtoNodeT,FtoF
  integer(pin), allocatable, dimension(:,:,:) ::   FVsNodes
  integer(pin), dimension(3,2):: Edge
  integer(pin) :: cell1,cell2
  integer(pin), dimension(2) :: vertices

  Faces = 3

  ! Edge_no = 0
  ! allocate(Face_array(3*amesh%Ncells,2))
  ! allocate(dummy_Edge(3*amesh%Ncells,4))
  ! dummy_Edge = 0


  write(*,*) 'starting allocation'
  !call flush()

  allocate(amesh%EToE(amesh%Ncells,Faces))
  !allocate(FVsNodes(amesh%Nnodes,amesh%Nnodes,2))  ! interger
  allocate(amesh%FVsNodes(amesh%Nnodes,amesh%Nnodes,2))  ! interger


  Edge(1,:) = (/ 1, 2/)
  Edge(2,:) = (/ 2, 3/)
  Edge(3,:) = (/ 3, 1/)


  write(*,*) 'create FVsNodes'
!  call flush()
!  FVsNodes = 0
  amesh%FVsNodes = 0

  do i = 1, amesh%Ncells
    do j = 1,Faces
       vertices = amesh%cell(i,Edge(j,:))    ! (Ncells, 3)
       do k = 1,2
          if (amesh%FVsNodes(vertices(3-k),vertices(k),1) == 0) then
            amesh%FVsNodes(vertices(3-k),vertices(k),1) = i
          elseif (amesh%FVsNodes(vertices(3-k),vertices(k),2) == 0) then
            amesh%FVsNodes(vertices(3-k),vertices(k),2) = i
          endif
       enddo
    enddo
  enddo

  do i =1,amesh%Ncells
      amesh%EToE(i,:) = i    ! null element for this matrix is identity
  enddo
  write(*,*) 'start EToE'
  do j = 1,amesh%nnodes-1
     do k = j+1,amesh%nnodes
            cell1 = amesh%FVsNodes(j,k,1)
            cell2 = amesh%FVsNodes(j,k,2)
            if ((cell1 /= 0).and.(cell2 /= 0 ))then
              do i = 1,3
                if (amesh%EToE(cell1,i) == cell1) then
                    amesh%EToE(cell1,i) = cell2
                    exit
                endif
              enddo
              do i = 1,3
                if (amesh%EToE(cell2,i) == cell2) then
                    amesh%EToE(cell2,i) = cell1
                    exit
                endif
              enddo
            endif
      enddo
  enddo



  ! do i = 1, amesh%Ncells
  !   do j = 1,Faces
  !      vertices = amesh%cell(i,Edge(j,:))    ! (Ncells, 3)
  !      do k = 1,2
  !         if (FVsNodes(vertices(3-k),vertices(k),1) == 0) then
  !              FVsNodes(vertices(3-k),vertices(k),1) = i
  !         elseif (FVsNodes(vertices(3-k),vertices(k),2) == 0) then
  !              FVsNodes(vertices(3-k),vertices(k),2) = i
  !         endif
  !      enddo
  !   enddo
  ! enddo

  ! do i =1,amesh%Ncells
  !     amesh%EToE(i,:) = i    ! null element for this matrix is identity
  ! enddo
  ! write(*,*) 'start EToE'
  ! do j = 1,amesh%nnodes-1
  !    do k = j+1,amesh%nnodes
  !           cell1 = FVsNodes(j,k,1)
  !           cell2 = FVsNodes(j,k,2)
  !           if ((cell1 /= 0).and.(cell2 /= 0 ))then
  !             do i = 1,3
  !               if (amesh%EToE(cell1,i) == cell1) then
  !                   amesh%EToE(cell1,i) = cell2
  !                   exit
  !               endif
  !             enddo
  !             do i = 1,3
  !               if (amesh%EToE(cell2,i) == cell2) then
  !                   amesh%EToE(cell2,i) = cell1
  !                   exit
  !               endif
  !             enddo
  !           endif
  !     enddo
  ! enddo



  ! ! find boundary of mesh
  ! do i = 1,amesh%Ncells
  !      amesh%EToE(i,:) = i   ! initial assumption : there is no overlap
  !      node_id = amesh%cell(i,:)
  !      noedges = 0
  !      neighbour_element = 0
  ! ! check node id against rest of elements
  !      do k = 1,amesh%Ncells
  !        if (i /= k) then
  !            if (any(node_id(1).eq.amesh%cell(k,:)).or.      &
  !                any(node_id(2).eq.amesh%cell(k,:)).or.      &
  !                any(node_id(3).eq.amesh%cell(k,:))) then
  ! ! compare all the nodes of the two elements
  !               call intersect(node_id,amesh%cell(k,:),common_node)
  !               if (size(common_node)==2) then
  !           ! two common nodes = face
  !                     noedges = noedges+1
  !                     neighbour_element(noedges) = k
  !               endif
  !             endif
  !        endif

  !      enddo
  !      if (noedges > 0) then
  !        do k = 1,noedges
  !           amesh%EToE(i,k) = neighbour_element(k)
  !        enddo
  !      endif
  ! ! if any of k index values in EToE = i value => there are less than 3 faces
  ! ! 3 faces = internal element
  ! ! 2 faces = boundary element
  ! ! 1 face = corner element
  ! !loop inside cells & find faces
  !        do j = 1,3
  !         new_face = .true.
  !         if (j < 3) then
  !          test_face = amesh%cell(i, (/ j, j+1 /))
  !         else
  !          test_face = amesh%cell(i ,(/ j, 1 /))
  !         endif
  ! ! check if face is already in list
  !         if (Edge_no > 0) then
  !           do k = 1,Edge_no
  !             if (test_face(1) == Face_array(k,1).and.test_face(2) == Face_array(k,2)          &
  !             .or.test_face(2) == Face_array(k,1).and.test_face(1) == Face_array(k,2)) then
  !                new_face =.false.
  !                dummy_Edge(k,2) = i
  !              endif
  !          enddo
  !         endif
  !         if(new_face) then
  !             Edge_no = Edge_no+1
  !             Face_array(Edge_no,:) = test_face
  !             dummy_Edge(Edge_no,1) = i
  !             dummy_Edge(Edge_no,3:4) = test_face
  !          endif
  !     enddo
  ! enddo
  ! amesh%Nedges = Edge_no
  ! allocate(amesh%Edge(Edge_no,4))
  ! amesh%Edge = 0
  ! amesh%Edge = dummy_Edge(1:Edge_no,:)
  ! deallocate(dummy_Edge)
  ! write(*,*) 'No of faces.......', amesh%Nedges

  return
end subroutine mesh_connectivity
!==============================================================================
! subroutine intersect(A,B,v)
!   !  Find the common values to both A and B
!   !
!   !
!   integer(pin),dimension(:),intent(in) :: A, B
!   integer(pin), allocatable,dimension(:),intent(out) :: v
!   integer(pin), allocatable,dimension(:) :: dummy
!   integer(pin),allocatable,dimension(:) :: haveit

!   integer(pin) :: i, inx

!   allocate(dummy(size(A)))
!   inx = 0
!   do  i = 1,size(A)
!          if (allocated(haveit)) deallocate(haveit)

!          call find_1dI(B,'==',A(i),haveit)
!          if ( haveit(1) > 0 ) then
!                 inx = inx +1
!                 dummy(inx) = A(i)
!          endif
!   enddo

!   allocate(v(inx))
!   v = dummy(1:inx)
!   deallocate(dummy,haveit)

!   return
!   end subroutine intersect
! !==============================================================================
!   !==============================================================================
! subroutine find_1dI(array,condt,target_value,idx)
!   implicit none
!   integer(pin),dimension(:),intent(in) :: array
!   integer(pin),allocatable,dimension(:), intent(out) :: idx
!   integer(pin),allocatable,dimension(:) :: loc
!   integer(pin) :: i,n,target_value,inx
!   character(2),intent(in) :: condt

!   n = size(array,1)
!   allocate(loc(n))
!   loc = 0
!   inx = 0

!   if (condt == '==') then
!          do i = 1,n
!                 if(array(i) == target_value) then
!                        inx = inx+1
!                        loc(inx) = i
!                 endif
!          enddo
!   else
!          do i = 1,n
!                 if(array(i) /= target_value) then
!                        inx = inx+1
!                        loc(inx) = i
!                 endif
!          enddo
!   endif


!   if (inx > 0) then  ! target_value is in array
!     allocate(idx(inx))
!     idx = loc(1:inx)
!   else   ! target_value is not in array
!     allocate(idx(1))
!     idx = -1 !loc(1)

!   endif
!   deallocate(loc)

!   return

!   end subroutine find_1dI
! !========================================================
!==============================================================================
! subroutine mesh_EToE(amesh)
! !  Create a matrix that defines neighbouring cells
! !  Dimension :  (No of cells , 3)
! !  Row relates to cells
! !  Columns give the cell numbers of the 3 cells with adjoining faces to current cell
! !  If a cell has less than 3 neighbours the index of the current cell is placed in the column
! !
! !
! type(mesh),intent (inout) :: amesh
! integer(pin), dimension(3) :: node_id,neighbour_element
! integer(pin),allocatable,dimension(:) :: common_node
! integer(pin) :: i,k,nofaces
!
!
! allocate(amesh%EToE(amesh%Ncells,3))
! ! find boundary of mesh
! do i = 1,amesh%Ncells
!      amesh%EToE(i,:) = i
!      node_id = amesh%cell(i,:)
!      nofaces = 0
!      neighbour_element = 0
! ! check node id against rest of array
!      do k = 1,amesh%Ncells
!        if (i /= k) then
!            if (any(node_id(1).eq.amesh%cell(k,:)).or.      &
!                any(node_id(2).eq.amesh%cell(k,:)).or.      &
!                any(node_id(3).eq.amesh%cell(k,:))) then
! ! compare all the nodes of the two elements
!               call intersect(node_id,amesh%cell(k,:),common_node)
!               if (size(common_node)==2) then
!           ! two common nodes = face
!                     nofaces = nofaces+1
!                     neighbour_element(nofaces) = k
!               endif
!             endif
!        endif
!
!      enddo
!      if (nofaces > 0) then
!        do k = 1,nofaces
!           amesh%EToE(i,k) = neighbour_element(k)
!        enddo
!      endif
! ! if any of k index values in EToE = i value => there are less than 3 faces
! ! 3 faces = internal element
! ! 2 faces = boundary element
! ! 1 face = edge element
!
! enddo
!
! return
! end subroutine mesh_EToE
!###############################################################################
subroutine read_vtk_mesh(amesh,mesh_file)


  type(mesh) :: amesh
  character(30),intent(in) :: mesh_file
  integer(pin),allocatable,dimension(:) :: layer_bc
  integer(pin) :: i, id,nsize,c1,c2,c3,element_type
  integer(pin) :: max_cells,t_cells,t_start,nnodes
  integer(pin) :: ipos_start,ipos_end
  logical :: tri_cells_start
  character(len=60) :: line

  open(13,file=trim(adjustl(mesh_file)),form = 'formatted', action = 'read')
  do i = 1,4
         read(13,*)
  enddo
  read(13,'(a)') line
  ipos_start = scan(line,"S",back=.true.)
  ipos_end = scan(line,"d",back=.false.)
!  write(*,*) ipos_start,ipos_end
  read (line(1+ipos_start:ipos_end),*) nnodes
!  amesh%Nnodes = nnodes-1  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
  amesh%Nnodes = nnodes  !
  write(*,*) 'No. of nodes....',amesh%Nnodes
  allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
  amesh%px = 0._pr
  amesh%py = 0._pr
  amesh%pz = 0._pr
!  read(13,*)  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
  do i = 1,amesh%Nnodes
     read(13,*) amesh%px(i),amesh%py(i),amesh%pz(i)
!     print*,amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  read(13,*)
  read(13,'(a)')  line
!  write(*,'(a)') line
  ipos_start = scan(line,"CELLS",back=.true.)
!  print*,'ipos ........',ipos_start
  read (line(1+ipos_start:),*) max_cells,nsize
!  read(13,'(5x,i6,i7)')amesh%Ncells,nsize
!  write(*,*) max_cells,nsize
!  write(*,*)'no. of cells',max_cells
!  stop
  t_cells = 0
  tri_cells_start = .false.
  do i = 1,max_cells
      read(13,*) id
      if (id == 3) then
         if (.not.tri_cells_start) then
           t_start = i
            tri_cells_start = .true.
         endif
         t_cells = t_cells+1
       endif
  enddo

!  write(*,*) 'first triangle at.....',t_start
!  write(*,*)'no. of triangular cells.....',t_cells
  rewind(13)
  do i = 1,6+amesh%Nnodes+t_start
!  do i = 1,6+amesh%Nnodes+t_start+1 !add one for the node you jumped
      read(13,*)
  enddo
!  write(*,*) 'current position......',6+amesh%Nnodes+t_start
  amesh%Ncells = t_cells
  allocate(amesh%cell(amesh%Ncells,3))
  write(*,*)  'no. of cells........ ',amesh%Ncells
  do i = 1,amesh%Ncells
   read(13,*) id,c1,c2,c3
   amesh%cell(i,:) = (/ c1+1, c2+1,c3+1 /)
!   amesh%cell(i,:) = (/ c1, c2,c3 /)  !this is because we have removed the first node at the origin

  enddo
  close(13)

return
end subroutine read_vtk_mesh
!###############################################################################
!###############################################################################
subroutine read_vtk_mesh_jupv2(amesh,mesh_file,adiff)
  use LAT_time

  type(mesh) :: amesh
  type(diff) :: adiff

  character(30),intent(in) :: mesh_file
  integer(pin),allocatable,dimension(:) :: layer_bc,list_diff
!  real(pr),allocatable,dimension(:) :: velocity,time
  integer(pin) :: i, id,nsize,c1,c2,c3,element_type
  integer(pin) :: max_cells,t_cells,t_start,nnodes
  integer(pin) :: ipos_start,ipos_end,counter,diff_switch
  real(pr) :: itime
  logical :: tri_cells_start
  character(len=60) :: line

  integer(pin) :: diff_level
  character(100) :: char2int



  open(13,file=trim(adjustl(mesh_file)),form = 'formatted', action = 'read')
  do i = 1,4
         read(13,*)
  enddo
  read(13,'(a)') line
  ipos_start = scan(line,"S",back=.true.)
  ipos_end = scan(line,"d",back=.false.)
!  write(*,*) ipos_start,ipos_end
  read (line(1+ipos_start:ipos_end),*) nnodes
!  amesh%Nnodes = nnodes-1  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
  amesh%Nnodes = nnodes  !
  write(*,*) 'No. of nodes....',amesh%Nnodes
  allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
  amesh%px = 0._pr
  amesh%py = 0._pr
  amesh%pz = 0._pr
!  read(13,*)  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
  do i = 1,amesh%Nnodes
     read(13,*) amesh%px(i),amesh%py(i),amesh%pz(i)
!     print*,amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  read(13,*)
  read(13,'(a)')  line
!  write(*,'(a)') line
  ipos_start = scan(line,"CELLS",back=.true.)
!  print*,'ipos ........',ipos_start
  read (line(1+ipos_start:),*) max_cells,nsize
!  read(13,'(5x,i6,i7)')amesh%Ncells,nsize
!  write(*,*) max_cells,nsize
!  write(*,*)'no. of cells',max_cells
!  stop
  t_cells = 0
  tri_cells_start = .false.
  do i = 1,max_cells
      read(13,*) id
      if (id == 3) then
         if (.not.tri_cells_start) then
           t_start = i
            tri_cells_start = .true.
         endif
         t_cells = t_cells+1
       endif
  enddo

!  write(*,*) 'first triangle at.....',t_start
!  write(*,*)'no. of triangular cells.....',t_cells
  rewind(13)
  do i = 1,6+amesh%Nnodes+t_start
!  do i = 1,6+amesh%Nnodes+t_start+1 !add one for the node you jumped
      read(13,*)
  enddo
!  write(*,*) 'current position......',6+amesh%Nnodes+t_start
  amesh%Ncells = t_cells
  allocate(amesh%cell(amesh%Ncells,3))
  write(*,*)  'no. of cells........ ',amesh%Ncells
  do i = 1,amesh%Ncells
   read(13,*) id,c1,c2,c3
   amesh%cell(i,:) = (/ c1+1, c2+1,c3+1 /)
!   amesh%cell(i,:) = (/ c1, c2,c3 /)  !this is because we have removed the first node at the origin
  enddo
  read(13,*)     ! blank line


! next read velocity
  do i = 1,3
   read(13,*)     ! skip header informatin on CellEntityIds
  enddo
  allocate(velocity(amesh%Ncells))
  velocity = 0._pr
  do i = 1,amesh%Ncells
    read(13,*) velocity(i)
  enddo
  read(13,*)     ! blank line

! set up time array
   allocate(time(amesh%Nnodes))
   allocate(kappa(amesh%Nnodes))  ! inverse of curvature 
   allocate(mode(amesh%Nnodes))    ! operator type 
  !  kappa = infinity !0._pr
  !  mode = -1 
   time = infinity


! set up for debugging solvers 
  allocate(sface(amesh%Nnodes))
  allocate(shead(amesh%Nnodes))
  allocate(splane(amesh%Nnodes))
  allocate(sedge(amesh%Nnodes))
  allocate(scplane(amesh%Nnodes))
  allocate(origidnode(amesh%Nnodes))

  sface = infinity
  shead = infinity
  splane = infinity
  sedge = infinity
  scplane = infinity
  !------ end of debugging---

   do i = 1,3
      read(13,*)     ! skip header information on initial time
   enddo

   do i = 1,amesh%Nnodes
      read(13,*) itime
      if (itime > -1.0) then
        time(i) = itime
        kappa(i) = 0._pr
        mode(i) = 0
      endif
   enddo


   if(command_argument_count().ne.1) then
     write(*,*)'Assuming no diffraction taking place'
     adiff%fast = .true.  ! turn off diffraction
   else
    call get_command_argument(1,char2int)
    read(char2int,*) diff_level
    if (diff_level== 0) then
        write(*,*) 'No diffraction points considered'
        adiff%fast = .true.  ! turn off diffraction
    elseif (diff_level== 1) then
        read(13,*)
        read(13,*)
        read(13,*)

        write(*,*) 'Defined diffraction points considered'
        adiff%fast = .false.  ! turn on diffraction
        allocate(list_diff(amesh%Nnodes))
        list_diff = 0
        counter = 0
        do i = 1,amesh%Nnodes
          read(13,*) diff_switch
          if (diff_switch > 0) then
              counter=counter+1
              list_diff(counter) = i
          endif
        enddo

        if (counter == 0) then
          write(*,*) 'Error: no diffraction points detected'
          write(*,*) 'Will run assuming no diffraction'
          adiff%fast = .true.  ! turn on diffraction
        endif
        adiff%Nnodes = counter
        allocate(adiff%nodes(adiff%Nnodes))
        adiff%nodes = list_diff(1:counter)
        deallocate(list_diff)
        write(*,*) 'no. of diffraction nodes:     ', adiff%Nnodes
        !HERE we need to
        !(1) read in the number of nodes that are nonzero
        !(2) test if number of nodes ==0  => set adiff%fast = .true. and flag it
        !(3) create an array with all non-zero diffraction nodes
    elseif (diff_level == 2) then
           write(*,*) 'All mesh points considered diffraction points'
           adiff%fast = .false.  ! turn on diffraction
           adiff%Nnodes = amesh%Nnodes
           allocate(adiff%nodes(adiff%Nnodes))
           do i = 1,adiff%Nnodes
             adiff%nodes(i) = i
           enddo
    endif
   endif
   !
   ! adiff%fast = .true.  ! turn off diffraction
   ! if (.not.adiff%fast) then
   !   write(*,*)'number of diffraction nodes.......',ibc
   !   adiff%Nnodes = ibc
   !   allocate(adiff%nodes(adiff%Nnodes))
   !   adiff%nodes = bc_nodes(1:ibc)
   ! else
   !   adiff%Nnodes = 1
   !   allocate(adiff%nodes(adiff%Nnodes))
   !   adiff%nodes = 1
   ! endif

  close(13)

return
end subroutine read_vtk_mesh_jupv2
!###############################################################################
!###############################################################################
subroutine read_vtk_mesh_jup(amesh,cell_id,mesh_file)


  type(mesh) :: amesh
  character(30),intent(in) :: mesh_file
  integer(pin),allocatable,dimension(:) :: layer_bc,cell_id
  integer(pin) :: i, id,nsize,c1,c2,c3,element_type
  integer(pin) :: max_cells,t_cells,t_start,nnodes
  integer(pin) :: ipos_start,ipos_end
  logical :: tri_cells_start
  character(len=60) :: line

  open(13,file=trim(adjustl(mesh_file)),form = 'formatted', action = 'read')
  do i = 1,4
         read(13,*)
  enddo
  read(13,'(a)') line
  ipos_start = scan(line,"S",back=.true.)
  ipos_end = scan(line,"d",back=.false.)
!  write(*,*) ipos_start,ipos_end
  read (line(1+ipos_start:ipos_end),*) nnodes
!  amesh%Nnodes = nnodes-1  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
  amesh%Nnodes = nnodes  !
  write(*,*) 'No. of nodes....',amesh%Nnodes
  allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
  amesh%px = 0._pr
  amesh%py = 0._pr
  amesh%pz = 0._pr
!  read(13,*)  !NOTE: REMOVING NODE AT ZERO AS IT IS NOT USED
  do i = 1,amesh%Nnodes
     read(13,*) amesh%px(i),amesh%py(i),amesh%pz(i)
!     print*,amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  read(13,*)
  read(13,'(a)')  line
!  write(*,'(a)') line
  ipos_start = scan(line,"CELLS",back=.true.)
!  print*,'ipos ........',ipos_start
  read (line(1+ipos_start:),*) max_cells,nsize
!  read(13,'(5x,i6,i7)')amesh%Ncells,nsize
!  write(*,*) max_cells,nsize
!  write(*,*)'no. of cells',max_cells
!  stop
  t_cells = 0
  tri_cells_start = .false.
  do i = 1,max_cells
      read(13,*) id
      if (id == 3) then
         if (.not.tri_cells_start) then
           t_start = i
            tri_cells_start = .true.
         endif
         t_cells = t_cells+1
       endif
  enddo

!  write(*,*) 'first triangle at.....',t_start
!  write(*,*)'no. of triangular cells.....',t_cells
  rewind(13)
  do i = 1,6+amesh%Nnodes+t_start
!  do i = 1,6+amesh%Nnodes+t_start+1 !add one for the node you jumped
      read(13,*)
  enddo
!  write(*,*) 'current position......',6+amesh%Nnodes+t_start
  amesh%Ncells = t_cells
  allocate(amesh%cell(amesh%Ncells,3))
  write(*,*)  'no. of cells........ ',amesh%Ncells
  do i = 1,amesh%Ncells
   read(13,*) id,c1,c2,c3
   amesh%cell(i,:) = (/ c1+1, c2+1,c3+1 /)
!   amesh%cell(i,:) = (/ c1, c2,c3 /)  !this is because we have removed the first node at the origin
  enddo

  read(13,*)     ! blank line
  read(13,*)      ! skip line : CELL_TYPES xxxx
  do i = 1,amesh%Ncells     ! skip reading in cell_types
        read(13,*)
  enddo
  read(13,*)     ! blank line

  do i = 1,3
    read(13,*)     ! skip header informatin on CellEntityIds
  enddo
  allocate(cell_id(amesh%Ncells))
  do i = 1,amesh%Ncells
    read(13,*) cell_id(i)
  enddo


  close(13)

return
end subroutine read_vtk_mesh_jup
!###############################################################################
subroutine dumpmeshvtk_jup(dev,amesh)

  type(mesh) :: amesh
  integer(pin) :: dev

  integer(pin) :: i,j

  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a8)') 'distance'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a25:)') 'DATASET UNSTRUCTURED_GRID'
  if( pr == 4) then
     write(dev,'(a7,i10,a6)') 'POINTS ',amesh%Nnodes,' float'
  else
     write(dev,'(a7,i10,a7)') 'POINTS ',amesh%Nnodes,' double'
  endif

  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a6,2i10)') 'CELLS ',amesh%Ncells,amesh%Ncells*4
  do i=1,amesh%Ncells
     write(dev,'(a2,$)') '3 '
     write(dev,*) (amesh%cell(i,j)-1,j=1,3)
  enddo

  write(dev,*)
  write(dev,'(a11,i10)') 'CELL_TYPES ',amesh%Ncells
  do i = 1,amesh%Ncells
       write(dev,*) 5
  enddo

  return
end subroutine dumpmeshvtk_jup
!###############################################################################!###############################################################################
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
subroutine idumpnodeattributevtk(dev,amesh,field,attname,init)
  type(mesh) :: amesh
  integer(pin), dimension(amesh%Nnodes) :: field
  integer(pin) :: dev
  character*(*) :: attname
  logical :: init

  integer(pin) :: i

  ! if (init) write(dev,'(a11,i10)') 'POINT_DATA ',amesh%Nnodes
  ! if( pr == 4) then
    write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' int 1'
  ! else
  !    write(dev,'(a8,a,a10)') 'SCALARS ',trim(attname),' double 1'
  ! endif
  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,amesh%Nnodes
     write(dev,*) field(i)
  enddo
  return
end subroutine idumpnodeattributevtk
!###############################################################################
subroutine rdumpnodeattributevtk(dev,amesh,field,attname,init)
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
end subroutine rdumpnodeattributevtk
!###############################################################################
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
!###############################################################################
end module LAT_mesh_util
