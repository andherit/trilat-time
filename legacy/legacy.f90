module legacy

contains

!###################################################
subroutine allvsall2d(amesh,distarray)

  type(mesh) :: amesh
  real(pr), dimension(:,:), allocatable :: distarray

  type(container), dimension(amesh%Nnodes) :: ntoc
  type(listn), pointer :: ntodo,cellcur,pcur,last
  real(pr) :: d01,d02,d12,r1,r2,dface,dedge,dtest
  integer(pin), dimension(2) :: otherninc
  integer(pin) :: i,j,k,toggle
  logical :: begin,hasbeenupdated

! intialisation of dist array
!     allocation
  allocate(distarray(amesh%Nnodes,amesh%Nnodes))
!     set to infinity
  distarray=infinity
!     trace to zero
  do i=1,amesh%Nnodes
     distarray(i,i)=0.
  enddo
! computing distance to neighbor nodes cell by cell
  do i=1,amesh%Ncells
!     distance 1-2
      call updclosedist(amesh,distarray,i,1,2)
!     distance 2-3
      call updclosedist(amesh,distarray,i,2,3)
!     distance 3-1
      call updclosedist(amesh,distarray,i,3,1)
  enddo
! computing the node-to-cell array ntoc
  call compntoc(amesh,ntoc)
! main loop on nodes - nodes is a starting point
  if (verbose) write(*,*) 'starting main loop on nodes'
  do k=1,amesh%Nnodes
  if (verbose) write(*,*) '##############################################################'
  if (verbose) write(*,*) 'distance from node #',k,'/',amesh%Nnodes
     call printperc(k,amesh%Nnodes)
! initializing  node todo list
     if (verbose) write(*,*) 'entering vicinode'
     call vicinode(amesh,k,ntoc(k)%ptr,ntodo)
! loop on the to do list
     if (verbose) call printlist(ntodo)
     if (verbose) write(*,*) 'entering in the ntodo list management loop'
     do while (associated(ntodo))
! searching the cells attached to the first element of the to do list
        if (verbose) write(*,*) 'state of ntodo :'
        if (verbose) call printlist(ntodo)
        if (verbose) write(*,*) 'propagating distance from node #',ntodo%idnode
        if (verbose) write(*,*) 'its own distance is :',distarray(k,ntodo%idnode)
        hasbeenupdated=.true.
        toggle=2
        do while(hasbeenupdated)
        hasbeenupdated=.false.
! switch the toggle
        toggle=3-toggle
        if (verbose) then
        if (toggle==1) write(*,*) 'seep forward on cellist'
        if (toggle==2) write(*,*) 'seep backward on cellist'
        endif
        if (toggle==1) then
           pcur=>ntoc(ntodo%idnode)%ptr
        else
           pcur=>last
        endif
        do while (associated(pcur))
!     find the two complemantory nodes in the current cell
            if (verbose) write(*,*) '    working on cell #',pcur%idnode
            call givencomp(amesh,otherninc,pcur,ntodo)
            if (verbose) write(*,*) '    complemantory nodes : ',otherninc
!     loop on the two complemantory nodes
            do i=1,2
               if (verbose) write(*,'(a4,4(i2.2))') 'code',k,ntodo%idnode,pcur%idnode,otherninc(i)
               if (verbose) write(*,*) '       working on complemantory nodes :',otherninc(i)
!      symmetry check
                if (verbose) write(*,*) '       symmetry check:'
                if (verbose) write(*,*) otherninc(i),'<',k,' ???'
                if (otherninc(i) < k) then
                   if (verbose) write(*,*) '       yes'
                   if (verbose) write(*,*) 'is distarray(k,otherninc(i)) > distarray(otherninc(i),k) ?'
                   if (verbose) write(*,*) distarray(k,otherninc(i)),' > ',distarray(otherninc(i),k)
                   if (distarray(k,otherninc(i)) > distarray(otherninc(i),k)) then
                      if (verbose) write(*,*) '       yes then take the complimentary'
                      distarray(k,otherninc(i))=distarray(otherninc(i),k)
                      call updatelist(otherninc(i),ntodo)
                   endif
                   cycle
                endif
!     edge propagation
               dedge=distarray(k,ntodo%idnode)+distarray(ntodo%idnode,otherninc(i))
               if (verbose) write(*,*) 'dedge : ',distarray(k,ntodo%idnode),'+',distarray(ntodo%idnode,otherninc(i))
               if (verbose) write(*,*) 'dedge : ',dedge
!     face propagation
               if (distarray(k,otherninc(3-i)) /= infinity) then
               if (verbose) write(*,*) 'computing dface because ',distarray(k,otherninc(3-i)),' is not infinity'
                  d12=distarray(ntodo%idnode,otherninc(3-i))
                  d01=distarray(otherninc(i),ntodo%idnode)
                  d02=distarray(otherninc(i),otherninc(3-i))
               if (verbose) write(*,*) 'd12,d01,d02 :'
               if (verbose) write(*,*) d12,d01,d02
                  r1=distarray(k,ntodo%idnode)
                  r2=distarray(k,otherninc(3-i))
               if (verbose) write(*,*) 'r1,r2 : ',r1,r2
                  dface=dcircle(d12,d01,d02,r1,r2)
               if (verbose) write(*,*) 'dface : ',dface
               else
                  dface=infinity
               endif
               dtest=min(dedge,dface)
!               dtest=dedge
               if (verbose) write(*,*) 'dtest=min(dedge,dface) : ',dtest
               if (verbose) write(*,*) 'dtest < distance between',k,' and ',otherninc(i)
               if (verbose) write(*,*) dtest,distarray(k,otherninc(i))
               if (dtest < distarray(k,otherninc(i))) then
                   if (verbose) write(*,*) 'better !'
! distance is better : if not in the list, add it
                  distarray(k,otherninc(i))=dtest
                  if (verbose) write(*,*) 'updatelist with :',otherninc(i)
                  call updatelist(otherninc(i),ntodo)
                  if (verbose) write(*,*) 'state of ntodo :'
                  if (verbose) call printlist(ntodo)
               endif
            enddo
! moving on the vicinity list
            if (toggle==1) then
               last=>pcur
               pcur=>pcur%next
            else
               pcur=>pcur%previous
            endif
        enddo
        enddo
! cancelling the cell vicinity list
! removing the first element of the to do list
        if (associated(ntodo%next)) then
           ntodo=>ntodo%next
           deallocate(ntodo%previous)
           nullify(ntodo%previous)
        else
           deallocate(ntodo)
           nullify(ntodo)
        endif
      enddo
  enddo
! deallocating ntoc
  call deallocntoc(ntoc,amesh%Nnodes)
end subroutine allvsall2d
!###################################################
!###############################################################################
subroutine timeonevsall2d(amesh,time,velocity)

! given a mesh (amesh) and velocities (velocity) defined on
! the mesh cells, timeonevsall2d computes the first arrival
! time on the mesh (time) from all the nodes in time array which
! have a finite time (i.e. not infinite)

  type(mesh) :: amesh
  real(pr), dimension(amesh%Nnodes) :: time
  real(pr), dimension(amesh%Ncells) :: velocity

  type(container), dimension(amesh%Nnodes) :: ntoc
  type(listn), pointer :: ntodo,cellcur,pcur,last,pmin
  real(pr) :: d13,d23,d12,d14,d24,t1,t2,t4,tface,tedge,ttest,tref,v1,v2,a
  integer(pin), dimension(2) :: otherninc
  integer(pin) :: i,j,toggle,c4,p4,nref,nbest
  logical :: begin,hasbeenupdated,logic_conic
  logical, dimension(amesh%Nnodes) :: inthelist
  integer(pin) :: nswp,mxswp
  type(localxy) :: ptA,ptB,ptC,ptP,ptN
  integer(pin) :: opposite_cell
  real(pr) :: v1ov2,cosPCB,cosPBM,tconic,dCP,dNP,dNC
  real(pr),parameter :: r2d = 180._pr/acos(-1._pr), pi = acos(-1._pr)
  real(pr),dimension(2) :: vec_a,vec_b
  real(pr) :: ang_PBM,ang_PCB

!  initialisation
   ptA%x=0._pr
   ptA%y=0._pr
   ptB%y=0._pr
   inthelist=.false.
!  mxswp=0
   nref=0
   nbest=0
! computing the node-to-cell array ntoc
  call compntoc(amesh,ntoc)
! if (verbose) write(*,*) '##############################################################'
! initializing  node todo list
  begin=.true.
  do i=1,amesh%Nnodes
     if (time(i) > .9_pr*infinity) cycle
     if (begin) then
        allocate(ntodo)
        nullify(ntodo%previous)
        nullify(ntodo%next)
        pcur=>ntodo
        begin=.false.
     else
        allocate(pcur%next)
        nullify(pcur%next%next)
        pcur%next%previous=>pcur
        pcur=>pcur%next
     endif
     pcur%idnode=i
     inthelist(i)=.true.
     ntodo%previous=>pcur
  enddo
! loop on the to do list
! if (verbose) call printlist(ntodo)
! if (verbose) write(*,*) 'entering in the ntodo list management loop'
  do while (associated(ntodo))
! searching the cells attached to the first element of the to do list
!       if (verbose) write(*,*) 'state of ntodo :'
!       if (verbose) call printlist(ntodo)
! searching the minimum time in the ntodo list
        call lookformin(ntodo,time,amesh%Nnodes,pmin)
!       if (verbose) write(*,*) 'propagating time from node #',pmin%idnode
!       if (verbose) write(*,*) 'its own time is :',time(pmin%idnode)
!    initilisation of the sweep loop
        hasbeenupdated=.true.
        toggle=2
!       nswp=0
!    sweep loop up to the stabilisation (hasbeenupdated=.false.)
        do while (hasbeenupdated)
!       nswp=nswp+1
        hasbeenupdated=.false.
!    toggle variable swap state between 1 and 2 at each sweep
        toggle=3-toggle
!       if (verbose) then
!       if (toggle==1) write(*,*) 'sweep forward on cellist'
!       if (toggle==2) write(*,*) 'sweep backward on cellist'
!       endif
!    pointer on the beginning or the end of the list of cells containing the node idnode
        if (toggle==1) then
           pcur=>ntoc(pmin%idnode)%ptr
        else
           pcur=>last
        endif
!    cell list sweeping
        do while (associated(pcur))
!     find the two complemantory nodes in the current cell
!           if (verbose) write(*,*) '    working on cell #',pcur%idnode
            call givencomp(amesh,otherninc,pcur,pmin)
!           if (verbose) write(*,*) '    complemantory nodes : ',otherninc
!========================== start of conic text <turned off >  ===========================
!     Is there a faster velocity close to the cell
!           opposite_cell=find_face_cell(pcur%idnode,otherninc(1),otherninc(2),amesh%Nedges,amesh%Edge)
!           if (opposite_cell == 0) then
              logic_conic=.false.
!           else
!              logic_conic=(velocity(opposite_cell) > velocity(pcur%idnode))
!           endif
!========================== end of conic text <turned off >  ===========================

!     loop on the two complemantory nodes
            do i=1,2
!              if (verbose) write(*,'(a4,3(i4.4))') 'code',pmin%idnode,pcur%idnode,otherninc(i)
!              if (verbose) write(*,*) '       working on complemantory nodes :',otherninc(i)
!     edge propagation
               tedge=time(pmin%idnode)+tonedge(amesh,pmin%idnode,otherninc(i),velocity(pcur%idnode))
!              if (verbose) write(*,*) 'tedge : ',time(pmin%idnode),'+',&
!                  tonedge(amesh,pmin%idnode,otherninc(i),velocity(pcur%idnode))
!              if (verbose) write(*,*) 'tedge : ',tedge
!     face propagation when time on 2 nodes are available
               if (time(otherninc(3-i)) /= infinity) then
!                 if (verbose) write(*,*) 'computing tface because ',time(otherninc(3-i)),&
!                    ' is not infinity'
!     edge lengths in the current cell
!                 if (verbose) then
!                    write(*,*) 'd12 distance between ',pmin%idnode,'and',otherninc(3-i)
!                    write(*,*) 'position of ',pmin%idnode,' : '
!                    write(*,*) amesh%px(pmin%idnode),amesh%py(pmin%idnode)
!                    write(*,*) 'position of ',otherninc(3-i),' : '
!                    write(*,*) amesh%px(otherninc(3-i)),amesh%py(otherninc(3-i))
!                 endif
                  d12=donedge(amesh,pmin%idnode,otherninc(3-i))
                  d13=donedge(amesh,otherninc(i),pmin%idnode)
                  d23=donedge(amesh,otherninc(i),otherninc(3-i))

!                 if (verbose) write(*,*) 'd12,d13,d23 :',d12,d13,d23
!     time associated to the 2 nodes for the lateration computation
                  t1=time(pmin%idnode)
                  t2=time(otherninc(3-i))
!                 if (verbose) write(*,*) 't1,t2 : ',t1,t2
!     multi lateration kernel
                  tface=tcircle(d12,d13,d23,t1,t2,velocity(pcur%idnode),ptC,ptP)
!                 if (verbose) write(*,*) 'tface : ',tface
               else
                  tface=infinity
               endif
! conic estimation
               tconic=infinity
               if (logic_conic) then
                  ptB%x = d12
                  v1ov2 = velocity(pcur%idnode)/velocity(opposite_cell)
                  cosPCB = (ptP%x-ptC%x)*(ptB%x-ptC%x)+(ptP%y-ptC%y)*(ptB%y-ptC%y)
                  cosPCB = cosPCB/sqrt((ptP%x-ptC%x)**2.+(ptP%y-ptC%y)**2.)
                  cosPCB = cosPCB/sqrt((ptB%x-ptC%x)**2.+ptC%y**2.)
                  ang_PCB = acos(cosPCB)
!                  ang_PCB = atan2(ptP%y-ptC%y, ptP%x-ptC%x) - atan2(ptB%y-ptC%y ,  ptB%x-ptC%x)
!                  vec_a = (/ ptP%x-ptC%x, ptP%y-ptC%y /)
!                  vec_b = (/ ptB%x-ptC%x, ptB%y-ptC%y /)
                  if (ang_PCB < acos(v1ov2)) then
!                  if (acos(cosPCB) < acos(v1ov2)) then
!                     cosPBM = ((ptP%x-ptB%x)*ptP%y-ptC%y*(ptB%x-ptC%x))         &
                     cosPBM = ((ptP%x-ptB%x)*(ptB%x-ptC%x)+(ptP%y-ptB%y)*(ptB%y-ptC%y))    &
                                 /sqrt((ptP%x-ptB%x)**2.+ptP%y**2.)                        &
                                 /sqrt((ptB%x-ptC%x)**2.+ptC%y**2.)
                     ang_PBM = acos(cosPBM)
!                     vec_a = (/ ptP%x-ptB%x, ptP%y-ptB%y /)
!                     vec_b = (/ ptB%x-ptC%x, ptB%y-ptC%y /)
!                     ang_PBM = atan2(ptP%y-ptB%y, ptP%x-ptB%x) - atan2(ptB%y-ptC%y ,  ptB%x-ptC%x)
                     if ( acos(v1ov2) < ang_PBM) then
!                        write(*,*) 'angles ....',acos(cosPCB)*r2d,acos(v1ov2)*r2d,acos(cosPBM)*r2d
!                     if ( v1ov2 < cosPBM) then
! A conic happens inside !
! computing the distance between ptC and ptP
                        dCP=sqrt((ptP%x-ptC%x)**2.+(ptP%y-ptC%y)**2.)
! computing the distance between ptP and the location of the critical point on BC (N)
                        dNP=dCP*sqrt((1._pr-cosPCB**2.)/(1._pr-v1ov2**2.))
! computing the distance between N and C, aka the path length in the opposite cell
                        dNC=dCP*cosPCB-dNP*v1ov2
! time through conic path
                        ptN%x=(d23-dNC)/d23*(ptC%x-ptB%x)+ptB%x
                        ptN%y=(d23-dNC)/d23*(ptC%y-ptB%y)+ptB%y
                        a=x_entering(ptN,ptP)
                        if ( a > 0._pr .and. a < d12) then
                           tconic=dNP/velocity(pcur%idnode)+dNC/velocity(opposite_cell)
!                           write(*,*) 'two contributions.....',dNP/velocity(pcur%idnode),dNC/velocity(opposite_cell)
!                           write(*,*) 'using tconic value in cell....',pcur%idnode
                        else
                           tconic=infinity
                        endif
! end conic estimation
                     endif
                  endif
               endif
!               ttest=min(tedge,tface)
               ttest=min(tedge,tface,tconic)

!     better result found - updating process
               if (ttest < time(otherninc(i))) then
!                 if (verbose) write(*,*) 'better !'
                  time(otherninc(i))=ttest
                  hasbeenupdated=.true.
!     time is better : if not in the list, add it
!                 if (verbose) write(*,*) 'updatelist with :',otherninc(i)
!                  call updatelist(otherninc(i),ntodo)
                  if (.not.inthelist(otherninc(i))) then
                      inthelist(otherninc(i))=.true.
                      allocate(ntodo%previous%next)
                      nullify(ntodo%previous%next%next)
                      ntodo%previous%next%previous=>ntodo%previous
                      ntodo%previous%next%idnode=otherninc(i)
                      ntodo%previous=>ntodo%previous%next
                  endif
!                 if (verbose) write(*,*) 'state of ntodo :'
!                 if (verbose) call printlist(ntodo)
               endif
            enddo
! moving on the vicinity list
            if (toggle==1) then
               last=>pcur
               pcur=>pcur%next
            else
               pcur=>pcur%previous
            endif
        enddo
        enddo
!        mxswp=max(mxswp,nswp)
! remove pmin element from the to do list
        inthelist(pmin%idnode)=.false.
        call removepminfromthelist(pmin,ntodo)
      enddo
! deallocating ntoc
  call deallocntoc(ntoc,amesh%Nnodes)
!  write(0,*) 'maximun sweep number : ',mxswp
!  write(*,*) 'nref,nbest',nref,nbest
end subroutine timeonevsall2d
!###############################################################################
subroutine timeonevsall2d_old(amesh,time,velocity)

   ! given a mesh (amesh) and velocities (velocity) defined on
   ! the mesh cells, timeonevsall2d computes the first arrival
   ! time on the mesh (time) from all the nodes in time array which
   ! have a finite time (i.e. not infinite)

     type(mesh) :: amesh
     real(pr), dimension(amesh%Nnodes) :: time
     real(pr), dimension(amesh%Ncells) :: velocity

     type(container), dimension(amesh%Nnodes) :: ntoc
     type(listn), pointer :: ntodo,cellcur,pcur,last,pmin
     real(pr) :: d13,d23,d12,d14,d24,t1,t2,t4,tface,tedge,ttest,tref,v1,v2

     integer(pin), dimension(2) :: otherninc
     integer(pin) :: i,j,toggle,c4,p4,nref,nbest
     logical :: begin,hasbeenupdated
     logical, dimension(amesh%Nnodes) :: inthelist
     integer(pin) :: nswp,mxswp

     type(localxy) :: ptA,ptB,ptC,ptP,ptN

   !  initialisation
      inthelist=.false.
   !  mxswp=0
      nref=0
      nbest=0
   ! computing the node-to-cell array ntoc
     call compntoc(amesh,ntoc)
   ! if (verbose) write(*,*) '##############################################################'
   ! initializing  node todo list
     begin=.true.
     do i=1,amesh%Nnodes
        if (time(i) > .9_pr*infinity) cycle
        if (begin) then
           allocate(ntodo)
           nullify(ntodo%previous)
           nullify(ntodo%next)
           pcur=>ntodo
           begin=.false.
        else
           allocate(pcur%next)
           nullify(pcur%next%next)
           pcur%next%previous=>pcur
           pcur=>pcur%next
        endif
        pcur%idnode=i
        inthelist(i)=.true.
        ntodo%previous=>pcur
     enddo
   ! loop on the to do list
   ! if (verbose) call printlist(ntodo)
   ! if (verbose) write(*,*) 'entering in the ntodo list management loop'
     do while (associated(ntodo))
   ! searching the cells attached to the first element of the to do list
   !       if (verbose) write(*,*) 'state of ntodo :'
   !       if (verbose) call printlist(ntodo)
   ! searching the minimum time in the ntodo list
           call lookformin(ntodo,time,amesh%Nnodes,pmin)
   !       if (verbose) write(*,*) 'propagating time from node #',pmin%idnode
   !       if (verbose) write(*,*) 'its own time is :',time(pmin%idnode)
   !    initilisation of the sweep loop
           hasbeenupdated=.true.
           toggle=2
   !       nswp=0
   !    sweep loop up to the stabilisation (hasbeenupdated=.false.)
           do while (hasbeenupdated)
   !       nswp=nswp+1
           hasbeenupdated=.false.
   !    toggle variable swap state between 1 and 2 at each sweep
           toggle=3-toggle
   !       if (verbose) then
   !       if (toggle==1) write(*,*) 'sweep forward on cellist'
   !       if (toggle==2) write(*,*) 'sweep backward on cellist'
   !       endif
   !    pointer on the beginning or the end of the list of cells containing the node idnode
           if (toggle==1) then
              pcur=>ntoc(pmin%idnode)%ptr
           else
              pcur=>last
           endif
   !    cell list sweeping
           do while (associated(pcur))
   !     find the two complemantory nodes in the current cell
   !           if (verbose) write(*,*) '    working on cell #',pcur%idnode
               call givencomp(amesh,otherninc,pcur,pmin)
   !           if (verbose) write(*,*) '    complemantory nodes : ',otherninc
   !     loop on the two complemantory nodes
               do i=1,2
   !              if (verbose) write(*,'(a4,3(i4.4))') 'code',pmin%idnode,pcur%idnode,otherninc(i)
   !              if (verbose) write(*,*) '       working on complemantory nodes :',otherninc(i)
   !     edge propagation
                  tedge=time(pmin%idnode)+tonedge(amesh,pmin%idnode,otherninc(i),velocity(pcur%idnode))
   !              if (verbose) write(*,*) 'tedge : ',time(pmin%idnode),'+',&
   !                  tonedge(amesh,pmin%idnode,otherninc(i),velocity(pcur%idnode))
   !              if (verbose) write(*,*) 'tedge : ',tedge
   !     face propagation when time on 2 nodes are available
                  if (time(otherninc(3-i)) /= infinity) then
   !                 if (verbose) write(*,*) 'computing tface because ',time(otherninc(3-i)),&
   !                    ' is not infinity'
   !     edge lengths in the current cell
   !                 if (verbose) then
   !                    write(*,*) 'd12 distance between ',pmin%idnode,'and',otherninc(3-i)
   !                    write(*,*) 'position of ',pmin%idnode,' : '
   !                    write(*,*) amesh%px(pmin%idnode),amesh%py(pmin%idnode)
   !                    write(*,*) 'position of ',otherninc(3-i),' : '
   !                    write(*,*) amesh%px(otherninc(3-i)),amesh%py(otherninc(3-i))
   !                 endif
                     d12=donedge(amesh,pmin%idnode,otherninc(3-i))
                     d13=donedge(amesh,otherninc(i),pmin%idnode)
                     d23=donedge(amesh,otherninc(i),otherninc(3-i))
   !                 if (verbose) write(*,*) 'd12,d13,d23 :',d12,d13,d23
   !     time associated to the 2 nodes for the lateration computation
                     t1=time(pmin%idnode)
                     t2=time(otherninc(3-i))
   !                 if (verbose) write(*,*) 't1,t2 : ',t1,t2
   !     multi lateration kernel
                     tface=tcircle(d12,d13,d23,t1,t2,velocity(pcur%idnode),ptC,ptP)

   !                 if (verbose) write(*,*) 'tface : ',tface
                  else
                     tface=infinity
                  endif
   !     result is the minimum between edge propagation and face propagation
   !              ttest=min(tedge,tface,tref)
                  ttest=min(tedge,tface)
                  if (abs(ttest-tref) < epsilon(tref)) nbest=nbest+1
   !              ttest=tedge
   !              if (verbose) write(*,*) 'ttest=min(tedge,tface) : ',ttest
   !              if (verbose) write(*,*) 'ttest < time between',pmin%idnode,' and ',otherninc(i)
   !              if (verbose) write(*,*) ttest,time(otherninc(i))
   !     better result found - updating process
                  if (ttest < time(otherninc(i))) then
   !                 if (verbose) write(*,*) 'better !'
                     time(otherninc(i))=ttest
                     hasbeenupdated=.true.
   !     time is better : if not in the list, add it
   !                 if (verbose) write(*,*) 'updatelist with :',otherninc(i)
   !                  call updatelist(otherninc(i),ntodo)
                     if (.not.inthelist(otherninc(i))) then
                         inthelist(otherninc(i))=.true.
                         allocate(ntodo%previous%next)
                         nullify(ntodo%previous%next%next)
                         ntodo%previous%next%previous=>ntodo%previous
                         ntodo%previous%next%idnode=otherninc(i)
                         ntodo%previous=>ntodo%previous%next
                     endif
   !                 if (verbose) write(*,*) 'state of ntodo :'
   !                 if (verbose) call printlist(ntodo)
                  endif
               enddo
   ! moving on the vicinity list
               if (toggle==1) then
                  last=>pcur
                  pcur=>pcur%next
               else
                  pcur=>pcur%previous
               endif
           enddo
           enddo
   !        mxswp=max(mxswp,nswp)
   ! remove pmin element from the to do list
           inthelist(pmin%idnode)=.false.
           call removepminfromthelist(pmin,ntodo)
         enddo
   ! deallocating ntoc
     call deallocntoc(ntoc,amesh%Nnodes)
   !  write(0,*) 'maximun sweep number : ',mxswp
     write(*,*) 'nref,nbest',nref,nbest

end subroutine timeonevsall2d_old
!###############################################################################
function find_face_cell(cell_no,n1,n2,n,Edge)
   integer(pin) ::find_face_cell,cell_no,n1,n2,n,i
   integer(pin) , dimension(n,4) :: Edge

   find_face_cell = -1
   do i = 1,n
       if(Edge(i,3)==n1.and.Edge(i,4)==n2.or.           &
          Edge(i,3)==n2.and.Edge(i,4)==n1) then
              if(Edge(i,1)==cell_no) then
                 find_face_cell = Edge(i,2)
              elseif(Edge(i,2)==cell_no) then
                 find_face_cell = Edge(i,1)
           endif
       endif
   enddo

   if (find_face_cell == -1) then
      write(0,*) 'ERROR in find_face_cell: no face found'
      stop
   endif
   return
end function find_face_cell
!###############################################################################
subroutine getfourth(amesh,clist,p1,p2,p4,c4,c1)

  type(mesh) :: amesh
  type(listn), pointer :: clist
  integer(pin) :: c4,p4,p1,p2,c1

  type(listn), pointer :: pcur
  integer(pin) :: i,score,idsum

  pcur=>clist
  p4=0
  do while (associated(pcur))
     if (pcur%idnode == c1) then
        pcur=>pcur%next
        cycle
     endif
     score=0
     idsum=0
     do i=1,3
        if (amesh%cell(pcur%idnode,i) == p1) score=score+1
        if (amesh%cell(pcur%idnode,i) == p2) score=score+1
        idsum=idsum+amesh%cell(pcur%idnode,i)
     enddo
     if (score /= 2) then
        pcur=>pcur%next
        cycle
     endif
     c4=pcur%idnode
     p4=idsum-p1-p2
     return
  enddo
  return
end subroutine getfourth
!###############################################################################
function trefract(d12,d13,d23,d14,d24,t4,v1,v2)

  real(pr) :: d12,d13,d23,d14,d24,t4,v1,v2,trefract

  real(pr) :: x4,y4,x3,y3,dl,tp1,tp1pdl,tp2,tp2mdl,x
  logical :: minotf
! nimp and infinity are global variables.

  dl=d12/float(nimp-1)
! 1st lateration : the location of v3
  x3=(d12**2._pr-d23**2._pr+d13**2._pr)/(2._pr*d12)
  y3=sqrt(max(d13**2._pr-x3**2._pr,0._pr))
! 2nd lateration : the location of v4
  x4=(d12**2._pr-d24**2._pr+d14**2._pr)/(2._pr*d12)
  y4=-sqrt(max(d14**2._pr-x4**2.,0._pr))
! a refraction exists in d12 if the time gradient along d12 is negative on v1 and
! positive on v2.
  tp1=d14/v2+d12/v1
  tp1pdl=sqrt((dl-x4)**2._pr+y4**2._pr)/v2+sqrt((dl-x3)**2._pr+y3**2._pr)/v1
  tp2=d24/v2+d23/v1
  tp2mdl=sqrt((d12-dl-x4)**2._pr+y4**2._pr)/v2+sqrt((d12-dl-x3)**2._pr+y3**2._pr)/v1
  if (.not.(tp1pdl < tp1 .and. tp2mdl < tp2)) then
     trefract=infinity
     return
  endif
! a minimum exists along d12
  minotf=.true.
  x=dl
  tp1=tp1pdl
  do while (minotf)
     x=x+dl
     tp2=sqrt((x-x4)**2._pr+y4**2._pr)/v2+sqrt((x-x3)**2._pr+y3**2._pr)/v1
     if (tp2 < tp1) then
        tp1=tp2
     else
        minotf=.false.
     endif
  enddo
  trefract=t4+tp1
  return
end function trefract
!###############################################################################
subroutine updclosedist(amesh,distarray,k,i,j)

   type(mesh) :: amesh
   real(pr), dimension(:,:), allocatable :: distarray
   integer(pin) :: i,j,k

   if (distarray(amesh%cell(k,i),amesh%cell(k,j)) == infinity) then
      distarray(amesh%cell(k,i),amesh%cell(k,j))=&
           sqrt((amesh%px(amesh%cell(k,i))-amesh%px(amesh%cell(k,j)))**2._pr+&
                (amesh%py(amesh%cell(k,i))-amesh%py(amesh%cell(k,j)))**2._pr+&
                (amesh%pz(amesh%cell(k,i))-amesh%pz(amesh%cell(k,j)))**2._pr)
! reciprocity
      distarray(amesh%cell(k,j),amesh%cell(k,i))=distarray(amesh%cell(k,i),amesh%cell(k,j))
   endif
end subroutine updclosedist
!###############################################################################
function dplane(d12,d13,d23,r1,r2)

  real(pr) :: dplane,d12,d13,d23,r1,r2
  real(pr) :: x,y,xc,yc
  real(pr) :: dr,xn,yn,dn,d13n,d23n

! position of v3
  x=(d12**2-d23**2._pr+d13**2._pr)/(2._pr*d12)
  y=sqrt(max(d13**2._pr-x**2._pr,0._pr))
! difference of distance at the base
  dr=abs(r1-r2)
! computation of the new base length
  dn=sqrt(d12**2._pr-dr**2._pr)
  if (r1 < r2) then
! case where v2 remains and V1 is changed
     xn=(d12**2._pr-dn**2._pr+dr**2._pr)/(2._pr*d12)
     yn=(sqrt(max(dr**2._pr-xn**2._pr,0._pr)))
     d13n=sqrt((x-xn)**2._pr+(y-yn)**2._pr)
! reposition of v3 in the new frame
     x=(dn**2._pr-d23**2._pr+d13n**2._pr)/(2._pr*dn)
     if (x > 0._pr .and. x < dn) then
        y=sqrt(max(d13n**2._pr-x**2._pr,0._pr))
        dplane=r2+y
     else
        dplane=infinity
     endif
  else
! case where v1 remains and V2 is changed
     xn=(d12**2._pr-dr**2._pr+dn**2._pr)/(2._pr*d12)
     yn=(sqrt(max(dn**2._pr-xn**2._pr,0._pr)))
     d23n=sqrt((x-xn)**2._pr+(y-yn)**2._pr)
! reposition of v3 in the new frame
     x=(dn**2._pr-d23n**2._pr+d13**2._pr)/(2._pr*dn)
     if (x > 0._pr .and. x < dn) then
        y=sqrt(max(d13**2._pr-x**2._pr,0._pr))
        dplane=r1+y
     else
        dplane=infinity
     endif
  endif
end function dplane
!###############################################################################
subroutine startdist(amesh,anode,nodelist,distarray)

! Compute the distance at the neightbor nodes list (nodelist) of a given node (anode)
! This routine is needed by onevsall routine after the initial vicinode call.
! The routine allvsall uses the updclosedist routine instead

  type(mesh) :: amesh
  type(listn), pointer :: nodelist
  integer(pin) :: anode
  real(pr), dimension(amesh%Nnodes) :: distarray

  type(listn), pointer :: pcur

  pcur=>nodelist
  do while (associated(pcur))
     distarray(pcur%idnode)=donedge(amesh,anode,pcur%idnode)
     pcur=>pcur%next
  enddo
  return
end subroutine startdist
!###############################################################################
subroutine startime(amesh,anode,celllist,nodelist,time,velocity)

  type(mesh) :: amesh
  type(listn), pointer :: nodelist,celllist
  integer(pin) :: anode
  real(pr), dimension(amesh%Nnodes) :: time
  real(pr), dimension(amesh%Ncells) :: velocity

  type(listn), pointer :: curcell,curnode
  logical :: begin
  integer(pin) :: i,j

!  write(*,*) 'startime : entering vicicell for node ',anode
   if (verbose) write(*,*) 'celllist%idnode',celllist%idnode
   begin=.true.
   curcell=>celllist
   do while (associated(curcell))
      do i=1,3
         if (amesh%cell(curcell%idnode,i).ne.anode) then
            if (begin) then
               allocate(nodelist)
               nullify(nodelist%previous)
               nullify(nodelist%next)
               curnode=>nodelist
               curnode%idnode=amesh%cell(curcell%idnode,i)
               begin=.false.
            else
               if (.not.isinit(amesh%cell(curcell%idnode,i),nodelist)) then
                  allocate(curnode%next)
                  curnode%next%previous=>curnode
                  nullify(curnode%next%next)
                  curnode=>curnode%next
                  curnode%idnode=amesh%cell(curcell%idnode,i)
               endif
           endif
           time(amesh%cell(curcell%idnode,i))=min(tonedge(amesh,anode,&
                amesh%cell(curcell%idnode,i),velocity(curcell%idnode)),&
                time(amesh%cell(curcell%idnode,i)))
        endif
     enddo
     curcell=>curcell%next
  enddo
end subroutine startime
!###############################################################################
subroutine vicinode(amesh,anode,celllist,nodelist)

! Compute the neightbor nodes list (nodelist) of a given node (anode)
! and its given neighbor cell list (celllist).

  type(mesh) :: amesh
  type(listn), pointer :: celllist,nodelist
  integer(pin) :: anode

  type(listn), pointer :: curcell,curnode
  logical :: begin
  integer(pin) :: i,j

  if (verbose) write(*,*) 'celllist%idnode',celllist%idnode
  begin=.true.
  curcell=>celllist
  do while (associated(curcell))
     do i=1,3
        if (amesh%cell(curcell%idnode,i).ne.anode) then
           if (begin) then
              allocate(nodelist)
              nullify(nodelist%previous)
              nullify(nodelist%next)
              curnode=>nodelist
              curnode%idnode=amesh%cell(curcell%idnode,i)
              begin=.false.
           else
              if (.not.isinit(amesh%cell(curcell%idnode,i),nodelist)) then
                 allocate(curnode%next)
                 curnode%next%previous=>curnode
                 nullify(curnode%next%next)
                 curnode=>curnode%next
                 curnode%idnode=amesh%cell(curcell%idnode,i)
              endif
           endif
        endif
     enddo
     curcell=>curcell%next
  enddo
end subroutine vicinode
!###############################################################################

end module legacy
