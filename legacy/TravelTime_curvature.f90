!########################################################
subroutine timeonevsall2dV2(amesh,adiff)
    use LAT_mesh
    use LAT_time
    use lists
    use distance
  ! given a mesh (amesh) and velocities (velocity) defined on
  ! the mesh cells, timeonevsall2d computes the first arrival
  ! time on the mesh (time) from all the nodes in time array which
  ! have a finite time (i.e. not infinite)
    type(diff) :: adiff
    type(mesh) :: amesh
    type(container), dimension(amesh%Nnodes) :: ntoc
    type(listn), pointer :: ntodo,cellcur,pcur,last,pmin
    real(pr) :: d13,d23,d12,d14,d24,t1,t2,t4,tface,tedge,ttest,tref,v1,v2,a,thead
    integer(pin), dimension(2) :: otherninc
    integer(pin) :: i,j,toggle,c4,p4,nref,nbest
    integer(pin) :: idiff,waitfordiff
    logical :: begin,hasbeenupdated,reloop,firstrun,diffoccur
    logical, dimension(amesh%Nnodes) :: inthelist
    integer(pin) :: nswp,mxswp
    integer(pin) :: opposite_cell
    real(pr),parameter :: r2d = 180._pr/acos(-1._pr), pi = acos(-1._pr)
  ! checksecondary : .true. = cell has been tested ; .false. = point still needs to be checked
  ! therefore set nodes on boundary as = .false. with all else positive
    real(pr), dimension(amesh%Nnodes) :: secondary_time
    real(pr) :: kedge,kface
    integer(pin) :: medge, mface
    logical, dimension(amesh%Nnodes) :: checksecondary
    logical :: logic_contrast,logic_diff,logic_conic
    integer(pin), parameter :: waitdiffthres=25
  
  !====================================
  ! computing the node-to-cell array ntoc
    call compntoc(amesh,ntoc)
  ! initialize the secondary distance array
  !   secondary_time = time   ! distarray = secondary sources
    checksecondary=.true.   !assume all points have been checked
  !  checksecondary=.false.
    inthelist=.false.
  ! initialization of ntodo from non infinite values inside distarray
    begin=.true.
    ! designate nodal points that will be checked as secondary sources
  !   if (.not.adiff%fast) then
  !      checksecondary(adiff%nodes)=.false.
  !   endif
    do i=1,amesh%Nnodes
     if (time(i) > .9_pr*infinity) cycle
     ! if (secondary_time(i) > .9_pr*infinity) cycle
  ! if a vertex has a null distance, it is not a secondary diffraction point
  ! add list of nodes to check not whole array
  !   if (.not.adiff%fast) checksecondary(i)=.true.   ! <sounds better>
  !     if (secondary_time(i) < epsilon(secondary_time(i))) checksecondary(i)=.true.
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
  ! initializing reloop to true
  ! The condition is changed or not at the end of the loop
  ! The fast first loop is mandatory with a null distance offset
  !   reloop=.true.
  !   firstrun=.true.
  !   diffoccur=.false.
  !    do while (reloop)
  ! loop on the to do list
     if (verbose) call printlist(ntodo)
     if (verbose) write(*,*) 'entering in the ntodo list management loop'
     do while (associated(ntodo))
  ! searching the cells attached to the first element of the to do list
  !       if (verbose) write(*,*) 'state of ntodo :'
  !       if (verbose) call printlist(ntodo)
  ! searching the minimum time in the ntodo list
        call lookformin(ntodo,time,amesh%Nnodes,pmin)
        ! call lookformin(ntodo,secondary_time,amesh%Nnodes,pmin)
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
          toggle=3-toggle   !    toggle variable swap state between 1 and 2 at each sweep
  !        if (verbose) then
  !          if (toggle==1) write(*,*) 'sweep forward on cellist'
  !          if (toggle==2) write(*,*) 'sweep backward on cellist'
  !        endif
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
  !    check for velocity contrast between cells where travel time is being calculated 
             logic_contrast = .false.
             logic_diff = .false. 
             logic_conic = .false.
  !           opposite_cell=find_face_cell(pcur%idnode,otherninc(1),otherninc(2),amesh%Nedges,amesh%Edge)
             opposite_cell=find_face_cell_v2(pcur%idnode,otherninc(1),otherninc(2),amesh%nNodes,amesh%FVsNodes)
             if (opposite_cell == 0) then  ! there is no opposite cell (i.e. an edge)
                 logic_contrast=.false.
             elseif(abs(velocity(pcur%idnode)-velocity(opposite_cell)) <= epsilon(1.)) then   ! there is no velocity contrast 
                 logic_contrast=.false.
             else
                 logic_contrast=.true.
                 if (velocity(opposite_cell) > velocity(pcur%idnode)) then 
                       ! apply conic 
                       logic_conic = .true.
  !                     write(*,*) 'applying conic', pcur%idnode
                 else
                       ! apply diffraction 
                       logic_diff = .true. 
                       ! write(*,*) 'applying diffraction', pcur%idnode
                 endif
              endif
  !     loop on the two complemantory nodes
              do i=1,2
  !              if (verbose) write(*,'(a4,3(i4.4))') 'code',pmin%idnode,pcur%idnode,otherninc(i)
  !              if (verbose) write(*,*) '       working on complemantory nodes :',otherninc(i)
  !     edge propagation
                 !  tedge=secondary_time(pmin%idnode)+tonedge(amesh,pmin%idnode,otherninc(i),velocity(pcur%idnode))
                 call calc_tonedge(amesh,pmin%idnode,otherninc(i),velocity(pcur%idnode),kedge,tedge,medge) 
                 ! tedge=secondary_time(pmin%idnode)+tedge
  
  !               if (verbose) write(*,*) 'tedge : ',time(pmin%idnode),'+',&
  !                  tonedge(amesh,pmin%idnode,otherninc(i),velocity(pcur%idnode))
  !               if (verbose) write(*,*) 'tedge : ',tedge
  !     face propagation when time on 2 nodes are available
                 if (time(otherninc(3-i)) /= infinity) then
  !                  if (secondary_time(otherninc(3-i)) /= infinity) then
  !               if (verbose) write(*,*) 'computing tface because ',time(otherninc(3-i)),&
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
                    ! t1=secondary_time(pmin%idnode)
                    ! t2=secondary_time(otherninc(3-i))
                    t1=time(pmin%idnode)
                    t2=time(otherninc(3-i))
  !                 if (verbose) write(*,*) 't1,t2 : ',t1,t2
  !     multi lateration kernel
                    thead = infinity 
                    tface = infinity
  !                   if  (logic_conic)  then ! use plane wave approximation at slower side of velocity interface 
  !                        thead=tplane(d12,d13,d23,t1,t2,velocity(pcur%idnode))
  ! !                   ! !    tface=tcircle_conic(d12,d13,d23,t1,t2,velocity(pcur%idnode))
  ! !                   ! ! elseif (logic_diff) then 
  ! ! !                   if (logic_diff) then 
  ! ! !                            tface=tcircle_diff(d12,d13,d23,t1,t2,velocity(pcur%idnode))
  ! ! !                   ! ! !     !tface=infinity
  !                   else
                       !  thead=tplane(d12,d13,d23,t1,t2,velocity(pcur%idnode))
                       !   tface=tcircle(d12,d13,d23,t1,t2,velocity(pcur%idnode))
                        if ( (kappa(pmin%idnode) /= infinity).and.(kappa(otherninc(3-i))/= infinity) ) then 
                         call calc_tcircle(d12,d13,d23,t1,t2,                   &
                                kappa(pmin%idnode),kappa(otherninc(3-i)),     &
                                velocity(pcur%idnode),kface,tface,mface)
                        else 
                             !  thead=tplane(d12,d13,d23,t1,t2,velocity(pcur%idnode))
                              thead=infinity
                        endif 
                    ! endif
  !                 if (verbose) write(*,*) 'tface : ',tface
                 else
                    thead=infinity 
                    tface=infinity
                 endif
  
                 ttest=min(tedge,thead,tface)
  !     better result found - updating process
  !              if (ttest < time(otherninc(i))) then
                 if (ttest < time(otherninc(i))) then
  !                  if (ttest < secondary_time(otherninc(i))) then
  !                 if (verbose) write(*,*) 'better !'
                    hasbeenupdated=.true.
  !     time is better : if not in the list, add it
                    time(otherninc(i))=ttest
                    ! secondary_time(otherninc(i))=ttest
                    if (abs(ttest-tedge) <= epsilon(ttest)) then 
                             mode(otherninc(i)) = medge 
                             kappa(otherninc(i)) = kedge 
                    elseif(abs(ttest-tface) <= epsilon(ttest)) then
                                mode(otherninc(i)) = mface
                                kappa(otherninc(i)) = kface
                    elseif(abs(ttest-thead) <= epsilon(ttest)) then 
                          mode(otherninc(i)) = 2
                          kappa(otherninc(i)) = infinity
                    endif
  
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
  ! diff case and no diffraction has been detected yet
  !         if (.not.firstrun.and..not.diffoccur) then
  !            if (secondary_time(pmin%idnode)+time(idiff) < time(pmin%idnode)) then
  ! !   diffraction detected
  !               diffoccur=.true.
  !            else
  ! !   no diffraction
  !               waitfordiff=waitfordiff-1   ! check to see if there is an improvement in cells near diff point
  !               if (waitfordiff == 0) then
  ! !   countdown expired for diffraction. hard exit on ntodo list to pass to another
  ! !                                      potential secondary source
  !                  call wipe(ntodo)
  !                  exit
  !               endif
  !            endif
  !         endif
  ! remove pmin element from the to do list
          inthelist(pmin%idnode)=.false.
          call removepminfromthelist(pmin,ntodo)
        enddo
  ! decision to reloop or not
  ! fast case
  !       if (adiff%fast) then
  !             time = secondary_time
  !             write(*,*) 'exiting based on fast protocol ....'
  ! ! first case of general exit : only one run is requested
  !             reloop=.false.
  !             cycle
  !       endif
  !       if (firstrun) then
  !          time = secondary_time
  !       else
  !          if (diffoccur) then
  !       ! a secondary diff has occured
  !                do i=1,amesh%Nnodes
  !                      time(i)=min(secondary_time(i)+time(idiff),time(i))
  !                enddo
  !          endif
  !       endif
  !       idiff=nextdiffid(time,amesh%Nnodes,checksecondary)
  ! ! second case of general exit : all the potential secondary source have been investigated.
  !       if (idiff == 0) then
  !          reloop=.false.
  !          cycle
  !       endif
  ! ! preparation for the next loop
  !       checksecondary(idiff)=.true.
  !       waitfordiff=waitdiffthres    ! reset counter for diffraction
  !       diffoccur=.false.            ! rest diffraction check
  ! ! At this point notodo is deallocated. reallocate it with idiff
  !       allocate(ntodo)
  !       ntodo%previous=>ntodo
  !       nullify(ntodo%next)
  !       ntodo%idnode=idiff
  !       inthelist=.false.
  !       inthelist(idiff)=.true.
  !    ! reinitialize distarray
  !       secondary_time=infinity
  !       secondary_time(idiff)=0._pr
  !       firstrun=.false.
  ! enddo ! end reloop section
  
  ! deallocating ntoc
    call deallocntoc(ntoc,amesh%Nnodes)
  !  write(0,*) 'maximun sweep number : ',mxswp
  !  write(*,*) 'nref,nbest',nref,nbest
  end subroutine timeonevsall2dV2
!###############################################################################
subroutine calc_tonedge(amesh,i,j,v,lkappa,ltedge,lmode)
    use LAT_mesh
    use LAT_time
 
    type(mesh) :: amesh
    integer(pin) :: i,j,lmode
    real(pr) :: ltedge,v,d13,lkappa
 
    d13=sqrt((amesh%px(i)-amesh%px(j))**2._pr+&
             (amesh%py(i)-amesh%py(j))**2._pr+&
             (amesh%pz(i)-amesh%pz(j))**2._pr)
    if (mode(i) == 0) then  ! starting a diffraction
       ltedge = time(i)+d13/v
       lmode = 1
       lkappa = d13
    elseif (mode(i) == 1) then ! continuing a diffraction 
       ltedge = time(i)+d13/v
       lmode = 1
       lkappa = kappa(i)+d13
    endif   
    return   
 end subroutine calc_tonedge
!###############################################################################
subroutine calc_tcircle(d12,d13,d23,t1,t2,r1,r2,v,r3,t3,lmode)
    ! computes the coordinates in 2D plane of three vertices V1,V2,V3 making the
    ! assumptions that V1 is the origin (0,0), V2 is on the x axis (0,d12). The
    ! first part gives the coordinates x,y of V3.
    ! Same set of equations are used to estimate the origin of a point distant by
    ! r1 and r2 from V1 and V2 respectively (xc,+/-yc). The function returns the
    ! distance between (xc,-yc) and V0 (x,y)
    ! The position of V3 and C (the virtual origin) in the local reference system is returned
    ! with po1 and po2 respectively
    
      real(pr) :: t3,d12,d13,d23,t1,t2,v
      real(pr) :: r3,t31,t32,d
      real(pr) :: x,y,xc,yc,a,r1,r2,yc_sqrt
      integer(pin) :: lmode
      type(localxy) :: po1,po2
    
      x=(d12**2._pr-d23**2._pr+d13**2._pr)/(2._pr*d12)
      y=sqrt(max(d13**2._pr-x**2._pr,0._pr))
      po1%x=x
      po1%y=y
    ! if (verbose) write(*,*) 'tcircle - x,y :',x,y
    !   r1=t1*v
    !   r2=t2*v
      xc=(d12**2._pr-r2**2._pr+r1**2._pr)/(2._pr*d12)
      yc_sqrt = r1**2._pr-xc**2._pr
      if (yc_sqrt < 0._pr) then 
         t3 = infinity
         return
      endif
      yc = sqrt(yc_sqrt)
    !   yc=sqrt(max(r1**2._pr-xc**2._pr,0._pr))
      po2%x = xc
      po2%y = -yc
    ! if (verbose) write(*,*) 'tcircle - xc,yc :',xc,yc
      a=x_entering(po1,po2)
      if (a < 0._pr .or. a > d12) then
       ! tcircle=infinity
         t3=infinity
    !    if (verbose) write(*,*) 'dcircle - a : ',a
    !    if (verbose) write(*,*) 'dcircle - a outside range'
         return
      endif
    !   tcircle=sqrt((x-xc)**2._pr+(y+yc)**2._pr)/v  ! remember ... distance to (xc,-yc)
      d=sqrt((x-xc)**2._pr+(y+yc)**2._pr)  ! remember ... distance to (xc,-yc)
 
      t31 = t1+(d-r1)/v
      t32 = t2+(d-r2)/v
      r3 = d
      lmode = 0 
      if (t31 <= t32) then
        t3 = t31 
       else
          t3 = t32
      endif
 
      return
    end subroutine calc_tcircle
 !###############################################################################