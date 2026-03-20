!###############################################################################
subroutine calc_tcplane(tri_ops,t3)
    type(ops) :: tri_ops
    real(pr) :: t3
 ! calculate the travel time to t3 using the plane wave approximation and 3 nodal points
 ! using two adject cells 
 ! convention is : first cell contains ending *a and second cell contains ending *b
 ! common values to both cells are : d13, t1     
    
 
    real(pr) :: v_mean,x2a,y2a,x2b,y2b
    real(pr) :: dx,dy,gamma,gamma_x,gamma_y ! ,beta
    real(pr) :: a,b,c,d
    real(pr) :: alpha,beta,vx_a,vx_b,x2,y2,adum
    real(pr), parameter :: pi=acos(-1._pr)
 
    real(pr) ::ax,ay,bx,by,m
 
 
    t3 = infinity
 
    x2a = (tri_ops%d13**2._pr-tri_ops%d23a**2._pr+tri_ops%d12a**2._pr)/(2._pr*tri_ops%d13)
    y2a = sqrt(max(tri_ops%d12a**2._pr-x2a**2._pr,0._pr))
 
    x2b = (tri_ops%d13**2._pr-tri_ops%d23b**2._pr+tri_ops%d12b**2._pr)/(2._pr*tri_ops%d13)
    y2b = -sqrt(max(tri_ops%d12b**2._pr-x2b**2._pr,0._pr))
 
    ax = abs(tri_ops%t2a-tri_ops%t1)/x2a
    ! ay = (tri_ops%t2a-tri_ops%t1)/y2a
 
    bx = abs(tri_ops%t2b-tri_ops%t1)/x2b
    ! by = (tri_ops%t2b-tri_ops%t1)/y2b
    ! m = (ay-by)/(ax-bx)
    ! write(*,*) 'Ray_a (deg),Ray_b(deg)', ay/ax*180._pr/pi,by/bx*180._pr/pi
     v_a = 2._pr/(ax+bx)
    
 ! !    ! x2 = (tri_ops%d12**2._pr-tri_ops%d23**2._pr+tri_ops%d13**2._pr)/(2._pr*tri_ops%d12)
 ! !    ! y2 = sqrt(max(tri_ops%d13**2._pr-tri_ops%x3**2._pr,0._pr))
 
 
    x2 = (tri_ops%d13**2._pr-tri_ops%d23a**2._pr+tri_ops%d12a**2._pr)/(2._pr*tri_ops%d13)
    y2 = sqrt(max(tri_ops%d12a**2._pr-x2**2._pr,0._pr))
 
    
   
    if (abs(tri_ops%t1-tri_ops%t2a) <= water_level(tri_ops%t2a)) then 
       beta = pi/2._pr 
    else  
       if (tri_ops%va > tri_ops%d12a/abs(tri_ops%t1-tri_ops%t2a) ) then 
          beta = 0._pr 
       else
          beta = acos(tri_ops%va/tri_ops%d12a*abs(tri_ops%t2a-tri_ops%t1))
       endif
    endif 
 
    alpha = acos(x2/tri_ops%d12a)
 
 
    adum=cos(alpha-beta)
    if (adum < epsilon(adum)) then
       ! t3 = tri_ops%t1
       return   
    else
       vx_a = tri_ops%va/cos(alpha-beta)
    endif
 
    write(*,*) 'v, v_xa, cos',tri_ops%va,vx_a,adum
 
    ! write(*,*) 'time',tri_ops%t1,tri_ops%t2a
    ! write(*,*) 'velocity', tri_ops%d12a/abs(tri_ops%t1-tri_ops%t2a),tri_ops%va
    ! write(*,*) 'alpha, beta',alpha,beta
 
 
 ! the other side 
    x2 = (tri_ops%d13**2._pr-tri_ops%d23b**2._pr+tri_ops%d12b**2._pr)/(2._pr*tri_ops%d13)
    y2 = sqrt(max(tri_ops%d12b**2._pr-x2**2._pr,0._pr))
 
    if (abs(tri_ops%t1-tri_ops%t2b) <= water_level(tri_ops%t2b)) then 
       beta = pi/2._pr 
    else  
       if (tri_ops%vb > tri_ops%d12b/abs(tri_ops%t1-tri_ops%t2b) ) then 
          beta = 0._pr 
       else
          beta = acos(tri_ops%vb/tri_ops%d12b*abs(tri_ops%t2b-tri_ops%t1))  
       endif
    endif 
 
 
    alpha = acos(x2/tri_ops%d12b)
 
    ! write(*,*) alpha, beta  
    adum=cos(alpha-beta)
    if (adum < epsilon(adum)) then
       ! t3 = tri_ops%t1
       return   
    else
       vx_b = tri_ops%vb/cos(alpha-beta)
    endif
 
    write(*,*) 'v, v_xb, cos',tri_ops%vb,vx_b,adum
 
    t3 = tri_ops%t1 + tri_ops%d13/(max(vx_a,vx_b)) 
 
    !  v_mean = (tri_ops%va+tri_ops%vb)/2._pr   !!!! check correct velocities !! 
 
    ! !  use trilateration to change coordinate system so that d13 is along x-axis 
    ! x2a = (tri_ops%d13**2._pr-tri_ops%d23a**2._pr+tri_ops%d12a**2._pr)/(2._pr*tri_ops%d13)
    ! y2a = sqrt(max(tri_ops%d12a**2._pr-x2a**2._pr,0._pr))
 
    ! x2b = (tri_ops%d13**2._pr-tri_ops%d23b**2._pr+tri_ops%d12b**2._pr)/(2._pr*tri_ops%d13)
    ! y2b = -sqrt(max(tri_ops%d12b**2._pr-x2b**2._pr,0._pr))
    ! ! write(*,*)  x2a, y2a, x2b, y2b
 
    ! alpha = (tri_ops%t2a - tri_ops%t1*(1._pr - x2a/tri_ops%d13) ) /y2a                   &
    !       + ( tri_ops%t1*(1._pr - x2b/tri_ops%d13) -tri_ops%t2b ) /y2b 
    
    ! ! beta = x2b/y2b/tri_ops%d13 - x2a/y2a/tri_ops%d13
    ! beta = x2a/y2a/tri_ops%d13 - x2b/y2b/tri_ops%d13
 
    ! a = 1._pr+(tri_ops%d13*beta/2._pr)**2._pr
    ! b = (alpha*beta*tri_ops%d13**2)/2._pr - 2._pr*tri_ops%t1 ! *tri_ops%d13
    ! c = tri_ops%t1**2._pr + (tri_ops%d13*alpha/2._pr)**2._pr -  (tri_ops%d13/v_mean)**2._pr
    
    ! d = b**2._pr - 4._pr*a*c
    ! if (d >= epsilon(d)) t3 = (-b+sqrt(d))/(2._pr*a)
 !    ! set a,b,c depending on situation 
 !     if ((abs(x2a) <= water_level(x2a)).and.(abs(x2b) <= water_level(x2b)) ) then 
 ! ! .t2a
 ! ! |\
 ! ! |  \
 ! ! |    \ 
 ! ! .t1_ _ _ .t3
 ! ! |    /   
 ! ! |  /
 ! ! |/
 ! ! . t2b
 !       ! write(*,*) 'option 1'
 
 !       a = 1._pr
 !       b = -2._pr*tri_ops%t1
 !       c = tri_ops%t1**2._pr+(tri_ops%d13/(y2a-y2b)*(tri_ops%t2a-tri_ops%t2b))**2._pr-(tri_ops%d13/v_mean)**2._pr
 
 !       elseif ( (abs(tri_ops%d13-x2a)<= water_level(x2a) ).and.( abs(tri_ops%d13-x2b)<= water_level(x2b) )) then  !! talk to Andre about conditions 
 ! !          .t2a
 ! !       /  |
 ! !     /    |
 ! !   /      |
 ! ! .t1_ _ _ .t3
 ! !   \      |
 ! !     \    |
 ! !       \  |
 ! !          . t2b   
 !     ! write(*,*) 'option 2b'
 !       a = 1._pr
 !       b = -2._pr*tri_ops%t1
 !       c = tri_ops%t1**2._pr+(tri_ops%d13/(y2a-y2b)*(tri_ops%t2a-tri_ops%t2b))**2._pr-(tri_ops%d13/v_mean)**2._pr
 
 
 !     elseif (abs(x2b)+abs(tri_ops%d13-x2a) >= abs(x2a)+abs(tri_ops%d13-x2b) ) then 
 ! !         .t2a
 ! !       /   \
 ! !     /      \
 ! !   /          \ 
 ! ! .t1_ _ _ _ _ _.t3
 ! !   \         /
 ! !     \     /
 ! !       . /
 ! !       t2b 
 
 !          dx = tri_ops%d13-x2a
 !          dy = y2b
 !          gamma_x = (tri_ops%t2b-tri_ops%t1)/x2b  - tri_ops%t2a/dx  
 !          gamma_y = (tri_ops%t2a-tri_ops%t1)/y2a + tri_ops%t2b/dy
 !          ! write(*,*) 'option 3'
 !          a = dx**2 + dy**2
 !          b = -2._pr*dx*dy*(dy*gamma_x+gamma_y*dx)
 !          c = (gamma_x**2._pr+gamma_y**2._pr-(2._pr/v_mean)**2._pr)*(dx*dy)**2._pr
 !     else 
 ! !       .t2a
 ! !     /     \
 ! !   /          \ 
 ! ! .t1_ _ _ _ _ _.t3
 ! !   \          /
 ! !     \       /
 ! !      \     /
 ! !        \  /
 ! !          . 
 ! !         t2b 
 
 
 !          dx = tri_ops%d13-x2b
 !          dy = y2a
 !          gamma_x = (tri_ops%t2a - tri_ops%t1)/x2a - tri_ops%t2b /dx
 !          gamma_y = (tri_ops%t2b - tri_ops%t1 )/y2b  + tri_ops%t2a / dy
 !          ! write(*,*) 'option 4'
 !          a = dx**2 + dy**2
 !          b = -2._pr*dx*dy*(dy*gamma_x+gamma_y*dx)
 !          c = (gamma_x**2._pr+gamma_y**2._pr-(2._pr/v_mean)**2._pr)*(dx*dy)**2._pr
 !     endif 
 
 
 !     if ((b**2._pr >= 4._pr*a*c)) then 
 !        t3 = (-b+sqrt(b**2._pr-4._pr*a*c))/(2._pr*a) 
 !     endif 
 return
 end subroutine calc_tcplane
 !###############################################################################



 Next attempt 



    v_mean = (tri_ops%va+tri_ops%vb)/2._pr  

   if (x2a+tri_ops%d13-x2b > tri_ops%d13-x2a+x2b) then
      dx = tri_ops%d13-x2b 
      eta = (tri_ops%t2a-tri_ops%t1)/x2a - tri_ops%t2b/dx
   else 
      dx = tri_ops%d13-x2a 
      eta = (tri_ops%t2b-tri_ops%t1)/x2a - tri_ops%t2a/dx
   endif

    alpha = (tri_ops%t2a-tri_ops%t1*(1._pr-x2a/tri_ops%d13))/y2a+(tri_ops%t1*(1._pr-x2b/tri_ops%d13)-tri_ops%t2b)/y2b
    beta = x2b/y2b/tri_ops%d13  - x2a/y2a/tri_ops%d13 

    a = beta**2._pr+1._pr/dx**2._pr 
    b = 2._pr*alpha*beta+2._pr*eta/dx
    c = alpha**2._pr+eta**2._pr-4._pr/v_mean**2
    write(*,*) 'y2a,y2b: ',x2a,x2b,tri_ops%d13,dx 
    write(*,*) 't1,t2a,t2b',tri_ops%t1,tri_ops%t2a,tri_ops%t2b
    write(*,*) 'eta parts ',(tri_ops%t2a-tri_ops%t1)/x2a,tri_ops%t2b/dx
    write(*,*) 'alpha, beta, eta: ', alpha, beta, eta 
    write(*,*) 'a,b,c: ',a,b,c
    d = b**2._pr - 4._pr*a*c
    if (d >= epsilon(d)) then 
      t3 = (-b+sqrt(d))/(2._pr*a)
      if (t3 < 0. )then 
        write(*,*) 'd,t3 ', d,t3
        stop
      endif
    endif


    !!! next attempt 


   call rotation(x2a,y2a,theta,x2a_r,y2a_r)
   call rotation(tri_ops%d13,0._pr,theta,x3_r,y3_r)
   ! call rotation(x2b,y2b,theta,x2b_r,y2b_r)
   ! write(*,*) 'check rotation',x2b_r,y2b_r,tri_ops%d12b

   ! x2b_r = tri_ops%d12b;   
   dx = x3_r-x2a_r
   dy = y3_r
   
   t3i = tri_ops%t1+(tri_ops%t2b-tri_ops%t1)*x3_r/tri_ops%d12b
   tai = tri_ops%t1+(tri_ops%t2b-tri_ops%t1)*x2a_r/tri_ops%d12b

   write(*,*) 'interpolated time  ',t3i,tai 

   alpha = (tri_ops%t2b-tri_ops%t1)/tri_ops%d12b - tri_ops%t2a/dx 
   beta = (tri_ops%t2a-tai)/y2a_r - t3i/dy
   v_mean = (tri_ops%va+tri_ops%vb)/2._pr  


   a = 1._pr/dx**2._pr+1._pr/dy**2._pr
   b = 2._pr*(alpha/dx+beta/dy)
   c = (alpha**2._pr + beta**2._pr) - 4._pr /v_mean**2._pr

   d = b**2._pr - 4._pr*a*c
   if (d >= epsilon(d)) then 
      t3 = (-b+sqrt(d))/2._pr/a
      write(*,*) 't3....',t3 
   endif