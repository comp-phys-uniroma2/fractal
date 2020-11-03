module runge_kutta
  use precision
  implicit none
  private
  public :: rk2, rk4, rkf45, dp45, adaptive45

contains  

 subroutine rk2(t0, y0, dt, f, y1)
  real(dp) :: t0, dt
  real(dp), dimension(:) :: y0, y1 
  interface
     subroutine f(t,y0,yp)
        use precision
        real(dp) :: t
        real(dp), dimension(:) :: y0, yp
     end subroutine f   
  end interface    
       
  integer :: n
  real(dp), dimension(:), allocatable :: k1, k2
  
  n=size(y0,1)
  allocate(k1(n))
  allocate(k2(n))

  call f(t0,y0,k1)
  call f(t0+dt/2.0d0,y0+k1*dt/2.0d0,k2)

  y1 = y0 + k2*dt

  deallocate(k1,k2)

 end subroutine rk2

 !----------------------------------------------------
 subroutine rk4(t0, y0, dt, f, y1)
  real(dp) :: t0, dt
  real(dp), dimension(:) :: y0, y1 
  interface
     subroutine f(t,y0,yp)
        use precision
        real(dp) :: t
        real(dp) :: y0(:), yp(:)
     end subroutine f   
  end interface    
       
  integer :: n
  real(dp), dimension(:), allocatable :: k1, k2, k3, k4
  
  n=size(y0,1)
  allocate(k1(n))
  allocate(k2(n))
  allocate(k3(n))
  allocate(k4(n))

  call f(t0,y0,k1)
  call f(t0+dt/2.0d0,y0+k1*dt/2.0d0,k2)
  call f(t0+dt/2.0d0,y0+k2*dt/2.0d0,k3)
  call f(t0+dt,y0+k3*dt,k4)

  y1 = y0 + dt*(k1+2.d0*(k2+k3)+k4)/6.d0

  deallocate(k1,k2,k3,k4)

 end subroutine rk4

 ! --------------------------------------------------
 subroutine adaptive45(f, t0, tout, y0, abserr, s)
  real(dp), intent(in) :: t0, tout, abserr
   real(dp), dimension(:), intent(inout) :: y0 
   real(dp), intent(inout) :: s
   interface
     subroutine f(t,y0,yp)
        use precision
        real(dp) :: t
        real(dp), dimension(:) :: y0, yp
     end subroutine f   
   end interface    

   real(dp), parameter :: smin = 0.0625d0
   real(dp) :: t, dt, tf, se, err, h
   real(dp), dimension(:), allocatable :: yh

   allocate(yh(size(y0)))
   t = t0
   h = s * (tout-t0)
   dt = h
 

   do while (t < tout-EPS)
    
     yh = y0
     call dp45(f, t, t+h, y0, err)
        
     print*,'t=',t,'s=',s,'err=',err
     ! estimate interval scaling     
     
     if (err > abserr) then
        s = s/2.d0
        h = s*dt
        y0 = yh
     else   
        se = sqrt(sqrt(abserr/(2.d0*err)*h/dt)) 
        if (se >= 2.d0 .and. s<1.d0) s = s*2.d0
        h = s*dt
        if (t+h > tout) then 
          y0 = yh
          call dp45(f,t,tout,y0,err)
          exit
        end if
        t = t + h
     endif
     
   enddo

 end subroutine adaptive45
 
 ! --------------------------------------------------
 subroutine rkf45(f, t0, tout, y0, abserr)
   real(dp), intent(in) :: t0, tout
   real(dp), dimension(:), intent(inout) :: y0 
   real(dp), intent(inout) :: abserr
   interface
     subroutine f(t,y0,yp)
        use precision
        real(dp) :: t
        real(dp), dimension(:) :: y0, yp
     end subroutine f   
   end interface    

   integer :: n
   real(dp), dimension(:,:), allocatable :: K
   real(dp), dimension(:), allocatable :: y, z 
   real(dp) :: dt, h

   
   n=size(y0,1)

   allocate(K(n,6))
   allocate(y(n))
   allocate(z(n))
   
   h = tout-t0
 
   call f(t0, y0, K(:,1))
   K(:,1) = K(:,1) * h 
   
   call f(t0 + 0.25_dp * h, y0(:) + 0.25_dp * K(:,1), K(:,2)) 
   K(:,2) = K(:,2) * h

   call f(t0 + 0.375_dp * h, y0(:) + 0.09375_dp * K(:,1) + 0.28125_dp * K(:,2), K(:,3)) 
   K(:,3) = K(:,3) * h

   call f(t0 + 12.0_dp/13.0_dp * h, y0(:) + &
        &(1932.0_dp*K(:,1)+7296.0_dp*K(:,3)-7200.0_dp*K(:,2))/2197.0_dp, K(:,4)) 

   K(:,4) = K(:,4) * h

   call f(t0 + h, y0(:) + 439.0_dp/216.0_dp*K(:,1)+ 3680.0_dp/513.0_dp*K(:,3) &
        & -8.0_dp*K(:,2)-845.0_dp/4104.0_dp*K(:,4), K(:,5)) 

   K(:,5) = K(:,5) * h

   call f(t0 + 0.5_dp * h, y0(:) - 8.0_dp/27.0_dp*K(:,1) - 3544.0_dp/2565.0_dp*K(:,3) &
        & -11.0_dp/40.0_dp*K(:,5)+2.0_dp*K(:,2) + 1859.0_dp/4104.0_dp*K(:,4), K(:,6)) 

   K(:,6) = K(:,6) * h

   ! order 4
   y = y0 + 25.0_dp/216.0_dp * K(:,1) + 1408.0_dp/2565.0_dp * K(:,3) &
        & + 2197.0_dp/4104.0_dp * K(:,4) - 0.2_dp * K(:,5)

   ! order 5
   z = y0 + 16.0_dp/135.0_dp * K(:,1) + 6656.0_dp/12825.0_dp * K(:,3) &
        & + 28561.0_dp/56430.0_dp * K(:,4) - 9.0_dp/50.0_dp * K(:,5) + 2.0_dp/55.0_dp * K(:,6)   
 

   abserr = maxval(abs(y-z))  

   y0 = z 

   deallocate(y)
   deallocate(z)
   deallocate(K)
   
 end subroutine rkf45

 ! --------------------------------------------------
 subroutine dp45(f, t0, tout, y0, abserr)
  real(dp), intent(in) :: t0, tout
  real(dp), dimension(:), intent(inout) :: y0 
  real(dp), intent(inout) :: abserr
  interface
     subroutine f(t,y0,yp)
        use precision
        real(dp) :: t
        real(dp) :: y0(:), yp(:)
     end subroutine f   
  end interface    

  real(dp) :: h
  real(dp), dimension(:,:), allocatable :: K
  !real(dp), dimension(:), allocatable :: z, y, yh
  integer :: n

  n=size(y0,1)
  !allocate(z(n))
  !allocate(y(n))
  !allocate(yh(n))
  allocate(K(n,7))

  h = (tout - t0)

  ! 5th order estimate
  call ode4(t0, y0, h, f, k, y0)

  ! 4th order
  !z = yh + (5179.d0/57600.d0*k(:,1)+7571.d0/16695.d0*k(:,3)+393.d0/640.d0*k(:,4)+&
  !           187.d0/2100.d0*k(:,6)+k(:,7)/40.d0 - 92097.d0/339200.d0*k(:,5))*h
 
  ! just estimate the error as  ( z - y ):
  abserr = maxval(abs( 71.d0/57600.d0*k(:,1)+22.d0/525.d0*k(:,6) &
      & +71.d0/1920.d0*k(:,4) &
      & -71.d0/16695.d0*k(:,3)-17253.d0/339200.d0*k(:,5)-1.d0/40.d0*k(:,7) ))

  !abserr = maxval(abs(z))

  deallocate(K)
  !deallocate(y,yh,z)

 end subroutine dp45
 ! --------------------------------------------------
 
 subroutine ode4(t0, y0, h, f, k, y1)
  real(dp), intent(in) :: t0, h
  real(dp), dimension(:), intent(in) :: y0
  real(dp), dimension(:,:), intent(out) :: K
  real(dp), dimension(:), intent(out) :: y1
  interface
     subroutine f(t,y0,yp)
        use precision
        real(dp) :: t
        real(dp) :: y0(:), yp(:)
     end subroutine f   
  end interface    

  real(dp), dimension(:), allocatable :: y
  allocate(y(size(y0)))
  
  call f(t0,y0,k(:,1))

  y= y0+k(:,1)*h/5.0d0
  call f(t0+h/5.0d0,y,k(:,2))

  y= y0+(3.d0/40.0d0 *k(:,1)+9.d0/40.0d0*k(:,2))*h
  call f(t0+3.d0*h/10.0d0,y,k(:,3))
  
  y= y0+(44.d0/45.d0*k(:,1)+160.d0/45.d0*k(:,3)-168.d0/45.d0*k(:,2))*h
  call f(t0+h*4.d0/5.d0,y,k(:,4))

  y= y0+(38744.d0/13122.d0*k(:,1)+128896.d0/13122.d0*k(:,3)-(152160.d0/13122.d0*k(:,2)+&
         3816.d0/13122.d0*k(:,4)))*h
  call f(t0+h*8.d0/9.d0,y,k(:,5))

  y= y0+(9017.d0/3168.d0*k(:,1)+46732.d0/5247.d0*k(:,3)+5194.d0/18656.d0*k(:,4) -&
        (34080.d0/3168.d0*k(:,2)+5103.d0/18656.d0*k(:,5)))*h
  call f(t0+h,y,k(:,6))

  y1 = y0 + (35.d0/384.d0*k(:,1)+500.d0/1113.d0*k(:,3)+125.d0/192.d0*k(:,4)+&
            11.d0/84.d0*k(:,6)-2187.d0/6784.d0*k(:,5))*h
  call f(t0+h,y1,k(:,7))
  
  deallocate(y)

 end subroutine ode4


end module runge_kutta

