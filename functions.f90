module functions
  use precision
  implicit none
  private

  public :: ff
  public :: pendolo
  public :: harmonic
  public :: pendolo4
  public :: pendolo8
  public :: pendolo16
  public :: pendolo32
  
  real(dp), public :: k1 
  real(dp), public :: k2
  real(dp), public :: Q
  real(dp), public :: A
  real(dp), public :: w  

  contains

  subroutine fsub(t,u,up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)     
    real(dp), intent(out) :: up(:)       

    up(1) = u(2)
    up(2) = k1*u(1) + k2*u(2)
      
  end subroutine fsub
      
  function f1(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))

    up(1) = k1*u(1) 

  end function f1

  
  function ff(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1) = u(2)
    up(2) = k1*u(1) + k2*u(2)

  end function ff

  function pendolo(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1) = u(2)
    up(2) = -sin(u(1)) - u(2)/Q + A*cos(w*t)

  end function pendolo

  function harmonic(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1) = u(2)
    up(2) = -u(1) - u(2)/Q + A*cos(w*t)

  end function harmonic

  function pendolo4(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1:4) = u(5:8)
    up(5:8) = -sin(u(1:4)) - u(5:8)/Q + A*cos(w*t)

  end function pendolo4

  
  function pendolo8(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1:8) = u(9:16)
    up(9:16) = -sin(u(1:8)) - u(9:16)/Q + A*cos(w*t)

  end function pendolo8

  function pendolo16(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1:16) = u(17:32)
    up(17:32) = -sin(u(1:16)) - u(17:32)/Q + A*cos(w*t)

  end function pendolo16
  
  function pendolo32(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1:32) = u(33:64)
    up(33:64) = -sin(u(1:32)) - u(33:64)/Q + A*cos(w*t)

  end function pendolo32


end module functions


