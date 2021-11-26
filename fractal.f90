program fractal 
 use precision
 use solvers 
 use functions 
 implicit none

 integer, parameter :: bufsize=4
 integer, parameter :: int1 = selected_int_kind(1)
 integer(int1), dimension(:,:), allocatable :: buf
 integer(int1) :: byte 
 real(dp) :: h, t, tfin, tout, relerr, abserr, s
 real(dp) :: Qstart, y0, v0, dy0, dv0
 real(dp) :: y0min, y0max, v0min, v0max
 real(dp) :: av_teta(0:bufsize-1), y(0:2*bufsize-1) 
 integer :: i, j, k, N, flag, Nstep, Ntr, NQ, nbuf, np
 integer :: pos, bytepos
 character(10) :: arg
 integer, external :: omp_get_num_threads, omp_get_num_procs

 if(iargc()<10) then
  print*,'pendulum tfinal Nsteps A w Q y0min y0max v0min v0max np'
  print*,'tfinal: Integration time ( ~10000)'
  print*,'Nstep:  dt=2Pi/w*1/Nstep ( ~100)'
  stop  
 endif

 call getarg(1,arg)
 read(arg,*) tfin

 call getarg(2,arg)
 read(arg,*) N
 
 call getarg(3,arg)
 read(arg,*) A

 call getarg(4,arg)
 read(arg,*) w
 
 call getarg(5,arg)
 read(arg,*) Qstart
 
 call getarg(6,arg)
 read(arg,*) y0min
 
 call getarg(7,arg)
 read(arg,*) y0max
 
 call getarg(8,arg)
 read(arg,*) v0min
 
 call getarg(9,arg)
 read(arg,*) v0max
 
 call getarg(10,arg)
 read(arg,*) np
 
 h = (2.d0 * pi/w)/N

 Nstep = int(tfin/h)
 Ntr = 3*Nstep/4
 NQ = 1 
 Q = Qstart

 if (mod(np,bufsize) /= 0) then
    print*, 'for optimal performance w must be multiple of',bufsize
    stop
 end if
  
 nbuf = np/8
 allocate(buf(nbuf, np))
 buf=0_int1

 dy0 = (y0max-y0min)/real(np,dp) 
 dv0 = (v0max-v0min)/real(np,dp) 

 !$omp parallel
 write(*,*) 'threads=',omp_get_num_threads(),'/',omp_get_num_procs()
 !$omp end parallel

 !$omp parallel do default(none) shared(buf, np)&
 !$omp& shared(y0min, v0min, dy0, dv0, Nstep, Ntr, h)&
 !$omp& private(i, j, k, bytepos, pos, byte)&
 !$omp& private(v0, y, t, av_teta, abserr)&
 !$omp& schedule(guided)
 rows:do i = 0, np-1
   pos=0
   byte=0_int1
   bytepos=8
   v0 = v0min + dv0*i
   cols:do j = 0, np-1, bufsize
     ! set initial conditions
     ! Equations are vectorized: bufsize eqs are solved in ||
     do k = 0, bufsize-1 
       y(k) = y0min + dy0*k + dy0*j
     end do
     y(bufsize:2*bufsize-1) = v0
     t=0.0_dp
     av_teta = 0.0_dp
     ! run time propagation
     do k = 1, Nstep
       if (k > Ntr) then
         where (y(0:bufsize-1)>Pi) y(0:bufsize-1)=y(0:bufsize-1)-2.0_dp*Pi     
         where (y(0:bufsize-1)<-Pi) y(0:bufsize-1)=y(0:bufsize-1)+2.0_dp*Pi     
         av_teta=av_teta+y(0:bufsize-1)
       end if  
       call dopri54(pendolo4,t,h,y,y,abserr)
       t=t+h
     end do
     ! write byte on buf bitmap array
     do k = 0, bufsize-1
        bytepos=bytepos-1
        if (av_teta(k)>0) then
           byte = ibset(byte, bytepos)
        end if
        if(bytepos==0) then
           bytepos=8
           pos = pos + 1
           buf(pos,i+1) = byte                         
           byte=0_int1
        end if
     end do

   end do cols
 end do rows
 !$omp end parallel do

 open(100, file='fractal.bmp')   
 write(100,'("P4",/,i0," ",i0)') np, np
 do i = np,1,-1
   write(100, '(10000000a1)', advance='no') buf(1:nbuf,i)
 end do
 close(100)

end program fractal 
