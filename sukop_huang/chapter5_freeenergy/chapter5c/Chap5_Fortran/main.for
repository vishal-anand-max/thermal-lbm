c In this porgram the Swift free-energy model for D2Q9 model is presented.

      program Swift_D2Q9
      implicit none
      include "head.inc"
      integer time, k, m, n

      integer obst(lx,ly),x,y
c 'obst' denotes the lattice node occupied by fluid (obst=0) or solid (obst=1)
c The 9-speed lattice is used here:
c              6   2   5
c                \ | /
c              3 - 0 - 1
c                / | \
c              7   4   8
c Array u_x and u_y are x- and y- components of the velocity, respectively.
c 'rho' denotes the density of the fluid.

      real*8  u_x(lx,ly),u_y(lx,ly),rho(lx,ly)

c Array 'ff' is the distribution function. Array 'p' is the pressure.
      real*8  ff(0:8,lx,ly),Fx(lx,ly),Fy(lx,ly),
     &   p(lx,ly)

      data xc/0.d0,1.d0,0.d0, -1.d0, 0.d0, 1.d0, -1.d0, -1.d0, 1.d0/,
     &   yc/0.d0,0.d0,1.d0, 0.d0, -1.d0, 1.d0, 1.d0, -1.d0, -1.d0/

      character filename*16, B2*2, C*3, D*2
c-----------------------------------------------------
c     Author: Haibo Huang, huanghb@ustc.edu.cn
c-----------------------------------------------------
c   lattice speed cc=1 lu/ts.
      cc = 1.d0
      c_squ = cc *cc / 3.d0

      Nwri= 2000

      write (6,*) '@@@  Swift D2Q9 starting ...           @@@'
      write (6,*) 'Computational domain lx = ,', lx, ' ly = ',ly

c Begin initialization-------
c Read parameters.
      call read_parametrs()
c Initialize the fluid nodes.
      call read_obstacles(obst)
c Initialize the macro variables and distribution functions.
      call init_density(obst,u_x ,u_y ,rho ,p, ff, Fx, Fy)
c Dump initial flow field
        call write_results(obst,rho,u_x,u_y,p, 0)

c Begin iteration ------

      do 100 time = 1, t_max
c Dump flow field data
        if ( mod(time, Nwri) .eq. 0) then
             write(*,*) time
             call write_results(obst,rho,u_x,u_y,p, time)
        end if
c Streaming step.
        call stream(obst,ff )

c Update the macro variables.
        call getuv(obst,u_x ,u_y, rho, ff )

c Collision step.
        call collision(tau,obst,u_x,u_y,rho ,ff, p )

  100 continue
c End iteration ------

      write (6,*) '------   end   ---------'
      end

c------------------------------------------------
      subroutine read_parametrs()
      implicit none
      include "head.inc"
      real*8 Re1, V1
! Parameters  in the van der Waals EOS.
      TT = 0.56d0
      a = 9.d0/49.d0
      b = 2.d0/21.d0

! Parameter to adjust the surface tension.
       kappa = 0.01d0

! The other parameters control the flow.
       t_max  = 40000
       tau = 1.0d0      ! Relaxation time constant
       nu = ( tau -0.5d0 )/3.d0  ! Kinematic viscosity
      end

c------------------------------------------------
      subroutine read_obstacles(obst)
      implicit none
      include "head.inc"
      integer  x,y,obst(lx,ly)
        do 10 y = 1, ly
          do 40 x = 1, lx
   40       obst(x,y) =  0    ! Set all nodes to be fluid nodes
   10   continue
      end

c------------------------------------------------
      subroutine init_density(obst,u_x,u_y,rho,p, ff, Fx, Fy)
      implicit none
      include "head.inc"
      integer i,j,x,y,k,n,obst(lx,ly)
      real*8  u_squ, xx, u_n(0:8),u_x(lx,ly),u_y(lx,ly),
     & rho(lx,ly),ff(0:8,lx,ly), p(lx,ly)
      real*8 Fx(lx,ly),Fy(lx,ly)

c Initialize the macro variables.
      do 10 y = 1, ly
       do 10 x = 1, lx
        u_x(x,y) = 0.d0
        u_y(x,y) = 0.d0
        Fx(x,y) = 0.d0
        Fy(x,y) = 0.d0
c        rho(x,y) = 2.3d0
   10 continue

      do 20 x = 1, lx
      do 30 y = 1, ly

c----------------------
c Initial a uniform density with small disturbance.
      call random_number (xx) ! Get a random number between [0,1]
         rho(x,y) = 3.8d0 + 0.01d0 * xx
c---------------------
c  Initialize a droplet or a bubble inside the computational domain
c      if(  sqrt(float(x-lx/2)**2+float(y-ly/2)**2) .lt. 30.d0)
c     &  rho(x,y) = 4.5d0

   30 continue
   20 continue

c  Initialize the distribution functions.
      call get_feq(rho,u_x,u_y,ff, p)

       end

c-----------------------------------------------------
      subroutine EOS(rho, tmp)
      include 'head.inc'
      real*8 rho, tmp
        tmp= rho*TT/(1.d0 -rho*b) -a *rho *rho
       return
      end

c-----------------------------------------------------
      subroutine  write_results(obst,rho,upx,upy,p, n)
      implicit none
      include "head.inc"
      integer  x,y,n,obst(lx,ly)
      real*8  rho1, rho(lx,ly), upx(lx,ly), upy(lx,ly)
     & , p(lx,ly)
      character filename*20,  D*7

      write(D,'(i7.7)') n
      filename='out/swift'//D//'.plt'

      open(41,file=filename)

      write(41,*) 'variables = x, y, u, v, rho1, p, obst'
      write(41,*) 'zone i=', lx, ', j=', ly, ', f=point'

c   Write results to file (an ASCII file)
      do 10 y = 1, ly
      do 10 x = 1, lx
          write(41,9) x, y,
     &       upx(x,y), upy(x,y), rho(x,y), p(x,y), obst(x,y)
   10 continue

   9  format(i4,i5, 3f15.8, f15.8, i4 )

      close(41)
      end
c-----------------------------------------------------




c FILE Getfeq.for
c============================================
