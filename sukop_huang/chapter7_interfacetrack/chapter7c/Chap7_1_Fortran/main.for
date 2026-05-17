      program d2q9lbm
      implicit none
      include "head.inc"
      integer t_max,time,k,i, j, x, y
      integer obstA(xA,yA)
c     a 9-speed lattice is used here, other geometries are possible
c              6   2   5
c                \ | /
c              3 - 0 - 1
c                / | \
c              7   4   8
c
c    the lattice nodes are numbered as follows:
      real*8  deltaA(8,xA,yA),UA(xA,yA),VA(xA,yA),rhoA(xA,yA)
      real*8  noteqA(0:8,xA,yA), ffA(0:8,xA,yA), gpA(0:8,xA,yA)
      real*8   rhA(xA,yA), up(xA,yA), vp(xA,yA)
      real*8  Force1(xA,yA)
      real*8  Force2(xA,yA),pA(xA,yA)
      character BT*10
      integer*4 now(3)
c-------------------------------
c Author:     Haibo Huang   email: huanghb@ustc.edu.cn
c-------------------------------
      data xv/0.d0,1.d0,0.d0, -1.d0, 0.d0, 1.d0, -1.d0, -1.d0, 1.d0/,
     &   yv/0.d0,0.d0,1.d0, 0.d0, -1.d0, 1.d0, 1.d0, -1.d0, -1.d0/

      data opp/3,4,1,2,7,8,5,6/  ! opposite direction of each velocity
! Opposite direction of each velocity in the D2Q9 model, opp(1)=3, opp(2)=4, ...
! This array is serve for the bounce back implementation.
c------------------------
      BD = 0

      pai = datan(1.d0)*4.d0
      c_squ = 1.d0 / 3.d0
      t_0 =  4.d0 / 9.d0
      t_1 =  1.d0 / 9.d0
      t_2 =  1.d0 / 36.d0

      t_k(0) = t_0
      do 1 k =1,4
        t_k(k) = t_1
    1 continue
      do 2 k =5,8
        t_k(k) = t_2
    2 continue
! t_k(k) are the constant coefficient $\omega_i$ in $f^{eq}$, for D2Q9,
! when k=0, $\omega_i=4/9$, when i=1,..4, $\omega_i=1/9$, ......

      con2 =  2.d0 * c_squ * c_squ
! Another constant will be used in calculation of $f^{eq}$.

      write (6,*) '@@@ He-Chen-Zhang multiphase LBM ...    @@@'
      write (6,*) '@@@ precompiled for lattice size xA = ',xA
      write (6,*) '@@@                              yA = ',yA

! Initialization
      call read_parametrs(yA,t_max)
      call read_obst(xA,yA,obstA)

      call init_density(xA,yA,obstA,UA,VA,rhA,
     &  rhoA,pA,ffA,gpA,Force1,Force2)
      call write_results(xA,yA,obstA,rhoA,UA,VA,pA,0)
!      Nwri = 0.5*int(Tn)  ! every 0.5 non-dimensional time

      call comp_rey(xA,yA,time,BT)

      open(39,file='./out/residue.dat')
      open(38,file='./out/mass.dat')

c---------------------------------
c     Begin iterations
c---------------------------------
      do 100 time = 1, t_max

      if(mod(time,Nwri) .eq. 0) then
       write(*,*) time
       call itime(now)
       write(*,"(i2.2, ':', i2.2, ':', i2.2)") now
       call write_results(xA,yA,obstA,rhoA,UA,VA,pA,time)
      endif
        if ( mod(time,20) .eq. 0)
     &  call track_interface(xA,yA,rhoA,time)
! this subroutine has demonstrate in the section "Capillary rise".
! Insert the subroutine to this fortran file is ok to compile.
c--------------
      call stream(xA,yA,obstA,ffA,gpA)
      call getuv(xA,yA,obstA,UA,VA,rhA, rhoA,ffA, pA,Force1, Force2)
      call geten(xA,yA,obstA,UA,VA,rhoA,gpA,pA, Force1,Force2)
      call collision(xA,yA,obstA,UA,VA,rhA,rhoA,pA,ffA,gpA,
     &  Force1,Force2)

  100 continue

      close(39)
      write (6,*) '@@@@@@@@@@@@@@@@    end     @@@@@@@@@@@@@@@@'
      stop
      end
c---------------------------------------------
      subroutine index_function(lx,ly,rh,rho)
      implicit none
      include "head.inc"
      integer lx,ly,x,y
      real*8 rh(lx,ly), rho(lx,ly)
      do 1 x = 1,lx
      do 1 y = 1,ly
       rho(x,y) = rho_l + ( rh(x,y)- psi_min0 )*const
! From the index function, we can get the density of the fluid at (x,y).
    1 continue
      end
c---------------------------------------------
