      program D3Q19LBM

      implicit none
      include "head.inc"
      integer t_max,time,k
c   USE DFPORT

      integer obst(lx,ly,lz)
c-------------------------------------
! Author:  Haibo Huang, huanghb@ustc.edu.cn
c-------------------------------------
c    the lattice nodes are numbered as follows:
      real*8  u_x(lx,ly,lz),u_y(lx,ly,lz),rho(lx,ly,lz)
      real*8  rh(lx,ly,lz)
      real*8  u_z(lx,ly,lz),  p(lx,ly,lz)
      real*8  ff(0:18,lx,ly,lz),gp(0:18,lx,ly,lz),Fx(lx,ly,lz)
     &  ,Fy(lx,ly,lz), Fz(lx,ly,lz)
! Velocity vectors of the D3Q19 velocity model, xc, yc, zc are
! the x, y, z components of the vector, respectively. The order of 
! the velocity vector is slightly different from that in Fig.1.1.
      data xc/0.d0, 1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 1.d0,
     &  -1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0, 0.d0, 0.d0,
     &  0.d0, 0.d0  /,
     &   yc/0.d0, 0.d0, 0.d0, 1.d0, -1.0d0, 0.d0, 0.d0, 1.d0, -1.d0,
     &  1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 1.d0,
     &  -1.d0, -1.d0/,
     &   zc/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, -1.d0, 0.d0, 0.d0,
     &  0.d0, 0.d0, 1.d0, 1.d0, -1.d0, -1.d0, 1.d0, -1.d0,
     &  1.d0, -1.d0/
      data ex/0, 1, -1, 0, 0, 0, 0, 1, 1,
     &  -1, -1, 1, -1, 1, -1, 0, 0,
     &  0, 0 /,
     &   ey/0, 0, 0, 1, -1, 0, 0, 1, -1,
     &  1, -1, 0, 0, 0, 0, 1, 1,
     &  -1, -1/,
     &   ez/0, 0, 0, 0, 0, 1, -1, 0, 0,
     &  0, 0, 1, 1, -1, -1, 1, -1,
     &  1, -1/
      data opp/2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15/
c-----------------

      BD = 0  ! It was not used at present code. It can be removed at present stage.

! It turn on the buoyancy force if BUOYANCY=1
      BUOYANCY = 0

! Lattice speed
      cc = 1.d0
! Square of sound speed in the LBM
        c_squ = cc *cc / 3.d0
c-------------
! Weighting coefficients in the equilibrium distribution function for D3Q19 model.
        t_0 =  1.d0 / 3.d0
        t_1 =  1.d0 / 18.d0
        t_2 =  1.d0 / 36.d0

      t_k(0) = t_0
      do 1 k =1,6
      t_k(k) = t_1
    1 continue
      do 2 k =7,18
      t_k(k) = t_2
    2 continue
c-------------

      write (6,*)
      write (6,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*'
      write (6,*) '@@@  LBM for He_Chen_ZhangLBM starting ...       @@@'
      write (6,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*'
      write (6,*) '@@@ precompiled for lattice size lx = ',lx
      write (6,*) '@@@                              ly = ',ly
      write (6,*) '@@@                              lz = ',lz
      write (6,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*'
c----------------------------
c     Begin initialization
c----------------------------
! Specify how many steps to dump flow data file (flow field can be viewed by Tecplot)
      Nwri = 100

! Read important parameters in the flow we want to investigate
      call read_parameters(t_max)
! Initial each lattice node to be a fluid node (obst=0) or a solid node (obst=1)
      call read_obstacles(obst)
! Initial the flow field and f_i, g_i
      call init_density(obst,u_x,u_y,u_z,rh, rho,p, ff
     &  ,gp, Fx,Fy,Fz)

! Begin main loop
c-----------------------------------------------
      do 100 time = 1, t_max
      if ( mod(time, 50).eq. 0)  write(*,*) time
! Each "Nwri" time step dump a result
        if ( mod(time, Nwri) .eq. 0 .or. time. eq. 1) then
        call write_results2(obst,rho,u_x,u_y,u_z,p,time)
        end if

! Streaming step
        call stream(obst,ff )
        call stream(obst,gp )

! Obtained the macroscopic varibales
      call getuv(obst,u_x,u_y,u_z,rh, rho,ff, p,
     &  Fx, Fy,Fz)
      call geten(obst,u_x,u_y,u_z, rho,gp, p, Fx, Fy, Fz)
! Collision step
      call collision(obst,u_x,u_y,u_z,rh,
     &  rho,p,ff,gp, Fx,Fy,Fz)
c-----------------------------------------------
! End of main loop

  100 continue
  101   call write_results2(obst,rho,u_x,u_y,u_z,p,time)

  999 continue
      write (6,*) '@@@@@@@@@@@@@@@@@@**    end     @@@@@@@@@@@@@@@@@@**'
      end


c-----------------------------------------------
      subroutine read_parameters(t_max)
      implicit none
      include "head.inc"
      real*8  Re1,  We
      integer  t_max

! This is the critical parameter to adjust surface tension independently.
      Kappa = 0.02d0
! Specify a characteristic velocity
      UU = 0.00d0
! Specify droplet's diameter
      DD = 30.d0

!  The \rho_h and \rho_l seems not able to choose randomly, it should
! be near to the values of  \psi_{max} and \psi_{min}.
! Usuually, the maxmum density ratio should not exceed the ratio between \psi_{max}
! and \psi_{min}.

      rho_h = 0.12d0
      rho_l = 0.04d0  ! Density ratio is 3
!  The \psi_{max} and \psi_{min} depend on the EOS we used.
      psi_max = 0.2516d0
      psi_min = 0.0241d0

      t_max  = 50000
! Relaxation time in the LBE for f_i.
      tau_f = 0.7d0
! Relaxation time in the LBE for g_i.
      tau_g = 0.7d0
! External force
      gforce = -0.0000d0
! Reynolds number
      Re1 = UU *DD/(1./3.)/tau_f
! Web number in the simulation
      We = rho_h *DD *UU*UU/(0.0011d0 *Kappa/0.1d0)
      write(*,'("tau = ", f13.7, 4X, "gforce=", f13.7)') tau_f,gforce
      write(*,'("UU=", 1f13.7, 10X, "Re1=",1f13.7, 10X, "We=", 1f13.7)')
     &   UU,Re1, We
       return
      end
c-----------------------------------------------

      subroutine read_obstacles(obst)
      implicit none
      include "head.inc"
      integer  x,y,z,obst(lx,ly,lz)
! Initially set all nodes to be fluid nodes
        do 5 z = 1, lz
         do 10 y = 1, ly
          do 40 x = 1, lx
            obst(x,y,z) =  0
   40     continue
   10   continue
   5  continue


      end
c-----------------------------------------------



      subroutine init_density(obst,u_x,u_y,u_z,rh, rho,p, ff
     &  ,gp, Fx,Fy,Fz)
      implicit none
      include "head.inc"

      integer i,j,x,y,z,k,n,obst(lx,ly,lz)
      real*8  u_squ,u_n(0:18),fequi(0:18),u_x(lx,ly,lz),u_y(lx,ly,lz),
     & rho(lx,ly,lz),ff(0:18,lx,ly,lz),u_z(lx,ly,lz),rh(lx,ly,lz)
     & ,p(lx,ly,lz), Fx(lx,ly,lz), Fy(lx,ly,lz),Fz(lx,ly,lz),
     & gp(0:18,lx,ly,lz)

      do 12 z = 1, lz
        do 11 y = 1, ly
          do 10 x = 1, lx
        u_x(x,y,z) = 0.d0
        u_y(x,y,z) = 0.d0
        u_z(x,y,z) = 0.d0
        rh(x,y,z) = psi_min  !psi_max
        rho(x,y,z) = rho_l   !rho_h
cccc ----------------- phase-seperate
c   call random_number (xx) ! pursu random
c          rh(x,y) = 0.1d0 + 0.01* xx
c          rho(x,y) = rh(x,y)
cccc ----------------- phase-seperate
        Fx(x,y,z) = 0.d0
        Fy(x,y,z) = 0.d0
        Fz(x,y,z) = 0.d0

   10     continue
   11   continue
   12 continue


      do 23 z = 1, lz
      do 22 y = 1, ly
      do 21 x = 1, lx

! Here a droplet is initialized.
      if(  sqrt(float(x-lx/2)**2 + float(y-ly/2)**2 +
     &   float(z-lz/2)**2) .lt. DD/2.d0 ) then

          rh(x,y,z) = psi_max !psi_min  !psi_max  ! !
          rho(x,y,z) = rho_h  !rho_l  !rho_h
          u_z(x,y,z) = 0.d0   !UU/2.d0
       endif

! Here another droplet is initialized. For the case of droplet collision, it will be used.
c---------------
      if(0 .eq. 1)  then
        if(  sqrt(float(x-lx/2)**2 + float(y-ly/2)**2 +
     &   float(z-62)**2) .lt. DD/2.d0 ) then !.or.
          rh(x,y,z) = psi_max !psi_min  !psi_max  ! !
          rho(x,y,z) = rho_h  !rho_l  !rho_h
          u_z(x,y,z) = -UU/2.d0
        endif
      endif
c---------------

   21 continue
   22 continue
   23 continue



      do 82 z = 1, lz
         do 81 y = 1, ly
          do 80 x = 1, lx
! Initially the hydrodynamic pressure is supposed to be equal to the thermodynamic one, which can be
! calculated from the C-S EOS (He, Chen, Zhang, JCP, 1999).
! In the C-S EOS: a=12RT=12c_s^2, b=4. The values are in lattice units
        p(x,y,z) =   rho(x,y,z)*rho(x,y,z)
     *      *c_squ*( 4.d0 - 2.d0*rho(x,y,z) )
     *           /(1.d0- rho(x,y,z))**3.d0
     *  - 12.d0* c_squ* rho(x,y,z)* rho(x,y,z) + rho(x,y,z)*c_squ

      if(obst(x,y,z) .eq. 1) then

      do 59 k=0, 18
        ff(k,x,y,z) = 0.d0
        gp(k,x,y,z) = 0.d0
   59 continue

      else
! Initial f_i is supposed to be the equilibrium value.
       call f_eq( rh(x,y,z),u_x(x,y,z),u_y(x,y,z),u_z(x,y,z),fequi)
      do 60 k = 0,18
          ff(k,x,y,z) = fequi(k)*rh(x,y,z)
   60  continue
! Initial g_i is supposed to be the equilibrium value.
      call g_eq(p(x,y,z),rho(x,y,z),u_x(x,y,z),u_y(x,y,z),
     &      u_z(x,y,z),fequi)
      do 61 k = 0,18
          gp(k,x,y,z) = fequi(k)
   61  continue
      endif

   80    continue
   81  continue
   82 continue
      end

