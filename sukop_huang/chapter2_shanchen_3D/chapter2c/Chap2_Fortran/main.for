!Please create a new folder 'out' in the working directory
! before running the code
      program D3Q19LBM
      implicit none
      include "head.inc"
      integer t_max,time,k
! This array defines which lattice positions are occupied by fluid nodes (obst=0)
! or solid nodes (obst=1)
      integer obst(lx,ly,lz)
! Velocity components
      real*8  u_x(lx,ly,lz),u_y(lx,ly,lz), u_z(lx,ly,lz)
! Pressure and density
      real*8  p(lx,ly,lz),rho(lx,ly,lz)
! The real fluid density
! which may differ from the velocity componets in the above; refer to SC model
      real*8  upx(lx,ly,lz),upy(lx,ly,lz),upz(lx,ly,lz)
!-----------------------------
! Author: Huanghb@ustc.edu.cn
!-----------------------------

! The force components: Fx, Fy, Fz for the interaction between fluid nodes.
! Sx, Sy, Sz are the interaction (components) between the fluid nodes and solid nodes
! ff is the distribution function
      real*8  ff(0:18,lx,ly,lz),Fx(lx,ly,lz),Fy(lx,ly,lz),
     &   Fz(lx,ly,lz), Sx(lx,ly,lz),Sy(lx,ly,lz), Sz(lx,ly,lz)

! TT0W is the value of T/T0;   RHW and RLW are the coexisting densities
! in the sepcified T/T0.
! For initialization, \rho_l (lower density)
! and \rho_h (higher density) are supposed to be known.
      real*8   TT0W(12), RHW(12), RLW(12)

! The below data define the D3Q19 velocity model, xc(ex), yc(ey), zc(ez)
! are the components of e_{ix}, e_{iy}, and e_{iz}, respectively.
      data xc/0.d0, 1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 1.d0,
     &  -1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0, 0.d0, 0.d0,
     &  0.d0, 0.d0  /,
     &   yc/0.d0, 0.d0, 0.d0, 1.d0, -1.0d0, 0.d0, 0.d0, 1.d0, -1.d0,
     &  1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 1.d0,
     &  -1.d0, -1.d0/,
     &   zc/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, -1.d0, 0.d0, 0.d0,
     &  0.d0, 0.d0, 1.d0, 1.d0, -1.d0, -1.d0, 1.d0, -1.d0,
     &  1.d0, -1.d0/
      data ex/0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 1, -1, 0, 0,
     &  0, 0 /,
     &   ey/0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0,  0, 0,  0, 1, 1,
     &  -1, -1/,
     &   ez/0, 0, 0, 0, 0, 1, -1, 0, 0,  0, 0, 1, 1, -1, -1, 1, -1,
     &  1, -1/

! This array gives the opposite direction for e_1, e_2, e_3, .....e_18
! It implements the simple bounce-back rule we use in the collision step
! for solid nodes (obst=1)
      data opp/2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15/

! C-S EOS
! RHW and RLW are the coexisting densities in the corresponding sepcified T/T0.
      data TT0W/0.975d0, 0.95d0, 0.925d0, 0.9d0, 0.875d0, 0.85d0,
     &    0.825d0,  0.8d0, 0.775d0, 0.75d0, 0.7d0, 0.65d0 /,
     &     RHW/ 0.16d0, 0.21d0,  0.23d0, 0.247d0, 0.265d0, 0.279d0,
     &    0.29d0, 0.314d0, 0.30d0,  0.33d0, 0.36d0, 0.38d0 /,
     &     RLW/0.08d0, 0.067d0, 0.05d0, 0.0405d0, 0.038d0, 0.032d0,
     &    0.025d0, 0.0245d0, 0.02d0, 0.015d0, 0.009d0, 0.006d0/
! Speeds and weighting factors
      cc = 1.d0
      c_squ = cc *cc / 3.d0
      t_0 =  1.d0 / 3.d0
      t_1 =  1.d0 / 18.d0
      t_2 =  1.d0 / 36.d0
! Weighting coefficient in the equilibrium distribution function
      t_k(0) = t_0
      do 1 k =1,6
       t_k(k) = t_1
    1 continue
      do 2 k =7,18
       t_k(k) = t_2
    2 continue

! Please specify which temperature
! and corresponding \rho_h, \rho_l in above 'data' are chosen.
! Initial T/T0, rho_h, and rho_l for the C-S EOS are listed in above 'data' section
      k = 5    ! important
      TT0 = TT0W(k)
      rho_h = RHW(k)
      rho_l = RLW(k)


c=======================================================================
c     Initialisation
c=======================================================================
      write (6,*) '@@@  3D LBM for single component multiphase @@@'
      write (6,*) '@@@ lattice size lx = ',lx
      write (6,*) '@@@              ly = ',ly
      write (6,*) '@@@              lz = ',lz
      call read_parameters(t_max)
      call read_obstacles(obst)

      call init_density(obst,u_x ,u_y ,rho ,ff )
      open(40,file='residue.dat')
c=======================================================================
! Begin iterations
      do 100 time = 1, t_max
        if ( mod(time, Nwri) .eq. 0 .or. time. eq. 1) then
        write(*,*) time
        call write_results2(obst,rho,p,upx,upy,upz,time)
        end if

        call stream(obst,ff )   ! streaming (propagation) step

! Obtain the macro variables
      call getuv(obst,u_x ,u_y, u_z, rho, ff )

! Calculate the actual velocity
      call calcu_upr(obst,u_x,u_y,u_z,Fx,Fy,Fz,
     &    Sx,Sy,Sz,rho, upx,upy,upz)

! Calculate the interaction force between fluid nodes,
! and the interaction force between solid and fluid nodes.
      call calcu_Fxy(obst,rho,Fx,Fy,Fz,Sx,Sy,Sz,p)

!     BGK model (a single relaxation parameter) is used
      call collision(tau,obst,u_x,u_y,u_z,rho ,ff ,Fx ,Fy ,Fz,
     &    Sx,Sy, Sz )    ! collision step ,

  100 continue
c===== End of the main loop
      close(40)
      write (6,*) '@@@@** end **@@@@'
      end
c-------------------------------------
      subroutine read_parameters(t_max)
      implicit none
      include "head.inc"
      integer  t_max
      real*8 visc

      open(1,file='./params.in')
! Initial radius of the droplet.
      read(1,*) RR
! \rho_w in calculation of fluid-wall interaction
      read(1,*) rho_w
! Relaxation parameter, which is related to viscosity
      read(1,*) tau
! Maximum iteration specified
      read(1,*) t_max
! Output data frequency (can be viewed with TECPLOT)
      read(1,*) Nwri
      close(1)
       visc =c_squ*(tau-0.5)
       write (*,'("kinematic viscosity=",f12.5, "lu^2/ts",
     &   2X, "tau=", f12.7)') visc, tau

      end

!---------------------------------------------------
! Initialize which nodes are wall node (obst=1) and
! which are fluid nodes (obst=0)
      subroutine read_obstacles(obst)
      implicit none
      include "head.inc"

      integer  x,y,z,obst(lx,ly,lz)
        do 11 z = 1, lz
         do 10 y = 1, ly
          do 40 x = 1, lx
           obst(x,y,z) =  0
           if(z .eq. 1)   obst(x, y,1) = 1
   40     continue
   10   continue
   11 continue

      end
!--------------------------------------------------
      subroutine init_density(obst,u_x,u_y,rho,ff)
      implicit none
      include "head.inc"

      integer i,j,x,y,z,k,n,obst(lx,ly,lz)
      real*8  u_squ,u_n(0:18),fequi(0:18),u_x(lx,ly,lz),u_y(lx,ly,lz),
     & rho(lx,ly,lz),ff(0:18,lx,ly,lz),u_z(lx,ly,lz)

      do 12 z = 1, lz
       do 11 y = 1, ly
        do 10 x = 1, lx
        u_x(x,y,z) = 0.d0
        u_y(x,y,z) = 0.d0
        u_z(x,y,z) = 0.d0
        rho(x,y,z) = rho_l
        if(real(x-lx/2)**2+real(y-ly/2)**2+real(z-5)**2< RR**2) then
           rho(x,y,z) = rho_h
         endif
   10   continue
   11  continue
   12 continue

      do 82 z = 1, lz
       do 81 y = 1, ly
        do 80 x = 1, lx
            u_squ  =   u_x(x,y,z)*u_x(x,y,z) + u_y(x,y,z)*u_y(x,y,z)
     &           +   u_z(x,y,z) *u_z(x,y,z)
         do 60 k = 0,18
            u_n(k)   = xc(k)*u_x(x,y,z) + yc(k)*u_y(x,y,z)
     &              + zc(k) *u_z(x,y,z)
            fequi(k) = t_k(k)* rho(x,y,z) * ( cc*u_n(k) / c_squ
     &               + (u_n(k)*cc) *(u_n(k)*cc) / (2.d0 * c_squ *c_squ)
     &               - u_squ / (2.d0 * c_squ)) + t_k(k) * rho(x,y,z)
            ff(k,x,y,z)= fequi(k)
   60    continue
   80   continue
   81  continue
   82 continue
      end
!---------------------------------------------------
