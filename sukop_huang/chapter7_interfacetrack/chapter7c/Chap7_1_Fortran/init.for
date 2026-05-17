      subroutine read_parametrs(ly,t_max)
      implicit none
      include "head.inc"
      real*8  Re1, V1
      integer  t_max,ly


      Kappa = 0.1d0
      t_max  = 15000
      tau_f = 0.515d0
      tau_g = 0.515d0
      gforce = 0.d0  !-0.00008d0

      psi_max = 0.251260d0
      psi_min = 0.02428d0
      psi_max0 = 0.230d0
      psi_min0 = 0.027d0
      rho_h = 0.251260d0
      rho_l = 0.02428d0

      open(1,file='./params.in')
        read(1,*) psi_max
        read(1,*) psi_min
        read(1,*) psi_max0
        read(1,*) psi_min0
        read(1,*) rho_h
        read(1,*) rho_l
        read(1,*)
        read(1,*) Re1
        read(1,*) V1
        read(1,*) t_max
        read(1,*) Nwri
        read(1,*) Kappa
      close(1)

      const=  (rho_h - rho_l)/(psi_max0 - psi_min0)

c   V1 = 0.10d0
      gforce = -V1*V1/h_ref    !
      tau_f= 1.d0/(Re1/V1/real(xA)*c_squ) +0.5d0  !
      tau_g = tau_f    !
      Tn = sqrt(-xA/gforce)! characteristic time

      RAYLEIGH =1
      BD =1
      write(*,'("tau = ", f13.7, 4X, "gforce=", f13.7)') tau_f,gforce
      write(*,'("V1=", 1f13.7, 10X, "Re1=",1f13.7)') V1,Re1
      write(*, '("time is normalized by:", f15.4 )') sqrt(-xA/gforce)

      end
c----------------------------------------
      subroutine read_obst(lx,ly,obst)
      implicit none
      include "head.inc"
      integer  x,y,lx,ly,obst(lx,ly)
       do 10 y = 1, ly
        do 40 x = 1, lx
          obst(x,y) = 0
   40   continue
   10  continue

      if (BD .eq. 1) then
      do 20 y = 1, ly, ly-1
       do 20 x = 1, lx
         obst(x,y) = 1
   20 continue
      endif

      return
      end
c--------------------------------
      subroutine init_density(lx,ly,obst,u_x,u_y,rh,rho,p,ff,gp,
     &  Force1,Force2)
      implicit none
      include "head.inc"
      integer lx,ly,i,j,x,y,k,n,obst(lx,ly)
      real*8  u_squ,u_n(0:8),fequi(0:8),u_x(lx,ly),u_y(lx,ly),
     & rho(lx,ly),ff(0:8,lx,ly), gp(0:8,lx,ly), p(lx,ly), gequi(0:8),
     &  Force1(lx,ly), Force2(lx,ly),rh(lx,ly), rhoh,rhol, xx

      do 10 y = 1, ly
       do 10 x = 1, lx
        u_x(x,y) = 0.d0
        u_y(x,y) = 0.d0
        rh(x,y) = psi_min  !psi_max
        Force1(x,y) = 0.d0
        Force2(x,y) = 0.d0
   10 continue

c----------------------------------
      if(RAYLEIGH .eq. 1) then
      do 21 x = 1, lx
      do 31 y = 1, ly
      if(y.gt. real(ly/2)+h_ref* 0.1d0 *dcos(real(x)/h_ref* 2.d0* pai))
     &  then
            rh(x,y) = psi_max !psi_min ! !
      endif

   31 continue
   21 continue
      endif

      call index_function(lx,ly, rh,rho)

      do y = 1, ly
      do x = 1, lx
        p(x,y) = rho(x,y)*rho(x,y) *c_squ*( 4.d0 - 2.d0*rho(x,y) )
     *           /(1.d0- rho(x,y))/(1.d0- rho(x,y))/(1.d0- rho(x,y))
     *          - 12.d0* c_squ* rho(x,y)* rho(x,y) + rho(x,y)*c_squ
      enddo
      enddo
c-----------------------------------
      do 80 y = 1, ly
      do 80 x = 1, lx

        call f_eq( rh(x,y),u_x(x,y),u_y(x,y),fequi)
        call g_eq(p(x,y),rho(x,y),fequi,gequi)
        do 60 k = 0,8
          ff(k,x,y) = fequi(k)*rh(x,y)
          gp(k,x,y) = gequi(k)
   60 continue

   80 continue
      end

