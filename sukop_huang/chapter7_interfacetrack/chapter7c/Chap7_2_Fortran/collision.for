      subroutine stream(obst,f)
      implicit none
      include "head.inc"
      integer  k, obst(lx,ly,lz)
      real*8 f(0:18,lx,ly,lz),f_hlp(0:18,lx,ly,lz)
      integer  x,y,z,x_e,x_w,y_n,y_s,z_n,z_s

      do 12 z = 1, lz
       do 11 y = 1, ly
        do 10 x = 1, lx
! In x, y, z directions, periodic boundary conditions are applied
          z_n = mod(z,lz) + 1
          y_n = mod(y,ly) + 1
          x_e = mod(x,lx) + 1

          z_s = lz - mod(lz + 1 - z, lz)
          y_s = ly - mod(ly + 1 - y, ly)
          x_w = lx - mod(lx + 1 - x, lx)
c.........density propagation
          f_hlp(1 ,x_e,y  ,z  ) = f(1,x,y,z)
          f_hlp(2 ,x_w,y  ,z  ) = f(2,x,y,z)
          f_hlp(3 ,x  ,y_n,z  ) = f(3,x,y,z)
          f_hlp(4 ,x  ,y_s,z  ) = f(4,x,y,z)
          f_hlp(5 ,x  ,y  ,z_n) = f(5,x,y,z)
          f_hlp(6 ,x  ,y  ,z_s) = f(6,x,y,z)
          f_hlp(7 ,x_e,y_n,z  ) = f(7,x,y,z)
          f_hlp(8 ,x_e,y_s,z  ) = f(8,x,y,z)
          f_hlp(9 ,x_w,y_n,z  ) = f(9,x,y,z)
          f_hlp(10,x_w,y_s,z  ) = f(10,x,y,z)
          f_hlp(11,x_e,y  ,z_n) = f(11,x,y,z)
          f_hlp(12,x_w,y  ,z_n) = f(12,x,y,z)
          f_hlp(13,x_e,y  ,z_s) = f(13,x,y,z)
          f_hlp(14,x_w,y  ,z_s) = f(14,x,y,z)
          f_hlp(15,x  ,y_n,z_n) = f(15,x,y,z)
          f_hlp(16,x  ,y_n,z_s) = f(16,x,y,z)
          f_hlp(17,x  ,y_s,z_n) = f(17,x,y,z)
          f_hlp(18,x  ,y_s,z_s) = f(18,x,y,z)
   10   continue
   11  continue
   12 continue
c
      do 22 z = 1, lz
       do 21 y = 1, ly
         do 20 x = 1, lx
          do k = 1, 18
          f(k,x,y,z) = f_hlp(k,x,y,z)
          enddo
   20   continue
   21  continue
   22 continue
      return
      end


c--------------------------------------------------
      subroutine getuv(obst,u_x,u_y,u_z,rh, rho,fp, p,
     &  Force1, Force2,Force3)
      implicit none
      include "head.inc"
      integer x,y,z,obst(lx,ly,lz),ip,jp,l,k,xp,yp, xn,yn,zp,zn
      real*8  temp(lx,ly,lz),
     &  u_x(lx,ly,lz),u_y(lx,ly,lz),rho(lx,ly,lz),rh(lx,ly,lz),
     & fp(0:18,lx,ly,lz),  Force1(lx,ly,lz), fai(lx,ly,lz),
     & Force2(lx,ly,lz), p(lx,ly,lz), prho(lx,ly,lz),
     & Force3(lx,ly,lz),u_z(lx,ly,lz)

      do 10 z = 1, lz
         do 10 y = 1, ly
          do 10 x = 1, lx
         if(obst(x,y,z) .eq. 0) then
! Macroscopic variable \phi
            rh(x,y,z) = fp(0,x,y,z) +fp(1,x,y,z) +fp(2,x,y,z)
     &                 +fp(3,x,y,z) +fp(4,x,y,z) +fp(5,x,y,z)
     &                 +fp(6,x,y,z) +fp(7,x,y,z) +fp(8,x,y,z)
     &                 +fp(9,x,y,z) +fp(10,x,y,z) +fp(11,x,y,z)
     &                 +fp(12,x,y,z) +fp(13,x,y,z) +fp(14,x,y,z)
     &       +fp(15,x,y,z) +fp(16,x,y,z) +fp(17,x,y,z) +fp(18,x,y,z)
            endif
   10 continue

! From \phi to get density of fluid \rho
      call index_function( rh,rho)
c-----------------------------------------------

      do 11 z = 1, lz
         do 11 y = 1, ly
          do 11 x = 1, lx
! To apply the EOS for thermodynamic pressure

      if(obst(x,y,z) .eq. 0) then
! \fai is (p_t- \rh* c_s^2), where p_t is the thermodynamic pressure (EOS is used)
! In the EOS: a=12RT=12c_s^2, b=4. The values are in lattice units
      fai(x,y,z) =
     &  rh(x,y,z)*rh(x,y,z) *c_squ*( 4.d0 - 2.d0*rh(x,y,z) )
     %           /(1.d0- rh(x,y,z))**3.d0
     %          - 12.d0* c_squ* rh(x,y,z)* rh(x,y,z)

!  prho means: p-\rho* c_s^2 , where p is the hydrodynamic pressure
      prho(x,y,z) = p(x,y,z) - rho(x,y,z) *c_squ

      endif
   11 continue

c---------------------------------------------------

      do 15 z = 1, lz
         do 15 y = 1, ly
          do 15 x = 1, lx
          xp = x+1
          yp = y+1
          zp = z+1
          xn = x-1
          yn = y-1
          zn = z-1
! Periodic boundary is applied here
      if (xp.gt.lx )    xp = 1
      if (xn.lt.1 )     xn = lx
      if (yp.gt.ly )    yp = 1
      if (yn.lt.1 )     yn = ly
      if (zp.gt.lz )    zp = 1
      if (zn.lt.1 )     zn = lz

! The value of Laplacian operator applied to the variable \rho
      temp(x,y,z) =
     *  ( ( rh(xp,y,z) + rh(xn,y,z) + rh(x,yp,z)+ rh(x,yn,z)
     *       +rh(x,y,zp) + rh(x,y,zn) )*2.d0/6.d0
     *     +( rh(xp,yp,z )+ rh(xn,yn,z )
     *     +  rh(xp,yn,z )+ rh(xn,yp,z )
     *     +  rh(xp,y ,zp)+ rh(xn,y ,zn)
     *     +  rh(xp,y ,zn)+ rh(xn,y ,zp)
     *     +  rh(x ,yp,zp)+ rh(x ,yn,zn)
     *     +  rh(x ,yp,zn)+ rh(x ,yn,zp)      )/6.d0
     *     - 24.d0* rh(x,y,z)/6.d0
     *     )

   15 continue

c----------------------------------------
      do 20 z = 1, lz
         do 20 y = 1, ly
          do 20 x = 1, lx
          xp = x+1
          yp = y+1
          xn = x-1
          yn = y-1
          zp = z+1
          zn = z-1
! In x, y, z directions, periodic boundary conditions are applied
      if (xp.gt.lx )    xp = 1
      if (xn.lt.1 )     xn = lx
      if (yp.gt.ly )    yp = 1
      if (yn.lt.1 )     yn = ly
      if (zp.gt.lz )    zp = 1
      if (zn.lt.1 )     zn = lz

      Force1(x,y,z) = Kappa *rh(x,y,z) *
     *  ( (temp(xp,y ,z )- temp(xn,y ,z ) )/6.d0
     *     +(temp(xp,yp,z )- temp(xn,yn,z ) )/12.d0
     *     +(temp(xp,yn,z )- temp(xn,yp,z ) )/12.d0
     *     +(temp(xp,y ,zp)- temp(xn,y ,zn) )/12.d0
     *     +(temp(xp,y ,zn)- temp(xn,y ,zp) )/12.d0
     *     )

      Force2(x,y,z) = Kappa *rh(x,y,z) *
     &  ( (temp(x ,yp,z )- temp(x ,yn,z ) )/6.d0
     &     +(temp(xp,yp,z )- temp(xn,yn,z ) )/12.d0
     &     +(temp(xn,yp,z )- temp(xp,yn,z ) )/12.d0
     &     +(temp(x ,yp,zp)- temp(x ,yn,zn) )/12.d0
     &     +(temp(x ,yp,zn)- temp(x ,yn,zp) )/12.d0  )

      Force3(x,y,z) = Kappa *rh(x,y,z) *
     &  ( (temp(x ,y ,zp)- temp(x ,y ,zn) )/6.d0
     &     +(temp(xp,y ,zp)- temp(xn,y ,zn) )/12.d0
     &     +(temp(xn,y ,zp)- temp(xp,y ,zn) )/12.d0
     &     +(temp(x ,yp,zp)- temp(x ,yn,zn) )/12.d0
     &     +(temp(x ,yn,zp)- temp(x ,yp,zn) )/12.d0  )

!  Buoyant force is added here
      if(BUOYANCY .eq. 1 ) then
      Force3(x,y,z) = Force3(x,y,z)
     &     + gforce *( rho(x,y,z)- rho_h )
      endif

   20 continue
c----------------------------------------

! To get the gradient (first derivative) of \fai and  p-\rho* c_s^2
c----------------------------------------------
      do 30 z = 1, lz
         do 31 y = 1, ly
          do 32 x = 1, lx
      xp = x+1
      yp = y+1
      zp = z+1
      xn = x-1
      yn = y-1
      zn   = z-1
! In x, y, z directions, periodic boundary conditions are applied
      if (xp.gt.lx )    xp = 1
      if (xn.lt.1 )     xn = lx
      if (yp.gt.ly )    yp = 1
      if (yn.lt.1 )     yn = ly
      if (zp.gt.lz )    zp = 1
      if (zn.lt.1 )     zn = lz

! Gradient of \fai in x direction
      dfai_x(x,y,z) =
     *  ( (fai(xp,y ,z )- fai(xn,y ,z ) )/6.d0
     *     +(fai(xp,yp,z )- fai(xn,yn,z ) )/12.d0
     *     +(fai(xp,yn,z )- fai(xn,yp,z ) )/12.d0
     *     +(fai(xp,y ,zp)- fai(xn,y ,zn) )/12.d0
     *     +(fai(xp,y ,zn)- fai(xn,y ,zp) )/12.d0
     *     )

! Gradient of \fai in y direction
      dfai_y(x,y,z) =
     &  ( (fai(x ,yp,z )- fai(x ,yn,z ) )/6.d0
     &     +(fai(xp,yp,z )- fai(xn,yn,z ) )/12.d0
     &     +(fai(xn,yp,z )- fai(xp,yn,z ) )/12.d0
     &     +(fai(x ,yp,zp)- fai(x ,yn,zn) )/12.d0
     &     +(fai(x ,yp,zn)- fai(x ,yn,zp) )/12.d0  )

! Gradient of \fai in z direction
      dfai_z(x,y,z) =
     &  ( (fai(x ,y ,zp)- fai(x ,y ,zn) )/6.d0
     &     +(fai(xp,y ,zp)- fai(xn,y ,zn) )/12.d0
     &     +(fai(xn,y ,zp)- fai(xp,y ,zn) )/12.d0
     &     +(fai(x ,yp,zp)- fai(x ,yn,zn) )/12.d0
     &     +(fai(x ,yn,zp)- fai(x ,yp,zn) )/12.d0  )
c----------------------------------
      dprho_x(x,y,z) =
     *  ( (prho(xp,y ,z )- prho(xn,y ,z ) )/6.d0
     *     +(prho(xp,yp,z )- prho(xn,yn,z ) )/12.d0
     *     +(prho(xp,yn,z )- prho(xn,yp,z ) )/12.d0
     *     +(prho(xp,y ,zp)- prho(xn,y ,zn) )/12.d0
     *     +(prho(xp,y ,zn)- prho(xn,y ,zp) )/12.d0
     *     )

      dprho_y(x,y,z) =
     &  ( (prho(x ,yp,z )- prho(x ,yn,z ) )/6.d0
     &     +(prho(xp,yp,z )- prho(xn,yn,z ) )/12.d0
     &     +(prho(xn,yp,z )- prho(xp,yn,z ) )/12.d0
     &     +(prho(x ,yp,zp)- prho(x ,yn,zn) )/12.d0
     &     +(prho(x ,yp,zn)- prho(x ,yn,zp) )/12.d0  )

      dprho_z(x,y,z) =
     &  ( (prho(x ,y ,zp)- prho(x ,y ,zn) )/6.d0
     &     +(prho(xp,y ,zp)- prho(xn,y ,zn) )/12.d0
     &     +(prho(xn,y ,zp)- prho(xp,y ,zn) )/12.d0
     &     +(prho(x ,yp,zp)- prho(x ,yn,zn) )/12.d0
     &     +(prho(x ,yn,zp)- prho(x ,yp,zn) )/12.d0  )
   32   continue
   31  continue
   30 continue
      end


c---------------------------------------------------------------------
      subroutine geten(obst,u_x,u_y,u_z,rho,gp, p,Fx,Fy, Fz)
      include "head.inc"
      integer x,y,obst(lx,ly,lz)
      real*8  u_x(lx,ly,lz),u_y(lx,ly,lz),rho(lx,ly,lz),
     & gp(0:18,lx,ly,lz),u_z(lx,ly,lz), Fx(lx,ly,lz),
     & Fy(lx,ly,lz), Fz(lx,ly,lz),p(lx,ly,lz)

       do 12 z = 1, lz
        do 11 y = 1, ly
         do 10 x = 1, lx
! Initialization in each step
        u_x(x,y,z) = 0.d0
        u_y(x,y,z) = 0.d0
        u_z(x,y,z) = 0.d0

       if(obst(x,y,z) .eq. 0) then


        u_x(x,y,z)=(gp(1,x,y,z)+ gp(7,x,y,z)+ gp(8,x,y,z) +
     &             gp(11,x,y,z) + gp(13,x,y,z)
     &           -(gp(2,x,y,z) + gp(9,x,y,z) + gp(10,x,y,z)+
     &           gp(12,x,y,z) + gp(14,x,y,z) ))
     &       + 0.5d0 * c_squ *( Fx(x,y,z) +0.d0 )

c------------------------------------------------------
        u_y(x,y,z) = (gp(3,x,y,z) + gp(7,x,y,z) + gp(9,x,y,z) +
     &             gp(15,x,y,z) + gp(16,x,y,z)
     &           -(gp(4,x,y,z) + gp(8,x,y,z) + gp(10,x,y,z) +
     &             gp(17,x,y,z) + gp(18,x,y,z) ))
     &       + 0.5d0 * c_squ *( Fy(x,y,z) +0.d0 )

c------------------------------------------------------
        u_z(x,y,z)= (gp(5,x,y,z) + gp(11,x,y,z) + gp(12,x,y,z)+
     &             gp(15,x,y,z) + gp(17,x,y,z)
     &           -(gp(6,x,y,z) + gp(13,x,y,z) + gp(14,x,y,z) +
     &             gp(16,x,y,z) + gp(18,x,y,z) ))
     &       + 0.5d0 * c_squ *( Fz(x,y,z) +0.d0 )

        u_x(x,y,z) = u_x(x,y,z)/rho(x,y,z)/c_squ
        u_y(x,y,z) = u_y(x,y,z)/rho(x,y,z)/c_squ
        u_z(x,y,z) = u_z(x,y,z)/rho(x,y,z)/c_squ

            p(x,y,z) = gp(0,x,y,z) +gp(1,x,y,z) +gp(2,x,y,z)
     &                +gp(3,x,y,z) +gp(4,x,y,z) +gp(5,x,y,z)
     &                +gp(6,x,y,z) +gp(7,x,y,z) +gp(8,x,y,z)
     &                +gp(9,x,y,z) +gp(10,x,y,z) +gp(11,x,y,z)
     &                +gp(12,x,y,z) +gp(13,x,y,z) +gp(14,x,y,z)
     &  +gp(15,x,y,z) +gp(16,x,y,z) +gp(17,x,y,z) +gp(18,x,y,z)

              p(x,y,z) = p(x,y,z) - 0.5d0 * ( u_x(x,y,z)* dprho_x(x,y,z)
     %          + u_y(x,y,z)* dprho_y(x,y,z)
     %                  + u_z(x,y,z)* dprho_z(x,y,z)  )

      endif

  10    continue
   11  continue
   12 continue
      end


c------------------------------------------------------
      subroutine collision(obst,ux,uy,uz,rh,
     &  rho,p,ff,gp, Fx,Fy,Fz)
      implicit none
      include "head.inc"
      integer  l,obst(lx,ly,lz)
      real*8 ux(lx,ly,lz),uy(lx,ly,lz),ff(0:18,lx,ly,lz),rho(lx,ly,lz)
      real*8 Fx(lx,ly,lz),Fy(lx,ly,lz), rh(lx,ly,lz)
      real*8 Fz(lx,ly,lz),uz(lx,ly,lz),p(lx,ly,lz),gp(0:18,lx,ly,lz)

      integer  x,y,z,k
      real*8  u_n(0:18),fequ(0:18),gequ(0:18),u_squ,temp(18),temp2(18)
      real*8  noteq(0:18),gneq(0:18)

      do 4 z = 1, lz
         do 5 y = 1, ly
          do 6 x = 1, lx
! Bounce back for solid nodes
           if(BD .eq. 1 .and. obst(x,y,z) .eq. 1) then
             do k=1, 18 
             temp(k)   = ff(k,x,y,z)
             temp2(k)   = gp(k,x,y,z)
             enddo
             do k=1, 18
             ff(opp(k),x,y,z) = temp(k)
             gp(opp(k),x,y,z) = temp2(k)
	     enddo		
        endif

        if(obst(x,y,z) .eq. 0) then

      call g_eq( p(x,y,z),rho(x,y,z),ux(x,y,z),uy(x,y,z),uz(x,y,z),gequ)
      call f_eq(rh(x,y,z),ux(x,y,z),uy(x,y,z),uz(x,y,z),fequ)

         do 60 k = 0,18
! Non-equilibrium distribution function
         noteq(k)= ff(k,x,y,z) - rh(x,y,z)*fequ(k)

             ff(k,x,y,z) = ff(k,x,y,z) -1.d0/tau_f * noteq(k)
! The source term in the LBE for f_i
     &       - (tau_f-0.5d0)/tau_f*
     &      ( (xc(k)-ux(x,y,z))* dfai_x(x,y,z)
     &       +(yc(k)-uy(x,y,z))* dfai_y(x,y,z)
     &           +(zc(k)-uz(x,y,z))* dfai_z(x,y,z) )*fequ(k)/c_squ
c----------------------------------------------------------------
! Non-equilibrium distribution function
         gneq(k)= gp(k,x,y,z) - gequ(k)

         gp(k,x,y,z) = gp(k,x,y,z) -1.d0/tau_g *gneq(k)
! The source term in the LBE for g_i
     &        + (tau_g-0.5d0)/tau_g *
     *      (
     &      (  (xc(k)-ux(x,y,z))*  Fx(x,y,z)
     &        +(yc(k)-uy(x,y,z))* (Fy(x,y,z) )
     &        +(zc(k)-uz(x,y,z))* (Fz(x,y,z) )
     &        )  *fequ(k)

     &          - (xc(k)-ux(x,y,z))*( fequ(k)-t_k(k) )*dprho_x(x,y,z)
     &          - (yc(k)-uy(x,y,z))*( fequ(k)-t_k(k) )*dprho_y(x,y,z)
     &          - (zc(k)-uz(x,y,z))*( fequ(k)-t_k(k) )*dprho_z(x,y,z)
     *      )
   60   continue

      endif

    6   continue
    5  continue
    4 continue

      return
      end
