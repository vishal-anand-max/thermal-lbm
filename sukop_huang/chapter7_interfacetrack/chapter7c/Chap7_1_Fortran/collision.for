
      subroutine stream(lx,ly,obst,fp,gp)
      implicit none
      integer  k,lx,ly,obst(lx,ly)
      real*8 fp(0:8,lx,ly), gp(0:8,lx,ly), fx(lx), fy(ly)
      real*8 gx(lx), gy(ly)
      integer  x,y,x_e,x_w,y_n,y_s,l,m,n,xi,yi
      real*8 f_hlp(0:8,lx,ly)

      do 10 y = 1, ly
         do 10 x = 1, lx
          y_n = mod(y,ly) + 1
          x_e = mod(x,lx) + 1
          y_s = ly - mod(ly + 1 - y, ly)
          x_w = lx - mod(lx + 1 - x, lx)

          f_hlp(1,x_e,y  ) = fp(1,x,y)
          f_hlp(2,x  ,y_n) = fp(2,x,y)
          f_hlp(3,x_w,y  ) = fp(3,x,y)
          f_hlp(4,x  ,y_s) = fp(4,x,y)
          f_hlp(5,x_e,y_n) = fp(5,x,y)
          f_hlp(6,x_w,y_n) = fp(6,x,y)
          f_hlp(7,x_w,y_s) = fp(7,x,y)
          f_hlp(8,x_e,y_s) = fp(8,x,y)

   10 continue

      do 20 y = 1, ly
         do 20 x = 1, lx
          do 20 k = 1, 8
            fp(k,x,y) = f_hlp(k,x,y)
   20 continue

c--------------
      do 30 y = 1, ly
         do 30 x = 1, lx
          y_n = mod(y,ly) + 1
          x_e = mod(x,lx) + 1
          y_s = ly - mod(ly + 1 - y, ly)
          x_w = lx - mod(lx + 1 - x, lx)

          f_hlp(1,x_e,y  ) = gp(1,x,y)
          f_hlp(2,x  ,y_n) = gp(2,x,y)
          f_hlp(3,x_w,y  ) = gp(3,x,y)
          f_hlp(4,x  ,y_s) = gp(4,x,y)
          f_hlp(5,x_e,y_n) = gp(5,x,y)
          f_hlp(6,x_w,y_n) = gp(6,x,y)
          f_hlp(7,x_w,y_s) = gp(7,x,y)
          f_hlp(8,x_e,y_s) = gp(8,x,y)

   30 continue

      do 40 y = 1, ly
         do 40 x = 1, lx
          do 40 k = 1, 8
            gp(k,x,y) = f_hlp(k,x,y)
   40 continue

      return
      end

c----------------------------------------------------------------
      subroutine getuv(lx,ly,obst,u_x,u_y,rh, rho,fp, p, Force1, Force2)
      implicit none
      include "head.inc"
      integer x,y,lx,ly,obst(lx,ly),ip,jp,l,k,xp,yp, xn,yn
      real*8  temp(lx,ly),
     &  u_x(lx,ly),u_y(lx,ly),rho(lx,ly),rh(lx,ly),
     & fp(0:8,lx,ly), Force1(lx,ly), fai(lx,ly),
     & Force2(lx,ly), p(lx,ly), prho(lx,ly)

      do 10 y = 1, ly
      do 10 x = 1, lx
      if(obst(x,y) .eq. 0) then
          rh(x,y) = fp(0,x,y) +fp(1,x,y) +fp(2,x,y) +fp(3,x,y)
     &         +fp(4,x,y) +fp(5,x,y) +fp(6,x,y) +fp(7,x,y) +fp(8,x,y)
      endif

   10 continue

      call index_function(lx,ly, rh,rho)
c-----------------------------------------------
      do 11 y = 1, ly
       do 11 x = 1, lx
        if(obst(x,y) .eq. 0) then
         fai(x,y) =
     &    rh(x,y)*rh(x,y) *c_squ*( 4.d0 - 2.d0*rh(x,y) )
     %       /(  (1.d0- rh(x,y))*(1.d0- rh(x,y))*(1.d0- rh(x,y))  )
     %  - 12.d0* c_squ* rh(x,y)* rh(x,y)

        prho(x,y) = p(x,y) - rho(x,y) *c_squ
       endif
   11 continue

ccccc=================================
      do 15 y = 1, ly
       do 15 x = 1, lx

        xp = x+1
        yp = y+1
        xn = x-1
        yn = y-1
        if (xp.gt.lx )  xp = 1
        if (xn.lt.1 )   xn = lx
        if (yp.gt.ly )  yp = 1
        if (yn.lt.1 )   yn = ly

        temp(x,y) =
     *  ( ( rh(xp,y) + rh(xn,y) + rh(x,yp)+ rh(x,yn) )*4.d0/6.d0
     *     +( rh(xp,yp)+ rh(xn,yn)
     *     +  rh(xp,yn)+ rh(xn,yp) )/6.d0
     *     - 20.d0* rh(x,y)/6.d0
     *     )
   15 continue

      do 14 x = 1, lx
       temp(x,1) = temp(x,2)
       temp(x,ly) = temp(x,ly-1)
   14 continue

c----------------------------------------
      do 20 y = 1, ly
      do 20 x = 1, lx
        xp = x+1
        yp = y+1
        xn = x-1
        yn = y-1
        if (xp.gt.lx )  xp = 1
        if (xn.lt.1 )   xn = lx
        if (yp.gt.ly )  yp = 1
        if (yn.lt.1 )   yn = ly

        Force1(x,y) = Kappa *rh(x,y) *
     *  ( (temp(xp,y)- temp(xn,y) )/3.d0
     *     +(temp(xp,yp)- temp(xn,yn) )/12.d0
     *     +(temp(xp,yn)- temp(xn,yp) )/12.d0
     *     )

        Force2(x,y) = Kappa *rh(x,y) *
     &  ( (temp(x,yp)- temp(x,yn) )/3.d0
     &     +(temp(xp,yp)- temp(xn,yn) )/12.d0
     &     +(temp(xn,yp)- temp(xp,yn) )/12.d0  )
     &     + rho(x,y)* gforce
c---------! the gravity force !!!!
   20 continue

      do 30 y = 1, ly
       do 30 x = 1, lx

        xp = x+1
        yp = y+1
        xn = x-1
        yn = y-1
        if (xp.gt.lx )  xp = 1
        if (xn.lt.1 )   xn = lx
        if (yp.gt.ly )  yp = 1
        if (yn.lt.1 )   yn = ly

        dfai_x(x,y) =
     *  ( (fai(xp,y)- fai(xn,y) )/3.d0
     *     +(fai(xp,yp)- fai(xn,yn) )/12.d0
     *     +(fai(xp,yn)- fai(xn,yp) )/12.d0
     *     )
        dfai_y(x,y) =
     *  ( (fai(x,yp)- fai(x,yn) )/3.d0
     *     +(fai(xp,yp)- fai(xn,yn) )/12.d0
     *     +(fai(xn,yp)- fai(xp,yn) )/12.d0
     *     )

       dprho_x(x,y) =
     *  ( (prho(xp,y)- prho(xn,y) )/3.d0
     *     +(prho(xp,yp)- prho(xn,yn) )/12.d0
     *     +(prho(xp,yn)- prho(xn,yp) )/12.d0
     *     )
        dprho_y(x,y) =
     *  ( (prho(x,yp)- prho(x,yn) )/3.d0
     *     +(prho(xp,yp)- prho(xn,yn) )/12.d0
     *     +(prho(xn,yp)- prho(xp,yn) )/12.d0
     *     )

   30 continue

! Derivatives in the nodes most near to the bottom and upper boundary, e.g., 
! dfai_x(x,2), dfai_y(x,2), ... are extrapolated from inner fluid nodes. 
      if(BD .eq. 1) then
        do 13 x = 1, lx

	dfai_x(x,2) = dfai_x(x,3)
	dfai_y(x,2) = dfai_y(x,3)
	dprho_x(x,2) = dprho_x(x,3)
	dprho_y(x,2) = dprho_y(x,3)

	dfai_x(x,ly-1) = dfai_x(x,ly-2)
	dfai_y(x,ly-1) = dfai_y(x,ly-2)
	dprho_x(x,ly-1) = dprho_x(x,ly-2)
	dprho_y(x,ly-1) = dprho_y(x,ly-2)

   13   continue
      endif

c---------------------------------------------------------------------
c It is noted that the calculations of  dfai_x, dfai_y, dprho_x, dprho_y
c for the lattice node in bottom and upper boundaries (walls) may be
c not necessary because in the collision step,  dfai_x, dfai_y, dprho_x, 
c dprho_y on the wall nodes are not used at all. ONLY bounce back is 
c necessary to implement for the wall nodes.

      end
c----------------------------------------------------------------
      subroutine geten(lx,ly,obst,u_x,u_y,rho,gp, p, Force1,Force2)
      implicit none
      include "head.inc"
      integer x,y,lx,ly,obst(lx,ly),ip,jp,k,l
      real*8  u_x(lx,ly),u_y(lx,ly),rho(lx,ly),
     & gp(0:8,lx,ly), p(lx,ly),Force1(lx,ly),
     & Force2(lx,ly),tmp

      do 10 y = 1, ly
       do 10 x = 1, lx
       if(obst(x,y) .eq. 0) then

        u_x(x,y)= (gp(1,x,y) + gp(5,x,y) + gp(8,x,y)
     &         -(gp(3,x,y) + gp(6,x,y) + gp(7,x,y)) )
     &       + 0.5d0 * c_squ *( Force1(x,y) +0.d0 )

        u_y(x,y)= (gp(2,x,y) + gp(5,x,y) + gp(6,x,y)
     &         -(gp(4,x,y) + gp(7,x,y) + gp(8,x,y)) )
     &         + 0.5d0 * c_squ *( Force2(x,y)  )

         u_x(x,y) = u_x(x,y)/rho(x,y)/c_squ
         u_y(x,y) = u_y(x,y)/rho(x,y)/c_squ

        p(x,y)=gp(0,x,y) +gp(1,x,y) +gp(2,x,y)
     &      +gp(3,x,y) +gp(4,x,y) +gp(5,x,y)
     &      +gp(6,x,y) +gp(7,x,y) +gp(8,x,y)
     %  + 0.5d0 * ( u_x(x,y)* (-dprho_x(x,y))
     %            + u_y(x,y)* (-dprho_y(x,y)) )
      endif
  10  continue

      end
c----------------------------------------------------------------
      subroutine collision(lx,ly,obst,u_x,u_y,rh,rho,p,fp,gp,
     &  Force1,Force2)
      implicit none
      include "head.inc"
      integer  l,lx,ly,obst(lx,ly),x,y,k,ip,jp
      real*8   u_x(lx,ly),u_y(lx,ly)
     &, fp(0:8,lx,ly),rho(lx,ly),Force1(lx,ly),Force2(lx,ly),
     &  gp(0:8,lx,ly), p(lx,ly),rh(lx,ly),
     &  gneq(0:8),fneq(0:8), gequ(0:8), temp(8),temp2(8)
      real*8  u_n(0:8),fequ(0:8),u_squ, div
c-----------------------------------
      do 5 y = 1, ly
      do 5 x = 1, lx

      if(obst(x,y) .eq. 0) then

      call f_eq(  rh(x,y), u_x(x,y), u_y(x,y), fequ)
      call g_eq(  p(x,y), rho(x,y), fequ, gequ)

        do 60 k = 0,8

          fneq(k)= fp(k,x,y) - rh(x,y)*fequ(k)

            fp(k,x,y) = fp(k,x,y) -1.d0/tau_f * fneq(k)
     &       + (tau_f-0.5d0)/tau_f*
     &       ( (xv(k)-u_x(x,y))* (-dfai_x(x,y))
     &        +(yv(k)-u_y(x,y))* (-dfai_y(x,y)) )*fequ(k)/c_squ
c----------------------------------------------------------------
         gneq(k)= gp(k,x,y) - gequ(k)

         gp(k,x,y) = gp(k,x,y) -1.d0/tau_g *gneq(k)
     &        + (tau_g-0.5d0)/tau_g *
     *      (
     &          (  (xv(k)-u_x(x,y))*  Force1(x,y)
     &        +(yv(k)-u_y(x,y))*  Force2(x,y)
     &        )  *fequ(k)
     &        + (xv(k)-u_x(x,y))*( fequ(k)-t_k(k) )*(-dprho_x(x,y))
     &        + (yv(k)-u_y(x,y))*( fequ(k)-t_k(k) )*(-dprho_y(x,y))

     *      )
   60 continue

      endif

c-----------------
      if( BD .eq. 1 .and. obst(x,y) .eq. 1) then
      do k = 1, 8
      temp(k)   = fp(k,x,y)
      temp2(k)   = gp(k,x,y)
      enddo

      do k = 1, 8
      fp(opp(k),x,y)= temp(k)
      gp(opp(k),x,y)= temp2(k)
      enddo

      endif
c----------------

   5  continue
c--------------------------------------
      end
c--------------------------------------------------
       subroutine comp_rey(lx,ly,time,BT)
! Here is output some parameter into a file to keep a record.
       implicit none
       include "head.inc"
       integer  time,lx,ly
       real*8  rey
       character BT*10

       open(44,file='./out/parameter.dat')
       write(44,'("He_Chen_Doolen multiphase model")')
       write(44,'("Nwri= ", i5)')  Nwri

       write(44,'("gforce: ", f15.8)')  gforce
       write(44,'("viscosity = ",f8.5,2X,"tau_f=",f17.8,2X,
     &   2X,"tau_g=",f17.8) ') c_squ*(tau_f-0.5), tau_f,tau_g
       write(44,'("@@@ mesh size= ",i5,"X",i5)')  lx,  ly
       write(44,'("character U = ",f10.5)')  sqrt(-gforce*h_ref)
       write(44,'("Re = ",f10.5)')
     &   h_ref*sqrt(-gforce*h_ref)/c_squ/(tau_f-0.5)
       write(44,'("Kappa = ",f10.5)')  Kappa

       write(44,'(2X,"rho_h=",f17.7,2X,
     &   2X,"rho_l=",f17.7) ') rho_h, rho_l
       write(44,'(2X,"psi_max=",f17.7,2X,
     &   2X,"psi_min=",f17.7) ') psi_max, psi_min
      close(44)
      return
      end

