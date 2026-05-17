       subroutine stream(obst,ff)
      implicit none
      include "head.inc"
      integer  k,obst(lx,ly)
      real*8 ff(0:8,lx,ly),f_hlp(0:8,lx,ly)
      integer  x,y,x_e,x_w,y_n,y_s,l,m,n

      do 10 y = 1, ly
         do 10 x = 1, lx

! Set periodic boundary conditions
          y_n = mod(y,ly) + 1
          x_e = mod(x,lx) + 1
          y_s = ly - mod(ly + 1 - y, ly)
          x_w = lx - mod(lx + 1 - x, lx)

          f_hlp(1,x_e,y  ) = ff(1,x,y)
          f_hlp(2,x  ,y_n) = ff(2,x,y)
          f_hlp(3,x_w,y  ) = ff(3,x,y)
          f_hlp(4,x  ,y_s) = ff(4,x,y)
          f_hlp(5,x_e,y_n) = ff(5,x,y)
          f_hlp(6,x_w,y_n) = ff(6,x,y)
          f_hlp(7,x_w,y_s) = ff(7,x,y)
          f_hlp(8,x_e,y_s) = ff(8,x,y)

   10 continue

      do 20 y = 1, ly
         do 20 x = 1, lx
           do 20 k = 1, 8
           ff(k,x,y) = f_hlp(k,x,y)
   20 continue

      end
c-------------------------------------------------------
      subroutine getuv(obst,u_x,u_y,rho,ff)
      include "head.inc"
      integer x,y,obst(lx,ly)
      real*8  u_x(lx,ly),u_y(lx,ly),rho(lx,ly),
     & ff(0:8,lx,ly)

      do 10 y = 1, ly
          do 10 x = 1, lx

      if(obst(x,y) .eq. 0) then
              rho(x,y) = ff(0,x,y) +ff(1,x,y) +ff(2,x,y)
     &                 + ff(3,x,y) +ff(4,x,y) +ff(5,x,y)
     &                 + ff(6,x,y) +ff(7,x,y) +ff(8,x,y)
      endif

      if(obst(x,y) .eq. 0) then
            u_x(x,y) = cc* (ff(1,x,y) + ff(5,x,y) + ff(8,x,y)
     &         -(ff(3,x,y) + ff(6,x,y) + ff(7,x,y))) /rho(x,y)

            u_y(x,y) = cc* (ff(2,x,y) + ff(5,x,y) + ff(6,x,y)
     &         -(ff(4,x,y) + ff(7,x,y) + ff(8,x,y))) /rho(x,y)
      else
      u_x(x,y) = 0.d0
      u_y(x,y) = 0.d0

      endif

  10  continue

      end

c-------------------------------------------------------

c The BGK collision
      subroutine collision(tauc,obst,u_x,u_y,rho,ff, pre)
      implicit none
      include "head.inc"
      integer  l,obst(lx,ly)
      real*8   u_x(lx,ly),u_y(lx,ly),ff(0:8,lx,ly),rho(lx,ly)
      real*8 Fx(lx,ly),Fy(lx,ly), pre(lx,ly)
      integer  x,y,k,ip,jp
      real*8  u_n(0:8),fequ(0:8, lx,ly),u_squ,temp,tauc,ux,uy

        call get_feq(rho, u_x, u_y, fequ, pre)


      do 5 y = 1, ly
       do 5 x = 1, lx
          do 10 k = 0,8
          ff(k,x,y) = fequ(k,x,y)
     &      + (1.d0-1.d0/tauc)*( ff(k,x,y) - fequ(k,x,y) )
   10 continue
   5  continue
      end
