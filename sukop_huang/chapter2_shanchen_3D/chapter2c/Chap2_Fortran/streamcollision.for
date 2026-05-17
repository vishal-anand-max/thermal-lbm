      subroutine stream(obst,ff)
      implicit none
      include "head.inc"
      integer  k,obst(lx,ly,lz)
      real*8 ff(0:18,lx,ly,lz),f_hlp(0:18,lx,ly,lz)
      integer  x,y,z,x_e,x_w,y_n,y_s,z_n,z_s

      do 12 z = 1, lz
       do 11 y = 1, ly
        do 10 x = 1, lx
!
          z_n = mod(z,lz) + 1
          y_n = mod(y,ly) + 1
          x_e = mod(x,lx) + 1

          z_s = lz - mod(lz + 1 - z, lz)
          y_s = ly - mod(ly + 1 - y, ly)
          x_w = lx - mod(lx + 1 - x, lx)

c......... Propagation
          f_hlp(1 ,x_e,y  ,z  ) = ff(1,x,y,z)
          f_hlp(2 ,x_w,y  ,z  ) = ff(2,x,y,z)
          f_hlp(3 ,x  ,y_n,z  ) = ff(3,x,y,z)
          f_hlp(4 ,x  ,y_s,z  ) = ff(4,x,y,z)
          f_hlp(5 ,x  ,y  ,z_n) = ff(5,x,y,z)
          f_hlp(6 ,x  ,y  ,z_s) = ff(6,x,y,z)
          f_hlp(7 ,x_e,y_n,z  ) = ff(7,x,y,z)
          f_hlp(8 ,x_e,y_s,z  ) = ff(8,x,y,z)
          f_hlp(9 ,x_w,y_n,z  ) = ff(9,x,y,z)
          f_hlp(10,x_w,y_s,z  ) = ff(10,x,y,z)
          f_hlp(11,x_e,y  ,z_n) = ff(11,x,y,z)
          f_hlp(12,x_w,y  ,z_n) = ff(12,x,y,z)
          f_hlp(13,x_e,y  ,z_s) = ff(13,x,y,z)
          f_hlp(14,x_w,y  ,z_s) = ff(14,x,y,z)
          f_hlp(15,x  ,y_n,z_n) = ff(15,x,y,z)
          f_hlp(16,x  ,y_n,z_s) = ff(16,x,y,z)
          f_hlp(17,x  ,y_s,z_n) = ff(17,x,y,z)
          f_hlp(18,x  ,y_s,z_s) = ff(18,x,y,z)
   10   continue
   11  continue
   12 continue
c---------------------Update distribution function
      do 22 z = 1, lz
        do 21 y = 1, ly
           do 20 x = 1, lx
            do k =1, 18
                ff(k,x,y,z) = f_hlp(k,x,y,z)
            enddo
   20     continue
   21   continue
   22 continue

      return
      end


c-----------------------------------------
      subroutine getuv(obst,u_x,u_y,u_z,rho,ff)
      include "head.inc"
      integer x,y,obst(lx,ly,lz)
      real*8  u_x(lx,ly,lz),u_y(lx,ly,lz),rho(lx,ly,lz),
     & ff(0:18,lx,ly,lz),u_z(lx,ly,lz)

      do 12 z = 1, lz
       do 11 y = 1, ly
        do 10 x = 1, lx
          rho(x,y,z) = 0.d0

      if(obst(x,y,z) .eq. 0 ) then
            do 5 k = 0 ,18
              rho(x,y,z) = rho(x,y,z) + ff(k,x,y,z)
    5    continue
c----------------------
       if(rho(x,y,z) .ne. 0.d0) then

        u_x(x,y,z)=(ff(1,x,y,z)+ ff(7,x,y,z)+ ff(8,x,y,z) +
     &             ff(11,x,y,z) + ff(13,x,y,z)
     &           -(ff(2,x,y,z) + ff(9,x,y,z) + ff(10,x,y,z)+
     &           ff(12,x,y,z) + ff(14,x,y,z) ))/rho(x,y,z)

        u_y(x,y,z) = (ff(3,x,y,z) + ff(7,x,y,z) + ff(9,x,y,z) +
     &             ff(15,x,y,z) + ff(16,x,y,z)
     &           -(ff(4,x,y,z) + ff(8,x,y,z) + ff(10,x,y,z) +
     &             ff(17,x,y,z) + ff(18,x,y,z) )) /rho(x,y,z)

        u_z(x,y,z)= (ff(5,x,y,z) + ff(11,x,y,z) + ff(12,x,y,z)+
     &             ff(15,x,y,z) + ff(17,x,y,z)
     &           -(ff(6,x,y,z) + ff(13,x,y,z) + ff(14,x,y,z) +
     &             ff(16,x,y,z) + ff(18,x,y,z) )) /rho(x,y,z)
       endif

      endif

  10    continue
   11  continue
   12 continue
      end
            subroutine calcu_upr(obst,u_x,u_y,u_z,
     &  Fx,Fy,Fz,Sx,Sy,Sz,rho,upx,upy,upz)
      implicit none
      include "head.inc"
      integer x,y,z ,obst(lx,ly,lz)
      real*8  u_x(lx,ly,lz),u_y(lx,ly,lz),rho(lx,ly,lz),
     % upx(lx,ly,lz), upy(lx,ly,lz), upz(lx,ly,lz),
     & u_z(lx,ly,lz), Fx(lx,ly,lz), Fy(lx,ly,lz), Fz(lx,ly,lz),
     &  Sx(lx,ly,lz), Sy(lx,ly,lz), Sz(lx,ly,lz)

      do  9 z = 1, lz
       do 10 y = 1, ly
        do 11 x = 1, lx

      if(obst(x,y,z) .eq. 0) then

      upx(x,y,z) = u_x(x,y,z) + (Fx(x,y,z)+Sx(x,y,z))/2.d0/ rho(x,y,z)
      upy(x,y,z) = u_y(x,y,z) + (Fy(x,y,z)+Sy(x,y,z))/2.d0/ rho(x,y,z)
      upz(x,y,z) = u_z(x,y,z) + (Fz(x,y,z)+Sz(x,y,z))/2.d0/ rho(x,y,z)

      else
      upx(x,y,z) = u_x(x,y,z)
      upy(x,y,z) = u_y(x,y,z)
      upz(x,y,z) = u_z(x,y,z)
      endif

  11    continue
  10   continue
   9  continue
      end
c-----------------------------------------
      subroutine collision(tauc,obst,u_x,u_y,u_z,
     &  rho,ff,Fx,Fy,Fz, Sx, Sy, Sz)
!
      implicit none
      include "head.inc"
      integer  l,obst(lx,ly,lz)
      real*8 u_x(lx,ly,lz),u_y(lx,ly,lz),ff(0:18,lx,ly,lz),rho(lx,ly,lz)
      real*8 Fx(lx,ly,lz),Fy(lx,ly,lz), Sx(lx,ly,lz), Sy(lx,ly,lz)
      real*8 Fz(lx,ly,lz),u_z(lx,ly,lz),Sz(lx,ly,lz),temp(18)

      integer  x,y,z,k
      real*8  u_n(0:18),fequ(0:18),fequ2(0:18),u_squ,tauc,ux,uy,uz

      do 4 z = 1, lz
       do 5 y = 1, ly
        do 6 x = 1, lx
         if(obst(x,y,z) .eq. 1) then
          do k =1, 18
           temp(k) =ff(k,x,y,z)
          enddo
          do k =1, 18
           ff(opp(k),x,y,z) = temp(k)
          enddo
         endif

         if(obst(x,y,z) .eq. 0) then
          ux = u_x(x,y,z) +tauc * ( Fx(x,y,z)+Sx(x,y,z) ) / rho(x,y,z)
          uy = u_y(x,y,z) +tauc * ( Fy(x,y,z)+Sy(x,y,z) ) / rho(x,y,z)
          uz = u_z(x,y,z) +tauc * ( Fz(x,y,z)+Sz(x,y,z) ) / rho(x,y,z)

            u_squ = ux * ux + uy * uy +uz *uz

        do 10 k = 0,18
c...........Equillibrium distribution function
            u_n(k)  = xc(k)*ux + yc(k)*uy + zc(k)*uz
            fequ(k) = t_k(k)* rho(x,y,z) * ( cc*u_n(k) / c_squ
     &               + (u_n(k)*cc) *(u_n(k)*cc) / (2.d0 * c_squ *c_squ)
     &               - u_squ / (2.d0 * c_squ)) + t_k(k) * rho(x,y,z)

c...........Collision step
      ff(k,x,y,z) = fequ(k) + (1.d0-1.d0/tauc)*( ff(k,x,y,z) - fequ(k) )
   10   continue

        endif
   6    continue
   5   continue
   4  continue
      end
c-----------------------------------------
      subroutine getf_equ(rh,u,v,w,f_equ)
      include 'head.inc'
      real*8 rh, u,v,w,u_squ, f_equ(0:18),u_n(0:18)

          u_squ =u*u +v*v +w*w

      do 10 i =0,18
        u_n(i) = u *xc(i) +v *yc(i)+ w *zc(i)
            f_equ(i) = t_k(i) * rh *( u_n(i)/c_squ
     &               + u_n(i) *u_n(i) / (2.d0 * c_squ *c_squ)
     &               - u_squ / (2.d0 * c_squ)) + t_k(i) * rh
  10  continue
      end


