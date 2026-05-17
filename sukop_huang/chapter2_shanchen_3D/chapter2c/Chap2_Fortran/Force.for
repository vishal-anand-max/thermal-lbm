       subroutine calcu_Fxy(obst,rho,Fx,Fy,Fz,Sx,Sy,Sz,p)
       implicit none
       include "head.inc"
       integer x,y,z,obst(lx,ly,lz),yn,yp,xn,xp,zp,zn, i,j,k
       real*8 Fx(lx,ly,lz),Fy(lx,ly,lz),Fz(lx,ly,lz)
     & ,psx(lx,ly,lz), sum_x, sum_y, sum_z, psx_w
     & ,rho(lx,ly,lz), Sx(lx,ly,lz), Sy(lx,ly,lz), Sz(lx,ly,lz)
     % , Fztemp, R,a,b, Tc, TT, alfa, omega, G1,p(lx,ly,lz)

!  Parameters in YUAN C-S EOS
        R = 1.0d0
        b = 4.d0
        a = 1.d0
        Tc = 0.3773d0*a/(b*R)
        TT= TT0 *Tc

        do 4 k = 1,lz
         do 5 j = 1,ly
          do 6 i = 1,lx
        if (obst(i,j,k ) .eq. 0 .and. rho(i,j,k).ne. 0.d0)  then

         if( (R*TT*
     &      (1.d0+(4.d0* rho(i,j,k)-2.d0* rho(i,j,k)* rho(i,j,k)
     &     )/(1.d0- rho(i,j,k))**3    )
     %       -a* rho(i,j,k) -1.d0/3.d0) .gt. 0.) then
           G1= 1.d0/3.d0
         else
           G1= -1.d0/3.d0
         endif

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        psx(i,j,k) = sqrt( 6.d0* rho(i,j,k) * ( R*TT*
     &      (1.d0+ (4.d0* rho(i,j,k)-2.d0*rho(i,j,k)*rho(i,j,k) )
     &  /(1.d0-rho(i,j,k))**3  )
     &       -a* rho(i,j,k) -1.d0/3.d0)
     &       /G1 )
c  Yuan C-S EOS
        p(i,j,k) = rho(i,j,k)/3.d0 + G1/6.d0 * psx(i,j,k) *psx(i,j,k)
       endif
    6 continue
    5 continue
    4 continue


      psx_w = sqrt( 6.d0* rho_w * ( R*TT*
     &      (1.d0+ (4.d0* rho_w-2.d0*rho_w * rho_w )
     &       /(1.d0- rho_w)**3  )
     &       -a* rho_w -1.d0/3.d0)
     &       /G1 )

      do 30 z = 1,lz
       do 20 y = 1,ly
        do 10 x = 1,lx
c.........interaction between neighbouring with periodic boundaries
            Fx(x,y,z) =0.d0
            Fy(x,y,z) =0.d0
            Fz(x,y,z) =0.d0

      if (obst(x,y,z) .eq. 0)  then

        sum_x = 0.d0
        sum_y = 0.d0
        sum_z = 0.d0

      do 11  k =1, 18
       xp=x+ex(k)
       yp=y+ey(k)
       zp=z+ez(k)
       if(xp .lt. 1)  xp = lx
       if(xp .gt. lx) xp =1
       if(yp .lt. 1)  yp = ly
       if(yp .gt. ly) yp =1
       if(zp .lt. 1)  zp = lz
       if(zp .gt. lz) zp =1

        if (obst(xp,yp,zp) .eq. 1)  then
! Interact with solid nodes (obst=1)
         sum_x = sum_x + t_k(k)*xc(k)
         sum_y = sum_y + t_k(k)*yc(k)
         sum_z = sum_z + t_k(k)*zc(k)
      else
! Interact with fluid nodes (obst=0)
         Fx(x,y,z)=Fx(x,y,z) +t_k(k)*xc(k)* psx(xp,yp,zp)
         Fy(x,y,z)=Fy(x,y,z) +t_k(k)*yc(k)* psx(xp,yp,zp)
         Fz(x,y,z)=Fz(x,y,z) +t_k(k)*zc(k)* psx(xp,yp,zp)
        endif

   11 continue
! Final wall-fluid interaction
         Sx(x,y,z) = -G1*sum_x *psx(x,y,z) *psx_w
         Sy(x,y,z) = -G1*sum_y *psx(x,y,z) *psx_w
         Sz(x,y,z) = -G1*sum_z *psx(x,y,z) *psx_w
! Final fluid-fluid interaction
         Fx (x,y,z)= -G1 *psx (x,y,z)* Fx(x,y,z)
         Fy (x,y,z)= -G1 *psx (x,y,z)* Fy(x,y,z)
         Fz (x,y,z)= -G1 *psx (x,y,z)* Fz(x,y,z)
      endif

   10 continue
   20 continue
   30 continue
      end

