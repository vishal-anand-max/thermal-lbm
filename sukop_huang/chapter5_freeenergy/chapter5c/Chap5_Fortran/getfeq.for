
      subroutine get_feq(rho, u_x, u_y, fequ, pre)
      implicit none
      include "head.inc"
      integer x,y,yn,yp,xn,xp
      real*8 Fx(lx,ly),Fy(lx,ly),temp(lx,ly),
     &  rho(lx,ly), u_x(lx,ly), u_y(lx,ly),R

      real*8 fequ(0:8,lx,ly),pre(lx,ly),u_n(0:8),rh,u_squ
      real*8 A1, A2, A0, B1, B2, D1, D2, C0, C1, C2, p, Gxx1, Gxx2,
     &  Gyy1, Gyy2, Gxy1, Gxy2
      integer k



       do 15 x = 1, lx
        do 15 y = 1, ly
! Periodic boundary condition
        xp = x+1
        yp = y+1
        xn = x-1
        yn = y-1
        if (xp.gt.lx )  xp = 1
        if (xn.lt.1 )   xn = lx
        if (yp.gt.ly )  yp = 1
        if (yn.lt.1 )   yn = ly
! Calculate Laplacian operator on density
      temp(x,y) =
     *  ( ( rho(xp,y) + rho(xn,y) + rho(x,yp)+ rho(x,yn) )*4.d0/6.d0
     *     +( rho(xp,yp)+ rho(xn,yn)
     *     +  rho(xp,yn)+ rho(xn,yp) )/6.d0
     *     - 20.d0* rho(x,y)/6.d0
     *     )

   15 continue

c----------------------------------------
      do 20 x = 1, lx
        do 20 y = 1, ly
! Periodic boundary condition
        xp = x+1
        yp = y+1
        xn = x-1
        yn = y-1
        if (xp.gt.lx )  xp = 1
        if (xn.lt.1 )   xn = lx
        if (yp.gt.ly )  yp = 1
        if (yn.lt.1 )   yn = ly

! Calculate density gradient. Either simple finite difference or homogeneous
! finite difference can be used.
       Fx(x,y) = (rho(xp,y)- rho(xn,y) )/2.d0
c     * ( (rho(xp,y)- rho(xn,y) )/3.d0
c     *     +(rho(xp,yp)- rho(xn,yn) )/12.d0
c     *     +(rho(xp,yn)- rho(xn,yp) )/12.d0
c     *     )

       Fy(x,y) =  (rho(x,yp)- rho(x,yn) )/2.d0
c     & ( (rho(x,yp)- rho(x,yn) )/3.d0
c     &     +(rho(xp,yp)- rho(xn,yn) )/12.d0
c     &     +(rho(xn,yp)- rho(xp,yn) )/12.d0  )


   20 continue

      do 30 x = 1, lx
        do 30 y = 1, ly

          rh = rho(x,y)

          u_squ = u_x(x,y)*u_x(x,y) + u_y(x,y)*u_y(x,y)
! Get thermodynamic pressure from the equation of state
          call EOS(rh, p)
          pre(x,y) = p

! Calculate the coefficients or terms in the equilibrium distribution function.

         A1 = 1.d0/3.d0 *(p- kappa*rh *temp(x,y) )
         A2 = A1/4.d0
         A0 = rho(x,y) - 5.d0/3.d0*(p- kappa*rh *temp(x,y))
         B1 = rho(x,y)/3.d0
         B2 = B1/4.d0
         D1 = rho(x,y)/2.d0
         D2 = D1/4.d0
         C1 = -rho(x,y)/6.d0
         C2 = C1/4.d0
         C0 = -rho(x,y)*2.d0/3.d0  !!!important
         Gxx1 = kappa/4.d0*(Fx(x,y)*Fx(x,y)- Fy(x,y)*Fy(x,y))
         Gxx2 = Gxx1/4.d0
         Gyy1 = -Gxx1
         Gyy2 = -Gxx2
         Gxy1 = kappa/2.d0 *Fx(x,y)*Fy(x,y)
         Gxy2 = Gxy1/4.d0

        do 65 k = 0, 8
            u_n(k)  = xc(k)*u_x(x,y) + yc(k)*u_y(x,y)
   65 continue

c   Calculate the equilibrium distribution functions.

            fequ(0,x,y) = A0 + C0* u_squ

          fequ(1,x,y) = A1 + B1*u_n(1) +C1*u_squ +D1*u_n(1)*u_n(1)
     &     + Gxx1
          fequ(2,x,y) = A1 + B1*u_n(2) +C1*u_squ +D1*u_n(2)*u_n(2)
     &     + Gyy1
          fequ(3,x,y) = A1 + B1*u_n(3) +C1*u_squ +D1*u_n(3)*u_n(3)
     &     + Gxx1
          fequ(4,x,y) = A1 + B1*u_n(4) +C1*u_squ +D1*u_n(4)*u_n(4)
     &     + Gyy1

          fequ(5,x,y) = A2 + B2*u_n(5) +C2*u_squ +D2*u_n(5)*u_n(5)
     &     + Gxx2 + Gyy2 + 2.d0 * Gxy2
          fequ(6,x,y) = A2 + B2*u_n(6) +C2*u_squ +D2*u_n(6)*u_n(6)
     &     + Gxx2 + Gyy2 - 2.d0 * Gxy2
          fequ(7,x,y) = A2 + B2*u_n(7) +C2*u_squ +D2*u_n(7)*u_n(7)
     &     + Gxx2 + Gyy2 + 2.d0 * Gxy2
          fequ(8,x,y) = A2 + B2*u_n(8) +C2*u_squ +D2*u_n(8)*u_n(8)
     &     + Gxx2 + Gyy2 - 2.d0 * Gxy2

   30 continue

      end
