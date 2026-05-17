! Linear map from index function \phi in JCP, 1999 (\rh) to density (\rho)
      subroutine index_function(rh,rho)
      implicit none
      include "head.inc"
      integer x,y,z
      real*8 rh(lx,ly,lz), rho(lx,ly,lz), psi_min0, psi_max0

!  The \psi_{max} and \psi_{min} depend on the EOS we used.

      psi_min0 = psi_min
      psi_max0 = psi_max

      do 1 x = 1,lx
      do 1 y = 1,ly
      do 1 z = 1,lz
      rho(x,y,z) = rho_l + ( rh(x,y,z)- psi_min0 )/(psi_max0 - psi_min0)
     &  * (rho_h - rho_l)
    1 continue
      end


c------------------------------------------------------
! Calculate the equilibrium distribution function of g_i, which is used to get pressure and velocity
      subroutine g_eq(p,rh,u,v,w,gequ)
      implicit none
      include "head.inc"
      real*8 u,v,w,gequ(0:18),u_n(0:18),u_squ,rh,p
      integer k
            u_squ = u*u + v*v +w*w
       do 5 k = 0, 18
            u_n(k)  = xc(k)*u + yc(k)*v +zc(k)*w

            gequ(k) = t_k(k)* ( p + rh* c_squ* (
     &                   u_n(k) / c_squ
     &               + u_n(k) *u_n(k) / (2.d0 * c_squ *c_squ)
     &               - u_squ / (2.d0 * c_squ)
     *              ) )
    5 continue
      end
cc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Calculate the equilibrium distribution function f_i, which is used to get index function

      subroutine f_eq(rh,u,v,w,fequ)
      implicit none
      include "head.inc"
      real*8 u,v,w,fequ(0:18),u_n(0:18),rh,u_squ
      integer k
            u_squ = u*u + v*v +w*w
        do 65 k = 0, 18
            u_n(k)  = xc(k)*u + yc(k)*v +zc(k)*w

            fequ(k) = t_k(k)* ( 1.0d0+  u_n(k) / c_squ
     &               + u_n(k) *u_n(k) / (2.d0 * c_squ *c_squ)
     &               - u_squ / (2.d0 * c_squ))
   65 continue
      end