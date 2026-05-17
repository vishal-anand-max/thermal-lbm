      subroutine g_eq(p,rho,fequ,gequ)
      implicit none
      include "head.inc"
      real*8 u_x,u_y,gequ(0:8),fequ(0:8),u_n(0:8),u_squ,rho,p, tmp1
      integer k	
      tmp1 = rho*c_squ
      do k= 0, 8	
      gequ(k) = (fequ(k)-t_k(k) )*tmp1 +t_k(k)*p
      enddo
      end
c-----------------------------------------
      subroutine f_eq(rho,u_x,u_y,fequ)
      implicit none
      include "head.inc"
      real*8 u_x,u_y,fequ(0:8),u_n(0:8),rho,u_squ, tmp1
      integer k
            u_squ = u_x * u_x + u_y * u_y

            tmp1 =  u_squ *1.5d0
          fequ(0) = t_k(0)* ( -tmp1 + 1.d0)
      do 10 k = 1, 8
            u_n(k) =  xv(k)* u_x + yv(k)* u_y
          fequ(k) = t_k(k)*  ( u_n(k) / c_squ
     &               + u_n(k) * u_n(k) / (con2)
     &               - tmp1 + 1.d0)
   10 continue
      end
c--------------------------------------
