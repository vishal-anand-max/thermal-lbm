
        integer xA,yA
        parameter(xA=128,yA=512)
        real*8 Rref(3), h_ref
        PARAMETER(h_ref=128.d0)

        common/AA/ Nwri,  BD, RAYLEIGH, opp(8)
        integer Nwri,  BD, RAYLEIGH, opp

        common/b/ error,vel,tau_g,tau_f,xv(0:8),yv(0:8),t_k(0:8)
        real*8 error,vel,tau_g, tau_f,xv,yv,t_k

        common/app/ limi,lam,gam,t_0,t_1,t_2,c_squ,visc,pai, con2
        real*8 limi,lam,gam,t_0,t_1,t_2,c_squ,visc,pai, con2

        common /shp/  rho_h, rho_l,psi_max,psi_min, const, Tn
        real*8 rho_h, rho_l,psi_max,psi_min, const, Tn

        common/psi/ psi_min0, psi_max0
        real*8  psi_min0, psi_max0
        common /velo/ Kappa, gforce,  We, R, Eo
        real*8  Kappa, gforce,  We, R , Eo

        common /nu/  dfai_x(xA,yA),dfai_y(xA,yA)
        real*8 dfai_x,dfai_y
        common /prho/ dprho_x(xA,yA), dprho_y(xA,yA)
        real*8 dprho_x,dprho_y


